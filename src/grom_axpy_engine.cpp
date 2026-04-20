// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <unordered_map>
#include <memory>
#include <vector>
#include <cstddef>
#include <algorithm>        // std::max_element, std::min
#include <cstring>          // std::memcpy
#include <stdexcept>        // std::runtime_error
//#include "h5_helpers.hpp"  // no longer used here
#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include <R_ext/BLAS.h>
}
#include "RPgenReader.h"    // your RPgenReader class

#include <fstream>
#include <cstdint>

using namespace Rcpp;
using std::size_t;

// ====================================================
// Raw binary matrix helpers (column-major)
// ====================================================

static inline void create_empty_raw_matrix(
    const std::string& path,
    size_t n_rows,
    size_t n_cols)
{
  std::uint64_t total_elems =
    static_cast<std::uint64_t>(n_rows) *
    static_cast<std::uint64_t>(n_cols);
  std::uint64_t total_bytes =
    total_elems * static_cast<std::uint64_t>(sizeof(double));

  std::ofstream f(path, std::ios::binary | std::ios::trunc);
  if (!f) {
    throw std::runtime_error("create_empty_raw_matrix: cannot create file");
  }

  if (total_bytes > 0) {
    // Pre-allocate file size by seeking and writing one byte
    f.seekp(static_cast<std::streamoff>(total_bytes - 1));
    char zero = 0;
    f.write(&zero, 1);
  }
}

struct RawWriter {
  std::string path;
  size_t n_rows;
  size_t n_cols;
  std::ofstream f;

  RawWriter(const std::string& path_,
            size_t n_rows_,
            size_t n_cols_)
    : path(path_), n_rows(n_rows_), n_cols(n_cols_),
      f(path_, std::ios::binary | std::ios::in | std::ios::out)
  {
    if (!f) {
      throw std::runtime_error("RawWriter: cannot open file for writing");
    }
  }

  // Write a contiguous chunk of rows [r0 .. r0 + n_chunk - 1] for column col_idx
  void write_rows_chunk(size_t col_idx,
                        size_t r0,
                        size_t n_chunk,
                        const double* src)
  {
    if (col_idx >= n_cols) return;
    if (n_chunk == 0) return;
    if (r0 >= n_rows) return;

    if (r0 + n_chunk > n_rows) {
      n_chunk = n_rows - r0;
    }
    if (n_chunk == 0) return;

    std::uint64_t offset_elems =
      static_cast<std::uint64_t>(col_idx) *
      static_cast<std::uint64_t>(n_rows) +
      static_cast<std::uint64_t>(r0);

    std::uint64_t offset_bytes =
      offset_elems * static_cast<std::uint64_t>(sizeof(double));

    f.seekp(static_cast<std::streamoff>(offset_bytes));
    f.write(reinterpret_cast<const char*>(src),
            static_cast<std::streamsize>(n_chunk * sizeof(double)));
    if (!f) {
      throw std::runtime_error("RawWriter: write failed");
    }
  }
};

// ----------------------
// Tileable column-major dense structure
// ----------------------
struct TileableDense {
  std::size_t n_rows;           // individuals (subset size)
  std::size_t n_cols;           // SNPs (length of variant_subset)
  std::vector<double> data;     // column-major: col j starts at data[j * n_rows]

  TileableDense() : n_rows(0), n_cols(0), data() {}

  TileableDense(std::size_t rows, std::size_t cols)
    : n_rows(rows), n_cols(cols), data(rows * cols) {}

  double* col_ptr(std::size_t j) {
    return data.data() + j * n_rows;
  }

  const double* col_ptr(std::size_t j) const {
    return data.data() + j * n_rows;
  }

  double& operator()(std::size_t i, std::size_t j) {
    return data[j * n_rows + i];
  }

  const double& operator()(std::size_t i, std::size_t j) const {
    return data[j * n_rows + i];
  }
};

// ----------------------
// Stream a SNP chunk from pgen into TileableDense using RPgenReader::ReadList
// ----------------------
TileableDense stream_pgen_snps_cpp(
    RPgenReader& reader,
    const Rcpp::IntegerVector& vsubset, // 1-based variant indices for this chunk
    bool meanimpute,
    bool is_round
) {
  // n_inds = subset size after Load()
  std::size_t n_inds  = static_cast<std::size_t>(reader.GetSubsetSize());
  std::size_t n_snps  = static_cast<std::size_t>(vsubset.size());

  // Convert IntegerVector -> std::vector<int> (1-based indices preserved)
  std::vector<int> vsubset_std;
  vsubset_std.reserve(n_snps);
  for (std::size_t j = 0; j < n_snps; ++j) {
    vsubset_std.push_back(vsubset[static_cast<int>(j)]);
  }

  // Buffer for ReadList: size = n_inds * n_snps, column-major in variant order
  std::vector<double> buf;
  buf.resize(n_inds * n_snps);

  reader.ReadList(buf, vsubset_std, meanimpute, is_round);

  // Wrap into TileableDense (copy; can be optimized later if needed)
  TileableDense Gc(n_inds, n_snps);
  std::copy(buf.begin(), buf.end(), Gc.data.begin());

  return Gc;
}


// -------------------------
// Structure of Arrays (SoA)
// -------------------------
struct SoA {
  size_t n_rows;
  std::unordered_map<size_t, std::unique_ptr<double[]>> col; // sparse cols

  explicit SoA(size_t n_rows_) : n_rows(n_rows_) {}

  double* ensure_alloc(size_t j) {
    auto& ptr = col[j];
    if (!ptr) ptr = std::make_unique<double[]>(n_rows); // value-init → 0
    return ptr.get();
  }
};


// ----------------------------------------------------
// AXPY kernel for a single gene *row chunk*
// ----------------------------------------------------
static inline void grex_axpy_for_snp(
    SoA& S,
    size_t r0,
    size_t n_chunk,
    size_t gene,
    double weight,
    const double* chunk_src
) {
  if (n_chunk == 0) return;

  double* y = S.ensure_alloc(gene) + r0;
  int n = static_cast<int>(n_chunk);
  int inc = 1;

  F77_CALL(daxpy)(
    &n,
    &weight,
    chunk_src, &inc,
    y,         &inc
  );
}


// SoA → dense R matrix [n_rows x n_cols], column-major.
// (kept for possible debugging; not used in main pipeline)
static inline Rcpp::NumericMatrix soa_to_dense(const SoA& S,
                                               size_t n_rows,
                                               size_t n_cols) {
  Rcpp::NumericMatrix M(n_rows, n_cols); // zero-initialized
  for (const auto& kv : S.col) {
    size_t j = kv.first;
    if (j >= n_cols) continue;
    const double* src = kv.second.get();
    double* dst_col = &M(0, j);
    std::memcpy(dst_col, src, n_rows * sizeof(double));
  }
  return M;
}

// Flush a batch of genes to RAW BINARY, then clear them from SoA
static inline void flush_gid_batch(SoA& S,
                                   RawWriter& writer,
                                   std::vector<size_t>& gid_toExport) {
  if (gid_toExport.empty()) return;

  for (size_t gid : gid_toExport) {
    auto it = S.col.find(gid);
    if (it != S.col.end()) {
      double* col_data = it->second.get();
      // Export full column [0 .. n_rows-1] for this gene index
      writer.write_rows_chunk(gid, 0, S.n_rows, col_data);
      S.col.erase(it);
    }
  }
  gid_toExport.clear();
}


// ----------------------------------------------------
// Unified per-SNP processing  (PARALLEL OVER ROW-CHUNKS)
// ----------------------------------------------------
static inline void process_one_snp(
    int snp,                         // global SNP index [0..K-1]
    const double* col,               // length n_inds, dosage vector
    size_t n_inds,
    size_t CHUNK,
    const int* ip,
    const int* idx,
    const double* w,
    const int* gp,
    size_t n_chunks_rows,
    SoA& S,
    std::vector<size_t>& gid_toExport,
    size_t exportChunk,
    RawWriter& writer,
    bool showWarnings
) {
  const size_t p_start = static_cast<size_t>(ip[snp]);
  const size_t p_end   = static_cast<size_t>(ip[snp + 1]);
  if (p_start == p_end) return;  // no genes for this SNP

  for (size_t p = p_start; p < p_end; ++p) {
    size_t gene   = static_cast<size_t>(idx[p]);
    double weight = w[p];
    int    flag   = gp[p];

    // --- PARALLELIZE OVER ROW-CHUNKS FOR THIS (snp, gene) ---
#ifdef _OPENMP
    double* y_base = S.ensure_alloc(gene);

    int n_chunks_int = static_cast<int>(n_chunks_rows);

    #pragma omp parallel for schedule(static)
    for (int chunk_id = 0; chunk_id < n_chunks_int; ++chunk_id) {
      size_t r0 = static_cast<size_t>(chunk_id) * CHUNK;
      if (r0 >= n_inds) continue;
      size_t n_chunk_rows = std::min(CHUNK, n_inds - r0);

      const double* chunk_src = col + r0;
      double*       y        = y_base + r0;

      // Optional debug
      if (showWarnings && chunk_id == 0) {
        Rcpp::Rcout << "snp " << snp
                    << " gene " << gene
                    << " using " << omp_get_num_threads()
                    << " threads\n";
      }

      int n = static_cast<int>(n_chunk_rows);
      int inc = 1;
      F77_CALL(daxpy)(
        &n,
        &weight,
        chunk_src, &inc,
        y,         &inc
      );
    }
#else
    // --- SERIAL FALLBACK (original behavior using helper) ---
    for (size_t chunk_id = 0; chunk_id < n_chunks_rows; ++chunk_id) {
      size_t r0 = chunk_id * CHUNK;
      if (r0 >= n_inds) break;
      size_t n_chunk_rows = std::min(CHUNK, n_inds - r0);

      const double* chunk_src = col + r0;
      grex_axpy_for_snp(S, r0, n_chunk_rows, gene, weight, chunk_src);
    }
#endif

    // Finished gene → mark for export (still serial, no races)
    if (flag == 0 || flag == 2) {
      gid_toExport.push_back(gene);
      if (gid_toExport.size() >= exportChunk) {
        flush_gid_batch(S, writer, gid_toExport);
      }
    }
  }
}


// ----------------------------------------------------
// Abstract SNP source interface
// ----------------------------------------------------
struct ISnpSource {
  virtual ~ISnpSource() = default;
  virtual size_t n_inds() const = 0;
  virtual int    K() const = 0;
  virtual const double* col_ptr(int snp) = 0;  // pointer valid until next re-chunk
};


// Dense matrix implementation
struct DenseSnpSource : public ISnpSource {
  Rcpp::NumericMatrix G;
  size_t n_rows;
  int K_;

  explicit DenseSnpSource(const Rcpp::NumericMatrix& G_)
    : G(G_),
      n_rows(static_cast<size_t>(G_.nrow())),
      K_(G_.ncol())
  {}

  size_t n_inds() const override { return n_rows; }
  int    K() const override { return K_; }

  const double* col_ptr(int snp) override {
    if (snp < 0 || snp >= K_) {
      throw std::runtime_error("DenseSnpSource::col_ptr: snp index out of range");
    }
    return &G(0, snp);  // column-major
  }
};


// PGEN implementation with chunked streaming via RPgenReader
struct PgenSnpSource : public ISnpSource {
  RPgenReader reader;        // default-constructed, then Load()'d
  std::string pgen_file_;
  int raw_sample_ct_;
  Rcpp::IntegerVector sample_subset_;
  std::size_t n_inds_;
  int K_;
  std::size_t snp_chunk_;
  bool meanimpute_;
  bool rounded_mean_;

  std::size_t   chunk_start;
  std::size_t   chunk_end;
  TileableDense Gc;

  PgenSnpSource(const std::string& pgen_file,
                int raw_sample_ct,
                const Rcpp::IntegerVector& sample_subset,
                bool meanimpute,
                bool rounded_mean,
                std::size_t snp_chunk,
                int K)
    : reader(),
      pgen_file_(pgen_file),
      raw_sample_ct_(raw_sample_ct),
      sample_subset_(sample_subset),
      n_inds_(0),
      K_(K),
      snp_chunk_(snp_chunk),
      meanimpute_(meanimpute),
      rounded_mean_(rounded_mean),
      chunk_start(0),
      chunk_end(0),
      Gc()
  {
    // Convert sample_subset (Rcpp) to std::vector<int> (1-based indices)
    std::vector<int> subset_std;
    subset_std.reserve(sample_subset_.size());
    for (int i = 0; i < sample_subset_.size(); ++i) {
      subset_std.push_back(sample_subset_[i]);
    }

    // rp == nullptr (we don't need pvar info here)
    reader.Load(pgen_file_, nullptr, raw_sample_ct_, subset_std);

    // After Load(), subset size is defined (either subset or all samples)
    n_inds_ = static_cast<std::size_t>(reader.GetSubsetSize());
  }

  std::size_t n_inds() const override { return n_inds_; }
  int         K() const override { return K_; }

  const double* col_ptr(int snp) override {
    if (snp < 0 || snp >= K_) {
      throw std::runtime_error("PgenSnpSource::col_ptr: snp index out of range");
    }
    std::size_t s = static_cast<std::size_t>(snp);

    // If snp is not in the current chunk, load the chunk that contains it
    if (s < chunk_start || s >= chunk_end) {
      chunk_start = (s / snp_chunk_) * snp_chunk_;
      chunk_end   = std::min(chunk_start + snp_chunk_,
                             static_cast<std::size_t>(K_));

      int len = static_cast<int>(chunk_end - chunk_start);
      Rcpp::IntegerVector vsubset(len);
      for (int j = 0; j < len; ++j) {
        vsubset[j] = static_cast<int>(chunk_start) + j + 1; // 1-based
      }

      // Stream this chunk from pgen into a tile
      Gc = stream_pgen_snps_cpp(
        reader,
        vsubset,
        meanimpute_,
        rounded_mean_
      );

      if (Gc.n_rows != n_inds_) {
        throw std::runtime_error("PgenSnpSource: inconsistent n_inds across chunks");
      }
    }

    std::size_t local = s - chunk_start;
    return Gc.col_ptr(local);
  }
};


// ----------------------------------------------------
// Main exported function (unified dense + pgen logic)
// Exports to a raw binary [n_inds x n_total_genes] matrix
// stored at grom_file in column-major order.
// ----------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix grom_axpy_engine(
    Rcpp::IntegerVector indptr,    // length K+1
    Rcpp::IntegerVector indices,   // length nnz
    Rcpp::NumericVector data,      // length nnz (weights)
    const std::string& pgen_file,  // path to .pgen ("" if using G_in)
    int raw_sample_ct,             // total samples in the pgen
    size_t snp_chunk,              // SNP chunk size (used by PGEN source)
    Rcpp::IntegerVector sample_subset, // 1-based samples; integer(0) → all
    bool meanimpute,               // mean-impute missing dosages? (PGEN)
    const std::string& grom_file,    // path to RAW BINARY file
    size_t CHUNK,                  // row-chunk size
    Rcpp::IntegerVector gene_pos,  // length nnz: 1,0,2,-1 flags
    size_t exportChunk,            // how many finished genes to batch before export
    size_t n_total_genes,          // NEW: total mg_id columns across ALL chromosomes
    bool create_new,               // NEW: TRUE only in first call
    Rcpp::Nullable<Rcpp::NumericMatrix> G_in = R_NilValue, // optional precomputed G
    bool showWarnings = false,
    bool rounded_mean = false           // after mean imputing it rounds to nearest integer when true
) {
  const int*    ip  = indptr.begin();
  const int*    idx = indices.begin();
  const double* w   = data.begin();
  const int*    gp  = gene_pos.begin();

  const int K = indptr.size() - 1;

  // Sanity check: local max idx must fit in global n_total_genes
  int max_local_gene = -1;
  if (indices.size() > 0) {
    max_local_gene = *std::max_element(indices.begin(), indices.end());
  }
  if (max_local_gene >= static_cast<int>(n_total_genes)) {
    Rcpp::stop("grex_axpy_all_snps_to_dense: local gene index exceeds n_total_genes.");
  }

  const size_t n_cols_dense = n_total_genes;

  // No SNPs: optionally create empty matrix only on first call
  if (K <= 0) {
    if (create_new) {
      create_empty_raw_matrix(grom_file, 0, n_cols_dense);
    }
    return Rcpp::NumericMatrix(0, n_cols_dense);
  }

  // Decide source type ONCE, then shared logic
  std::unique_ptr<ISnpSource> src;

  if (G_in.isNotNull()) {
    Rcpp::NumericMatrix G = G_in.get();
    if (G.ncol() != K) {
      Rcpp::stop("grex_axpy_all_snps_to_dense: K mismatch between indptr and G.");
    }
    src = std::make_unique<DenseSnpSource>(G);
  } else {
    if (pgen_file.empty()) {
      Rcpp::stop("grex_axpy_all_snps_to_dense: either G_in or pgen_file must be provided.");
    }
    src = std::make_unique<PgenSnpSource>(
      pgen_file,
      raw_sample_ct,
      sample_subset,
      meanimpute,
      rounded_mean,
      snp_chunk,
      K
    );
  }

  size_t n_inds = src->n_inds();
  SoA S(n_inds);

  // Create RAW BINARY matrix [n_inds x n_total_genes] only in first call
  if (create_new) {
    create_empty_raw_matrix(grom_file, n_inds, n_cols_dense);
  }

  // Open existing file for in-place writes
  RawWriter writer(grom_file, n_inds, n_cols_dense);

  std::vector<size_t> gid_toExport;
  gid_toExport.reserve(exportChunk);

  const size_t n_chunks_rows = (n_inds + CHUNK - 1) / CHUNK;

  // Unified per-SNP loop
  const int K_src = src->K();
  if (K_src != K) {
    Rcpp::stop("grex_axpy_all_snps_to_dense: K mismatch between indptr and SNP source.");
  }

  for (int snp = 0; snp < K; ++snp) {
    const double* col = src->col_ptr(snp);
    process_one_snp(
      snp, col,
      n_inds, CHUNK,
      ip, idx, w, gp,
      n_chunks_rows,
      S, gid_toExport, exportChunk,
      writer,
      showWarnings
    );
  }

  // Flush any leftover finished genes
  flush_gid_batch(S, writer, gid_toExport);

  // Safety net: export any remaining columns still in SoA
  if (!S.col.empty()) {
    for (const auto& kv : S.col) {
      size_t gid = kv.first;
      double* col_data = kv.second.get();
      writer.write_rows_chunk(gid, 0, S.n_rows, col_data);
    }
    S.col.clear();
  }

  // We don't return the big matrix to R (it lives on disk)
  return Rcpp::NumericMatrix(0, 0);
}


// [[Rcpp::export]]
Rcpp::IntegerVector grex_omp_info() {
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  int num_procs   = omp_get_num_procs();
  return Rcpp::IntegerVector::create(
    _OPENMP,       // OpenMP version (e.g. 201511)
    max_threads,   // max threads OpenMP plans to use
    num_procs      // number of HW cores it sees
  );
#else
  // -1 means "no OpenMP in this compilation"
  return Rcpp::IntegerVector::create(-1, -1, -1);
#endif
}
