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
#include <chrono>
#include <unordered_set>
//#include "h5_helpers.hpp"  // no longer used here
#ifdef _OPENMP
#include <omp.h>
#endif

#include "RPgenReader.h"    // your RPgenReader class

#include <fstream>
#include <cstdint>

using namespace Rcpp;
using std::size_t;

static inline double wall_seconds_now() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  using clock = std::chrono::steady_clock;
  static const clock::time_point t0 = clock::now();
  return std::chrono::duration<double>(clock::now() - t0).count();
#endif
}

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

  TileableDense Gc(n_inds, n_snps);
  reader.ReadList(Gc.data, vsubset_std, meanimpute, is_round);

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


struct SnpChunkView {
  size_t global_snp_start;
  size_t n_rows;
  size_t n_cols;
  const double* base;
  size_t ld;

  const double* col_ptr(size_t local_snp) const {
    return base + local_snp * ld;
  }
};

struct WaveOp {
  size_t local_snp;
  size_t gene;
  double weight;
  int flag;
  double* y;
};

struct EngineTiming {
  double pgen_load = 0.0;
  double wave_prepare = 0.0;
  double parallel_compute = 0.0;
  double output_flush = 0.0;
  size_t pgen_chunks = 0;
  size_t compute_waves = 0;
  size_t operations = 0;
  size_t max_wave_ops = 0;
  size_t max_live_gene_columns = 0;
};

static inline int available_omp_workers(int requested_threads) {
#ifdef _OPENMP
  const int max_threads = std::max(1, omp_get_max_threads());
  if (requested_threads <= 0) {
    return max_threads;
  }
  return std::max(1, std::min(requested_threads, max_threads));
#else
  return 1;
#endif
}

static inline size_t choose_row_tile_count(size_t n_rows, int workers) {
  if (n_rows == 0) return 0;
  const size_t target = static_cast<size_t>(std::max(1, workers)) * 4u;
  return std::max<size_t>(1, std::min(n_rows, target));
}

static inline void flush_if_needed(SoA& S,
                                   RawWriter& writer,
                                   std::vector<size_t>& gid_toExport,
                                   size_t exportChunk) {
  if (gid_toExport.size() >= exportChunk) {
    flush_gid_batch(S, writer, gid_toExport);
  }
}

static void process_decoded_snp_chunk(
    const SnpChunkView& chunk,
    const int* ip,
    const int* idx,
    const double* w,
    const int* gp,
    SoA& S,
    RawWriter& writer,
    std::vector<size_t>& gid_toExport,
    size_t exportChunk,
    int requested_threads,
    EngineTiming& timing)
{
  if (chunk.n_rows == 0 || chunk.n_cols == 0) return;

  const int workers = available_omp_workers(requested_threads);
  const size_t row_tiles = choose_row_tile_count(chunk.n_rows, workers);
  const size_t max_wave_unique_genes = std::max<size_t>(1, exportChunk);
  const size_t max_wave_ops = std::max<size_t>(max_wave_unique_genes, static_cast<size_t>(workers) * 4u);

  size_t snp = chunk.global_snp_start;
  const size_t snp_end = chunk.global_snp_start + chunk.n_cols;
  size_t edge = static_cast<size_t>(ip[snp]);

  while (snp < snp_end) {
    while (snp < snp_end && edge >= static_cast<size_t>(ip[snp + 1])) {
      ++snp;
      if (snp < snp_end) edge = static_cast<size_t>(ip[snp]);
    }
    if (snp >= snp_end) break;

    const double prep_start = wall_seconds_now();
    std::vector<WaveOp> ops;
    ops.reserve(max_wave_ops);
    std::unordered_set<size_t> wave_genes;
    wave_genes.reserve(max_wave_unique_genes * 2u);

    while (snp < snp_end) {
      const size_t edge_end = static_cast<size_t>(ip[snp + 1]);
      if (edge >= edge_end) {
        ++snp;
        if (snp < snp_end) edge = static_cast<size_t>(ip[snp]);
        continue;
      }

      const size_t gene = static_cast<size_t>(idx[edge]);
      const bool new_wave_gene = wave_genes.find(gene) == wave_genes.end();
      if (!ops.empty() &&
          (ops.size() >= max_wave_ops ||
           (new_wave_gene && wave_genes.size() >= max_wave_unique_genes))) {
        break;
      }

      wave_genes.insert(gene);
      WaveOp op;
      op.local_snp = snp - chunk.global_snp_start;
      op.gene = gene;
      op.weight = w[edge];
      op.flag = gp[edge];
      op.y = S.ensure_alloc(gene);
      ops.push_back(op);
      ++edge;
    }

    timing.wave_prepare += wall_seconds_now() - prep_start;
    if (ops.empty()) continue;

    const double compute_start = wall_seconds_now();
    const int row_tiles_int = static_cast<int>(row_tiles);
    const int thread_count = std::max(1, std::min(workers, row_tiles_int));
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(thread_count)
#endif
    for (int tile = 0; tile < row_tiles_int; ++tile) {
      const size_t r0 = (chunk.n_rows * static_cast<size_t>(tile)) / row_tiles;
      const size_t r1 = (chunk.n_rows * static_cast<size_t>(tile + 1)) / row_tiles;
      for (const WaveOp& op : ops) {
        const double* x = chunk.col_ptr(op.local_snp);
        double* y = op.y;
        const double weight = op.weight;
        for (size_t r = r0; r < r1; ++r) {
          y[r] += weight * x[r];
        }
      }
    }
    timing.parallel_compute += wall_seconds_now() - compute_start;

    const double flush_start = wall_seconds_now();
    for (const WaveOp& op : ops) {
      if (op.flag == 0 || op.flag == 2) {
        gid_toExport.push_back(op.gene);
        flush_if_needed(S, writer, gid_toExport, exportChunk);
      }
    }
    timing.output_flush += wall_seconds_now() - flush_start;

    ++timing.compute_waves;
    timing.operations += ops.size();
    timing.max_wave_ops = std::max(timing.max_wave_ops, ops.size());
    timing.max_live_gene_columns = std::max(timing.max_live_gene_columns, S.col.size());
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
    bool rounded_mean = false,          // after mean imputing it rounds to nearest integer when true
    int threads = 0                     // <= 0 uses omp_get_max_threads()
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

  size_t n_inds = 0;
  bool using_dense = false;
  Rcpp::NumericMatrix G_dense;
  RPgenReader pgen_reader;

  if (G_in.isNotNull()) {
    G_dense = G_in.get();
    if (G_dense.ncol() != K) {
      Rcpp::stop("grex_axpy_all_snps_to_dense: K mismatch between indptr and G.");
    }
    using_dense = true;
    n_inds = static_cast<size_t>(G_dense.nrow());
  } else {
    if (pgen_file.empty()) {
      Rcpp::stop("grex_axpy_all_snps_to_dense: either G_in or pgen_file must be provided.");
    }

    std::vector<int> subset_std;
    subset_std.reserve(sample_subset.size());
    for (int i = 0; i < sample_subset.size(); ++i) {
      subset_std.push_back(sample_subset[i]);
    }
    pgen_reader.Load(pgen_file, nullptr, raw_sample_ct, subset_std);
    n_inds = static_cast<size_t>(pgen_reader.GetSubsetSize());
  }

  SoA S(n_inds);

  // Create RAW BINARY matrix [n_inds x n_total_genes] only in first call
  if (create_new) {
    create_empty_raw_matrix(grom_file, n_inds, n_cols_dense);
  }

  // Open existing file for in-place writes
  RawWriter writer(grom_file, n_inds, n_cols_dense);

  std::vector<size_t> gid_toExport;
  gid_toExport.reserve(exportChunk);

  if (exportChunk == 0) {
    Rcpp::stop("exportChunk must be positive.");
  }
  if (snp_chunk == 0) {
    Rcpp::stop("snp_chunk must be positive.");
  }

  EngineTiming timing;
  const size_t snp_chunk_safe = std::max<size_t>(1, snp_chunk);

  for (size_t chunk_start = 0; chunk_start < static_cast<size_t>(K); chunk_start += snp_chunk_safe) {
    const size_t chunk_end = std::min(static_cast<size_t>(K), chunk_start + snp_chunk_safe);
    const size_t n_chunk_snps = chunk_end - chunk_start;
    SnpChunkView chunk;
    chunk.global_snp_start = chunk_start;
    chunk.n_rows = n_inds;
    chunk.n_cols = n_chunk_snps;
    chunk.ld = n_inds;

    TileableDense pgen_chunk;
    if (using_dense) {
      chunk.base = &G_dense(0, static_cast<int>(chunk_start));
    } else {
      Rcpp::IntegerVector vsubset(static_cast<int>(n_chunk_snps));
      for (size_t j = 0; j < n_chunk_snps; ++j) {
        vsubset[static_cast<int>(j)] = static_cast<int>(chunk_start + j + 1);
      }
      const double load_start = wall_seconds_now();
      pgen_chunk = stream_pgen_snps_cpp(
        pgen_reader,
        vsubset,
        meanimpute,
        rounded_mean
      );
      timing.pgen_load += wall_seconds_now() - load_start;
      if (pgen_chunk.n_rows != n_inds) {
        throw std::runtime_error("grom_axpy_engine: inconsistent n_inds across PGEN chunks");
      }
      chunk.base = pgen_chunk.data.data();
    }

    process_decoded_snp_chunk(
      chunk, ip, idx, w, gp, S, writer, gid_toExport, exportChunk, threads, timing
    );
    ++timing.pgen_chunks;
  }

  // Flush any leftover finished genes
  {
    const double flush_start = wall_seconds_now();
    flush_gid_batch(S, writer, gid_toExport);
    timing.output_flush += wall_seconds_now() - flush_start;
  }

  // Safety net: export any remaining columns still in SoA
  if (!S.col.empty()) {
    const double flush_start = wall_seconds_now();
    for (const auto& kv : S.col) {
      size_t gid = kv.first;
      double* col_data = kv.second.get();
      writer.write_rows_chunk(gid, 0, S.n_rows, col_data);
    }
    S.col.clear();
    timing.output_flush += wall_seconds_now() - flush_start;
  }

  Rcpp::Rcout << "grom_axpy_engine timing: "
              << "chunks=" << timing.pgen_chunks
              << " waves=" << timing.compute_waves
              << " ops=" << timing.operations
              << " pgen_load=" << timing.pgen_load
              << "s wave_prepare=" << timing.wave_prepare
              << "s parallel_compute=" << timing.parallel_compute
              << "s output_flush=" << timing.output_flush
              << "s max_wave_ops=" << timing.max_wave_ops
              << " max_live_gene_columns=" << timing.max_live_gene_columns
              << " omp_workers=" << available_omp_workers(threads)
              << "\n";

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
