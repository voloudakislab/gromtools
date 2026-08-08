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
#include <exception>
#include <sstream>
#include <limits>
#include <cerrno>
#include <climits>
#include <string>
//#include "h5_helpers.hpp"  // no longer used here
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_WIN32)
#define GROMTOOLS_HAS_POSIX_PWRITE 0
#else
#define GROMTOOLS_HAS_POSIX_PWRITE 1
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

#if GROMTOOLS_HAS_POSIX_PWRITE
using RawFileOffset = off_t;
#else
using RawFileOffset = long long;
#endif

#include "RPgenReader.h"    // your RPgenReader class

#include <fstream>
#include <cstdint>

using namespace Rcpp;
using std::size_t;

static size_t checked_multiply_size_t(size_t lhs, size_t rhs, const char* context) {
  if (lhs != 0 && rhs > std::numeric_limits<size_t>::max() / lhs) {
    throw std::runtime_error(std::string(context) + ": size overflow");
  }
  return lhs * rhs;
}

static size_t checked_add_size_t(size_t lhs, size_t rhs, const char* context) {
  if (rhs > std::numeric_limits<size_t>::max() - lhs) {
    throw std::runtime_error(std::string(context) + ": size overflow");
  }
  return lhs + rhs;
}

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
  const size_t total_elems = checked_multiply_size_t(
    n_rows, n_cols, "create_empty_raw_matrix elements"
  );
  const size_t total_bytes = checked_multiply_size_t(
    total_elems, sizeof(double), "create_empty_raw_matrix bytes"
  );

  std::ofstream f(path, std::ios::binary | std::ios::trunc);
  if (!f) {
    throw std::runtime_error("create_empty_raw_matrix: cannot create file");
  }

  if (total_bytes > 0) {
    // Pre-allocate file size by seeking and writing one byte
    if (total_bytes - 1 > static_cast<size_t>(std::numeric_limits<std::streamoff>::max())) {
      throw std::runtime_error("create_empty_raw_matrix: file is too large for stream offset");
    }
    f.seekp(static_cast<std::streamoff>(total_bytes - 1));
    char zero = 0;
    f.write(&zero, 1);
  }
}

struct RawWriter {
  std::string path;
  size_t n_rows;
  size_t n_cols;
#if GROMTOOLS_HAS_POSIX_PWRITE
  int fd;
#else
  std::ofstream f;
#endif

  RawWriter(const std::string& path_,
            size_t n_rows_,
            size_t n_cols_)
    : path(path_), n_rows(n_rows_), n_cols(n_cols_)
#if GROMTOOLS_HAS_POSIX_PWRITE
      , fd(-1)
#else
      , f(path_, std::ios::binary | std::ios::in | std::ios::out)
#endif
  {
#if GROMTOOLS_HAS_POSIX_PWRITE
#ifdef O_CLOEXEC
    fd = ::open(path.c_str(), O_RDWR | O_CLOEXEC);
#else
    fd = ::open(path.c_str(), O_RDWR);
#endif
    if (fd < 0) {
      const int saved_errno = errno;
      throw std::runtime_error("RawWriter: cannot open file for positioned writes: " +
                               std::string(std::strerror(saved_errno)));
    }
    try {
      validate_file_size();
    } catch (...) {
      ::close(fd);
      fd = -1;
      throw;
    }
#else
    if (!f) {
      throw std::runtime_error("RawWriter: cannot open file for writing");
    }
    validate_file_size();
#endif
  }

  RawWriter(const RawWriter&) = delete;
  RawWriter& operator=(const RawWriter&) = delete;
  RawWriter(RawWriter&&) = delete;
  RawWriter& operator=(RawWriter&&) = delete;

  ~RawWriter() {
#if GROMTOOLS_HAS_POSIX_PWRITE
    if (fd >= 0) {
      ::close(fd);
    }
#endif
  }

  bool supports_parallel_positioned_writes() const {
    return GROMTOOLS_HAS_POSIX_PWRITE != 0;
  }

  size_t expected_file_bytes() const {
    return checked_multiply_size_t(
      checked_multiply_size_t(n_rows, n_cols, "RawWriter expected file elements"),
      sizeof(double),
      "RawWriter expected file bytes"
    );
  }

  void validate_file_size() const {
    const size_t expected = expected_file_bytes();
#if GROMTOOLS_HAS_POSIX_PWRITE
    struct stat st;
    if (::fstat(fd, &st) != 0) {
      const int saved_errno = errno;
      throw std::runtime_error("RawWriter: fstat failed: " +
                               std::string(std::strerror(saved_errno)));
    }
    if (st.st_size < 0) {
      throw std::runtime_error("RawWriter: negative file size");
    }
    if (static_cast<uintmax_t>(st.st_size) != static_cast<uintmax_t>(expected)) {
      throw std::runtime_error("RawWriter: file size does not match expected matrix dimensions");
    }
#else
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) {
      throw std::runtime_error("RawWriter: cannot stat file");
    }
    const std::streampos observed = in.tellg();
    if (observed < 0 ||
        static_cast<uintmax_t>(observed) != static_cast<uintmax_t>(expected)) {
      throw std::runtime_error("RawWriter: file size does not match expected matrix dimensions");
    }
#endif
  }

  size_t column_offset_bytes(size_t col_idx, size_t r0 = 0) const {
    if (col_idx >= n_cols) {
      throw std::runtime_error("RawWriter: column index out of bounds");
    }
    if (r0 > n_rows) {
      throw std::runtime_error("RawWriter: row offset out of bounds");
    }
    return checked_multiply_size_t(
      checked_add_size_t(
        checked_multiply_size_t(col_idx, n_rows, "RawWriter offset elements"),
        r0,
        "RawWriter offset elements"
      ),
      sizeof(double),
      "RawWriter offset bytes"
    );
  }

  size_t column_byte_count(size_t n_chunk) const {
    return checked_multiply_size_t(n_chunk, sizeof(double), "RawWriter write bytes");
  }

  RawFileOffset checked_file_offset(size_t value, const char* context) const {
#if GROMTOOLS_HAS_POSIX_PWRITE
    if (value > static_cast<size_t>(std::numeric_limits<off_t>::max())) {
      throw std::runtime_error(std::string(context) + ": offset exceeds off_t range");
    }
    return static_cast<off_t>(value);
#else
    (void)context;
    return static_cast<off_t>(value);
#endif
  }

  void write_at(const void* src, size_t bytes, size_t offset_bytes, const std::string& context) {
    if (bytes == 0) return;
    if (!src) {
      throw std::runtime_error("RawWriter: null write buffer");
    }
    checked_add_size_t(offset_bytes, bytes, "RawWriter write extent");
    const size_t expected = expected_file_bytes();
    if (offset_bytes > expected || bytes > expected - offset_bytes) {
      throw std::runtime_error("RawWriter: write extent exceeds output file dimensions");
    }
#if GROMTOOLS_HAS_POSIX_PWRITE
    RawFileOffset offset = checked_file_offset(offset_bytes, "RawWriter write offset");
    const char* ptr = static_cast<const char*>(src);
    size_t remaining = bytes;
    while (remaining > 0) {
      size_t chunk = remaining;
      if (chunk > static_cast<size_t>(SSIZE_MAX)) {
        chunk = static_cast<size_t>(SSIZE_MAX);
      }
      const ssize_t written = ::pwrite(fd, ptr, chunk, offset);
      if (written < 0) {
        const int saved_errno = errno;
        if (saved_errno == EINTR) {
          continue;
        }
        throw std::runtime_error("RawWriter: pwrite failed for " + context +
                                 ": " + std::string(std::strerror(saved_errno)));
      }
      if (written == 0) {
        throw std::runtime_error("RawWriter: pwrite returned zero bytes for " + context);
      }
      const size_t written_size = static_cast<size_t>(written);
      ptr += written_size;
      remaining -= written_size;
      offset = checked_file_offset(
        checked_add_size_t(static_cast<size_t>(offset), written_size, "RawWriter pwrite offset"),
        "RawWriter pwrite offset"
      );
    }
#else
    if (offset_bytes > static_cast<size_t>(std::numeric_limits<std::streamoff>::max())) {
      throw std::runtime_error("RawWriter: offset is too large for stream");
    }
    if (bytes > static_cast<size_t>(std::numeric_limits<std::streamsize>::max())) {
      throw std::runtime_error("RawWriter: write size is too large for stream");
    }
    f.seekp(static_cast<std::streamoff>(offset_bytes));
    f.write(reinterpret_cast<const char*>(src), static_cast<std::streamsize>(bytes));
    if (!f) {
      throw std::runtime_error("RawWriter: write failed for " + context);
    }
#endif
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

    if (n_chunk > n_rows - r0) {
      n_chunk = n_rows - r0;
    }
    if (n_chunk == 0) return;

    write_at(
      src,
      column_byte_count(n_chunk),
      column_offset_bytes(col_idx, r0),
      "column " + std::to_string(col_idx)
    );
  }
};

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

struct EngineTiming {
  double reader_init = 0.0;
  double pgen_load = 0.0;
  double wave_prepare = 0.0;
  double parallel_compute = 0.0;
  double output_flush = 0.0;
  double output_write_prepare = 0.0;
  double output_write_positioned = 0.0;
  size_t pgen_chunks = 0;
  size_t decoded_variants = 0;
  size_t decoded_dosages = 0;
  size_t compute_waves = 0;
  size_t operations = 0;
  size_t max_wave_ops = 0;
  size_t max_live_gene_columns = 0;
  size_t decoded_buffer_bytes = 0;
  size_t reader_workspace_bytes = 0;
  size_t write_batches = 0;
  size_t gene_columns_written = 0;
  size_t output_bytes_written = 0;
  size_t max_active_write_threads = 0;
  size_t duplicate_write_gids = 0;
};

struct WriteJob {
  size_t gid;
  const double* data;
  size_t offset_bytes;
  size_t bytes;
};

// Flush a batch of genes to RAW BINARY, then clear them from SoA
static inline void flush_gid_batch(SoA& S,
                                   RawWriter& writer,
                                   std::vector<size_t>& gid_toExport,
                                   int resolved_write_threads,
                                   EngineTiming& timing) {
  if (gid_toExport.empty()) return;

  const double flush_start = wall_seconds_now();
  const double prepare_start = wall_seconds_now();
  std::vector<WriteJob> jobs;
  jobs.reserve(gid_toExport.size());
  std::unordered_set<size_t> seen;
  seen.reserve(gid_toExport.size() * 2u);

  for (size_t gid : gid_toExport) {
    if (!seen.insert(gid).second) {
      ++timing.duplicate_write_gids;
      continue;
    }
    auto it = S.col.find(gid);
    if (it != S.col.end()) {
      const double* col_data = it->second.get();
      if (!col_data) {
        throw std::runtime_error("flush_gid_batch: null SoA column");
      }
      jobs.push_back(WriteJob{
        gid,
        col_data,
        writer.column_offset_bytes(gid, 0),
        writer.column_byte_count(S.n_rows)
      });
    }
  }

  for (const WriteJob& job : jobs) {
    checked_add_size_t(job.offset_bytes, job.bytes, "flush_gid_batch write extent");
  }
  timing.output_write_prepare += wall_seconds_now() - prepare_start;

  if (!jobs.empty()) {
    if (jobs.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("flush_gid_batch: too many write jobs");
    }
    size_t active_writers_size = std::min(
      static_cast<size_t>(std::max(1, resolved_write_threads)),
      jobs.size()
    );
    if (!writer.supports_parallel_positioned_writes() || jobs.size() <= 1) {
      active_writers_size = 1;
    }
    if (active_writers_size > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("flush_gid_batch: too many writer threads");
    }
    const int active_writers = static_cast<int>(active_writers_size);

    std::exception_ptr first_error;
    const double write_start = wall_seconds_now();
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(active_writers)
#endif
    for (int job_idx = 0; job_idx < static_cast<int>(jobs.size()); ++job_idx) {
      const WriteJob& job = jobs[static_cast<size_t>(job_idx)];
      try {
        writer.write_at(
          job.data,
          job.bytes,
          job.offset_bytes,
          "gene column " + std::to_string(job.gid)
        );
      } catch (const std::exception& e) {
        std::ostringstream msg;
        msg << "RAW writer failed for gid " << job.gid
            << " at offset " << job.offset_bytes
            << " for " << job.bytes << " bytes: " << e.what();
#ifdef _OPENMP
        #pragma omp critical(grom_raw_writer_error)
#endif
        {
          if (!first_error) {
            first_error = std::make_exception_ptr(std::runtime_error(msg.str()));
          }
        }
      } catch (...) {
        std::ostringstream msg;
        msg << "RAW writer failed for gid " << job.gid
            << " at offset " << job.offset_bytes
            << " for " << job.bytes << " bytes: unknown exception";
#ifdef _OPENMP
        #pragma omp critical(grom_raw_writer_error)
#endif
        {
          if (!first_error) {
            first_error = std::make_exception_ptr(std::runtime_error(msg.str()));
          }
        }
      }
    }
    timing.output_write_positioned += wall_seconds_now() - write_start;

    if (first_error) {
      std::rethrow_exception(first_error);
    }

    size_t batch_bytes = 0;
    for (const WriteJob& job : jobs) {
      batch_bytes = checked_add_size_t(batch_bytes, job.bytes, "flush_gid_batch bytes written");
      S.col.erase(job.gid);
    }
    timing.write_batches = checked_add_size_t(timing.write_batches, 1, "flush_gid_batch write batches");
    timing.gene_columns_written = checked_add_size_t(
      timing.gene_columns_written,
      jobs.size(),
      "flush_gid_batch gene columns written"
    );
    timing.output_bytes_written = checked_add_size_t(
      timing.output_bytes_written,
      batch_bytes,
      "flush_gid_batch output bytes written"
    );
    timing.max_active_write_threads = std::max(
      timing.max_active_write_threads,
      active_writers_size
    );
  }

  gid_toExport.clear();
  timing.output_flush += wall_seconds_now() - flush_start;
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

struct PgenLoadRange {
  size_t destination_offset;
  size_t destination_size;
  size_t first_variant_1based;
  size_t last_variant_1based;
  size_t variant_count;
};

static void parallel_load_pgen_chunk(
    std::vector<std::unique_ptr<RPgenReader>>& readers,
    double* destination,
    size_t destination_size,
    size_t n_inds,
    size_t global_variant_start_1based,
    size_t n_chunk_snps,
    int resolved_pgen_threads,
    bool meanimpute,
    bool rounded_mean)
{
  if (n_chunk_snps == 0) return;
  if (n_inds == 0) {
    throw std::runtime_error("parallel_load_pgen_chunk: n_inds must be positive");
  }
  if (!destination) {
    throw std::runtime_error("parallel_load_pgen_chunk: destination is null");
  }
  if (readers.empty()) {
    throw std::runtime_error("parallel_load_pgen_chunk: reader pool is empty");
  }
  if (resolved_pgen_threads <= 0) {
    throw std::runtime_error("parallel_load_pgen_chunk: resolved_pgen_threads must be positive");
  }
  const size_t active_loaders_size = std::max<size_t>(
    1,
    std::min(static_cast<size_t>(resolved_pgen_threads), n_chunk_snps)
  );
  if (active_loaders_size > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("parallel_load_pgen_chunk: too many active loaders");
  }
  const int active_loaders = static_cast<int>(active_loaders_size);
  if (static_cast<size_t>(active_loaders) > readers.size()) {
    throw std::runtime_error("parallel_load_pgen_chunk: insufficient readers for active loaders");
  }
  const size_t active_destination_size = checked_multiply_size_t(
    n_chunk_snps,
    n_inds,
    "parallel_load_pgen_chunk destination"
  );
  if (destination_size != active_destination_size) {
    throw std::runtime_error("parallel_load_pgen_chunk: destination size is " +
                             std::to_string(destination_size) +
                             "; expected exactly " +
                             std::to_string(active_destination_size));
  }

  std::vector<PgenLoadRange> ranges(static_cast<size_t>(active_loaders));
  for (int loader = 0; loader < active_loaders; ++loader) {
    const size_t range_start = checked_multiply_size_t(
      n_chunk_snps,
      static_cast<size_t>(loader),
      "parallel_load_pgen_chunk range start"
    ) / static_cast<size_t>(active_loaders);
    const size_t range_end = checked_multiply_size_t(
      n_chunk_snps,
      static_cast<size_t>(loader + 1),
      "parallel_load_pgen_chunk range end"
    ) / static_cast<size_t>(active_loaders);
    if (range_end < range_start) {
      throw std::runtime_error("parallel_load_pgen_chunk: invalid range partition");
    }
    const size_t range_count = range_end - range_start;
    const size_t destination_offset = checked_multiply_size_t(
      range_start,
      n_inds,
      "parallel_load_pgen_chunk destination offset"
    );
    const size_t range_destination_size = checked_multiply_size_t(
      range_count,
      n_inds,
      "parallel_load_pgen_chunk range destination"
    );
    if (checked_add_size_t(destination_offset,
                           range_destination_size,
                           "parallel_load_pgen_chunk range extent") > destination_size) {
      throw std::runtime_error("parallel_load_pgen_chunk: range exceeds destination");
    }
    const size_t first_variant_1based = checked_add_size_t(
      global_variant_start_1based,
      range_start,
      "parallel_load_pgen_chunk first variant"
    );
    const size_t last_variant_1based = checked_add_size_t(
      first_variant_1based,
      range_count - 1,
      "parallel_load_pgen_chunk last variant"
    );
    ranges[static_cast<size_t>(loader)] = PgenLoadRange{
      destination_offset,
      range_destination_size,
      first_variant_1based,
      last_variant_1based,
      range_count
    };
  }

  std::exception_ptr first_error;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(active_loaders)
#endif
  for (int loader = 0; loader < active_loaders; ++loader) {
    const PgenLoadRange& range = ranges[static_cast<size_t>(loader)];
    if (range.variant_count == 0) continue;

    try {
      readers[static_cast<size_t>(loader)]->ReadRangeInto(
        destination + range.destination_offset,
        range.destination_size,
        range.first_variant_1based,
        range.variant_count,
        meanimpute,
        rounded_mean
      );
    } catch (const std::exception& e) {
      std::ostringstream msg;
      msg << "PGEN loader " << loader
          << " failed while reading variants "
          << range.first_variant_1based
          << ".."
          << range.last_variant_1based
          << ": " << e.what();
#ifdef _OPENMP
      #pragma omp critical(grom_pgen_loader_error)
#endif
      {
        if (!first_error) {
          first_error = std::make_exception_ptr(std::runtime_error(msg.str()));
        }
      }
    } catch (...) {
      std::ostringstream msg;
      msg << "PGEN loader " << loader
          << " failed while reading variants "
          << range.first_variant_1based
          << ".."
          << range.last_variant_1based
          << ": unknown exception";
#ifdef _OPENMP
      #pragma omp critical(grom_pgen_loader_error)
#endif
      {
        if (!first_error) {
          first_error = std::make_exception_ptr(std::runtime_error(msg.str()));
        }
      }
    }
  }

  if (first_error) {
    std::rethrow_exception(first_error);
  }
}

static inline void flush_if_needed(SoA& S,
                                   RawWriter& writer,
                                   std::vector<size_t>& gid_toExport,
                                   size_t exportChunk,
                                   int resolved_write_threads,
                                   EngineTiming& timing) {
  if (gid_toExport.size() >= exportChunk) {
    flush_gid_batch(S, writer, gid_toExport, resolved_write_threads, timing);
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
    int resolved_write_threads,
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

    for (const WaveOp& op : ops) {
      if (op.flag == 0 || op.flag == 2) {
        gid_toExport.push_back(op.gene);
        flush_if_needed(S, writer, gid_toExport, exportChunk, resolved_write_threads, timing);
      }
    }

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
    int threads = 0,                    // <= 0 uses omp_get_max_threads()
    int pgen_threads = 0,               // <= 0 inherits resolved compute threads
    int write_threads = 1               // 1 serial positioned writes; <= 0 inherits compute threads
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

  if (exportChunk == 0) {
    Rcpp::stop("exportChunk must be positive.");
  }
  if (snp_chunk == 0) {
    Rcpp::stop("snp_chunk must be positive.");
  }
  const size_t snp_chunk_safe = snp_chunk;
  const size_t max_chunk_variants = std::min(
    snp_chunk_safe,
    static_cast<size_t>(K)
  );
  if (max_chunk_variants == 0) {
    Rcpp::stop("maximum PGEN chunk variant count must be positive.");
  }

  size_t n_inds = 0;
  bool using_dense = false;
  Rcpp::NumericMatrix G_dense;
  std::vector<std::unique_ptr<RPgenReader>> pgen_readers;
  std::vector<int> subset_std;
  const int resolved_compute_threads = available_omp_workers(threads);
  const int requested_write_threads = write_threads;
  const int resolved_write_threads = available_omp_workers(
    (write_threads <= 0) ? resolved_compute_threads : write_threads
  );
  int requested_pgen_threads = pgen_threads;
  int resolved_pgen_threads = 0;
  EngineTiming timing;

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

    subset_std.reserve(sample_subset.size());
    for (int i = 0; i < sample_subset.size(); ++i) {
      subset_std.push_back(sample_subset[i]);
    }

    if (requested_pgen_threads <= 0) {
      requested_pgen_threads = resolved_compute_threads;
    }
    resolved_pgen_threads = available_omp_workers(requested_pgen_threads);
    resolved_pgen_threads = std::max(
      1,
      std::min(resolved_pgen_threads, static_cast<int>(max_chunk_variants))
    );
    const double reader_init_start = wall_seconds_now();
    pgen_readers.reserve(static_cast<size_t>(resolved_pgen_threads));
    for (int reader_idx = 0; reader_idx < resolved_pgen_threads; ++reader_idx) {
      std::unique_ptr<RPgenReader> reader(new RPgenReader());
      reader->Load(pgen_file, nullptr, raw_sample_ct, subset_std);
      if (reader_idx == 0) {
        n_inds = static_cast<size_t>(reader->GetSubsetSize());
        if (n_inds == 0) {
          throw std::runtime_error("grom_axpy_engine: PGEN reader returned zero samples");
        }
        if (static_cast<size_t>(K) > static_cast<size_t>(reader->GetVariantCt())) {
          throw std::runtime_error("grom_axpy_engine: requested variants exceed PGEN variant count");
        }
      } else if (static_cast<size_t>(reader->GetSubsetSize()) != n_inds) {
        throw std::runtime_error("grom_axpy_engine: inconsistent subset size across PGEN readers");
      }
      timing.reader_workspace_bytes = checked_add_size_t(
        timing.reader_workspace_bytes,
        reader->GetEstimatedWorkspaceBytes(meanimpute),
        "grom_axpy_engine reader workspace"
      );
      pgen_readers.push_back(std::move(reader));
    }
    timing.reader_init = wall_seconds_now() - reader_init_start;
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

  std::vector<double> pgen_buffer;
  if (!using_dense) {
    const size_t decoded_buffer_elems = checked_multiply_size_t(
      n_inds,
      max_chunk_variants,
      "grom_axpy_engine decoded buffer"
    );
    pgen_buffer.resize(decoded_buffer_elems);
    timing.decoded_buffer_bytes = checked_multiply_size_t(
      decoded_buffer_elems,
      sizeof(double),
      "grom_axpy_engine decoded buffer bytes"
    );
  }

  for (size_t chunk_start = 0; chunk_start < static_cast<size_t>(K); chunk_start += snp_chunk_safe) {
    const size_t chunk_end = std::min(static_cast<size_t>(K), chunk_start + snp_chunk_safe);
    const size_t n_chunk_snps = chunk_end - chunk_start;
    SnpChunkView chunk;
    chunk.global_snp_start = chunk_start;
    chunk.n_rows = n_inds;
    chunk.n_cols = n_chunk_snps;
    chunk.ld = n_inds;

    if (using_dense) {
      chunk.base = &G_dense(0, static_cast<int>(chunk_start));
    } else {
      const double load_start = wall_seconds_now();
      const size_t active_destination_size = checked_multiply_size_t(
        n_chunk_snps,
        n_inds,
        "grom_axpy_engine active decoded chunk"
      );
      parallel_load_pgen_chunk(
        pgen_readers,
        pgen_buffer.data(),
        active_destination_size,
        n_inds,
        chunk_start + 1,
        n_chunk_snps,
        resolved_pgen_threads,
        meanimpute,
        rounded_mean
      );
      timing.pgen_load += wall_seconds_now() - load_start;
      timing.decoded_variants = checked_add_size_t(
        timing.decoded_variants,
        n_chunk_snps,
        "grom_axpy_engine decoded variants"
      );
      timing.decoded_dosages = checked_add_size_t(
        timing.decoded_dosages,
        active_destination_size,
        "grom_axpy_engine decoded dosages"
      );
      chunk.base = pgen_buffer.data();
    }

    process_decoded_snp_chunk(
      chunk, ip, idx, w, gp, S, writer, gid_toExport, exportChunk,
      threads, resolved_write_threads, timing
    );
    ++timing.pgen_chunks;
  }

  // Flush any leftover finished genes
  flush_gid_batch(S, writer, gid_toExport, resolved_write_threads, timing);

  // Safety net: export any remaining columns still in SoA
  if (!S.col.empty()) {
    std::vector<size_t> remaining_gids;
    remaining_gids.reserve(S.col.size());
    for (const auto& kv : S.col) {
      remaining_gids.push_back(kv.first);
    }
    flush_gid_batch(S, writer, remaining_gids, resolved_write_threads, timing);
  }

  Rcpp::Rcout << "grom_axpy_engine timing: "
              << "chunks=" << timing.pgen_chunks
              << " waves=" << timing.compute_waves
              << " ops=" << timing.operations
              << " compute_threads=" << resolved_compute_threads
              << " requested_pgen_threads=" << (using_dense ? 0 : pgen_threads)
              << " resolved_pgen_threads=" << resolved_pgen_threads
              << " requested_write_threads=" << requested_write_threads
              << " resolved_write_threads=" << resolved_write_threads
              << " write_batches=" << timing.write_batches
              << " gene_columns_written=" << timing.gene_columns_written
              << " output_bytes_written=" << timing.output_bytes_written
              << " max_active_write_threads=" << timing.max_active_write_threads
              << " duplicate_write_gids=" << timing.duplicate_write_gids
              << " reader_init=" << timing.reader_init
              << "s decoded_variants=" << timing.decoded_variants
              << " decoded_dosages=" << timing.decoded_dosages
              << " pgen_load=" << timing.pgen_load
              << "s pgen_dosages_per_sec="
              << ((timing.pgen_load > 0.0) ? (static_cast<double>(timing.decoded_dosages) / timing.pgen_load) : 0.0)
              << " wave_prepare=" << timing.wave_prepare
              << "s parallel_compute=" << timing.parallel_compute
              << "s output_flush=" << timing.output_flush
              << "s output_write_prepare=" << timing.output_write_prepare
              << "s output_write_positioned=" << timing.output_write_positioned
              << "s output_bytes_per_sec="
              << ((timing.output_write_positioned > 0.0) ? (static_cast<double>(timing.output_bytes_written) / timing.output_write_positioned) : 0.0)
              << " max_wave_ops=" << timing.max_wave_ops
              << " max_live_gene_columns=" << timing.max_live_gene_columns
              << " decoded_buffer_bytes=" << timing.decoded_buffer_bytes
              << " reader_workspace_bytes=" << timing.reader_workspace_bytes
              << " omp_workers=" << resolved_compute_threads
              << "\n";

  // We don't return the big matrix to R (it lives on disk)
  return Rcpp::NumericMatrix(0, 0);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix grom_decode_pgen_range(
    const std::string& pgen_file,
    int raw_sample_ct,
    size_t first_variant_1based,
    size_t variant_count,
    Rcpp::IntegerVector sample_subset,
    bool meanimpute = true,
    bool rounded_mean = false,
    int pgen_threads = 1)
{
  if (pgen_file.empty()) {
    Rcpp::stop("pgen_file must not be empty.");
  }
  if (variant_count == 0) {
    return Rcpp::NumericMatrix(0, 0);
  }
  if (first_variant_1based == 0) {
    Rcpp::stop("first_variant_1based must be positive.");
  }

  std::vector<int> subset_std;
  subset_std.reserve(sample_subset.size());
  for (int i = 0; i < sample_subset.size(); ++i) {
    subset_std.push_back(sample_subset[i]);
  }

  int resolved_pgen_threads = available_omp_workers(pgen_threads);
  resolved_pgen_threads = static_cast<int>(std::max<size_t>(
    1,
    std::min(static_cast<size_t>(resolved_pgen_threads), variant_count)
  ));
  std::vector<std::unique_ptr<RPgenReader>> readers;
  readers.reserve(static_cast<size_t>(resolved_pgen_threads));
  size_t n_inds = 0;
  for (int reader_idx = 0; reader_idx < resolved_pgen_threads; ++reader_idx) {
    std::unique_ptr<RPgenReader> reader(new RPgenReader());
    reader->Load(pgen_file, nullptr, raw_sample_ct, subset_std);
    if (reader_idx == 0) {
      n_inds = static_cast<size_t>(reader->GetSubsetSize());
      if (n_inds == 0) {
        throw std::runtime_error("grom_decode_pgen_range: PGEN reader returned zero samples");
      }
      const size_t raw_variant_ct = static_cast<size_t>(reader->GetVariantCt());
      const size_t first_variant_0based = first_variant_1based - 1;
      if (first_variant_0based >= raw_variant_ct ||
          variant_count > raw_variant_ct - first_variant_0based) {
        throw std::runtime_error("grom_decode_pgen_range: variant range out of bounds");
      }
    } else if (static_cast<size_t>(reader->GetSubsetSize()) != n_inds) {
      throw std::runtime_error("grom_decode_pgen_range: inconsistent subset size across PGEN readers");
    }
    readers.push_back(std::move(reader));
  }

  if (n_inds > static_cast<size_t>(std::numeric_limits<int>::max()) ||
      variant_count > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("grom_decode_pgen_range: dimensions exceed R matrix limits");
  }
  const size_t destination_size = checked_multiply_size_t(
    n_inds,
    variant_count,
    "grom_decode_pgen_range destination"
  );
  Rcpp::NumericMatrix out(static_cast<int>(n_inds), static_cast<int>(variant_count));
  if (static_cast<size_t>(out.size()) != destination_size) {
    throw std::runtime_error("grom_decode_pgen_range: allocated matrix size mismatch");
  }
  parallel_load_pgen_chunk(
    readers,
    out.begin(),
    destination_size,
    n_inds,
    first_variant_1based,
    variant_count,
    resolved_pgen_threads,
    meanimpute,
    rounded_mean
  );
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix grom_decode_pgen_range_readlist(
    const std::string& pgen_file,
    int raw_sample_ct,
    size_t first_variant_1based,
    size_t variant_count,
    Rcpp::IntegerVector sample_subset,
    bool meanimpute = true,
    bool rounded_mean = false)
{
  if (pgen_file.empty()) {
    Rcpp::stop("pgen_file must not be empty.");
  }
  if (variant_count == 0) {
    return Rcpp::NumericMatrix(0, 0);
  }
  if (first_variant_1based == 0) {
    Rcpp::stop("first_variant_1based must be positive.");
  }

  std::vector<int> subset_std;
  subset_std.reserve(sample_subset.size());
  for (int i = 0; i < sample_subset.size(); ++i) {
    subset_std.push_back(sample_subset[i]);
  }

  RPgenReader reader;
  reader.Load(pgen_file, nullptr, raw_sample_ct, subset_std);
  const size_t n_inds = static_cast<size_t>(reader.GetSubsetSize());
  if (n_inds == 0) {
    throw std::runtime_error("grom_decode_pgen_range_readlist: PGEN reader returned zero samples");
  }
  const size_t raw_variant_ct = static_cast<size_t>(reader.GetVariantCt());
  const size_t first_variant_0based = first_variant_1based - 1;
  if (first_variant_0based >= raw_variant_ct ||
      variant_count > raw_variant_ct - first_variant_0based) {
    throw std::runtime_error("grom_decode_pgen_range_readlist: variant range out of bounds");
  }
  if (n_inds > static_cast<size_t>(std::numeric_limits<int>::max()) ||
      variant_count > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("grom_decode_pgen_range_readlist: dimensions exceed R matrix limits");
  }

  std::vector<int> variants;
  variants.reserve(variant_count);
  for (size_t j = 0; j < variant_count; ++j) {
    const size_t variant_1based = checked_add_size_t(
      first_variant_1based,
      j,
      "grom_decode_pgen_range_readlist variant index"
    );
    if (variant_1based > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("grom_decode_pgen_range_readlist: variant index exceeds int range");
    }
    variants.push_back(static_cast<int>(variant_1based));
  }

  const size_t destination_size = checked_multiply_size_t(
    n_inds,
    variant_count,
    "grom_decode_pgen_range_readlist destination"
  );
  std::vector<double> decoded(destination_size);
  reader.ReadList(decoded, variants, meanimpute, rounded_mean);

  Rcpp::NumericMatrix out(static_cast<int>(n_inds), static_cast<int>(variant_count));
  if (static_cast<size_t>(out.size()) != destination_size) {
    throw std::runtime_error("grom_decode_pgen_range_readlist: allocated matrix size mismatch");
  }
  std::copy(decoded.begin(), decoded.end(), out.begin());
  return out;
}

// [[Rcpp::export]]
Rcpp::List grom_test_readlist_buffer_contract(
    const std::string& pgen_file,
    int raw_sample_ct,
    size_t first_variant_1based,
    size_t variant_count,
    Rcpp::IntegerVector sample_subset,
    bool meanimpute = true,
    bool rounded_mean = false)
{
  if (variant_count == 0) {
    Rcpp::stop("variant_count must be positive for this buffer contract test.");
  }

  std::vector<int> subset_std;
  subset_std.reserve(sample_subset.size());
  for (int i = 0; i < sample_subset.size(); ++i) {
    subset_std.push_back(sample_subset[i]);
  }

  RPgenReader reader;
  reader.Load(pgen_file, nullptr, raw_sample_ct, subset_std);
  const size_t n_inds = static_cast<size_t>(reader.GetSubsetSize());
  const size_t required_size = checked_multiply_size_t(
    n_inds,
    variant_count,
    "grom_test_readlist_buffer_contract destination"
  );
  if (required_size == 0) {
    Rcpp::stop("required_size must be positive for this buffer contract test.");
  }

  std::vector<int> variants;
  variants.reserve(variant_count);
  for (size_t j = 0; j < variant_count; ++j) {
    const size_t variant_1based = checked_add_size_t(
      first_variant_1based,
      j,
      "grom_test_readlist_buffer_contract variant index"
    );
    if (variant_1based > static_cast<size_t>(std::numeric_limits<int>::max())) {
      Rcpp::stop("variant index exceeds int range.");
    }
    variants.push_back(static_cast<int>(variant_1based));
  }

  std::vector<double> reference(required_size);
  reader.ReadList(reference, variants, meanimpute, rounded_mean);

  std::vector<double> exact(required_size);
  reader.ReadList(exact, variants, meanimpute, rounded_mean);
  const bool exact_matches = std::equal(reference.begin(), reference.end(), exact.begin());

  const double sentinel = -987654321.0;
  std::vector<double> larger(required_size + 3, sentinel);
  reader.ReadList(larger, variants, meanimpute, rounded_mean);
  const bool larger_prefix_matches = std::equal(reference.begin(), reference.end(), larger.begin());
  const bool larger_tail_unchanged =
    larger[required_size] == sentinel &&
    larger[required_size + 1] == sentinel &&
    larger[required_size + 2] == sentinel;

  std::vector<int> empty_variants;
  std::vector<double> empty_buffer;
  bool empty_ok = true;
  try {
    reader.ReadList(empty_buffer, empty_variants, meanimpute, rounded_mean);
  } catch (...) {
    empty_ok = false;
  }

  std::vector<double> capacity_only;
  capacity_only.reserve(required_size);
  capacity_only.resize(required_size - 1);
  bool undersized_rejected = false;
  std::string undersized_error;
  try {
    reader.ReadList(capacity_only, variants, meanimpute, rounded_mean);
  } catch (const std::exception& e) {
    undersized_rejected = true;
    undersized_error = e.what();
  }

  return Rcpp::List::create(
    Rcpp::Named("required_size") = static_cast<double>(required_size),
    Rcpp::Named("capacity_only_size") = static_cast<double>(capacity_only.size()),
    Rcpp::Named("capacity_only_capacity") = static_cast<double>(capacity_only.capacity()),
    Rcpp::Named("undersized_rejected") = undersized_rejected,
    Rcpp::Named("undersized_error") = undersized_error,
    Rcpp::Named("exact_matches") = exact_matches,
    Rcpp::Named("larger_prefix_matches") = larger_prefix_matches,
    Rcpp::Named("larger_tail_unchanged") = larger_tail_unchanged,
    Rcpp::Named("empty_ok") = empty_ok
  );
}

// [[Rcpp::export]]
Rcpp::List grom_test_pgen_workspace_estimate(
    const std::string& pgen_file,
    int raw_sample_ct,
    Rcpp::IntegerVector sample_subset,
    size_t variant_1based,
    bool rounded_mean = false)
{
  std::vector<int> subset_std;
  subset_std.reserve(sample_subset.size());
  for (int i = 0; i < sample_subset.size(); ++i) {
    subset_std.push_back(sample_subset[i]);
  }

  RPgenReader reader;
  reader.Load(pgen_file, nullptr, raw_sample_ct, subset_std);
  const size_t n_inds = static_cast<size_t>(reader.GetSubsetSize());
  const size_t before_false = reader.GetEstimatedWorkspaceBytes(false);
  const size_t before_true = reader.GetEstimatedWorkspaceBytes(true);
  const size_t expected_scratch_bytes = checked_add_size_t(
    checked_multiply_size_t(n_inds, sizeof(double), "workspace estimate double scratch"),
    checked_multiply_size_t(n_inds, sizeof(unsigned char), "workspace estimate missing scratch"),
    "workspace estimate scratch"
  );

  std::vector<double> decoded(n_inds);
  reader.ReadRangeInto(
    decoded.data(),
    decoded.size(),
    variant_1based,
    1,
    true,
    rounded_mean
  );
  const size_t after_false = reader.GetEstimatedWorkspaceBytes(false);
  const size_t after_true = reader.GetEstimatedWorkspaceBytes(true);

  return Rcpp::List::create(
    Rcpp::Named("n_inds") = static_cast<double>(n_inds),
    Rcpp::Named("expected_scratch_bytes") = static_cast<double>(expected_scratch_bytes),
    Rcpp::Named("before_false") = static_cast<double>(before_false),
    Rcpp::Named("before_true") = static_cast<double>(before_true),
    Rcpp::Named("after_false") = static_cast<double>(after_false),
    Rcpp::Named("after_true") = static_cast<double>(after_true),
    Rcpp::Named("before_delta") = static_cast<double>(before_true - before_false),
    Rcpp::Named("after_delta") = static_cast<double>(after_true - after_false),
    Rcpp::Named("true_estimate_non_decreasing") = after_true >= before_true
  );
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
