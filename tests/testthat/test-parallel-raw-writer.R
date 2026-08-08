sha_file <- function(path) {
  unname(strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1]][1])
}

run_dense_writer_case <- function(write_threads, exportChunk = 2L, suffix = "") {
  G <- matrix(
    c(
      0, 1, 2, 1,
      2, 0, 1, 0,
      1, 1, 0, 2,
      0, 2, 2, 1,
      1, 0, 1, 2,
      2, 1, 0, 0
    ),
    nrow = 6,
    ncol = 4
  )
  indptr <- as.integer(c(0, 2, 4, 6, 7))
  indices <- as.integer(c(0, 7, 3, 3, 1, 7, 5))
  data <- c(1.5, -2.0, 0.0, 3.0, -1.0, 2.5, 4.0)
  gene_pos <- as.integer(c(1, 1, 0, 0, 2, 0, 5))
  out <- tempfile(paste0("dense-writer-", suffix), fileext = ".grom")
  grom_axpy_engine(
    indptr = indptr,
    indices = indices,
    data = data,
    pgen_file = "",
    raw_sample_ct = 0L,
    snp_chunk = 2L,
    sample_subset = integer(),
    meanimpute = TRUE,
    grom_file = out,
    CHUNK = 2L,
    gene_pos = gene_pos,
    exportChunk = exportChunk,
    n_total_genes = 8L,
    create_new = TRUE,
    G_in = G,
    showWarnings = FALSE,
    rounded_mean = FALSE,
    threads = 3L,
    pgen_threads = 0L,
    write_threads = write_threads
  )
  out
}

run_pgen_writer_case <- function(write_threads, exportChunk = 3L, suffix = "") {
  pgen_dir <- system.file("extdata", "synthetic_chromosomes", package = "gromtools")
  pgen_file <- file.path(pgen_dir, "synth.chr1.pgen")
  psam_file <- file.path(pgen_dir, "synth.chr1.psam")
  raw_sample_ct <- nrow(data.table::fread(psam_file))
  K <- 19L
  out <- tempfile(paste0("pgen-writer-", suffix), fileext = ".grom")
  grom_axpy_engine(
    indptr = as.integer(0:K),
    indices = as.integer(0:(K - 1L)),
    data = rep(1, K),
    pgen_file = pgen_file,
    raw_sample_ct = raw_sample_ct,
    snp_chunk = 7L,
    sample_subset = as.integer(c(3L, 4L, 8L, 9L, 10L, 20L, 40L, 80L, 100L)),
    meanimpute = TRUE,
    grom_file = out,
    CHUNK = 4L,
    gene_pos = rep(0L, K),
    exportChunk = exportChunk,
    n_total_genes = K,
    create_new = TRUE,
    G_in = NULL,
    showWarnings = FALSE,
    rounded_mean = FALSE,
    threads = 2L,
    pgen_threads = 2L,
    write_threads = write_threads
  )
  out
}

test_that("parallel RAW writer is byte-identical on dense input", {
  ref <- run_dense_writer_case(1L, exportChunk = 1L, suffix = "ref")
  ref_raw <- readBin(ref, "raw", n = file.info(ref)$size)
  ref_hash <- sha_file(ref)

  for (write_threads in c(2L, 3L, 4L, 8L, 16L, 32L, 0L)) {
    got <- run_dense_writer_case(write_threads, exportChunk = 1L, suffix = paste0("wt", write_threads))
    expect_equal(file.info(got)$size, file.info(ref)$size)
    expect_identical(sha_file(got), ref_hash)
    expect_identical(readBin(got, "raw", n = file.info(got)$size), ref_raw)
  }

  ref2 <- run_dense_writer_case(1L, exportChunk = 5L, suffix = "ref-large")
  for (write_threads in c(2L, 4L, 16L)) {
    got <- run_dense_writer_case(write_threads, exportChunk = 5L, suffix = paste0("large", write_threads))
    expect_identical(sha_file(got), sha_file(ref2))
  }
})

test_that("parallel RAW writer is byte-identical on PGEN input", {
  ref <- run_pgen_writer_case(1L, exportChunk = 3L, suffix = "ref")
  ref_raw <- readBin(ref, "raw", n = file.info(ref)$size)
  ref_hash <- sha_file(ref)

  for (write_threads in c(2L, 3L, 4L, 8L, 16L, 32L, 0L)) {
    got <- run_pgen_writer_case(write_threads, exportChunk = 3L, suffix = paste0("wt", write_threads))
    expect_equal(file.info(got)$size, file.info(ref)$size)
    expect_identical(sha_file(got), ref_hash)
    expect_identical(readBin(got, "raw", n = file.info(got)$size), ref_raw)
  }
})

test_that("parallel RAW writer handles create_new false continuation", {
  G1 <- matrix(1:12, nrow = 3, ncol = 4)
  G2 <- matrix(13:24, nrow = 3, ncol = 4)
  run_split <- function(write_threads) {
    out <- tempfile("split-writer-", fileext = ".grom")
    grom_axpy_engine(
      indptr = as.integer(c(0, 1, 1, 1, 1)),
      indices = as.integer(0L),
      data = 2,
      pgen_file = "",
      raw_sample_ct = 0L,
      snp_chunk = 2L,
      sample_subset = integer(),
      meanimpute = TRUE,
      grom_file = out,
      CHUNK = 2L,
      gene_pos = as.integer(0L),
      exportChunk = 1L,
      n_total_genes = 6L,
      create_new = TRUE,
      G_in = G1,
      showWarnings = FALSE,
      rounded_mean = FALSE,
      threads = 2L,
      pgen_threads = 0L,
      write_threads = write_threads
    )
    grom_axpy_engine(
      indptr = as.integer(c(0, 0, 0, 0, 1)),
      indices = as.integer(5L),
      data = -1,
      pgen_file = "",
      raw_sample_ct = 0L,
      snp_chunk = 2L,
      sample_subset = integer(),
      meanimpute = TRUE,
      grom_file = out,
      CHUNK = 2L,
      gene_pos = as.integer(0L),
      exportChunk = 1L,
      n_total_genes = 6L,
      create_new = FALSE,
      G_in = G2,
      showWarnings = FALSE,
      rounded_mean = FALSE,
      threads = 2L,
      pgen_threads = 0L,
      write_threads = write_threads
    )
    out
  }

  ref <- run_split(1L)
  for (write_threads in c(2L, 8L, 16L)) {
    got <- run_split(write_threads)
    expect_equal(file.info(got)$size, file.info(ref)$size)
    expect_identical(sha_file(got), sha_file(ref))
  }
})

test_that("parallel RAW writer is deterministic across repeated runs", {
  hashes <- character(20)
  for (i in seq_along(hashes)) {
    hashes[i] <- sha_file(run_pgen_writer_case(8L, exportChunk = 3L, suffix = paste0("stress", i)))
  }
  expect_length(unique(hashes), 1L)
})
