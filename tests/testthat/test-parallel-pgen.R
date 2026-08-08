test_that("parallel PGEN decoder matches one-loader decoding exactly", {
  pgen_dir <- system.file("extdata", "synthetic_chromosomes", package = "gromtools")
  pgen_file <- file.path(pgen_dir, "synth.chr1.pgen")
  psam_file <- file.path(pgen_dir, "synth.chr1.psam")
  raw_sample_ct <- nrow(data.table::fread(psam_file))

  scenarios <- list(
    list(first = 1L, count = 1L, subset = integer(0), mean = FALSE, round = FALSE),
    list(first = 3L, count = 7L, subset = integer(0), mean = TRUE, round = TRUE),
    list(first = 5L, count = 11L, subset = 1:20, mean = FALSE, round = FALSE),
    list(
      first = 9L,
      count = 17L,
      subset = c(1L, 2L, 5L, 7L, 13L, 21L, 34L, 55L, 89L),
      mean = TRUE,
      round = FALSE
    )
  )

  for (sc in scenarios) {
    ref <- grom_decode_pgen_range_readlist(
      pgen_file, raw_sample_ct, sc$first, sc$count,
      as.integer(sc$subset), sc$mean, sc$round
    )
    for (pgen_threads in c(1L, 2L, 3L, 4L, 8L, 32L)) {
      got <- grom_decode_pgen_range(
        pgen_file, raw_sample_ct, sc$first, sc$count,
        as.integer(sc$subset), sc$mean, sc$round, pgen_threads
      )
      expect_identical(dim(got), dim(ref))
      expect_identical(as.vector(got), as.vector(ref))
    }
  }

  expect_error(
    grom_decode_pgen_range(
      pgen_file, raw_sample_ct, 1L, 3L, c(5L, 1L, 3L), TRUE, FALSE, 2L
    ),
    "strictly increasing"
  )
})

test_that("parallel PGEN decoder matches ReadList with missing genotypes", {
  plink2 <- Sys.which("plink2")
  skip_if(!nzchar(plink2), "plink2 is required to create the missing-genotype PGEN fixture")

  fixture_dir <- file.path(tempdir(), paste0("gromtools-missing-pgen-", Sys.getpid()))
  dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
  vcf <- file.path(fixture_dir, "missing.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=1>",
    paste(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
      "s1", "s2", "s3", "s4", "s5", "s6",
      sep = "\t"
    ),
    "1\t1001\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.\t0/0\t0/1",
    "1\t1002\tm2\tC\tT\t.\tPASS\t.\tGT\t./.\t0/0\t0/1\t1/1\t./.\t0/0",
    "1\t1003\tm3\tG\tA\t.\tPASS\t.\tGT\t1/1\t./.\t0/0\t0/1\t1/1\t./.",
    "1\t1004\tm4\tT\tC\t.\tPASS\t.\tGT\t0/1\t1/1\t./.\t0/0\t0/1\t1/1",
    "1\t1005\tm5\tA\tC\t.\tPASS\t.\tGT\t0/0\t./.\t./.\t0/1\t1/1\t0/0"
  ), vcf)

  out_prefix <- file.path(fixture_dir, "missing")
  status <- system2(
    plink2,
    c("--vcf", vcf, "--make-pgen", "--out", out_prefix),
    stdout = TRUE,
    stderr = TRUE
  )
  status_code <- attr(status, "status")
  if (is.null(status_code)) status_code <- 0L
  expect_equal(status_code, 0L)
  if (!file.exists(paste0(out_prefix, ".pgen"))) {
    fail(paste(status, collapse = "\n"))
  }

  pgen_file <- paste0(out_prefix, ".pgen")
  scenarios <- list(
    list(subset = integer(0), mean = FALSE, round = FALSE),
    list(subset = integer(0), mean = TRUE, round = FALSE),
    list(subset = integer(0), mean = TRUE, round = TRUE),
    list(subset = c(1L, 3L, 6L), mean = FALSE, round = FALSE),
    list(subset = c(1L, 3L, 6L), mean = TRUE, round = TRUE)
  )

  for (sc in scenarios) {
    ref <- grom_decode_pgen_range_readlist(
      pgen_file, 6L, 1L, 5L, as.integer(sc$subset), sc$mean, sc$round
    )
    for (pgen_threads in c(1L, 2L, 3L, 4L, 8L)) {
      got <- grom_decode_pgen_range(
        pgen_file, 6L, 1L, 5L, as.integer(sc$subset), sc$mean, sc$round, pgen_threads
      )
      expect_identical(as.vector(got), as.vector(ref))
    }
  }
})

test_that("parallel PGEN engine output is byte-identical to one-loader output", {
  pgen_dir <- system.file("extdata", "synthetic_chromosomes", package = "gromtools")
  pgen_file <- file.path(pgen_dir, "synth.chr1.pgen")
  psam_file <- file.path(pgen_dir, "synth.chr1.psam")
  raw_sample_ct <- nrow(data.table::fread(psam_file))

  run_identity <- function(pgen_threads) {
    K <- 19L
    out <- tempfile(fileext = ".grom")
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
      exportChunk = 3L,
      n_total_genes = K,
      create_new = TRUE,
      G_in = NULL,
      showWarnings = FALSE,
      rounded_mean = FALSE,
      threads = 2L,
      pgen_threads = pgen_threads
    )
    readBin(out, "raw", n = file.info(out)$size)
  }

  ref <- run_identity(1L)
  for (pgen_threads in c(2L, 4L, 8L, 32L)) {
    expect_identical(run_identity(pgen_threads), ref)
  }
})
