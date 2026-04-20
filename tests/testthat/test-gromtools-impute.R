test_that("grom_impute reproduces the documented example subset", {
  pgen_dir <- system.file("extdata", "synthetic_chromosomes", package = "gromtools")
  db_dir <- system.file("extdata", "synth_small_variant_weights_db", package = "gromtools")

  expect_true(nzchar(pgen_dir))
  expect_true(nzchar(db_dir))

  weights_table <- read_db_dir(db_dir)

  out_dir <- file.path(
    tempdir(),
    paste0("gromtools-impute-test-", Sys.getpid(), "-", as.integer(Sys.time()))
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  grom_pfx <- file.path(out_dir, "synth_example")

  grom_impute(
    weights_table = weights_table,
    grom_pfx = grom_pfx,
    pgen_dir = pgen_dir
  )

  expect_true(file.exists(paste0(grom_pfx, ".grom")))
  expect_true(file.exists(paste0(grom_pfx, ".gid")))
  expect_true(file.exists(paste0(grom_pfx, ".sid")))

  grom_mat <- grom_read(
    grom_pfx = grom_pfx,
    models = c("AMR_subclass_IN_SST", "AMR_subclass_VLMC"),
    genes = c("ENSG00000003987", "ENSG00000053900"),
    samples = c("sample0001", "sample0014", "sample0034", "sample0039", "sample0043")
  )

  expected <- matrix(
    c(
      0.3191238, 0.4542425, 0.4138578, 0.2676161, 0.3721308,
      0.2737647, -0.4481538, -0.2212396, -0.3217697, 0.2609762,
      -0.5058091, -0.9249124, -0.5259653, -0.9044481, -0.6188168
    ),
    nrow = 5,
    ncol = 3
  )

  rownames(expected) <- c(
    "sample0001", "sample0014", "sample0034", "sample0039", "sample0043"
  )
  colnames(expected) <- c(
    "AMR_subclass_IN_SST_ENSG00000003987",
    "AMR_subclass_IN_SST_ENSG00000053900",
    "AMR_subclass_VLMC_ENSG00000053900"
  )

  expect_equal(dim(grom_mat), c(5L, 3L))
  expect_equal(rownames(grom_mat), rownames(expected))
  expect_equal(colnames(grom_mat), colnames(expected))
  expect_equal(unname(grom_mat), unname(expected), tolerance = 1e-7)
})
