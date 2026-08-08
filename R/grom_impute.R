#' Impute genetically regulated omics values
#'
#' Build sparse weight triplets and generate `.grom`, `.gid`, and `.sid`
#' outputs from model weights and PLINK2 genotype inputs.
#'
#' @param weights_table A `data.frame` or `data.table` of model weights.
#' @param grom_pfx Output prefix for the generated grom files.
#' @param pgen_dir Path to a directory containing PLINK2 `.pgen`, `.pvar`, and
#'   `.psam` files.
#' @param snp_chunk Integer chunk size used while streaming SNPs. Specifies how
#'   many genotype records are loaded at a time and can be used to control peak
#'   and average memory usage.
#' @param sample_subset Integer vector specifying which samples to load and use
#'   for imputation.
#' @param CHUNK Integer chunk size passed to the low-level imputation engine.
#'   Specifies the tile size used during computation and can help maintain
#'   efficient CPU-cache utilization.
#' @param exportChunk Integer chunk size used when exporting results. Specifies
#'   how many output columns are written at a time and can be used to control
#'   peak and average memory usage.
#' @param meanimpute Logical indicating whether missing dosages should be mean
#'   imputed.
#' @param is_round Logical indicating whether mean imputation should use rounded
#'   means.
#' @param threads Optional positive integer limiting the number of OpenMP worker
#'   threads used by the low-level imputation engine. `NULL` or non-positive
#'   values preserve the engine default.
#' @param pgen_threads Optional positive integer limiting the number of
#'   independent PGEN loader threads. `NULL` or non-positive values inherit the
#'   resolved compute-thread count.
#'
#' @return Invisibly returns the result of the imputation pipeline after writing
#'   output files to disk.
#'
#' @examples
#' library(gromtools)
#'
#' pgen_dir <- system.file(
#'   "extdata",
#'   "synthetic_chromosomes",
#'   package = "gromtools"
#' )
#' db_directory <- system.file(
#'   "extdata",
#'   "synth_small_variant_weights_db",
#'   package = "gromtools"
#' )
#' model_weights_table <- read_db_dir(db_dir = db_directory)
#'
#' out_dir <- file.path(tempdir(), "tmp_grom_run")
#' dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#' grom_pfx <- file.path(out_dir, "synth_example")
#'
#' grom_impute(
#'   weights_table = model_weights_table,
#'   grom_pfx = grom_pfx,
#'   pgen_dir = pgen_dir
#' )
#'
#' @details
#' Chromosome-to-file mapping uses the `chromosomes` column in
#' `weights_table`. The basenames of the `.pgen` and `.pvar` files in
#' `pgen_dir` must contain the corresponding chromosome tag. For example:
#'
#' \preformatted{
#' weights_table$chromosomes:
#'   chr1
#'   chr2
#'   chr10
#'
#' Matching filenames in pgen_dir:
#'   cohort.chr1.pgen
#'   cohort.chr1.pvar
#'   cohort.chr2.pgen
#'   cohort.chr2.pvar
#'   cohort.chr10.pgen
#'   cohort.chr10.pvar
#' }
#'
#' Matching is unique by chromosome tag, so `chr1` does not match `chr10`.
#' Each chromosome tag must resolve to exactly one `.pgen` and one `.pvar`
#' file. The `.psam` file is discovered separately and is not chromosome-tag
#' matched.
#'
#' To confirm the resolved mapping after a run, inspect the generated files
#' `meta/chr_to_pgen_map.tsv` and `meta/chr_to_pvar_map.tsv` under the output
#' directory. These record which chromosome tag was matched to which genotype
#' file basename. On completion, `grom_impute()` also reports the output
#' directory path in the console message stream.
grom_impute <- function(
  weights_table,
  grom_pfx,
  pgen_dir,
  snp_chunk     = 1000L,
  sample_subset = integer(0),
  CHUNK         = 1000L,
  exportChunk   = 256L,
  meanimpute    = TRUE,
  is_round      = TRUE,
  threads       = NULL,
  pgen_threads  = NULL
) {
  if (!data.table::is.data.table(weights_table)) {
    if (is.data.frame(weights_table)) {
      weights_table <- data.table::as.data.table(weights_table)
    } else {
      stop("weights_table must be a data.frame or data.table.")
    }
  }

  out_dir <- dirname(grom_pfx)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  build_csc_triplets(
    grom_pfx = grom_pfx,
    pgen_dir = pgen_dir,
    variant_weights = weights_table
  )

  impute_grom(
    pgen_dir = pgen_dir,
    grom_pfx = grom_pfx,
    variant_weights = weights_table,
    snp_chunk = snp_chunk,
    sample_subset = sample_subset,
    CHUNK = CHUNK,
    exportChunk = exportChunk,
    meanimpute = meanimpute,
    is_round = is_round,
    threads = threads,
    pgen_threads = pgen_threads
  )
}
