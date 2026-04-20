#' Read selected values from a grom output
#'
#' Stream selected model, gene, and sample combinations from an existing grom
#' output on disk.
#'
#' @param grom_pfx Path prefix to an existing grom output without the file
#'   extension.
#' @param models Optional character vector of model IDs to retain. `NULL`
#'   selects all models.
#' @param genes Optional character vector of gene IDs to retain. `NULL` selects
#'   all genes.
#' @param samples Optional character vector of sample IDs to retain. `NULL`
#'   selects all samples.
#'
#' @return A numeric matrix with samples in rows and selected model-gene pairs
#'   in columns.
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
#' grom_mat <- grom_read(
#'   grom_pfx = grom_pfx,
#'   models = c("AMR_subclass_IN_SST", "AMR_subclass_VLMC"),
#'   genes = c("ENSG00000003987", "ENSG00000053900"),
#'   samples = c(
#'     "sample0001", "sample0014", "sample0034", "sample0039", "sample0043"
#'   )
#' )
#'
#' print(grom_mat)
grom_read <- function(grom_pfx, models = NULL, genes = NULL, samples = NULL) {
  assert_string_vector <- function(x, arg) {
    if (is.null(x)) {
      return(invisible(NULL))
    }
    if (data.table::is.data.table(x) || is.data.frame(x) || !is.atomic(x) || !is.null(dim(x))) {
      stop(arg, " must be a character vector or NULL.")
    }
    if (!is.character(x)) {
      stop(arg, " must be a character vector or NULL.")
    }
    invisible(NULL)
  }

  if (!is.character(grom_pfx) || length(grom_pfx) != 1L || is.na(grom_pfx) || !nzchar(grom_pfx)) {
    stop("grom_pfx must be a single non-empty character string.")
  }

  assert_string_vector(models, "models")
  assert_string_vector(genes, "genes")
  assert_string_vector(samples, "samples")

  required_files <- paste0(grom_pfx, c(".grom", ".gid", ".sid"))
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files)) {
    stop(
      "The following grom output files do not exist: ",
      paste(missing_files, collapse = ", ")
    )
  }

  read_grom_gid <- function(grom_prefix) {
    gid <- fread(paste0(grom_prefix, ".gid"))
    if (!"mg_id" %in% names(gid)) {
      gid[, mg_id := .I - 1L]
    }
    gid
  }

  read_grom_sid <- function(grom_prefix) {
    sid <- fread(paste0(grom_prefix, ".sid"))
    if ("#IID" %in% names(sid) && !"IID" %in% names(sid)) {
      setnames(sid, "#IID", "IID")
    }
    if (!"GROM_SID" %in% names(sid)) {
      sid[, GROM_SID := .I]
    }
    sid
  }

  gid <- read_grom_gid(grom_pfx)
  sid <- read_grom_sid(grom_pfx)

  gid_use <- copy(gid)
  if (!is.null(models)) {
    gid_use <- gid_use[model_ID %in% unique(models)]
  }
  if (!is.null(genes)) {
    gid_use <- gid_use[gene %in% unique(genes)]
  }

  keep_use <- NULL
  if (!is.null(samples)) {
    keep_idx <- match(unique(samples), sid$IID)
    keep_idx <- keep_idx[!is.na(keep_idx)]
    keep_use <- sid$IID[keep_idx]
  }

  read_grom(
    prefix = grom_pfx,
    extract = gid_use,
    keep = keep_use
  )
}
