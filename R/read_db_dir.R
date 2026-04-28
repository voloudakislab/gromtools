#' Read model database files from a directory
#'
#' Load all SQLite model databases in a directory into a single weights table.
#'
#' @param db_dir Path to a directory containing `.db` model files.
#' @param extra_cols Optional character vector of additional columns to read from
#'   each database when available.
#'
#' @return A `data.table` combining the contents of all `.db` files in
#'   `db_dir`.
#'
#' @examples
#' library(gromtools)
#'
#' db_directory <- system.file(
#'   "extdata",
#'   "synth_small_variant_weights_db",
#'   package = "gromtools"
#' )
#'
#' model_weights_table <- read_db_dir(db_dir = db_directory)
#' head(model_weights_table)
read_db_dir <- function(db_dir, extra_cols = NULL) {
  if (!is.character(db_dir) || length(db_dir) != 1L || is.na(db_dir) || !nzchar(db_dir)) {
    stop("db_dir must be a single non-empty character string.")
  }

  if (dir.exists(db_dir)) {
    db_files <- list.files(
      db_dir,
      pattern = "[.]db$",
      full.names = TRUE
    )
    db_files <- sort(db_files)
  } else if (file.exists(db_dir) && grepl("[.]db$", db_dir, ignore.case = TRUE)) {
    db_files <- normalizePath(db_dir, winslash = "/", mustWork = TRUE)
  } else {
    stop("db_dir does not exist or is not a .db file: ", db_dir)
  }

  if (!length(db_files)) {
    stop("No .db files found in db_dir: ", db_dir)
  }

  tables <- lapply(
    db_files,
    function(db_path) sqlite_read_model_db(db_path, extra_cols = extra_cols)
  )

  data.table::rbindlist(tables, use.names = TRUE, fill = TRUE)
}
