#' Read a model SQLite database into a table-like R object
#'
#' @param db_path Path to a SQLite model database.
#' @param extra_cols Optional character vector of additional columns to look for
#'   in the `weights` and `extra` tables.
#' @return A `data.table` with required columns plus any discovered
#'   `extra_cols`.
#' @export
sqlite_read_model_db <- function(db_path, extra_cols = NULL) {
  data.table::as.data.table(sqlite_read_model_db_cpp(db_path, extra_cols))
}
