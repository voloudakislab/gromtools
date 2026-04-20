# Internal helper used by read_db_dir().
sqlite_read_model_db <- function(db_path, extra_cols = NULL) {
  data.table::as.data.table(sqlite_read_model_db_cpp(db_path, extra_cols))
}
