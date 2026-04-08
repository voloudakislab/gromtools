// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include "sqlite3.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

enum class ColType {
  text,
  numeric
};

struct ColMeta {
  std::string name;
  std::string decl_type;
  ColType col_type;
};

struct QueryCol {
  std::string expr;
  std::string alias;
  ColType col_type;
};

struct ResultCol {
  std::string name;
  ColType col_type;
  std::vector<std::string> text_values;
  std::vector<double> num_values;
  std::vector<int> is_na;
};

struct SqliteDbCloser {
  void operator()(sqlite3* db) const {
    if (db != nullptr) {
      sqlite3_close(db);
    }
  }
};

struct SqliteStmtFinalizer {
  void operator()(sqlite3_stmt* stmt) const {
    if (stmt != nullptr) {
      sqlite3_finalize(stmt);
    }
  }
};

std::string basename_cpp(const std::string& path) {
  const std::string::size_type pos = path.find_last_of("/\\");
  std::string out = (pos == std::string::npos) ? path : path.substr(pos + 1);
  if (out.size() >= 3 && out.substr(out.size() - 3) == ".db") {
    out.resize(out.size() - 3);
  }
  return out;
}

std::string quote_ident(const std::string& ident) {
  std::string out = "\"";
  for (char ch : ident) {
    if (ch == '"') {
      out += "\"\"";
    } else {
      out.push_back(ch);
    }
  }
  out.push_back('"');
  return out;
}

std::string make_unique_alias(
    const std::string& base,
    const std::string& suffix,
    std::set<std::string>& used_aliases) {
  std::string alias = base;
  if (used_aliases.count(alias) > 0) {
    alias = base + suffix;
  }
  int i = 2;
  while (used_aliases.count(alias) > 0) {
    alias = base + suffix + "_" + std::to_string(i);
    ++i;
  }
  used_aliases.insert(alias);
  return alias;
}

std::string join_strings(const std::vector<std::string>& values, const std::string& sep) {
  std::string out;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i > 0) {
      out += sep;
    }
    out += values[i];
  }
  return out;
}

bool starts_with_chr(const std::string& value) {
  return value.size() >= 3 && value[0] == 'c' && value[1] == 'h' && value[2] == 'r';
}

std::string normalize_chromosome(const std::string& value) {
  if (value.size() >= 3) {
    const char c0 = static_cast<char>(std::tolower(static_cast<unsigned char>(value[0])));
    const char c1 = static_cast<char>(std::tolower(static_cast<unsigned char>(value[1])));
    const char c2 = static_cast<char>(std::tolower(static_cast<unsigned char>(value[2])));
    if (c0 == 'c' && c1 == 'h' && c2 == 'r') {
      return std::string("chr") + value.substr(3);
    }
  }
  return std::string("chr") + value;
}

std::pair<std::string, double> parse_varid_chr_pos(const std::string& varid, const std::string& db_path) {
  const std::string::size_type first_us = varid.find('_');
  if (first_us == std::string::npos) {
    Rcpp::stop(
      "Cannot derive chromosome/position from weights.varID in '%s': '%s'",
      db_path.c_str(),
      varid.c_str()
    );
  }

  const std::string chr = varid.substr(0, first_us);
  if (!starts_with_chr(chr)) {
    Rcpp::stop(
      "weights.varID does not begin with a 'chr' prefix in '%s': '%s'",
      db_path.c_str(),
      varid.c_str()
    );
  }

  const std::string::size_type second_us = varid.find('_', first_us + 1);
  if (second_us == std::string::npos) {
    Rcpp::stop(
      "Cannot derive position from weights.varID in '%s': '%s'",
      db_path.c_str(),
      varid.c_str()
    );
  }

  const std::string pos_str = varid.substr(first_us + 1, second_us - first_us - 1);
  char* endptr = nullptr;
  const double pos = std::strtod(pos_str.c_str(), &endptr);
  if (endptr == pos_str.c_str() || *endptr != '\0') {
    Rcpp::stop(
      "weights.varID has a non-numeric position component in '%s': '%s'",
      db_path.c_str(),
      varid.c_str()
    );
  }

  return std::make_pair(normalize_chromosome(chr), pos);
}

ColType infer_col_type(const std::string& decl_type) {
  std::string upper = decl_type;
  std::transform(
    upper.begin(),
    upper.end(),
    upper.begin(),
    [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); }
  );

  if (upper.find("INT") != std::string::npos ||
      upper.find("REAL") != std::string::npos ||
      upper.find("FLOA") != std::string::npos ||
      upper.find("DOUB") != std::string::npos ||
      upper.find("NUM") != std::string::npos) {
    return ColType::numeric;
  }
  return ColType::text;
}

bool table_exists(sqlite3* db, const std::string& table_name) {
  const char* sql =
    "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?1 LIMIT 1";
  sqlite3_stmt* raw_stmt = nullptr;
  if (sqlite3_prepare_v2(db, sql, -1, &raw_stmt, nullptr) != SQLITE_OK) {
    Rcpp::stop("Failed to prepare sqlite_master lookup: %s", sqlite3_errmsg(db));
  }
  std::unique_ptr<sqlite3_stmt, SqliteStmtFinalizer> stmt(raw_stmt);

  if (sqlite3_bind_text(stmt.get(), 1, table_name.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
    Rcpp::stop("Failed to bind table name '%s': %s", table_name.c_str(), sqlite3_errmsg(db));
  }

  const int rc = sqlite3_step(stmt.get());
  if (rc == SQLITE_ROW) {
    return true;
  }
  if (rc == SQLITE_DONE) {
    return false;
  }
  Rcpp::stop("Failed to check for table '%s': %s", table_name.c_str(), sqlite3_errmsg(db));
  return false;
}

std::unordered_map<std::string, ColMeta> get_table_columns(sqlite3* db, const std::string& table_name) {
  std::string sql = "PRAGMA table_info(" + quote_ident(table_name) + ")";
  sqlite3_stmt* raw_stmt = nullptr;
  if (sqlite3_prepare_v2(db, sql.c_str(), -1, &raw_stmt, nullptr) != SQLITE_OK) {
    Rcpp::stop("Failed to inspect table '%s': %s", table_name.c_str(), sqlite3_errmsg(db));
  }
  std::unique_ptr<sqlite3_stmt, SqliteStmtFinalizer> stmt(raw_stmt);

  std::unordered_map<std::string, ColMeta> cols;
  while (true) {
    const int rc = sqlite3_step(stmt.get());
    if (rc == SQLITE_DONE) {
      break;
    }
    if (rc != SQLITE_ROW) {
      Rcpp::stop("Failed to read schema for table '%s': %s", table_name.c_str(), sqlite3_errmsg(db));
    }

    const unsigned char* name_ptr = sqlite3_column_text(stmt.get(), 1);
    const unsigned char* type_ptr = sqlite3_column_text(stmt.get(), 2);
    const std::string name = name_ptr == nullptr ? "" : reinterpret_cast<const char*>(name_ptr);
    const std::string decl_type = type_ptr == nullptr ? "" : reinterpret_cast<const char*>(type_ptr);

    cols[name] = ColMeta{name, decl_type, infer_col_type(decl_type)};
  }

  return cols;
}

void require_columns(
    const std::unordered_map<std::string, ColMeta>& cols,
    const std::string& table_name,
    const std::vector<std::string>& required_cols) {
  std::vector<std::string> missing;
  for (const std::string& col : required_cols) {
    if (cols.find(col) == cols.end()) {
      missing.push_back(col);
    }
  }
  if (!missing.empty()) {
    Rcpp::stop(
      "Missing required column(s) in table '%s': %s",
      table_name.c_str(),
      join_strings(missing, ", ").c_str()
    );
  }
}

Rcpp::List build_result_frame(const std::vector<ResultCol>& cols, int n_rows) {
  Rcpp::List out(cols.size());
  Rcpp::CharacterVector out_names(cols.size());

  for (std::size_t i = 0; i < cols.size(); ++i) {
    out_names[i] = cols[i].name;
    if (cols[i].col_type == ColType::text) {
      Rcpp::CharacterVector vec(n_rows);
      for (int j = 0; j < n_rows; ++j) {
        if (cols[i].is_na[j]) {
          vec[j] = NA_STRING;
        } else {
          vec[j] = cols[i].text_values[j];
        }
      }
      out[i] = vec;
    } else {
      Rcpp::NumericVector vec(n_rows);
      for (int j = 0; j < n_rows; ++j) {
        if (cols[i].is_na[j]) {
          vec[j] = NA_REAL;
        } else {
          vec[j] = cols[i].num_values[j];
        }
      }
      out[i] = vec;
    }
  }

  out.attr("names") = out_names;
  out.attr("class") = "data.frame";
  out.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n_rows);
  return out;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List sqlite_read_model_db_cpp(
    const std::string& db_path,
    Rcpp::Nullable<Rcpp::CharacterVector> extra_cols = R_NilValue) {
  sqlite3* raw_db = nullptr;
  const int open_rc = sqlite3_open_v2(
    db_path.c_str(),
    &raw_db,
    SQLITE_OPEN_READONLY,
    nullptr
  );
  std::unique_ptr<sqlite3, SqliteDbCloser> db(raw_db);

  if (open_rc != SQLITE_OK || raw_db == nullptr) {
    const char* msg = raw_db == nullptr ? "unknown sqlite open failure" : sqlite3_errmsg(raw_db);
    Rcpp::stop("Failed to open SQLite database '%s': %s", db_path.c_str(), msg);
  }

  if (!table_exists(db.get(), "weights")) {
    Rcpp::stop("Required table 'weights' was not found in '%s'.", db_path.c_str());
  }
  if (!table_exists(db.get(), "extra")) {
    Rcpp::stop("Required table 'extra' was not found in '%s'.", db_path.c_str());
  }

  const auto weights_cols = get_table_columns(db.get(), "weights");
  const auto extra_table_cols = get_table_columns(db.get(), "extra");

  require_columns(weights_cols, "weights", {"gene", "rsid", "weight", "ref_allele", "eff_allele"});
  require_columns(extra_table_cols, "extra", {"gene"});

  const bool have_extra_chr = extra_table_cols.find("chr") != extra_table_cols.end();
  const bool have_weights_pos = weights_cols.find("pos") != weights_cols.end();
  const bool have_weights_varid = weights_cols.find("varID") != weights_cols.end();

  if ((!have_extra_chr || !have_weights_pos) && !have_weights_varid) {
    Rcpp::stop(
      "Database '%s' needs weights.varID when extra.chr or weights.pos is missing.",
      db_path.c_str()
    );
  }

  std::vector<QueryCol> query_cols = {
    {"w." + quote_ident("gene"), "gene", ColType::text},
    {"w." + quote_ident("rsid"), "rsid", ColType::text}
  };

  if (have_extra_chr) {
    query_cols.push_back(QueryCol{
      "e." + quote_ident("chr"),
      "chromosomes",
      extra_table_cols.at("chr").col_type
    });
  } else {
    query_cols.push_back(QueryCol{
      "w." + quote_ident("varID"),
      "chromosomes",
      ColType::text
    });
  }

  if (have_weights_pos) {
    query_cols.push_back(QueryCol{
      "w." + quote_ident("pos"),
      "position",
      weights_cols.at("pos").col_type
    });
  } else {
    query_cols.push_back(QueryCol{
      "w." + quote_ident("varID"),
      "position",
      ColType::text
    });
  }

  query_cols.push_back(QueryCol{"w." + quote_ident("ref_allele"), "ref_allele", ColType::text});
  query_cols.push_back(QueryCol{"w." + quote_ident("eff_allele"), "eff_allele", ColType::text});
  query_cols.push_back(QueryCol{"w." + quote_ident("weight"), "weight", weights_cols.at("weight").col_type});

  std::set<std::string> used_aliases = {
    "gene", "rsid", "chromosomes", "position", "ref_allele", "eff_allele", "weight", "model_ID"
  };
  std::vector<std::string> missing_extra_cols;

  if (extra_cols.isNotNull()) {
    Rcpp::CharacterVector extra_cols_vec(extra_cols);
    std::set<std::string> seen_requested;

    for (R_xlen_t i = 0; i < extra_cols_vec.size(); ++i) {
      if (extra_cols_vec[i] == NA_STRING) {
        continue;
      }

      const std::string requested = Rcpp::as<std::string>(extra_cols_vec[i]);
      if (requested.empty()) {
        continue;
      }
      if (!seen_requested.insert(requested).second) {
        continue;
      }

      const bool in_weights = weights_cols.find(requested) != weights_cols.end();
      const bool in_extra = extra_table_cols.find(requested) != extra_table_cols.end();

      if (!in_weights && !in_extra) {
        missing_extra_cols.push_back(requested);
        continue;
      }

      if (in_weights && in_extra) {
        query_cols.push_back(QueryCol{
          "w." + quote_ident(requested),
          make_unique_alias(requested + "_weights", "", used_aliases),
          weights_cols.at(requested).col_type
        });
        query_cols.push_back(QueryCol{
          "e." + quote_ident(requested),
          make_unique_alias(requested + "_extra", "", used_aliases),
          extra_table_cols.at(requested).col_type
        });
      } else if (in_weights) {
        query_cols.push_back(QueryCol{
          "w." + quote_ident(requested),
          make_unique_alias(requested, "_weights", used_aliases),
          weights_cols.at(requested).col_type
        });
      } else {
        query_cols.push_back(QueryCol{
          "e." + quote_ident(requested),
          make_unique_alias(requested, "_extra", used_aliases),
          extra_table_cols.at(requested).col_type
        });
      }
    }
  }

  if (!missing_extra_cols.empty()) {
    Rcpp::warning(
      "Requested extra_cols not found in either table for '%s': %s",
      db_path.c_str(),
      join_strings(missing_extra_cols, ", ").c_str()
    );
  }

  std::string sql = "SELECT ";
  for (std::size_t i = 0; i < query_cols.size(); ++i) {
    if (i > 0) {
      sql += ", ";
    }
    sql += query_cols[i].expr + " AS " + quote_ident(query_cols[i].alias);
  }
  sql += " FROM " + quote_ident("weights") + " w";
  sql += " LEFT JOIN " + quote_ident("extra") + " e";
  sql += " ON w." + quote_ident("gene") + " = e." + quote_ident("gene");

  sqlite3_stmt* raw_stmt = nullptr;
  if (sqlite3_prepare_v2(db.get(), sql.c_str(), -1, &raw_stmt, nullptr) != SQLITE_OK) {
    Rcpp::stop("Failed to prepare SQLite query for '%s': %s", db_path.c_str(), sqlite3_errmsg(db.get()));
  }
  std::unique_ptr<sqlite3_stmt, SqliteStmtFinalizer> stmt(raw_stmt);

  std::vector<ResultCol> result_cols;
  result_cols.reserve(query_cols.size() + 1);
  result_cols.push_back(ResultCol{"model_ID", ColType::text, {}, {}, {}});
  for (const QueryCol& col : query_cols) {
    const ColType out_type = (col.alias == "position") ? ColType::numeric : col.col_type;
    result_cols.push_back(ResultCol{col.alias, out_type, {}, {}, {}});
  }

  const std::string model_id = basename_cpp(db_path);
  int n_rows = 0;

  while (true) {
    const int rc = sqlite3_step(stmt.get());
    if (rc == SQLITE_DONE) {
      break;
    }
    if (rc != SQLITE_ROW) {
      Rcpp::stop("Failed while reading '%s': %s", db_path.c_str(), sqlite3_errmsg(db.get()));
    }

    result_cols[0].text_values.push_back(model_id);
    result_cols[0].is_na.push_back(0);

    std::vector<int> col_is_null(query_cols.size(), 0);
    std::vector<std::string> col_text(query_cols.size());
    std::vector<double> col_num(query_cols.size(), NA_REAL);

    for (int i = 0; i < static_cast<int>(query_cols.size()); ++i) {
      const int sqlite_type = sqlite3_column_type(stmt.get(), i);
      if (sqlite_type == SQLITE_NULL) {
        col_is_null[static_cast<std::size_t>(i)] = 1;
        continue;
      }

      if (query_cols[static_cast<std::size_t>(i)].col_type == ColType::text) {
        const unsigned char* value = sqlite3_column_text(stmt.get(), i);
        col_text[static_cast<std::size_t>(i)] =
          value == nullptr ? "" : reinterpret_cast<const char*>(value);
      } else {
        col_num[static_cast<std::size_t>(i)] = sqlite3_column_double(stmt.get(), i);
      }
    }

    int out_idx = 1;
    for (int i = 0; i < static_cast<int>(query_cols.size()); ++i) {
      const QueryCol& qcol = query_cols[static_cast<std::size_t>(i)];
      ResultCol& out_col = result_cols[static_cast<std::size_t>(out_idx)];
      ++out_idx;

      if (qcol.alias == "chromosomes" && !have_extra_chr) {
        if (col_is_null[static_cast<std::size_t>(i)]) {
          out_col.is_na.push_back(1);
          out_col.text_values.emplace_back();
        } else {
          const auto chr_pos = parse_varid_chr_pos(col_text[static_cast<std::size_t>(i)], db_path);
          out_col.is_na.push_back(0);
          out_col.text_values.push_back(chr_pos.first);
        }
        continue;
      }

      if (qcol.alias == "chromosomes" && have_extra_chr) {
        if (col_is_null[static_cast<std::size_t>(i)]) {
          out_col.is_na.push_back(1);
          out_col.text_values.emplace_back();
        } else {
          out_col.is_na.push_back(0);
          out_col.text_values.push_back(
            normalize_chromosome(col_text[static_cast<std::size_t>(i)])
          );
        }
        continue;
      }

      if (qcol.alias == "position" && !have_weights_pos) {
        if (col_is_null[static_cast<std::size_t>(i)]) {
          out_col.is_na.push_back(1);
          out_col.num_values.push_back(NA_REAL);
        } else {
          const auto chr_pos = parse_varid_chr_pos(col_text[static_cast<std::size_t>(i)], db_path);
          out_col.is_na.push_back(0);
          out_col.num_values.push_back(chr_pos.second);
        }
        continue;
      }

      if (col_is_null[static_cast<std::size_t>(i)]) {
        out_col.is_na.push_back(1);
        if (out_col.col_type == ColType::text) {
          out_col.text_values.emplace_back();
        } else {
          out_col.num_values.push_back(NA_REAL);
        }
      } else {
        out_col.is_na.push_back(0);
        if (out_col.col_type == ColType::text) {
          out_col.text_values.push_back(col_text[static_cast<std::size_t>(i)]);
        } else {
          out_col.num_values.push_back(col_num[static_cast<std::size_t>(i)]);
        }
      }
    }

    ++n_rows;
  }

  return build_result_frame(result_cols, n_rows);
}
