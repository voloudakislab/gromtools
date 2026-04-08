# ================================ build_csc_triplets() ================================ #
build_csc_triplets <- function(
  grom_pfx,
  pgen_dir,
  variant_weights
) {
  # Define main output directories
  out_dir <- dirname(grom_pfx)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tools_dir = file.path(out_dir, "tools")
  csc_root = file.path(out_dir, "tools/csc_triplets")
  extras_dir = file.path(tools_dir, "extra")
  grom_dir <- out_dir
  
  # Define grom's file name (for .gid file)
  fname <- basename(grom_pfx)
  grom_file_name  <- tools::file_path_sans_ext(fname) 



  ## ---------- setup dirs ----------
  dir.create(tools_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(csc_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(extras_dir, recursive = TRUE, showWarnings = FALSE)

  # Extract chromosome tags (used to source pvar file)
  CHROMOSOME_TAGS <- variant_weights$chromosomes
  all_pvar_files <- list.files(pgen_dir, full.names = TRUE, pattern = "[.]pvar$")

  # do a mapping sanity check, just like in impute_grom()
  chr_to_pvar <- map_chrom_to_pvar(
    CHROMS         = variant_weights$chromosomes,
    all_pvar_files = all_pvar_files
  )

  # For debugging
  meta_dir <- file.path(out_dir, "meta")
  dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(chr_to_pvar,
                   file.path(meta_dir, "chr_to_pvar_map.tsv"),
                   sep = "\t")

  ## ---------- load + normalize weights ----------
  has_ancestry <- "ancestry" %in% names(variant_weights)
  cols_to_keep <- c("model_ID", "gene", "rsid", "weight")
  if (has_ancestry) {
    cols_to_keep <- c("ancestry", cols_to_keep)
  }
  pre_binary <- variant_weights[, ..cols_to_keep]
  pre_binary[, model_ID := gsub(
    "230719_MegaAnalysis_|230721_MegaAnalysis_|230802_MegaAnalysis_|_prediXcan_noprior_alpha0.5_window1e6_filtered",
    "",
    model_ID
  )]

  ## ---------- optional ancestry metadata ----------
  dtg <- pre_binary
  anc <- NULL
  if (has_ancestry) {
    anc <- dtg$ancestry[1L]
  }

  # dedup weights rows we care about
  dtg <- unique(dtg[, .(model_ID, gene, rsid, weight)])

  # build (model,gene) row space ONCE per ancestry
  model_index <- data.table(model_ID = sort(unique(dtg$model_ID)))
  model_index[, model_id := seq_len(.N) - 1L]

  gene_index <- data.table(gene = sort(unique(dtg$gene)))
  gene_index[, gene_id := seq_len(.N) - 1L]

  mg_index <- unique(dtg[, .(model_ID, gene)])
  mg_index <- merge(mg_index, model_index, by = "model_ID", sort = FALSE)
  mg_index <- merge(mg_index, gene_index,  by = "gene",     sort = FALSE)
  setorder(mg_index, model_id, gene_id)
  mg_index[, mg_id := seq_len(.N) - 1L]
  R <- nrow(mg_index)

  # Write lookup tables once
  fwrite(model_index, file.path(tools_dir, "model_index.tsv.gz"))
  fwrite(gene_index,  file.path(tools_dir, "gene_index.tsv.gz"))
  fwrite(
    mg_index[, .(mg_id, model_id, gene_id, model_ID, gene)],
    file.path(grom_dir, paste0(grom_file_name, ".gid")), 
    sep = "\t"
  )

  # add pvar mapping check


  ## ---------- chromosome loop ----------
  for (i in seq_len(nrow(chr_to_pvar))) {
    CHR <- chr_to_pvar$chromosome[i]
    pvar_path <- all_pvar_files[basename(all_pvar_files) == chr_to_pvar$pvar_file[i]]

    chr_tag_clean <- gsub("_", "", CHR)

    # Read pvar (keeps 0-based snp_id in pvar order)
    pvar <- read_pvar(pvar_path)
    K    <- nrow(pvar)

    # subdirs per chromosome
    extras_dir_chr <- file.path(extras_dir, chr_tag_clean)
    csc_root_chr   <- file.path(csc_root, chr_tag_clean)
    dir.create(extras_dir_chr, recursive = TRUE, showWarnings = FALSE)
    dir.create(csc_root_chr, recursive = TRUE, showWarnings = FALSE)

    # write snp index for this chromosome in EXACT pvar order
    fwrite(
      pvar[, .(snp_id, rsid = ID, CHROM, POS, REF, ALT)],
      file.path(extras_dir_chr, "snp_index.tsv.gz")
    )

    # map weights → (mg_id, snp_id) for SNPs present on this chromosome
    dt_chr <- merge(
      dtg,
      pvar[, .(rsid = ID, snp_id)],
      by = "rsid",
      all = FALSE,
      sort = FALSE
    )

    if (nrow(dt_chr) == 0L) {
      # still emit empty CSC for completeness
      indptr   <- integer(K + 1L)
      indices  <- integer(0L)
      csc_data <- numeric(0L)
    } else {
      dt_chr <- merge(
        dt_chr,
        mg_index[, .(model_ID, gene, mg_id)],
        by = c("model_ID", "gene"),
        all.x = TRUE,
        sort = FALSE
      )

      # collapse duplicates and order CSC-friendly
      edges <- dt_chr[, .(weight = sum(weight)), by = .(snp_id, mg_id)]
      setorder(edges, snp_id, mg_id)
      setcolorder(edges, c("snp_id", "mg_id", "weight"))

      # export human-readable edges for this chr
      fwrite(edges, file.path(extras_dir_chr, "edges.tsv.gz"))

      # build CSC triplets
      counts <- edges[, .N, by = snp_id][order(snp_id)]
      counts_full <- integer(K)
      counts_full[counts$snp_id + 1L] <- counts$N

      indptr   <- c(0L, as.integer(cumsum(counts_full)))   # length K+1
      indices  <- as.integer(edges$mg_id)                  # nnz
      csc_data <- as.numeric(edges$weight)                 # nnz
    }

    # binary dumps
    con <- file(file.path(csc_root_chr, "indptr.u32.bin"), "wb")
    writeBin(as.integer(indptr), con, size = 4, endian = "little")
    close(con)

    con <- file(file.path(csc_root_chr, "indices.u32.bin"), "wb")
    writeBin(as.integer(indices), con, size = 4, endian = "little")
    close(con)

    con <- file(file.path(csc_root_chr, "data.f32.bin"), "wb")
    writeBin(as.numeric(csc_data), con, size = 4, endian = "little")
    close(con)

    # manifest for this chromosome
    man <- list(
      format     = "CSC",
      rows       = "mg_id",
      columns    = "snp_id (pvar order)",
      dtype      = list(
        indptr  = "uint32_le",
        indices = "uint32_le",
        data    = "float32_le"
      ),
      dims       = list(
        n_rows = R,
        n_snps = K,
        nnz    = length(indices)
      ),
      files      = list(
        indptr  = "indptr.u32.bin",
        indices = "indices.u32.bin",
        data    = "data.f32.bin"
      ),
      chromosome = as.character(chr_tag_clean),
      pvar_path  = pvar_path
    )
    if (!is.null(anc)) {
      man$ancestry <- anc
    }

    write_json(
      man,
      file.path(csc_root_chr, "manifest.json"),
      auto_unbox = TRUE,
      pretty     = TRUE
    )
  }

  invisible(TRUE)
}

# Optional helper (still unused but kept from your original)
safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# Helper: read a .pvar and return pvar table with exact order and 0-based snp_id
read_pvar <- function(pvar_path) {
  pvar <- fread(pvar_path, sep = "\t", header = TRUE, skip = "#CHROM")
  setnames(pvar, old = intersect("#CHROM", colnames(pvar)), new = "CHROM")
  if (!"CHROM" %in% names(pvar)) setnames(pvar, old = "CHROM", new = "CHROM")
  if (!"POS"   %in% names(pvar)) stop("PVAR missing POS column")
  if (!"ID"    %in% names(pvar)) stop("PVAR missing ID column (rsid)")
  if (!"REF"   %in% names(pvar)) pvar[, REF := NA_character_]
  if (!"ALT"   %in% names(pvar)) pvar[, ALT := NA_character_]

  pvar[, snp_id := .I - 1L]  # 0-based, EXACT file order
  pvar[, .(snp_id, CHROM, POS, ID, REF, ALT)]
}



map_chrom_to_pvar <- function(CHROMS, all_pvar_files) {
  # ---- global sanity checks ----
  if (length(CHROMS) == 0L) {
    stop(
      "Input CHROMS is empty. There are no chromosome values to map.\n",
      "Check the source table of variants weights.",
      call. = FALSE
    )
  }

  if (length(all_pvar_files) == 0L) {
    stop(
      "Input all_pvar_files is empty. No .pvar files were found to map against.",
      call. = FALSE
    )
  }

  # drop duplicates, but keep original values for error messages
  CHROMS <- unique(CHROMS)

  # check for NA / empty chromosome codes
  bad_na <- CHROMS[is.na(CHROMS) | CHROMS == ""]
  if (length(bad_na) > 0L) {
    stop(
      "Found NA or empty chromosome values in the table of variants weights:\n",
      paste(bad_na, collapse = ", "),
      "\nPlease fix the 'chromosome' column.",
      call. = FALSE
    )
  }

  # precompute basenames once
  pvar_basenames <- basename(all_pvar_files)

  # ---- per-chromosome mapping ----
  out <- setNames(character(length(CHROMS)), CHROMS)

  for (CHR in CHROMS) {
    # If CHR looks like a chromosome code (e.g. "1", "22", "chr1", "CHR22"),
    # enforce "no extra trailing digit" (so "chr1" != "chr10").
    if (grepl("^(chr|CHR)?[0-9]+$", CHR)) {
      pat  <- paste0("\\Q", CHR, "\\E(?![0-9])")  # escape CHR, then negative lookahead
      hits <- grepl(pat, pvar_basenames, perl = TRUE)
    } else {
      # Otherwise treat CHR as a plain literal substring (e.g. "ukb18")
      hits <- grepl(CHR, pvar_basenames, fixed = TRUE)
    }

    n_hits <- sum(hits)

    if (n_hits == 0L) {
      stop(
        "Chromosome value '", CHR,
        "' from the table of variants weights ",
        "could not be mapped to any of the .pvar files below:\n",
        paste(all_pvar_files, collapse = "\n"),
        call. = FALSE
      )
    }

    if (n_hits > 1L) {
      stop(
        "Chromosome value '", CHR,
        "' from the table of variants weights ",
        "matched ", n_hits, " .pvar files, but a unique mapping is required.\n",
        "Matching files (up to first 10):\n",
        paste(head(all_pvar_files[hits], 10L), collapse = "\n"),
        call. = FALSE
      )
    }

    # store only the basename for the hit
    out[CHR] <- pvar_basenames[hits]
  }

  # ---- convert named vector to table ----
  data.table::data.table(
    chromosome = names(out),
    pvar_file  = unname(out)
  )
}






# ================================ impute_grom() ================================ #

impute_grom <- function(
  pgen_dir,
  grom_pfx,
  variant_weights,
  snp_chunk     = 1000L,
  sample_subset = integer(0),
  CHUNK         = 1000L,
  exportChunk   = 256L,
  meanimpute    = TRUE,
  is_round      = TRUE
) {
  message("### impute_grom() started at: ", Sys.time())
  # Define main output directories
  out_dir <- dirname(grom_pfx)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tools_dir = file.path(out_dir, "tools")
  csc_root = file.path(out_dir, "tools/csc_triplets")
  extras_dir = file.path(tools_dir, "extra")
  grom_dir <- out_dir
  
  # Define grom's file name (for .gid file)
  fname <- basename(grom_pfx)
  grom_file_name  <- tools::file_path_sans_ext(fname) 
  
  # Extract Chromosome tags (used to source pgen files)
  CHROMOSOME_TAGS <- unique(variant_weights$chromosomes)

  # Directories
  dir.create(grom_dir, showWarnings = FALSE, recursive = TRUE)

  meta_dir <- file.path(out_dir, "meta")
  dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)

  # PGEN / PSAM / mg_index
  all_pgen_files <- list.files(pgen_dir, full.names = TRUE, pattern = "[.]pgen$")
  if (length(all_pgen_files) == 0L) {
    stop("No .pgen files found in pgen_dir: ", pgen_dir)
  }

  psam_file <- list.files(pgen_dir, full.names = TRUE, pattern = "[.]psam$")
  if (length(psam_file) == 0L) {
    stop("No .psam file found in pgen_dir: ", pgen_dir)
  }
  psam_file <- psam_file[1L]

  mg_index_file <- list.files(grom_dir, full.names = TRUE, pattern = ".gid")
  if (length(mg_index_file) == 0L) {
    stop("No mg_index* file found in tools_dir: ", tools_dir)
  }
  mg_index_file <- mg_index_file[1L]

  # grom file path
  grom_file <- if (tools::file_ext(grom_file_name) == "grom") {
    file.path(grom_dir, grom_file_name)
  } else {
    file.path(grom_dir, paste0(grom_file_name, ".grom"))
  }

  # Dimensions
  n_total_genes <- nrow(data.table::fread(mg_index_file))
  psam <- data.table::fread(psam_file)
  n_inds <- nrow(psam)

  if (length(sample_subset) > 0L) {
    sample_subset <- unique(as.integer(sample_subset))
    sample_subset <- sample_subset[sample_subset >= 1L & sample_subset <= n_inds]
    if (length(sample_subset) == 0L) {
      stop("sample_subset produced 0 valid sample indices.")
    }
    sample_subset <- sort(sample_subset)
    psam_out <- psam[sample_subset]
  } else {
    psam_out <- psam
  }

  # sid file
  sid_file <- sub("[.]grom$", ".sid", grom_file_name)
  if (!grepl("[.]sid$", sid_file)) {
    sid_file <- paste0(grom_file_name, ".sid")
  }
  data.table::fwrite(psam_out, file.path(grom_dir, sid_file), sep = "\t")

  # Map chromosomes → pgen
  options(warning.length = 8170L)
  chr_to_pgen <- map_chrom_to_pgen(
    CHROMS         = CHROMOSOME_TAGS,
    all_pgen_files = all_pgen_files
  )
  data.table::fwrite(chr_to_pgen, file.path(meta_dir, "chr_to_pgen_map.tsv"), sep = "\t")
  message("All chromosome tags were uniquely mapped to a .pgen file.")

  # Main loop
  start <- Sys.time()
  for (i in seq_len(nrow(chr_to_pgen))) {
    chr_start <- Sys.time()

    CHR <- chr_to_pgen$chromosome[i]
    pgen_file_chr <- file.path(pgen_dir, chr_to_pgen$pgen_file[i])

    # CSC triplets
    path_indptr  <- file.path(csc_root, CHR, "indptr.u32.bin")
    path_indices <- file.path(csc_root, CHR, "indices.u32.bin")
    path_weights <- file.path(csc_root, CHR, "data.f32.bin")

    if (!file.exists(path_indptr)) {
      stop("CSC 'indptr.u32.bin' for chromosome '", CHR, "' not found in ", csc_root)
    }
    if (!file.exists(path_indices)) {
      stop("CSC 'indices.u32.bin' for chromosome '", CHR, "' not found in ", csc_root)
    }
    if (!file.exists(path_weights)) {
      stop("CSC 'data.f32.bin' for chromosome '", CHR, "' not found in ", csc_root)
    }

    # u32 → integer
    n_indptr  <- n_from_file(path_indptr)
    n_indices <- n_from_file(path_indices)
    n_weights <- n_from_file(path_weights)

    indptr_chr <- readBin(path_indptr,
                          what   = "integer",
                          n      = n_indptr,
                          size   = 4L,
                          endian = "little")

    indices_chr <- readBin(path_indices,
                           what   = "integer",
                           n      = n_indices,
                           size   = 4L,
                           endian = "little")

    # f32 → numeric
    weights_chr <- readBin(path_weights,
                           what   = "numeric",
                           n      = n_weights,
                           size   = 4L,
                           endian = "little")

    gene_pos_chr <- flag_start_end(indices_chr)

    message("Beginning calculations for chromosome-tag: ", CHR, "...")

    grom_axpy_engine(
      pgen_file     = pgen_file_chr,
      grom_file     = grom_file,
      gene_pos      = gene_pos_chr,
      n_total_genes = n_total_genes,
      indptr        = indptr_chr,
      indices       = indices_chr,
      data          = weights_chr,
      raw_sample_ct = n_inds,
      snp_chunk     = snp_chunk,
      sample_subset = sample_subset,
      CHUNK         = CHUNK,
      exportChunk   = exportChunk,
      meanimpute    = meanimpute,
      create_new    = (i == 1L),
      showWarnings  = FALSE,
      rounded_mean  = is_round
    )

    chr_end <- Sys.time()
    message("Total time for chromosome-tag ", CHR, ": ", chr_end - chr_start)
  }
  end <- Sys.time()

  # Sanity check of mapping
  check_chrom_pgen_mapping(chr_to_pgen, CHROMOSOME_TAGS, all_pgen_files, pgen_dir)

  total_time <- end - start
  message("impute_grom() completed. Total time: ", total_time)

}



flag_start_end <- function(indices) {
  n   <- length(indices)
  out <- integer(n)  # initialized to 0

  # group positions by index value
  groups <- split(seq_len(n), indices)

  for (pos in groups) {
    if (length(pos) == 1L) {
      # only occurrence → both start and end
      out[pos] <- 2L
    } else {
      # mark all occurrences of this gene as -1 (middle by default)
      out[pos] <- -1L
      # first occurrence → start
      out[pos[1]] <- 1L
      # last occurrence → end
      out[pos[length(pos)]] <- 0L
    }
  }
  out
}

## Helper to get length = file_size / 4 bytes
n_from_file <- function(path) {
  sz <- file.info(path)$size
  stopifnot(sz %% 4L == 0L)
  as.integer(sz / 4L)
}

map_chrom_to_pgen <- function(CHROMS, all_pgen_files) {
  # ---- global sanity checks ----
  if (length(CHROMS) == 0L) {
    stop(
      "Input CHROMS is empty. There are no chromosome values to map.\n",
      "Check the source table of variants weights.",
      call. = FALSE
    )
  }

  if (length(all_pgen_files) == 0L) {
    stop(
      "Input all_pgen_files is empty. No .pgen files were found to map against.",
      call. = FALSE
    )
  }

  # drop duplicates, but keep original values for error messages
  CHROMS <- unique(CHROMS)

  # check for NA / empty chromosome codes
  bad_na <- CHROMS[is.na(CHROMS) | CHROMS == ""]
  if (length(bad_na) > 0L) {
    stop(
      "Found NA or empty chromosome values in the table of variants weights:\n",
      paste(bad_na, collapse = ", "),
      "\nPlease fix the 'chromosome' column.",
      call. = FALSE
    )
  }

  # precompute basenames once
  pgen_basenames <- basename(all_pgen_files)

  # ---- per-chromosome mapping ----
  out <- setNames(character(length(CHROMS)), CHROMS)

  for (CHR in CHROMS) {
    # If CHR looks like a chromosome code (e.g. "1", "22", "chr1", "CHR22"),
    # enforce "no extra trailing digit" (so "chr1" != "chr10").
    if (grepl("^(chr|CHR)?[0-9]+$", CHR)) {
      pat  <- paste0("\\Q", CHR, "\\E(?![0-9])")  # escape CHR, then negative lookahead
      hits <- grepl(pat, pgen_basenames, perl = TRUE)
    } else {
      # Otherwise treat CHR as a plain literal substring (e.g. "ukb18")
      hits <- grepl(CHR, pgen_basenames, fixed = TRUE)
    }

    n_hits <- sum(hits)

    if (n_hits == 0L) {
      stop(
        "Chromosome value '", CHR,
        "' from the table of variants weights ",
        "could not be mapped to any of the .pgen files below:\n",
        paste(all_pgen_files, collapse = "\n"),
        call. = FALSE
      )
    }

    if (n_hits > 1L) {
      stop(
        "Chromosome value '", CHR,
        "' from the table of variants weights ",
        "matched ", n_hits, " .pgen files, but a unique mapping is required.\n",
        "Matching files (up to first 10):\n",
        paste(head(all_pgen_files[hits], 10L), collapse = "\n"),
        call. = FALSE
      )
    }

    # store only the basename for the hit
    out[CHR] <- pgen_basenames[hits]
  }

  # ---- convert named vector to table ----
  data.table::data.table(
    chromosome = names(out),
    pgen_file  = unname(out)
  )
}


check_chrom_pgen_mapping <- function(chr_to_pgen,
                                     CHROMS,
                                     all_pgen_files,
                                     pgen_dir) {
  # 1) Report chromosome tags (no dependency on a specific weights filename)
  ch_unique <- unique(CHROMS)
  message(
    "Chromosome tags detected in weights table:\n",
    paste(ch_unique, collapse = ", ")
  )

  # 2) Check that all CHROMS are present in the mapping table (first column)
  mapped_chrs   <- chr_to_pgen[["chromosome"]]
  missing_chrom <- setdiff(ch_unique, mapped_chrs)

  if (length(missing_chrom) == 0L) {
    message("All chromosome tags of column 'chromosome' were mapped successfully.")
  } else {
    message(
      "The following chromosome tags could not be mapped:\n",
      paste(missing_chrom, collapse = ", ")
    )
  }

  # 3) Check that all .pgen basenames are present in the mapping table (second column)
  pgen_basenames <- basename(all_pgen_files)
  mapped_pgen    <- chr_to_pgen[["pgen_file"]]

  missing_idx <- which(!(pgen_basenames %in% mapped_pgen))

  if (length(missing_idx) == 0L) {
    message(
      "All .pgen files from directory '", pgen_dir,
      "' were successfully mapped."
    )
  } else {
    message(
      "The following .pgen files could not be found in the mapping table:\n",
      paste(all_pgen_files[missing_idx], collapse = "\n")
    )
  }

  # return missing items invisibly if you want to inspect programmatically
  invisible(list(
    missing_chrom_tags = missing_chrom,
    missing_pgen_files = if (length(missing_idx)) all_pgen_files[missing_idx] else character(0)
  ))
}









# ============================= read_grom_annotated() ============================= #


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

grom_streamer <- function(path, n_rows, n_cols, col_indices, samples = NULL) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)

  col_indices <- sort(unique(as.integer(col_indices)))
  col_indices <- col_indices[col_indices >= 0L & col_indices < n_cols]

  if (is.null(samples)) {
    row_idx <- seq_len(n_rows)
  } else {
    row_idx <- unique(as.integer(samples))
    row_idx <- row_idx[row_idx >= 1L & row_idx <= n_rows]
  }

  out <- matrix(NA_real_, nrow = length(row_idx), ncol = length(col_indices))
  use_targeted_rows <- !is.null(samples)

  for (j_out in seq_along(col_indices)) {
    j0 <- col_indices[[j_out]]
    column_offset <- j0 * n_rows * 8

    if (use_targeted_rows) {
      values <- numeric(length(row_idx))
      for (i in seq_along(row_idx)) {
        value_offset <- column_offset + (row_idx[[i]] - 1L) * 8
        seek(con, where = value_offset, origin = "start")
        values[[i]] <- readBin(con, what = "double", n = 1L, size = 8, endian = "little")
      }
      out[, j_out] <- values
    } else {
      seek(con, where = column_offset, origin = "start")
      full_col <- readBin(con, what = "double", n = n_rows, size = 8, endian = "little")
      out[, j_out] <- full_col[row_idx]
    }
  }

  out
}

read_grom <- function(prefix, extract = NULL, keep = NULL) {
  gid <- read_grom_gid(prefix)
  sid <- read_grom_sid(prefix)

  if (is.null(extract)) {
    gid_use <- copy(gid)
  } else {
    gid_use <- as.data.table(extract)
    if (!"mg_id" %in% names(gid_use)) {
      stop("extract must include mg_id in this tutorial version.")
    }
    gid_use <- unique(gid_use, by = "mg_id")
    setorder(gid_use, mg_id)
  }

  if (is.null(keep)) {
    sid_use <- copy(sid)
    sample_rows <- NULL
  } else {
    keep_idx <- match(keep, sid$IID)
    keep_idx <- keep_idx[!is.na(keep_idx)]
    sid_use <- sid[keep_idx]
    sample_rows <- sid_use$GROM_SID
  }

  mat <- grom_streamer(
    path = paste0(prefix, ".grom"),
    n_rows = nrow(sid),
    n_cols = max(gid$mg_id) + 1L,
    col_indices = gid_use$mg_id,
    samples = sample_rows
  )

  rownames(mat) <- sid_use$IID
  colnames(mat) <- paste(gid_use$model_ID, gid_use$gene, sep = "_")
  mat
}

gromtools_read <- function(grom_pfx, models = NULL, genes = NULL, samples = NULL) {
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

read_db_dir <- function(db_dir, extra_cols = NULL) {
  if (!is.character(db_dir) || length(db_dir) != 1L || is.na(db_dir) || !nzchar(db_dir)) {
    stop("db_dir must be a single non-empty character string.")
  }
  if (!dir.exists(db_dir)) {
    stop("db_dir does not exist: ", db_dir)
  }

  db_files <- list.files(
    db_dir,
    pattern = "[.]db$",
    full.names = TRUE
  )
  db_files <- sort(db_files)

  if (!length(db_files)) {
    stop("No .db files found in db_dir: ", db_dir)
  }

  tables <- lapply(
    db_files,
    function(db_path) sqlite_read_model_db(db_path, extra_cols = extra_cols)
  )

  data.table::rbindlist(tables, use.names = TRUE, fill = TRUE)
}


# ================================ gromtools_impute() ================================ #

gromtools_impute <- function(
  weights_table,
  grom_pfx,
  pgen_dir,
  snp_chunk     = 1000L,
  sample_subset = integer(0),
  CHUNK         = 1000L,
  exportChunk   = 256L,
  meanimpute    = TRUE,
  is_round      = TRUE
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
    is_round = is_round
  )
}






