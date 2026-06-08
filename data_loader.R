# ============================================================
# data_loader.R  v5
# ============================================================
# Decisions recorded:
#   Pck012, 013, 017  → EXCLUDED (splicing/skip per user)
#   Pck009            → RNA-seq + scRNA-seq (5 cell types)
#   Pck010            → Both 8h + 24h timepoints
#   Pck011            → mean per-patient logFC + Fisher combined p-value
#   Pck014            → TMT phosphoproteomics: GroupA-GroupB as log2FC, FDR
#   Pck015            → Top-down proteomics: logFC + adj.P.Val
#   Pck016            → cyt vs. NT comparison only
#   Pck018            → CT-NT vs NoCT-NT cytokine effect + Student T-test p
#   Pck019            → Cytokine effect + Parp12 KD effect (2 comparisons)
#   Pck020            → CT replicates mean as log2FC + one-sample t-test
#   Pck021            → Proteomics + Metabolomics (CT_Eth vs NCT_Eth t-test)
#   Pck022            → ln(FC) converted to log2FC; all 3 comparisons
#   Pck024            → 6h: IL-1β, IFNγ, IFNγ+IL-1β; 18h: IFNγ+IL-1β
#   Pck025            → Discovery + Validation RNA-seq + Proteomics (sex-stratified)
# ============================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
})

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !is.na(a[[1]]) && nzchar(a[[1]])) a else b
}

# ---- SC treatment extractor -------------------------------------------------
# Pulls the cytokine name out of verbose SC comparison labels such as:
#   "IL-1β (S2) vs. Control (S1)"          → "IL-1β"
#   "IFNγ + IL-1β (S4) vs. Control(S1) 6h" → "IFNγ + IL-1β"
#   "IFNγ + IL-1β + NMMA (S3) vs. ..."     → "IFNγ + IL-1β + NMMA"
extract_sc_treatment <- function(comp) {
  if (is.na(comp) || !nzchar(trimws(comp))) return(NA_character_)
  s <- trimws(comp)
  # Strip "(SN) vs. ..." suffix — keep everything before the first "(SN)"
  s <- sub("\\s*\\(S\\d+\\).*$", "", s, perl = TRUE)
  # Also strip any remaining "vs. ..." in case there was no (SN) label
  s <- sub("\\s+[Vv]s\\..*$", "", s, perl = TRUE)
  # Strip trailing time labels like "6h", "18h"
  s <- trimws(sub("\\s+\\d+h\\s*$", "", s, perl = TRUE))
  trimws(s)
}
extract_sc_treatment_v <- Vectorize(extract_sc_treatment)

# ---- Path helpers --------------------------------------------------------
resolve_data_dir <- function() {
  candidates <- c(
    Sys.getenv("MULTIOMICS_DATA_DIR"),
    file.path(getwd(), "Files"),
    file.path(getwd(), "..", "Files"),
    "/Users/sark224/Library/CloudStorage/OneDrive-PNNL/Documents/Projects/In-vitroCT treated exp-Shinnyapp/Files"
  )
  for (p in candidates)
    if (nzchar(p) && dir.exists(p) &&
        (file.exists(file.path(p, "Table 1-Metadata-Aap.xlsx")) ||
         file.exists(file.path(p, "Table 1-Metadata.xlsx"))))
      return(normalizePath(p))
  stop("Cannot find data directory. Set MULTIOMICS_DATA_DIR env var.")
}

# Pick the best available metadata file (new flat format preferred)
resolve_metadata_file <- function(data_dir) {
  new_file <- file.path(data_dir, "Table 1-Metadata-Aap.xlsx")
  old_file <- file.path(data_dir, "Table 1-Metadata.xlsx")
  if (file.exists(new_file)) return(new_file)
  if (file.exists(old_file)) return(old_file)
  stop("No metadata file found in ", data_dir)
}
pkg_path <- function(data_dir, id)
  file.path(data_dir, "Pck_Analysis", paste0(id, ".xlsx"))

# ---- Column-name normaliser ---------------------------------------------
normalize_colnames <- function(nms) {
  lnms <- tolower(trimws(as.character(nms)))
  find_col <- function(pats) {
    for (p in pats) { idx <- grep(p, lnms, perl = TRUE); if (length(idx)) return(idx[1]) }
    NA_integer_
  }
  list(
    gene = find_col(c(
      "external_gene_name",
      "gene.?symbols?$",        # gene_symbol, gene symbol, gene-symbol  ← check before names
      "gene.?names?$",          # gene_name, gene name, gene-name
      "genesymbol", "^symbol$", "^gene$", "^genes$", "^name$",
      "^column1$", "^protein$", "^proteins$", "majority.?protein",
      "metabolite", "featurename", "refseq.?accession"
    )),
    lfc  = find_col(c(
      "log2foldchange", "log2fc", "log2_fc", "avg[._]?log2[._]?fc",
      "log2.?fold", "log2.?\\[fold", "logfc[._]", "log2fc[._]",
      "log2_odds_ratio", "median.?logfc", "log2.?fold.?change",
      "log2.?\\(?fc",           # Log2(FC), log2(fc), log2 fc
      "^logfc$"                 # plain logFC with no suffix (e.g. Pck003 RNA-seq HI)
    )),
    pval = find_col(c(
      "^pvalue$", "^p_val$", "^p.value$", "^p-value$", "^pvalues$",
      "^pval$", "p[._]?value[._]?1$", "^anova.?p", "^t.test.?p",
      "student.*t.*p", "\\bp-?value\\b",
      "pvalue[._]",             # Pvalue_2h, Pvalue_8h etc. (Pck003 RNA-seq Endo)
      "p_val[._]"               # p_val_2h etc.
    )),
    padj = find_col(c(
      "padj", "adj[._]?p[._]?val", "p[._]?adj", "^fdr$", "^qvalues?$",
      "q.value", "q-value", "adjusted.?p", "adj.p.val",
      "^fdr[._]",               # FDR_2h, FDR_8h etc. (Pck003 RNA-seq Endo / HI)
      "p.val.adj"               # p_val_adj (Pck024 Seurat output)
    ))
  )
}

# ---- Standard output constructor ----------------------------------------
make_layer <- function(gene_name, log2FC, pvalue = NA_real_, padj = NA_real_,
                       comparison = "cytokine_vs_ctrl",
                       pkg_id, layer_name, omics_type,
                       model = NA_character_, treatment = NA_character_,
                       time_h = NA_character_, pmid = NA_character_,
                       repository = NA_character_,
                       is_sc = FALSE, cell_type = NA_character_,
                       has_stats = TRUE, replicate_averaged = FALSE,
                       is_genomic = FALSE) {
  n <- max(lengths(list(gene_name, log2FC, pvalue, padj, cell_type)))
  df <- data.frame(
    gene_name          = rep_len(as.character(gene_name), n),
    log2FC             = rep_len(suppressWarnings(as.numeric(log2FC)), n),
    pvalue             = rep_len(suppressWarnings(as.numeric(pvalue)), n),
    padj               = rep_len(suppressWarnings(as.numeric(padj)), n),
    comparison         = rep_len(as.character(comparison), n),
    pkg_id             = pkg_id, layer_name = layer_name, omics_type = omics_type,
    model              = model, treatment = treatment, time_h = time_h,
    pmid               = pmid, repository  = repository,
    is_sc              = is_sc,
    cell_type          = rep_len(as.character(cell_type), n),
    has_stats          = has_stats, replicate_averaged = replicate_averaged,
    is_genomic         = is_genomic,
    stringsAsFactors   = FALSE
  )
  df <- df[!is.na(df$gene_name) & nzchar(df$gene_name) & df$gene_name != "NA", ]
  df <- df[!is.na(df$log2FC), ]
  df
}

# ---- Auto-detecting header row -------------------------------------------
find_header_row <- function(raw, max_scan = 25) {
  n <- min(max_scan, nrow(raw) - 1L)
  for (i in seq_len(n)) {
    v <- as.character(unlist(raw[i, ]))
    v <- v[!is.na(v) & v != "NA" & nzchar(trimws(v))]
    if (length(v) < 2) next
    if (all(!is.na(suppressWarnings(as.numeric(v))))) next
    if (i + 1 <= nrow(raw)) {
      nv <- suppressWarnings(as.numeric(as.character(unlist(raw[i + 1, ]))))
      if (sum(!is.na(nv)) >= 1) return(i)
    }
  }
  1L
}

# ---- Generic single-block sheet parser ----------------------------------
parse_sheet <- function(path, sheet, skip = NULL, ln_to_log2 = FALSE) {
  raw_full <- suppressMessages(
    read_excel(path, sheet = sheet, col_names = FALSE, .name_repair = "minimal")
  )
  if (is.null(skip)) skip <- find_header_row(raw_full) - 1L
  df <- suppressMessages(
    read_excel(path, sheet = sheet, skip = skip, .name_repair = "minimal")
  )
  names(df) <- make.unique(as.character(names(df)))
  df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  nm <- normalize_colnames(names(df))
  if (is.na(nm$lfc)) return(NULL)
  gene_v <- if (!is.na(nm$gene)) as.character(df[[nm$gene]]) else NA_character_
  lfc_v  <- suppressWarnings(as.numeric(df[[nm$lfc]]))
  pval_v <- if (!is.na(nm$pval)) suppressWarnings(as.numeric(df[[nm$pval]])) else NA_real_
  padj_v <- if (!is.na(nm$padj)) suppressWarnings(as.numeric(df[[nm$padj]])) else NA_real_
  if (ln_to_log2) lfc_v <- lfc_v / log(2)
  list(gene = gene_v, lfc = lfc_v, pval = pval_v, padj = padj_v)
}

# ---- Wide-format parser: multiple comparison blocks side by side --------
parse_wide <- function(path, sheet, group_row_i, header_row_i, data_start_i,
                       ln_to_log2 = FALSE, first_block_only = FALSE) {
  raw <- suppressMessages(
    read_excel(path, sheet = sheet, col_names = FALSE, .name_repair = "minimal")
  )
  if (data_start_i > nrow(raw)) return(NULL)
  nc <- ncol(raw)

  grp <- as.character(unlist(raw[group_row_i, ]))
  grp[is.na(grp) | grp == "NA" | trimws(grp) == ""] <- NA_character_
  hdr <- as.character(unlist(raw[header_row_i, ]))
  hdr[is.na(hdr) | hdr == "NA" | trimws(hdr) == ""] <- NA_character_

  block_starts <- which(!is.na(grp))
  if (length(block_starts) == 0) block_starts <- 1L
  if (first_block_only) block_starts <- block_starts[1]
  block_ends   <- c(block_starts[-1] - 1L, nc)

  valid_idx <- vapply(seq_along(block_starts), function(i) {
    !is.na(normalize_colnames(hdr[block_starts[i]:block_ends[i]])$lfc)
  }, logical(1))
  block_starts  <- block_starts[valid_idx]
  block_ends    <- block_ends[valid_idx]
  group_labels  <- grp[block_starts]
  if (length(block_starts) == 0) return(NULL)

  data_rows <- raw[data_start_i:nrow(raw), , drop = FALSE]

  results <- lapply(seq_along(block_starts), function(i) {
    cols  <- block_starts[i]:block_ends[i]
    b_hdr <- hdr[cols]
    nm    <- normalize_colnames(b_hdr)
    if (is.na(nm$lfc)) return(NULL)
    if (!is.na(nm$gene)) {
      gene_v <- as.character(unlist(data_rows[, cols[nm$gene]]))
    } else if (i > 1) {
      # Try to borrow gene from block 1; if block 1 also has no named gene col,
      # fall back to block 1's first column (gene names without a header label).
      fc  <- block_starts[1]:block_ends[1]
      fn  <- normalize_colnames(hdr[fc])
      if (!is.na(fn$gene)) {
        gene_v <- as.character(unlist(data_rows[, fc[fn$gene]]))
      } else {
        gene_v <- as.character(unlist(data_rows[, block_starts[1]]))
      }
    } else {
      # Block 1 has no named gene column — use its first column as gene names.
      gene_v <- as.character(unlist(data_rows[, cols[1]]))
    }
    lfc_v  <- suppressWarnings(as.numeric(unlist(data_rows[, cols[nm$lfc]])))
    pval_v <- if (!is.na(nm$pval)) suppressWarnings(as.numeric(unlist(data_rows[, cols[nm$pval]]))) else rep(NA_real_, nrow(data_rows))
    padj_v <- if (!is.na(nm$padj)) suppressWarnings(as.numeric(unlist(data_rows[, cols[nm$padj]]))) else rep(NA_real_, nrow(data_rows))
    if (ln_to_log2) lfc_v <- lfc_v / log(2)
    comp <- trimws(group_labels[i]) %||% paste0("comparison_", i)
    data.frame(gene_name = gene_v, log2FC = lfc_v, pvalue = pval_v,
               padj = padj_v, comparison = comp, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out[!is.na(out$gene_name) & nzchar(out$gene_name) & !is.na(out$log2FC), ]
}

# ---- Replicate-averaging helper -----------------------------------------
avg_reps <- function(mat) {
  if (is.vector(mat)) return(mat)
  rowMeans(mat, na.rm = TRUE)
}

# ---- Pivot multi-timepoint columns (in-row) to long format --------------
pivot_timepoints <- function(df, gene_col_idx, lfc_prefix, padj_prefix,
                              timepoints, pval_prefix = NULL) {
  gene_v <- as.character(df[[gene_col_idx]])
  do.call(rbind, lapply(timepoints, function(tp) {
    lc  <- paste0(lfc_prefix,  tp)
    pc  <- paste0(padj_prefix, tp)
    pvc <- if (!is.null(pval_prefix)) paste0(pval_prefix, tp) else NULL
    if (!lc %in% names(df)) return(NULL)
    data.frame(
      gene_name  = gene_v,
      log2FC     = suppressWarnings(as.numeric(df[[lc]])),
      pvalue     = if (!is.null(pvc) && pvc %in% names(df)) suppressWarnings(as.numeric(df[[pvc]])) else NA_real_,
      padj       = if (pc %in% names(df)) suppressWarnings(as.numeric(df[[pc]])) else NA_real_,
      comparison = tp, stringsAsFactors = FALSE
    )
  }))
}

# ==========================================================================
# PACKAGE LOADERS
# ==========================================================================

load_pck001 <- function(path, meta) {
  layers <- list()
  d <- parse_sheet(path, "RNAseq", skip = 2)
  if (!is.null(d))
    layers$rnaseq <- make_layer(d$gene, d$lfc, d$pval, d$padj,
      pkg_id="Pck001", layer_name="RNA-seq", omics_type="RNA-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  d <- parse_sheet(path, "UMI-4C", skip = 3)
  if (!is.null(d))
    layers$umic <- make_layer(d$gene, d$lfc, d$pval, d$padj,
      pkg_id="Pck001", layer_name="UMI-4C", omics_type="UMI-4C",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_genomic=TRUE)
  d <- parse_sheet(path, "ATAC-seq", skip = 3)
  if (!is.null(d))
    layers$atac <- make_layer(d$gene, d$lfc, d$pval, d$padj,
      pkg_id="Pck001", layer_name="ATAC-seq", omics_type="ATAC-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_genomic=TRUE)
  d <- parse_sheet(path, "H3K27ac", skip = 3)
  if (!is.null(d))
    layers$h3k <- make_layer(d$gene, d$lfc, d$pval, d$padj,
      pkg_id="Pck001", layer_name="H3K27ac", omics_type="ChIP-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_genomic=TRUE)
  Filter(Negate(is.null), layers)
}

load_pck002 <- function(path, meta) {
  layers <- list()
  for (info in list(
    list(sh="RNA-seq",    sk=3, lt="RNA-seq",     ot="RNA-seq",     geo=FALSE),
    list(sh="ATAC-seq",   sk=3, lt="ATAC-seq",    ot="ATAC-seq",    geo=TRUE),
    list(sh="H3K27ac",    sk=3, lt="H3K27ac",     ot="ChIP-seq",    geo=TRUE),
    list(sh="Proteomics", sk=2, lt="Proteomics",  ot="Proteomics",  geo=FALSE),
    list(sh="Lipidomics", sk=5, lt="Lipidomics",  ot="Lipidomics",  geo=FALSE),
    list(sh="Metabolomics",sk=2,lt="Metabolomics",ot="Metabolomics",geo=FALSE)
  )) {
    if (!info$sh %in% excel_sheets(path)) next
    d <- parse_sheet(path, info$sh, skip = info$sk)
    if (!is.null(d))
      layers[[info$sh]] <- make_layer(d$gene, d$lfc, d$pval, d$padj,
        pkg_id="Pck002", layer_name=info$lt, omics_type=info$ot,
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository, is_genomic=info$geo)
  }
  Filter(Negate(is.null), layers)
}

load_pck003 <- function(path, meta) {
  # Pck003 has two cell models:
  #   RNA-seq HI   -> Human pancreatic islets
  #   RNA-seq Endo -> EndoC-betaH1
  # All other omics (ATAC-seq, Proteomics, Lipidomics, Metabolomics) are EndoC-betaH1
  # (confirmed from sheet headers; meta$model only captures the first metadata row)
  HI_MODEL   <- "Human pancreatic islets"
  ENDO_MODEL <- "EndoC-βH1"

  layers <- list()
  model_map <- c("RNA-seq HI" = HI_MODEL, "RNA-seq Endo" = ENDO_MODEL)
  for (sh in c("RNA-seq HI", "RNA-seq Endo")) {
    if (!sh %in% excel_sheets(path)) next
    d <- parse_wide(path, sh, 4, 5, 6)
    if (!is.null(d) && nrow(d) > 0)
      layers[[sh]] <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
        comparison=d$comparison, pkg_id="Pck003", layer_name=sh, omics_type="RNA-seq",
        model=model_map[[sh]], treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository)
  }
  d <- parse_sheet(path, "ATAC-seq")
  if (!is.null(d))
    layers$atac <- make_layer(d$gene, d$lfc, d$pval, d$padj,
      pkg_id="Pck003", layer_name="ATAC-seq", omics_type="ATAC-seq",
      model=ENDO_MODEL, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_genomic=TRUE)
  d <- parse_wide(path, "Proteomics",  4, 5, 6)
  if (!is.null(d) && nrow(d) > 0)
    layers$prot <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck003", layer_name="Proteomics", omics_type="Proteomics",
      model=ENDO_MODEL, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  d <- parse_wide(path, "Lipidomics",  3, 4, 5)
  if (!is.null(d) && nrow(d) > 0)
    layers$lipid <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck003", layer_name="Lipidomics", omics_type="Lipidomics",
      model=ENDO_MODEL, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  d <- parse_wide(path, "Metabolomics", 3, 4, 5)
  if (!is.null(d) && nrow(d) > 0)
    layers$metab <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck003", layer_name="Metabolomics", omics_type="Metabolomics",
      model=ENDO_MODEL, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  Filter(Negate(is.null), layers)
}

load_pck004 <- function(path, meta) {
  # Data 1-5 = scRNA-seq wide-format: ~33 side-by-side comparison×celltype blocks per sheet
  # Data 6   = ATAC-seq — excluded (incompatible structure)
  #
  # Two block styles (both anchored on the "stressor" column at index si):
  #   Alpha-style: gene(col 1 shared), pval(si-5), lfc(si-4), pct.1(si-3),
  #                pct.2(si-2), padj(si-1), stressor(si), celltype(si+1)
  #   Beta-style:  pval(si-6), lfc(si-5), pct.1(si-4), pct.2(si-3),
  #                padj(si-2), gene(si-1), stressor(si), celltype(si+1)
  #   Detection: if col si-1 is mostly character → beta-style, else alpha-style.
  #
  # Stressor codes → normalised cytokine treatment (BFA/TG = non-cytokine → NA)
  stressor_norm <- c(
    BFA      = NA_character_,          # Brefeldin A  — ER stressor, non-cytokine
    TG       = NA_character_,          # Thapsigargin — ER stressor, non-cytokine
    CM       = "IFNγ + IL-1β + TNFα", # cytokine mix (also appears as "Cm")
    Cm       = "IFNγ + IL-1β + TNFα",
    IFNG     = "IFNγ",
    IL1B     = "IL-1β",
    Il1B     = "IL-1β",               # case variant
    IL1BIFNG = "IFNγ + IL-1β",
    Il1BIFNG = "IFNγ + IL-1β",        # case variant
    TNFA     = "TNFα"
  )

  all_sheets <- grep("^Data", excel_sheets(path), value = TRUE)
  sheets     <- setdiff(all_sheets, "Data 6")
  layers     <- list()

  for (sh in sheets) {
    tryCatch({
      df <- suppressMessages(read_excel(path, sheet=sh, skip=4, .name_repair="minimal"))
      names(df) <- make.unique(as.character(names(df)))
      df <- df[rowSums(!is.na(df)) > 0, , drop=FALSE]
      if (nrow(df) == 0) next
      nc  <- ncol(df)
      hdr <- tolower(names(df))

      # Find every "stressor" column (exact match after make.unique adds suffixes)
      st_cols <- which(hdr == "stressor" | grepl("^stressor\\.\\d+$", hdr))
      if (length(st_cols) == 0) next

      for (si in st_cols) {
        if (si < 2 || si + 1 > nc) next

        stressor_v <- as.character(df[[si]])
        ct_v       <- as.character(df[[si + 1L]])

        # Detect block style from col si-1
        col_prev   <- suppressWarnings(as.numeric(as.character(df[[si - 1L]])))
        frac_num   <- sum(!is.na(col_prev)) / max(nrow(df), 1)
        beta_style <- frac_num < 0.3   # si-1 is gene names → beta-style

        if (beta_style) {
          if (si < 6) next
          gene_v <- as.character(df[[si - 1L]])
          lfc_v  <- suppressWarnings(as.numeric(df[[si - 5L]]))
          pval_v <- suppressWarnings(as.numeric(df[[si - 6L]]))
          padj_v <- suppressWarnings(as.numeric(df[[si - 2L]]))
        } else {
          # Alpha-style: all alpha blocks share gene at col 1
          if (si < 5) next
          gene_v <- as.character(df[[1L]])
          lfc_v  <- suppressWarnings(as.numeric(df[[si - 4L]]))
          pval_v <- suppressWarnings(as.numeric(df[[si - 5L]]))
          padj_v <- suppressWarnings(as.numeric(df[[si - 1L]]))
        }

        # Map stressor code → normalised cytokine name
        treat_v <- stressor_norm[stressor_v]
        treat_v[is.na(names(treat_v))] <- NA_character_

        # Drop PP/poly cells (pancreatic polypeptide — excluded from analysis)
        # and non-cytokine BFA/TG rows (treatment NA after mapping)
        keep <- !is.na(gene_v) & nzchar(gene_v) & gene_v != "NA" &
                !is.na(lfc_v)  & !is.na(stressor_v) & nzchar(stressor_v) &
                !is.na(as.character(treat_v)) &   # drops BFA, TG
                !tolower(ct_v) %in% c("pp","poly") # drops PP/poly cell types
        if (!any(keep)) next

        # Tag with patient/sheet so individual donors are traceable
        patient_label <- paste0("P", gsub("[^0-9]", "", sh))  # "Data 1" → "P1"

        layer_key <- paste0(sh, "_si", si)
        layers[[layer_key]] <- make_layer(
          gene_v[keep], lfc_v[keep], pval_v[keep], padj_v[keep],
          comparison = paste0(stressor_v[keep], " (", patient_label, ")"),
          pkg_id     = "Pck004", layer_name = "scRNA-seq", omics_type = "scRNA-seq",
          model      = meta$model,
          treatment  = as.character(treat_v[keep]),
          time_h     = meta$time_h,
          pmid       = meta$pmid, repository = meta$repository,
          is_sc      = TRUE, cell_type = ct_v[keep])
      }
    }, error = function(e) NULL)
  }
  Filter(Negate(is.null), layers)
}

load_pck005 <- function(path, meta) {
  sheets <- intersect(c("Global_Protein","Redox_Cysteine","Phospho_STY","Acetyl_Lysine"),
                      excel_sheets(path))
  tps    <- c("4h","8h","24h")

  # Map each sheet to a human-readable omics_type
  omics_map <- c(
    Global_Protein  = "Proteomics",
    Redox_Cysteine  = "Redox Proteomics",
    Phospho_STY     = "Phosphoproteomics",
    Acetyl_Lysine   = "Acetylomics"
  )

  layers <- list()
  for (sh in sheets) {
    tryCatch({
      df <- suppressMessages(read_excel(path, sheet=sh, .name_repair="minimal"))
      nm <- normalize_colnames(names(df))
      g_i <- if (!is.na(nm$gene)) nm$gene else 1L

      # gene_name = "PROTEIN_SITEID" for PTM sheets (e.g. "STAT1_T701"),
      # plain protein name for Global_Protein. This lets Gene Search do prefix
      # matching (STAT1 → finds STAT1, STAT1_T701, STAT1_S727 etc.)
      gene_v <- as.character(df[[g_i]])

      # For PTM sheets, fuse site ID into gene_name as "PROTEIN_SITE"
      site_col <- grep("site_id", tolower(names(df)), value=FALSE)
      if (length(site_col)) {
        site_v <- as.character(df[[site_col[1]]])
        gene_v <- paste0(gene_v, "_", site_v)
      }

      long <- pivot_timepoints(df, g_i, "logFC_", "adj.P.Val_", tps)
      if (is.null(long) || nrow(long) == 0) next   # skip this sheet, not the whole function

      # Re-attach gene names (pivot_timepoints already fills gene_name from df,
      # but we want our fused "PROTEIN_SITE" version here)
      long$gene_name <- rep(gene_v, times = length(tps))[seq_len(nrow(long))]

      long <- long[!is.na(long$gene_name) & nzchar(long$gene_name) & !is.na(long$log2FC), ]
      if (nrow(long) == 0) next

      otype <- omics_map[sh]
      if (is.na(otype)) otype <- "Proteomics"

      layers[[sh]] <- make_layer(long$gene_name, long$log2FC, long$pvalue, long$padj,
        comparison=long$comparison, pkg_id="Pck005", layer_name=sh, omics_type=otype,
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository)
    }, error=function(e) NULL)
  }
  Filter(Negate(is.null), layers)
}

load_pck006 <- function(path, meta) {
  layers <- list()
  d <- parse_wide(path, "Proteomics", 7, 8, 9)
  if (!is.null(d) && nrow(d)>0)
    layers$prot <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck006", layer_name="Proteomics", omics_type="Proteomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  d <- parse_sheet(path, "Lipidomics",  skip=5)
  if (!is.null(d))
    layers$lipid <- make_layer(d$gene, d$lfc, d$pval, d$padj,
      pkg_id="Pck006", layer_name="Lipidomics", omics_type="Lipidomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  d <- parse_sheet(path, "Metabolomics", skip=2)
  if (!is.null(d))
    layers$metab <- make_layer(d$gene, d$lfc, d$pval, d$padj,
      pkg_id="Pck006", layer_name="Metabolomics", omics_type="Metabolomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  Filter(Negate(is.null), layers)
}

load_pck007 <- function(path, meta) {
  sheets <- intersect(
    c("RNA-seq-IFNalpha","RNA-seq-IFNgamma","RNA-seq-IL1beta","RNA-seq-TNFalpha"),
    excel_sheets(path))
  layers <- list()
  for (sh in sheets) {
    sk <- if (sh == "RNA-seq-IFNalpha") 3L else 0L
    d  <- parse_sheet(path, sh, skip=sk)
    if (!is.null(d)) {
      cyt <- sub("RNA-seq-","",sh)
      layers[[sh]] <- make_layer(d$gene, d$lfc, d$pval, d$padj,
        comparison=paste0(cyt,"_vs_ctrl"),
        pkg_id="Pck007", layer_name=paste("RNA-seq",cyt), omics_type="RNA-seq",
        model=meta$model, treatment=cyt, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository)
    }
  }
  Filter(Negate(is.null), layers)
}

load_pck008 <- function(path, meta) {
  d <- parse_wide(path, "RNA-seq (Stat)", 3, 4, 5)
  if (is.null(d)||nrow(d)==0) return(list())
  list(rnaseq=make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
    comparison=d$comparison,
    pkg_id="Pck008", layer_name="RNA-seq", omics_type="RNA-seq",
    model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
    pmid=meta$pmid, repository=meta$repository))
}

load_pck009 <- function(path, meta) {
  layers <- list()
  d <- parse_sheet(path, "RNA-seq", skip=0)
  if (!is.null(d))
    layers$rnaseq <- make_layer(d$gene, d$lfc, d$pval, d$padj,
      pkg_id="Pck009", layer_name="RNA-seq", omics_type="RNA-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  d <- parse_wide(path, "scRNA-seq", 4, 5, 6)
  if (!is.null(d) && nrow(d)>0) {
    ct_v <- gsub("DGE_MAST_","",d$comparison)
    layers$scrna <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison,
      pkg_id="Pck009", layer_name="scRNA-seq", omics_type="scRNA-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_sc=TRUE, cell_type=ct_v)
  }
  Filter(Negate(is.null), layers)
}

load_pck010 <- function(path, meta) {
  sheets <- setdiff(excel_sheets(path), c("ReadMe","ReadMe "))
  layers <- list()
  for (sh in sheets) {
    d <- parse_sheet(path, sh, skip=0)
    if (!is.null(d))
      layers[[sh]] <- make_layer(d$gene, d$lfc, d$pval, d$padj,
        comparison=paste0("IFNa_",sh),
        pkg_id="Pck010", layer_name=sh, omics_type="RNA-seq",
        model=meta$model, treatment=meta$treatment, time_h=sh,
        pmid=meta$pmid, repository=meta$repository)
  }
  Filter(Negate(is.null), layers)
}

load_pck011 <- function(path, meta) {
  # Structure: row 9 = header; .up block cols 1-17, .down block cols 20-36 (1-based)
  # logFC1-5 at cols 8,10,12,14,16 (.up) / 27,29,31,33,35 (.down)
  # pval1-5  at cols 9,11,13,15,17 (.up) / 28,30,32,34,36 (.down)
  fisher_pval <- function(pvec) {
    # Fisher's combined probability: chi2 = -2*sum(ln(pi)), df = 2k
    pv <- suppressWarnings(as.numeric(pvec))
    pv <- pv[!is.na(pv) & pv > 0 & pv <= 1]
    if (length(pv) == 0) return(NA_real_)
    chi2 <- -2 * sum(log(pv))
    pchisq(chi2, df = 2 * length(pv), lower.tail = FALSE)
  }
  mean_lfc <- function(vec) {
    v <- suppressWarnings(as.numeric(vec))
    v <- v[!is.na(v) & is.finite(v)]
    if (length(v) == 0) return(NA_real_)
    mean(v)
  }
  tryCatch({
    df <- suppressMessages(read_excel(path, "RNA-seq", skip=8, .name_repair="minimal"))
    nc <- ncol(df)
    # Parse .up block (cols 1-17, gene at col 3)
    up_gene  <- as.character(df[[3]])
    up_lfc   <- vapply(seq_len(nrow(df)), function(i)
      mean_lfc(df[i, intersect(c(8,10,12,14,16), seq_len(nc))]), numeric(1))
    up_pval  <- vapply(seq_len(nrow(df)), function(i)
      fisher_pval(df[i, intersect(c(9,11,13,15,17), seq_len(nc))]), numeric(1))
    up_keep  <- !is.na(up_gene) & nzchar(up_gene) & up_gene != "NA" & !is.na(up_lfc)
    # Parse .down block (cols 20-36, gene at col 22)
    if (nc >= 22) {
      dn_gene  <- as.character(df[[22]])
      dn_lfc   <- vapply(seq_len(nrow(df)), function(i)
        mean_lfc(df[i, intersect(c(27,29,31,33,35), seq_len(nc))]), numeric(1))
      dn_pval  <- vapply(seq_len(nrow(df)), function(i)
        fisher_pval(df[i, intersect(c(28,30,32,34,36), seq_len(nc))]), numeric(1))
      dn_keep  <- !is.na(dn_gene) & nzchar(dn_gene) & dn_gene != "NA" & !is.na(dn_lfc)
      gene_v <- c(up_gene[up_keep], dn_gene[dn_keep])
      lfc_v  <- c(up_lfc[up_keep],  dn_lfc[dn_keep])
      pval_v <- c(up_pval[up_keep], dn_pval[dn_keep])
    } else {
      gene_v <- up_gene[up_keep]
      lfc_v  <- up_lfc[up_keep]
      pval_v <- up_pval[up_keep]
    }
    if (length(gene_v) == 0) return(list())
    list(rnaseq=make_layer(gene_v, lfc_v, pval_v, NA_real_,
      pkg_id="Pck011", layer_name="RNA-seq", omics_type="RNA-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE))
  }, error=function(e) list())
}

load_pck012 <- function(path, meta) list()   # splicing — excluded
load_pck013 <- function(path, meta) list()   # skipped per user decision

load_pck014 <- function(path, meta) {
  # MaxQuant TMT phosphoproteomics — human islets IFNγ+IL-1β 24h
  # Row 4 = column headers (patient IDs P1+..P8-).
  # Key aggregate cols (1-based): gene=41, position=37, groupA_mean=19,
  # groupB_mean=21, raw_pval=25, FDR=26, localization_prob=27
  tryCatch({
    df <- suppressMessages(read_excel(path, "Phosphoproteomics", skip=3, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    nc <- ncol(df)
    if (nc < 41) return(list())
    gene_v  <- as.character(df[[41]])
    pos_v   <- as.character(df[[37]])
    grpA    <- suppressWarnings(as.numeric(df[[19]]))
    grpB    <- suppressWarnings(as.numeric(df[[21]]))
    pval_v  <- suppressWarnings(as.numeric(df[[25]]))
    padj_v  <- suppressWarnings(as.numeric(df[[26]]))
    # FC = GroupA (cytokine) − GroupB (control) in log2 space (TMT log-normalized)
    lfc_v   <- grpA - grpB
    # Site-tagged gene name: GENE_POSITION (e.g. YLPM1_634)
    gene_site <- ifelse(!is.na(gene_v) & !is.na(pos_v) & nzchar(gene_v) & nzchar(pos_v),
                        paste0(gene_v, "_", pos_v), gene_v)
    keep <- !is.na(gene_site) & nzchar(gene_site) & gene_site != "NA" & !is.na(lfc_v)
    if (!any(keep)) return(list())
    list(phospho=make_layer(gene_site[keep], lfc_v[keep], pval_v[keep], padj_v[keep],
      comparison="cytokine_vs_ctrl",
      pkg_id="Pck014", layer_name="Phosphoproteomics", omics_type="Phosphoproteomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository))
  }, error=function(e) list())
}

load_pck015 <- function(path, meta) {
  # Top-down proteomics — human islets IFNγ+IL-1β 24h
  # Row 3 = column headers; Gene=col2, logFC=col13, P.Value=col14, adj.P.Val=col15
  tryCatch({
    df <- suppressMessages(read_excel(path, "Top-down Proteomics", skip=2, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    nc <- ncol(df)
    if (nc < 15) return(list())
    gene_v  <- as.character(df[[2]])
    lfc_v   <- suppressWarnings(as.numeric(df[[13]]))
    pval_v  <- suppressWarnings(as.numeric(df[[14]]))
    padj_v  <- suppressWarnings(as.numeric(df[[15]]))
    keep <- !is.na(gene_v) & nzchar(gene_v) & gene_v != "NA" & !is.na(lfc_v)
    if (!any(keep)) return(list())
    list(prot=make_layer(gene_v[keep], lfc_v[keep], pval_v[keep], padj_v[keep],
      comparison="cytokine_vs_ctrl",
      pkg_id="Pck015", layer_name="Top-down Proteomics", omics_type="Proteomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository))
  }, error=function(e) list())
}

load_pck016 <- function(path, meta) {
  # DE_analysis has 4 side-by-side comparisons:
  # cols 1=Gene_symbol, 8-10=Drug vs NT, 11-13=cyt vs NT (our target),
  # 14-16=cyt+drug vs NT, 17-19=cyt+drug vs cyt
  # We extract only "cyt vs. NT" (cols 1, 11, 12, 13)
  tryCatch({
    df <- suppressMessages(read_excel(path, "DE_analysis", skip=1, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    if (ncol(df) < 13) return(list())
    gene_v  <- as.character(df[[1]])
    lfc_v   <- suppressWarnings(as.numeric(df[[11]]))
    pval_v  <- suppressWarnings(as.numeric(df[[12]]))
    padj_v  <- suppressWarnings(as.numeric(df[[13]]))
    keep    <- !is.na(gene_v) & nzchar(gene_v) & gene_v != "NA" & !is.na(lfc_v)
    if (!any(keep)) return(list())
    list(rnaseq=make_layer(gene_v[keep], lfc_v[keep], pval_v[keep], padj_v[keep],
      comparison="cyt_vs_NT",
      pkg_id="Pck016", layer_name="RNA-seq", omics_type="RNA-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository))
  }, error=function(e) list())
}

load_pck017 <- function(path, meta) list()   # splicing — excluded

load_pck018 <- function(path, meta) {
  tryCatch({
    df <- suppressMessages(read_excel(path, "Proteomics", skip=5, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    gene_v <- as.character(df[[4]])
    mat_ctrl <- suppressWarnings(sapply(df[, 11:13], as.numeric))
    mat_ct   <- suppressWarnings(sapply(df[, 17:19], as.numeric))
    pval_v   <- suppressWarnings(as.numeric(df[[23]]))
    lfc_v    <- avg_reps(mat_ct) - avg_reps(mat_ctrl)
    keep     <- !is.na(gene_v) & !is.na(lfc_v)
    list(prot=make_layer(gene_v[keep], lfc_v[keep], pval_v[keep], NA_real_,
      comparison="CT-NT vs NoCT-NT",
      pkg_id="Pck018", layer_name="Proteomics", omics_type="Proteomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE))
  }, error=function(e) list())
}

load_pck019 <- function(path, meta) {
  # Row 5 = col headers; after skip=4: gene=col3,
  # CT1-NonTarget siRNA replicates = cols 12-15 (log2FC vs NoCT_NT baseline)
  # CT1-PARP12 siRNA replicates    = cols 16-19
  # T-test NoCT_NT vs CT_NT        = col 21
  # T-test CT_NT vs CT_si (PARP12 KD in CT context) = col 22
  tryCatch({
    df <- suppressMessages(read_excel(path, "Proteomics", skip=4, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    gene_v      <- as.character(df[[3]])
    ct_nt_mat   <- suppressWarnings(sapply(df[, 12:15], as.numeric))
    ct_parp_mat <- suppressWarnings(sapply(df[, 16:19], as.numeric))
    pval_cyt    <- suppressWarnings(as.numeric(df[[21]]))   # NoCT_NT vs CT_NT
    pval_kd     <- suppressWarnings(as.numeric(df[[22]]))   # CT_NT vs CT_si
    lfc_cyt     <- avg_reps(ct_nt_mat)                      # cytokine effect vs NoCT_NT
    lfc_kd      <- avg_reps(ct_parp_mat) - avg_reps(ct_nt_mat)  # net KD effect in CT
    keep_c      <- !is.na(gene_v) & !is.na(lfc_cyt)
    keep_k      <- !is.na(gene_v) & !is.na(lfc_kd)
    layers <- list()
    if (any(keep_c))
      layers$cyt <- make_layer(gene_v[keep_c], lfc_cyt[keep_c], pval_cyt[keep_c], NA_real_,
        comparison="Cytokine_vs_Control",
        pkg_id="Pck019", layer_name="Proteomics (CT-NT)", omics_type="Proteomics",
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE)
    if (any(keep_k))
      layers$kd <- make_layer(gene_v[keep_k], lfc_kd[keep_k], pval_kd[keep_k], NA_real_,
        comparison="Parp12_KD_vs_ctrl_in_CT",
        pkg_id="Pck019", layer_name="Proteomics (PARP12 KD)", omics_type="Proteomics",
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE)
    Filter(Negate(is.null), layers)
  }, error=function(e) list())
}

load_pck020 <- function(path, meta) {
  # Cols 8-12 are log2(Cytokine_treated / WT) per replicate (already log2FC vs WT).
  # log2FC = mean of 5 replicates; p-value = one-sample t-test against 0.
  tryCatch({
    df <- suppressMessages(read_excel(path, "Proteomics", skip=7, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    gene_v  <- as.character(df[[4]])
    mat_ct  <- suppressWarnings(sapply(df[, 8:12], as.numeric))
    if (is.vector(mat_ct)) mat_ct <- matrix(mat_ct, ncol=1)
    lfc_v <- rowMeans(mat_ct, na.rm=TRUE)
    # One-sample t-test: H0: mean log2FC == 0 for each protein
    pval_v <- vapply(seq_len(nrow(mat_ct)), function(i) {
      x <- mat_ct[i, !is.na(mat_ct[i, ])]
      if (length(x) < 2) return(NA_real_)
      tryCatch(t.test(x, mu=0)$p.value, error=function(e) NA_real_)
    }, numeric(1))
    keep <- !is.na(gene_v) & nzchar(gene_v) & gene_v != "NA" & !is.na(lfc_v)
    list(prot=make_layer(gene_v[keep], lfc_v[keep], pval_v[keep], NA_real_,
      comparison="Cytokine_vs_WildType",
      pkg_id="Pck020", layer_name="Proteomics", omics_type="Proteomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE))
  }, error=function(e) list())
}

load_pck021 <- function(path, meta) {
  layers <- list()
  # --- Proteomics ---
  tryCatch({
    df <- suppressMessages(read_excel(path, "Proteomics", skip=6, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    gene_v   <- as.character(df[[3]])
    pval_v   <- suppressWarnings(as.numeric(df[[5]]))
    mat_ctrl <- suppressWarnings(sapply(df[, 6:8],  as.numeric))
    mat_ct   <- suppressWarnings(sapply(df[, 9:11], as.numeric))
    # Ctrl = ethanol vehicle; CT = cytokine + ethanol. Values are already log2FC vs ref.
    lfc_v    <- avg_reps(mat_ct) - avg_reps(mat_ctrl)
    keep     <- !is.na(gene_v) & nzchar(gene_v) & gene_v != "NA" & !is.na(lfc_v)
    if (any(keep))
      layers$prot <- make_layer(gene_v[keep], lfc_v[keep], pval_v[keep], NA_real_,
        comparison="Cytokine_vs_Control",
        pkg_id="Pck021", layer_name="Proteomics", omics_type="Proteomics",
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE)
  }, error=function(e) NULL)
  # --- Metabolomics: CT_Eth (cols 2-5) vs NCT_Eth (cols 10-13), peak area intensities ---
  tryCatch({
    df_m <- suppressMessages(read_excel(path, "Metabolomics", skip=3, .name_repair="minimal"))
    names(df_m) <- make.unique(as.character(names(df_m)))
    metab_v  <- as.character(df_m[[1]])
    ct_mat   <- suppressWarnings(sapply(df_m[, 2:5],  as.numeric))
    nct_mat  <- suppressWarnings(sapply(df_m[, 10:13], as.numeric))
    if (is.vector(ct_mat))  ct_mat  <- matrix(ct_mat,  ncol=1)
    if (is.vector(nct_mat)) nct_mat <- matrix(nct_mat, ncol=1)
    # log2 transform (add 1 to handle zeros)
    lct  <- log2(ct_mat  + 1)
    lnct <- log2(nct_mat + 1)
    lfc_v <- rowMeans(lct, na.rm=TRUE) - rowMeans(lnct, na.rm=TRUE)
    # Welch t-test per metabolite
    pval_v <- vapply(seq_len(nrow(ct_mat)), function(i) {
      a <- lct[i, !is.na(lct[i, ])]
      b <- lnct[i, !is.na(lnct[i, ])]
      if (length(a) < 2 || length(b) < 2) return(NA_real_)
      tryCatch(t.test(a, b)$p.value, error=function(e) NA_real_)
    }, numeric(1))
    keep_m <- !is.na(metab_v) & nzchar(metab_v) & metab_v != "NA" & !is.na(lfc_v)
    if (any(keep_m))
      layers$metab <- make_layer(metab_v[keep_m], lfc_v[keep_m], pval_v[keep_m], NA_real_,
        comparison="CT_Eth_vs_NCT_Eth",
        pkg_id="Pck021", layer_name="Metabolomics", omics_type="Metabolomics",
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE)
  }, error=function(e) NULL)
  Filter(Negate(is.null), layers)
}

load_pck022 <- function(path, meta) {
  sc_sheets <- setdiff(excel_sheets(path), c("ReadMe"))
  layers <- list()
  for (sh in sc_sheets) {
    d <- parse_wide(path, sh, 2, 3, 4, ln_to_log2=TRUE)
    if (is.null(d)||nrow(d)==0) next
    ct      <- gsub("_markers","",sh,ignore.case=TRUE)
    treat_v <- extract_sc_treatment_v(d$comparison)   # e.g. "IL-1β (S2) vs. Control" → "IL-1β"
    layers[[sh]] <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck022", layer_name="scRNA-seq", omics_type="scRNA-seq",
      model=meta$model, treatment=as.character(treat_v), time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_sc=TRUE, cell_type=ct)
  }
  Filter(Negate(is.null), layers)
}

load_pck023 <- function(path, meta) {
  sc_sheets <- setdiff(excel_sheets(path), c("ReadMe"))
  layers <- list()
  for (sh in sc_sheets) {
    d <- parse_wide(path, sh, 2, 3, 4, ln_to_log2=FALSE)
    if (is.null(d)||nrow(d)==0) next
    ct      <- gsub("_markers","",sh,ignore.case=TRUE)
    # Keep only pure cytokine-vs-Control comparisons (drop NMMA rows)
    keep_comp <- grepl("^IFN|^IL|^TNF", extract_sc_treatment_v(d$comparison)) &
                 grepl("[Vv]s\\.?.*[Cc]ontrol", d$comparison, perl=TRUE)
    d <- d[keep_comp, , drop=FALSE]
    if (nrow(d) == 0) next
    treat_v <- extract_sc_treatment_v(d$comparison)
    layers[[sh]] <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck023", layer_name="scRNA-seq", omics_type="scRNA-seq",
      model=meta$model, treatment=as.character(treat_v), time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_sc=TRUE, cell_type=ct)
  }
  Filter(Negate(is.null), layers)
}

load_pck024 <- function(path, meta) {
  # Each sheet has 8 comparison blocks. We keep S2, S3, S4 (all 6h) and S6 (18h IFNγ+IL-1β).
  # The comparison label has the subject sample in parens BEFORE "vs",
  # so matching \\(S[2346]\\).*vs.*[Cc]ontrol avoids NMMA comparisons where S6 is the reference.
  sc_sheets <- setdiff(excel_sheets(path), c("ReadMe"))
  layers <- list()
  for (sh in sc_sheets) {
    d <- parse_wide(path, sh, 2, 3, 4, ln_to_log2=FALSE)
    if (is.null(d) || nrow(d) == 0) next
    # Keep only pure cytokine-vs-Control (S2/S3/S4/S6); drops NMMA and non-ctrl comparisons
    d <- d[grepl("\\(S[2346]\\).*[Vv]s\\.?.*[Cc]ontrol", d$comparison, perl=TRUE) &
           !grepl("NMMA", d$comparison, ignore.case=TRUE), ]
    if (nrow(d) == 0) next
    ct      <- gsub("_markers","", sh, ignore.case=TRUE)
    treat_v <- extract_sc_treatment_v(d$comparison)
    layers[[sh]] <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck024", layer_name="scRNA-seq", omics_type="scRNA-seq",
      model=meta$model, treatment=as.character(treat_v), time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_sc=TRUE, cell_type=ct)
  }
  Filter(Negate(is.null), layers)
}
load_pck025 <- function(path, meta) {
  layers <- list()
  for (sh in c("Discovery RNA-seq","Validation RNA-seq")) {
    if (!sh %in% excel_sheets(path)) next
    d <- parse_wide(path, sh, 2, 3, 4)
    if (!is.null(d)&&nrow(d)>0)
      layers[[sh]] <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
        comparison=d$comparison, pkg_id="Pck025", layer_name=sh, omics_type="RNA-seq",
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository)
  }
  if ("Proteomics" %in% excel_sheets(path)) {
    d <- parse_wide(path, "Proteomics", 2, 3, 4)
    if (!is.null(d)&&nrow(d)>0)
      layers$prot <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
        comparison=d$comparison, pkg_id="Pck025", layer_name="Proteomics", omics_type="Proteomics",
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository)
  }
  Filter(Negate(is.null), layers)
}

# ==========================================================================
# MASTER LOADER
# ==========================================================================

load_all_data <- function() {
  data_dir <- resolve_data_dir()
  cat("Data directory:", data_dir, "\n")

  # skip=1: row 1 is the merged title cell; row 2 has the real column headers
  meta_file <- resolve_metadata_file(data_dir)
  cat("Metadata file:", basename(meta_file), "\n")
  meta_raw <- suppressMessages(
    read_excel(meta_file, skip=1, .name_repair="minimal")
  )
  meta_raw <- as.data.frame(meta_raw)
  names(meta_raw) <- tolower(gsub("[[:space:]]+","_",trimws(names(meta_raw))))
  # Drop fully empty rows (new flat file has no blanks in col1, but guard anyway)
  meta_raw <- meta_raw[!is.na(meta_raw[[1]]) & nzchar(trimws(meta_raw[[1]])), , drop=FALSE]

  # ---- Metadata helpers ---------------------------------------------------
  # get_meta(): returns a scalar list from the FIRST row for a package.
  # All existing loaders use this and are backward-compatible.
  get_meta <- function(id) {
    rows <- meta_raw[grepl(paste0("^", id, "$"), trimws(meta_raw[[1]])), , drop=FALSE]
    if (nrow(rows) == 0) rows <- meta_raw[1, , drop=FALSE]
    row  <- rows[1, , drop=FALSE]
    cols <- names(row)
    safe <- function(nm) {
      idx <- grep(nm, cols, ignore.case=TRUE)
      if (length(idx)) as.character(row[[idx[1]]]) else NA_character_
    }
    list(
      model      = safe("model|cell.?type|islet"),
      treatment  = safe("treatment|cytokine|condition"),
      time_h     = safe("time|duration|hour"),
      pmid       = safe("pmid|pubmed"),
      repository = safe("repository|geo|accession|gse"),
      rows       = rows   # all rows for this package — for loaders that need per-omics lookup
    )
  }

  # meta_for_omics(): pick the metadata row whose omics type best matches `omics`.
  # Falls back to first row if no match. Used by multi-omics / multi-treatment loaders.
  meta_for_omics <- function(meta, omics) {
    rows <- meta$rows
    if (is.null(rows) || nrow(rows) == 0) return(meta)
    omics_col <- grep("omics|data.?type", names(rows), ignore.case=TRUE, value=FALSE)
    if (length(omics_col)) {
      idx <- grep(omics, rows[[omics_col[1]]], ignore.case=TRUE)
      if (length(idx)) rows <- rows[idx[1], , drop=FALSE]
      else             rows <- rows[1,       , drop=FALSE]
    } else {
      rows <- rows[1, , drop=FALSE]
    }
    cols <- names(rows)
    safe <- function(nm) {
      i <- grep(nm, cols, ignore.case=TRUE)
      if (length(i)) as.character(rows[[i[1]]]) else NA_character_
    }
    list(
      model      = safe("model|cell.?type|islet"),
      treatment  = safe("treatment|cytokine|condition"),
      time_h     = safe("time|duration|hour"),
      pmid       = safe("pmid|pubmed"),
      repository = safe("repository|geo|accession|gse"),
      rows       = meta$rows
    )
  }

  loaders <- list(
    Pck001=load_pck001, Pck002=load_pck002, Pck003=load_pck003,
    Pck004=load_pck004, Pck005=load_pck005, Pck006=load_pck006,
    Pck007=load_pck007, Pck008=load_pck008, Pck009=load_pck009,
    Pck010=load_pck010, Pck011=load_pck011, Pck012=load_pck012,
    Pck013=load_pck013, Pck014=load_pck014, Pck015=load_pck015,
    Pck016=load_pck016, Pck017=load_pck017, Pck018=load_pck018,
    Pck019=load_pck019, Pck020=load_pck020, Pck021=load_pck021,
    Pck022=load_pck022, Pck023=load_pck023, Pck024=load_pck024,
    Pck025=load_pck025
  )

  all_packages <- list()
  for (id in names(loaders)) {
    cat(sprintf("  Loading %s ...", id))
    path   <- pkg_path(data_dir, id)
    meta   <- get_meta(id)
    layers <- tryCatch(
      loaders[[id]](path, meta),
      error = function(e) {
        message(sprintf("\n  WARNING: %s failed — %s", id, conditionMessage(e)))
        list()
      }
    )
    if (length(layers) > 0) {
      all_packages[[id]] <- list(pkg_id=id, meta=meta, layers=layers)
      cat(sprintf(" %d layer(s)\n", length(layers)))
    } else {
      cat(" skipped\n")
    }
  }

  flat_list <- lapply(all_packages, function(pkg) do.call(rbind, pkg$layers))
  flat_all  <- do.call(rbind, Filter(Negate(is.null), flat_list))

  list(
    packages  = all_packages,
    metadata  = meta_raw,
    flat      = flat_all,
    flat_bulk = flat_all[!flat_all$is_sc, ],
    flat_sc   = flat_all[flat_all$is_sc,  ]
  )
}
