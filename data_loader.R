# ============================================================
# data_loader.R  v3
# ============================================================
# Decisions recorded:
#   Pck012, 017  → EXCLUDED (splicing ΔPSI data)
#   Pck022       → ln(FC) silently converted to log2FC
#   Pck013       → RNA-seq block only (miRNA block skipped)
#   Pck011,018,019,020,021 → replicate-averaged + caveat flag
#   Pck014,015   → skipped (no gene IDs / no stats)
# ============================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
})

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !is.na(a[[1]]) && nzchar(a[[1]])) a else b
}

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
        file.exists(file.path(p, "Table 1-Metadata.xlsx")))
      return(normalizePath(p))
  stop("Cannot find data directory. Set MULTIOMICS_DATA_DIR env var.")
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
      fc  <- block_starts[1]:block_ends[1]
      fn  <- normalize_colnames(hdr[fc])
      gene_v <- if (!is.na(fn$gene)) as.character(unlist(data_rows[, fc[fn$gene]])) else NA_character_
    } else {
      gene_v <- NA_character_
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
  layers <- list()
  for (sh in c("RNA-seq HI", "RNA-seq Endo")) {
    if (!sh %in% excel_sheets(path)) next
    d <- parse_wide(path, sh, 4, 5, 6)
    if (!is.null(d) && nrow(d) > 0)
      layers[[sh]] <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
        comparison=d$comparison, pkg_id="Pck003", layer_name=sh, omics_type="RNA-seq",
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository)
  }
  d <- parse_sheet(path, "ATAC-seq")
  if (!is.null(d))
    layers$atac <- make_layer(d$gene, d$lfc, d$pval, d$padj,
      pkg_id="Pck003", layer_name="ATAC-seq", omics_type="ATAC-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_genomic=TRUE)
  d <- parse_wide(path, "Proteomics",  4, 5, 6)
  if (!is.null(d) && nrow(d) > 0)
    layers$prot <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck003", layer_name="Proteomics", omics_type="Proteomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  d <- parse_wide(path, "Lipidomics",  3, 4, 5)
  if (!is.null(d) && nrow(d) > 0)
    layers$lipid <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck003", layer_name="Lipidomics", omics_type="Lipidomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  d <- parse_wide(path, "Metabolomics", 3, 4, 5)
  if (!is.null(d) && nrow(d) > 0)
    layers$metab <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck003", layer_name="Metabolomics", omics_type="Metabolomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository)
  Filter(Negate(is.null), layers)
}

load_pck004 <- function(path, meta) {
  sheets <- grep("^Data", excel_sheets(path), value = TRUE)
  layers <- list()
  for (sh in sheets) {
    tryCatch({
      df <- suppressMessages(read_excel(path, sheet=sh, skip=4, .name_repair="minimal"))
      names(df) <- make.unique(as.character(names(df)))
      df <- df[rowSums(!is.na(df)) > 0, , drop=FALSE]
      if (nrow(df) == 0) return(NULL)
      nm   <- normalize_colnames(names(df))
      g_i  <- if (!is.na(nm$gene)) nm$gene else 1L
      l_i  <- nm$lfc; p_i <- nm$pval; pa_i <- nm$padj
      if (is.na(l_i)) return(NULL)
      ct_col <- grep("celltype",  tolower(names(df)), value=FALSE)[1]
      st_col <- grep("stressor",  tolower(names(df)), value=FALSE)[1]
      gene_v <- as.character(df[[g_i]])
      lfc_v  <- suppressWarnings(as.numeric(df[[l_i]]))
      pval_v <- if (!is.na(p_i))  suppressWarnings(as.numeric(df[[p_i]]))  else NA_real_
      padj_v <- if (!is.na(pa_i)) suppressWarnings(as.numeric(df[[pa_i]])) else NA_real_
      ct_v   <- if (!is.na(ct_col)) as.character(df[[ct_col]]) else NA_character_
      comp_v <- if (!is.na(st_col)) as.character(df[[st_col]]) else "cytokine_vs_ctrl"
      layers[[sh]] <- make_layer(gene_v, lfc_v, pval_v, padj_v,
        comparison=comp_v, pkg_id="Pck004", layer_name="scRNA-seq", omics_type="scRNA-seq",
        model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
        pmid=meta$pmid, repository=meta$repository, is_sc=TRUE, cell_type=ct_v)
    }, error=function(e) NULL)
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
    ct_v <- sub("\\s.*","",ct_v)
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
  tryCatch({
    df <- suppressMessages(read_excel(path, "RNA-seq", skip=8, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    gene_v <- as.character(df[[3]])
    lfc_v  <- suppressWarnings(as.numeric(df[[6]]))
    pval_v <- suppressWarnings(as.numeric(df[[9]]))
    keep   <- !is.na(gene_v) & !grepl("[._]", gene_v) & !is.na(lfc_v) & is.finite(lfc_v)
    if (!any(keep)) return(list())
    list(rnaseq=make_layer(gene_v[keep], lfc_v[keep], pval_v[keep], NA_real_,
      pkg_id="Pck011", layer_name="RNA-seq", omics_type="RNA-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE))
  }, error=function(e) list())
}

load_pck012 <- function(path, meta) list()   # splicing — excluded

load_pck013 <- function(path, meta) {
  d <- parse_wide(path, "RNAseq", 3, 4, 5, first_block_only=TRUE)
  if (is.null(d)||nrow(d)==0) return(list())
  # "Fold Change" is raw FC → log2 transform
  d$log2FC <- suppressWarnings(log2(abs(d$log2FC)) * sign(d$log2FC))
  d <- d[is.finite(d$log2FC), ]
  if (nrow(d)==0) return(list())
  list(rnaseq=make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
    pkg_id="Pck013", layer_name="RNA-seq", omics_type="RNA-seq",
    model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
    pmid=meta$pmid, repository=meta$repository))
}

load_pck014 <- function(path, meta) list()   # phospho matrix, no gene IDs
load_pck015 <- function(path, meta) list()   # top-down catalog, no stats

load_pck016 <- function(path, meta) {
  d <- parse_sheet(path, "DE_analysis", skip=1)
  if (is.null(d)) return(list())
  list(rnaseq=make_layer(d$gene, d$lfc, d$pval, d$padj,
    pkg_id="Pck016", layer_name="RNA-seq (miRNA study)", omics_type="RNA-seq",
    model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
    pmid=meta$pmid, repository=meta$repository))
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
  tryCatch({
    raw <- suppressMessages(
      read_excel(path, "Proteomics", col_names=FALSE, .name_repair="minimal")
    )
    cond_labels <- as.character(unlist(raw[4, ]))
    cond_labels[is.na(cond_labels)|cond_labels=="NA"|trimws(cond_labels)==""] <- NA_character_
    df   <- suppressMessages(read_excel(path, "Proteomics", skip=4, .name_repair="minimal"))
    gene_v     <- as.character(df[[3]])
    block_cols <- which(!is.na(cond_labels))
    block_ends <- c(block_cols[-1]-1L, ncol(raw))
    # Take first 4 condition blocks (cols 4-19)
    main <- block_cols[block_cols >= 4 & block_cols <= 16]
    results <- lapply(main, function(s) {
      e   <- min(block_ends[which(block_cols==s)], s+3)
      mat <- suppressWarnings(sapply(df[, s:e], as.numeric))
      if (is.vector(mat)) mat <- matrix(mat, ncol=1)
      data.frame(gene_name=gene_v, log2FC=avg_reps(mat),
                 pvalue=NA_real_, padj=NA_real_,
                 comparison=trimws(cond_labels[s]), stringsAsFactors=FALSE)
    })
    out <- do.call(rbind, results)
    out <- out[!is.na(out$gene_name)&!is.na(out$log2FC), ]
    if (nrow(out)==0) return(list())
    list(prot=make_layer(out$gene_name, out$log2FC, out$pvalue, out$padj,
      comparison=out$comparison,
      pkg_id="Pck019", layer_name="Proteomics", omics_type="Proteomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, has_stats=FALSE, replicate_averaged=TRUE))
  }, error=function(e) list())
}

load_pck020 <- function(path, meta) {
  tryCatch({
    df <- suppressMessages(read_excel(path, "Proteomics", skip=7, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    gene_v  <- as.character(df[[4]])
    pval_v  <- suppressWarnings(as.numeric(df[[7]]))
    mat_wt  <- suppressWarnings(sapply(df[, 8:12],  as.numeric))
    mat_ct  <- suppressWarnings(sapply(df[, 13:17], as.numeric))
    lfc_v   <- avg_reps(mat_ct) - avg_reps(mat_wt)
    keep    <- !is.na(gene_v) & !is.na(lfc_v)
    list(prot=make_layer(gene_v[keep], lfc_v[keep], pval_v[keep], NA_real_,
      comparison="Cytokine vs WildType",
      pkg_id="Pck020", layer_name="Proteomics", omics_type="Proteomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE))
  }, error=function(e) list())
}

load_pck021 <- function(path, meta) {
  layers <- list()
  tryCatch({
    df <- suppressMessages(read_excel(path, "Proteomics", skip=6, .name_repair="minimal"))
    names(df) <- make.unique(as.character(names(df)))
    gene_v   <- as.character(df[[3]])
    pval_v   <- suppressWarnings(as.numeric(df[[5]]))
    mat_ctrl <- suppressWarnings(sapply(df[, 6:8],  as.numeric))
    mat_ct   <- suppressWarnings(sapply(df[, 9:11], as.numeric))
    lfc_v    <- avg_reps(mat_ct) - avg_reps(mat_ctrl)
    keep     <- !is.na(gene_v) & !is.na(lfc_v)
    layers$prot <- make_layer(gene_v[keep], lfc_v[keep], pval_v[keep], NA_real_,
      comparison="Cytokine vs Control",
      pkg_id="Pck021", layer_name="Proteomics", omics_type="Proteomics",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, replicate_averaged=TRUE)
  }, error=function(e) NULL)
  # Pck021 Metabolomics: sheet contains raw intensity columns only (no fold change
  # or statistical comparisons). Excluded — cannot derive meaningful log2FC or p-values.
  Filter(Negate(is.null), layers)
}

load_pck022 <- function(path, meta) {
  sc_sheets <- setdiff(excel_sheets(path), c("ReadMe"))
  layers <- list()
  for (sh in sc_sheets) {
    d <- parse_wide(path, sh, 2, 3, 4, ln_to_log2=TRUE)
    if (is.null(d)||nrow(d)==0) next
    ct <- gsub("_markers","",sh,ignore.case=TRUE)
    layers[[sh]] <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck022", layer_name="scRNA-seq", omics_type="scRNA-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
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
    ct <- gsub("_markers","",sh,ignore.case=TRUE)
    layers[[sh]] <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck023", layer_name="scRNA-seq", omics_type="scRNA-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
      pmid=meta$pmid, repository=meta$repository, is_sc=TRUE, cell_type=ct)
  }
  Filter(Negate(is.null), layers)
}

load_pck024 <- function(path, meta) {
  # Each sheet has two comparison blocks side by side (e.g. IL-1β vs ctrl | IFNγ vs ctrl).
  # Use parse_wide(2,3,4) so both comparisons are captured per cell-type sheet.
  sc_sheets <- setdiff(excel_sheets(path), c("ReadMe"))
  layers <- list()
  for (sh in sc_sheets) {
    d <- parse_wide(path, sh, 2, 3, 4, ln_to_log2=FALSE)
    if (is.null(d)||nrow(d)==0) next
    ct <- gsub("_markers","",sh,ignore.case=TRUE)
    layers[[sh]] <- make_layer(d$gene_name, d$log2FC, d$pvalue, d$padj,
      comparison=d$comparison, pkg_id="Pck024", layer_name="scRNA-seq", omics_type="scRNA-seq",
      model=meta$model, treatment=meta$treatment, time_h=meta$time_h,
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

  meta_raw <- suppressMessages(
    read_excel(file.path(data_dir, "Table 1-Metadata.xlsx"), .name_repair="minimal")
  )
  meta_raw <- as.data.frame(meta_raw)
  names(meta_raw) <- tolower(gsub("[[:space:]]+","_",trimws(names(meta_raw))))

  get_meta <- function(id) {
    row <- meta_raw[grepl(id, meta_raw[[1]], fixed=TRUE), , drop=FALSE]
    if (nrow(row)==0) row <- meta_raw[1,,drop=FALSE]
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
      repository = safe("repository|geo|accession|gse")
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
