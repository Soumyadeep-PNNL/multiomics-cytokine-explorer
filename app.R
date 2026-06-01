# ============================================================
# Multiomics In Vitro Cytokine Explorer — v6
# Tabs: Overview | Bulk Gene Search | Single Cell
# ============================================================
# GIT REMINDER — run after every session:
#   cd "/Users/sark224/Library/CloudStorage/OneDrive-PNNL/Documents/Projects/In-vitroCT treated exp-Shinnyapp"
#   git add App/app.R App/data_loader.R
#   git commit -m "describe changes here"
#   git push
# ============================================================

local({
  sd <- tryCatch(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)),
                 error=function(e)"")
  if (!nzchar(sd))
    sd <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile)), error=function(e)"")
  if (nzchar(sd) && file.exists(file.path(sd,"data_loader.R"))) setwd(sd)
})

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(dplyr); library(tidyr)
  library(ggplot2); library(plotly); library(DT)
  library(ggrepel); library(scales)
})

source("data_loader.R")

# ── Load / build cache ───────────────────────────────────────
CACHE_FILE <- "data_cache.rds"
cache_ok   <- file.exists(CACHE_FILE) && file.info(CACHE_FILE)$size > 0
if (cache_ok) {
  DATA <- tryCatch({
    cat("Loading from cache...\n")
    readRDS(CACHE_FILE)
  }, error = function(e) {
    message("Cache corrupt — rebuilding: ", conditionMessage(e))
    NULL
  })
  if (is.null(DATA)) cache_ok <- FALSE
}
if (!cache_ok) {
  cat("Parsing all packages (first run — may take 1–2 min)...\n")
  DATA <- load_all_data()
  tryCatch(saveRDS(DATA, CACHE_FILE, compress="xz"), error=function(e) NULL)
  cat("Cache saved.\n")
}

FLAT      <- DATA$flat
FLAT_BULK <- DATA$flat_bulk
FLAT_SC   <- DATA$flat_sc
META_RAW  <- DATA$metadata  # raw, multi-row (for display in Overview)

# ── Model abbreviations (defined BEFORE the normalisation loop below) ────────
abbrev_model <- function(m) {
  if (is.na(m) || !nzchar(m)) return("")
  lu <- c(
    "Human pancreatic islets$"            = "HPI",
    "EndoC"                               = "EndoC",
    "β-TC-6|bTC.?6|beta.TC.6"       = "βTC-6",
    "β-TC-3|bTC.?3|beta.TC.3"       = "βTC-3",
    "MIN6"                                = "MIN6",
    "Mouse pancreatic islets"             = "MPI",
    "iPSC|HEL115"                         = "iPSC",
    "Adprhl2|Parp12|Pla2g6|palmitate"    = "MIN6"
  )
  for (pat in names(lu))
    if (grepl(pat, m, ignore.case = TRUE, perl = TRUE)) return(lu[[pat]])
  substr(m, 1, 8)
}
abbrev_model_v <- Vectorize(abbrev_model)

# ── Treatment normaliser — strips non-cytokine additives; returns NA for pure
#    non-cytokine treatments (BFA, TG, NMMA, DFMO, E2, etc.) ─────────────────
normalize_treatment <- function(t) {
  if (is.na(t) || !nzchar(trimws(t))) return(NA_character_)
  pats <- c(
    "IFN.?[αa]|IFN.?alpha|INF.?[αa]" = "IFNα",
    "IFN.?[γg]|IFN.?gamma"                  = "IFNγ",
    "IL.?1.?[βb]|IL.?1\\s*beta"             = "IL-1β",
    "TNF.?[αa]|TNF.?alpha"                   = "TNFα"
  )
  found <- character(0)
  for (p in names(pats))
    if (grepl(p, t, ignore.case=TRUE, perl=TRUE))
      found <- c(found, pats[[p]])
  if (length(found)==0) return(NA_character_)
  ord <- c("IFNα","IFNγ","IL-1β","TNFα")
  paste(ord[ord %in% found], collapse=" + ")
}

# ── Per-row time point extractor ──────────────────────────────────────────────
# Priority: (1) number+h pattern in comparison, (2) single number in time_h meta
derive_time_h <- function(time_h_meta, comparison) {
  extract_h <- function(s) {
    m <- regexpr("\\b(\\d+)\\s*[hH]\\b", s, perl=TRUE)
    if (m[1] > 0) {
      num <- regmatches(s, regexpr("\\d+", substring(s, m[1])))
      if (length(num) > 0 && nzchar(num)) return(paste0(num, "h"))
    }
    NA_character_
  }
  if (!is.na(comparison) && nzchar(comparison) &&
      !comparison %in% c("cytokine_vs_ctrl","")) {
    r <- extract_h(comparison)
    if (!is.na(r)) return(r)
  }
  if (!is.na(time_h_meta) && nzchar(trimws(time_h_meta))) {
    nums <- regmatches(time_h_meta, gregexpr("\\d+", time_h_meta))[[1]]
    if (length(nums) == 1) return(paste0(nums, "h"))
  }
  NA_character_
}

# ── X-axis label builder ──────────────────────────────────────────────────────
make_ds_label <- function(treat_c, model_a, time_c, fmt="treat_time_model") {
  t <- ifelse(is.na(treat_c)|!nzchar(treat_c), "?",  treat_c)
  m <- ifelse(is.na(model_a)|!nzchar(model_a), "?",  model_a)
  d <- ifelse(is.na(time_c) |!nzchar(time_c),  "?h", time_c)
  switch(fmt,
    treat_time_model = paste0(t,"\n",d,"\n",m),
    treat_model_time = paste0(t,"\n",m,"\n",d),
    model_treat_time = paste0(m,"\n",t,"\n",d),
    treat_time       = paste0(t," \xb7 ",d),
    paste0(t,"\n",d,"\n",m)
  )
}


# ── Normalise model + treatment strings in the actual data ───────
# Filters are driven from FLAT_BULK so skipped packages never appear.
for (flt in c("FLAT", "FLAT_BULK", "FLAT_SC")) {
  d <- get(flt)
  if (!is.null(d) && nrow(d) > 0) {
    d$model_abbrev    <- abbrev_model_v(d$model)
    d$treatment_clean <- vapply(d$treatment, normalize_treatment, character(1))
    d$time_h_clean    <- mapply(derive_time_h, d$time_h, d$comparison,
                                SIMPLIFY=TRUE, USE.NAMES=FALSE)
    assign(flt, d)
  }
}
# Keep only cytokine-treated rows in FLAT_BULK (excludes BFA, TG, NMMA, E2, etc.)
if (!is.null(FLAT_BULK) && nrow(FLAT_BULK) > 0)
  FLAT_BULK <- FLAT_BULK[!is.na(FLAT_BULK$treatment_clean), , drop=FALSE]

# ── Omics type groupings ─────────────────────────────────────
PTM_TYPES  <- c("Phosphoproteomics","Redox Proteomics","Acetylomics")
BULK_TYPES <- c("RNA-seq","Proteomics","Lipidomics","Metabolomics",
                "ATAC-seq","ChIP-seq","UMI-4C")

OMICS_COLORS <- c(
  "RNA-seq"           = "#2196F3",
  "Proteomics"        = "#F44336",
  "Phosphoproteomics" = "#FF5722",
  "Redox Proteomics"  = "#FF9800",
  "Acetylomics"       = "#FFC107",
  "Lipidomics"        = "#00BCD4",
  "Metabolomics"      = "#795548",
  "ATAC-seq"          = "#4CAF50",
  "ChIP-seq"          = "#8BC34A",
  "UMI-4C"            = "#607D8B",
  "scRNA-seq"         = "#9C27B0"
)

# Mild pastel backgrounds per package (Material Design 100-level)
PKG_BG <- c(
  Pck001="#BBDEFB", Pck002="#C8E6C9", Pck003="#FFE0B2", Pck004="#F8BBD0",
  Pck005="#E1BEE7", Pck006="#B3E5FC", Pck007="#DCEDC8", Pck008="#FFF9C4",
  Pck009="#F0F4C3", Pck010="#B2EBF2", Pck011="#FFCCBC", Pck012="#D7CCC8",
  Pck013="#CFD8DC", Pck014="#BBDEFB", Pck015="#C8E6C9", Pck016="#FFE0B2",
  Pck017="#F8BBD0", Pck018="#E1BEE7", Pck019="#B3E5FC", Pck020="#DCEDC8",
  Pck021="#FFF9C4", Pck022="#F0F4C3", Pck023="#B2EBF2", Pck024="#FFCCBC",
  Pck025="#D7CCC8"
)




# ── Significance helper ──────────────────────────────────────
add_sig <- function(df, pval_thr, padj_thr, lfc_thr) {
  if (is.null(df)||nrow(df)==0) return(df)
  use_adj <- !is.na(df$padj)
  p_eff   <- ifelse(use_adj, df$padj, df$pvalue)
  sig     <- !is.na(p_eff) & p_eff < ifelse(use_adj,padj_thr,pval_thr) &
             !is.na(df$log2FC) & abs(df$log2FC) >= lfc_thr
  df$sig_status <- dplyr::case_when(!sig~"NS", df$log2FC>0~"Up", TRUE~"Down")
  df
}

# ── Gene match (exact + PTM prefix) ──────────────────────────
match_genes <- function(df, genes) {
  gp        <- paste(vapply(genes,function(g) gsub("([.\\^$*+?|(){}\\[\\]])","\\\\\\1",g),""),collapse="|")
  pat_exact <- paste0("(?i)^(",gp,")$")
  pat_ptm   <- paste0("(?i)^(",gp,")_")
  df[grepl(pat_exact,df$gene_name,perl=TRUE)|grepl(pat_ptm,df$gene_name,perl=TRUE),]
}

# ── Safe empty datatable (avoids DT's "length zero" warning) ─
empty_dt <- function(msg="No results") {
  datatable(data.frame(message=msg), rownames=FALSE,
            options=list(dom="t", ordering=FALSE))
}

# ── Bulk dot plot ─────────────────────────────────────────────
make_bulk_dotplot <- function(df, xlabel_fmt="treat_time_model") {
  if (is.null(df)||nrow(df)==0)
    return(ggplot()+labs(title="No results — enter gene names and click Search")+theme_minimal())

  df$neglog10p <- pmin(-log10(pmax(df$pvalue,1e-300)),20)

  # Use cleaned columns when available, fall back to raw
  tc <- if("treatment_clean" %in% names(df) && !all(is.na(df$treatment_clean)))
          df$treatment_clean else df$treatment
  ma <- if("model_abbrev" %in% names(df) && !all(is.na(df$model_abbrev)))
          df$model_abbrev else abbrev_model_v(df$model)
  th <- if("time_h_clean" %in% names(df) && !all(is.na(df$time_h_clean)))
          df$time_h_clean else df$time_h

  df$ds_label <- make_ds_label(tc, ma, th, xlabel_fmt)

  # ── Hierarchical sort: treatment → time (numeric h) → model → pkg ─────
  df$sort_treat <- tc
  df$sort_time  <- suppressWarnings(as.numeric(gsub("h$","",th)))
  df$sort_model <- ma
  df <- df %>%
    arrange(omics_type, sort_treat, sort_time, sort_model, pkg_id)
  df$ds_label <- factor(df$ds_label, levels=unique(df$ds_label))
  df$gene_name <- droplevels(df$gene_name)

  n_genes <- nlevels(df$gene_name)

  bg <- df %>%
    distinct(ds_label, pkg_id, omics_type) %>%
    arrange(omics_type, as.integer(ds_label)) %>%
    group_by(omics_type) %>%
    mutate(x_pos = row_number()) %>%
    ungroup()

  pck_in_plot <- sort(unique(bg$pkg_id))
  pck_palette <- PKG_BG[pck_in_plot]

  df$tip <- paste0("<b>",df$gene_name,"</b>\n",
                   df$pkg_id," · ",df$omics_type,"\n",
                   "Treatment: ",tc,"\n",
                   "Model: ",ma,"  Time: ",th,"\n",
                   "log2FC: ",round(df$log2FC,2),
                   ifelse(!is.na(df$padj),
                          paste0("  padj: ",signif(df$padj,3)),
                          ifelse(!is.na(df$pvalue),
                                 paste0("  p: ",signif(df$pvalue,3)),"")))

  ggplot(df, aes(x=ds_label, y=gene_name)) +
    geom_rect(data=bg,
              aes(xmin=x_pos-0.5, xmax=x_pos+0.5,
                  ymin=0.5,       ymax=n_genes+0.5,
                  fill=pkg_id),
              alpha=0.28, inherit.aes=FALSE) +
    geom_point(aes(size=neglog10p, color=log2FC, text=tip), alpha=0.90) +
    facet_grid(.~omics_type, scales="free_x", space="free_x") +
    scale_color_gradient2(low="#1565C0", mid="white", high="#C62828",
                          midpoint=0, name="log2FC", na.value="grey70") +
    scale_size_continuous(range=c(2,9), name="-log10(p)") +
    scale_fill_manual(values=pck_palette, name="Study",
                      guide=guide_legend(ncol=1, override.aes=list(alpha=0.6,size=4))) +
    labs(x=NULL, y=NULL) +
    theme_bw(base_size=12) +
    theme(axis.text.x      = element_text(angle=55, hjust=1, size=7.5, lineheight=0.85),
          axis.text.y      = element_text(size=10),
          strip.text        = element_text(face="bold", size=10, color="#2c3e50"),
          strip.background  = element_rect(fill="#ecf0f1", color="#bdc3c7"),
          panel.grid.major  = element_line(color="grey92"),
          panel.spacing     = unit(0.5,"lines"),
          legend.position   = "right",
          legend.box        = "vertical",
          legend.box.just   = "top",
          legend.margin     = margin(0,0,4,0),
          legend.key.size   = unit(0.65,"lines"),
          legend.text       = element_text(size=8),
          legend.title      = element_text(size=9, face="bold"),
          plot.margin       = margin(4,8,4,4))
}

# ── Bulk bar chart ────────────────────────────────────────────
make_bulk_bar <- function(df) {
  if (is.null(df)||nrow(df)==0)
    return(ggplot()+labs(title="No results")+theme_minimal())

  tc <- if("treatment_clean" %in% names(df) && !all(is.na(df$treatment_clean))) df$treatment_clean else df$treatment
  ma <- if("model_abbrev"    %in% names(df) && !all(is.na(df$model_abbrev)))    df$model_abbrev    else abbrev_model_v(df$model)
  th <- if("time_h_clean"    %in% names(df) && !all(is.na(df$time_h_clean)))    df$time_h_clean    else df$time_h
  df$ds_label <- make_ds_label(tc, ma, th)
  df2 <- as.data.frame(df) %>%
    mutate(label_om = paste0("[",omics_type,"] ", ds_label)) %>%
    group_by(gene_name, label_om, omics_type) %>%
    summarise(log2FC=mean(log2FC,na.rm=TRUE), .groups="drop") %>%
    arrange(omics_type, log2FC)
  df2$label_om <- factor(df2$label_om, levels=unique(df2$label_om))

  ggplot(df2, aes(x=log2FC, y=label_om, fill=omics_type, text=gene_name)) +
    geom_col() +
    geom_vline(xintercept=0, color="black", linewidth=0.4) +
    facet_wrap(~gene_name, scales="free_x") +
    scale_fill_manual(values=OMICS_COLORS, na.value="grey50") +
    labs(x="log2FC", y=NULL, fill="Omics") +
    theme_bw(base_size=11) +
    theme(strip.text=element_text(face="bold"), axis.text.y=element_text(size=7))
}

# ── Single-cell dot plot ──────────────────────────────────────
make_sc_dotplot <- function(df) {
  if (is.null(df)||nrow(df)==0)
    return(ggplot()+labs(title="No results — enter gene names and click Search")+theme_minimal())

  df$neglog10p <- pmin(-log10(pmax(df$pvalue,1e-300)),20)
  df$ct_label  <- gsub("_markers","",df$cell_type,ignore.case=TRUE)
  df$gene_name <- droplevels(df$gene_name)

  df$tip <- paste0("<b>",df$gene_name,"</b>\n",
                   "Cell: ",df$ct_label,"\n",
                   "Treatment: ",df$comparison,"\n",
                   "log2FC: ",round(df$log2FC,2),
                   ifelse(!is.na(df$padj),paste0("\npadj: ",signif(df$padj,3)),
                          ifelse(!is.na(df$pvalue),paste0("\np: ",signif(df$pvalue,3)),"")))

  ggplot(df, aes(x=ct_label, y=gene_name,
                 size=neglog10p, color=log2FC, text=tip)) +
    geom_point(alpha=0.85) +
    facet_wrap(~comparison, scales="free_x") +
    scale_color_gradient2(low="#1565C0", mid="white", high="#C62828",
                          midpoint=0, name="log2FC", na.value="grey70") +
    scale_size_continuous(range=c(2,9), name="-log10(p)") +
    labs(x="Cell type", y=NULL) +
    theme_bw(base_size=12) +
    theme(axis.text.x   = element_text(angle=50, hjust=1, size=9),
          axis.text.y   = element_text(size=10),
          strip.text    = element_text(face="bold", size=10, color="#2c3e50"),
          strip.background = element_rect(fill="#ecf0f1", color="#bdc3c7"),
          panel.grid.major = element_line(color="grey92"),
          panel.spacing    = unit(0.4,"lines"))
}

# ===========================================================================
# UI
# ===========================================================================
ui <- page_navbar(
  title = "Multiomics Cytokine Explorer",
  theme = bs_theme(bootswatch="flatly", base_font=font_google("Inter")),
  bg    = "#2c3e50",

  header = tags$head(tags$style(HTML("
    .nav-tabs .nav-link.active {
      color:#2c3e50!important; background-color:#fff!important;
      border-color:#dee2e6 #dee2e6 #fff!important;
    }
    .nav-tabs .nav-link       { color:#495057!important; }
    .nav-tabs .nav-link:hover { color:#2c3e50!important; }
    .sidebar { overflow-y:auto; }
    .ptm-section { background:#fff8e1; border-left:4px solid #FFC107;
                   padding:12px 16px; border-radius:4px; margin-top:18px; }
    .ptm-section h5 { color:#795548; margin-bottom:6px; font-weight:700; }
  "))),

  # ── TAB 1: Overview ─────────────────────────────────────────
  nav_panel("Overview",
    layout_sidebar(
      sidebar = sidebar(width=280,
        h5("About"),
        tags$p(style="font-size:13px",
          "25 data packages from cytokine-treated in vitro models. ",
          "Use the tabs to search genes across bulk or single-cell data."),
        hr(),
        h6("Excluded / partial packages"),
        tags$small(
          tags$b("Not shown in plots:"),tags$br(),
          tags$b("Pck012, Pck017")," — alternative splicing (ΔPSI, no FC)",tags$br(),
          tags$b("Pck014")," — phospho-intensity matrix (no gene IDs)",tags$br(),
          tags$b("Pck015")," — top-down proteomics catalog (no stats)",tags$br(),
          tags$b("Pck021 Metabolomics")," — raw intensities only",tags$br(),tags$br(),
          tags$b("Replicate-averaged (exploratory):"),tags$br(),
          "Pck011, Pck018–021 Proteomics — treat as exploratory."
        )
      ),
      card(card_header("Package metadata — 25 data packages"),
           DTOutput("meta_table"))
    )
  ),

  # ── TAB 2: Bulk Gene Search ──────────────────────────────────
  nav_panel("Bulk Gene Search",
    layout_sidebar(
      sidebar = sidebar(width=300,
        h5("Search genes / proteins"),
        textAreaInput("gs_genes","Gene names (comma or newline separated)",
                      placeholder="STAT1, IRF1, MX1, IFIT1", rows=4),
        hr(),
        h6("Filter by experiment"),
        selectizeInput("gs_model","In vitro model",   choices=NULL, multiple=TRUE),
        selectizeInput("gs_treat","Treatment",        choices=NULL, multiple=TRUE),
        selectizeInput("gs_time", "Duration",         choices=NULL, multiple=TRUE),
        selectizeInput("gs_pkgs", "Data packages",    choices=NULL, multiple=TRUE),
        hr(),
        h6("X-axis label format"),
        selectInput("gs_xlabel", NULL,
          choices = c(
            "Treatment · Time · Model" = "treat_time_model",
            "Treatment · Model · Time" = "treat_model_time",
            "Model · Treatment · Time" = "model_treat_time",
            "Treatment · Time (compact)" = "treat_time"
          ), selected="treat_time_model"),
        hr(),
        h6("Data types"),
        uiOutput("gs_omics_ui"),
        hr(),
        h6("Significance thresholds"),
        sliderInput("gs_pval","p-value  ≤", 0,1,0.05,step=0.01),
        sliderInput("gs_lfc", "|log2FC| ≥",0,5,0.5, step=0.25),
        actionButton("gs_search","Search",class="btn-primary w-100 mt-2"),
        hr(),
        downloadButton("dl_gs","Download all results",
                       class="btn-sm btn-outline-primary")
      ),

      tagList(
        navset_tab(
          nav_panel("Dot Plot",      plotlyOutput("gs_dotplot", height="560px")),
          nav_panel("Bar Chart",     plotlyOutput("gs_bar",     height="520px")),
          nav_panel("Summary Table", DTOutput("gs_table"))
        ),
        div(class="ptm-section",
          h5("⚗ PTM Results — Phosphoproteomics / Redox Proteomics / Acetylomics"),
          tags$p(tags$small(style="color:#888",
            "PTM data is shown here separately. ",
            "gene_name format: PROTEIN_SITE (e.g. STAT1_T701). ",
            "Search by protein name — all matching sites are returned.")),
          DTOutput("gs_ptm_table")
        )
      )
    )
  ),

  # ── TAB 3: Single Cell ───────────────────────────────────────
  nav_panel("Single Cell",
    layout_sidebar(
      sidebar = sidebar(width=300,
        h5("Search genes in scRNA / snRNA"),
        textAreaInput("gsc_genes","Gene names (comma or newline separated)",
                      placeholder="STAT1, IRF1, IFIT1", rows=4),
        hr(),
        h6("Filter by experiment"),
        selectizeInput("gsc_model","In vitro model",     choices=NULL, multiple=TRUE),
        selectizeInput("gsc_ct",  "Cell types",          choices=NULL, multiple=TRUE),
        selectizeInput("gsc_comp","Treatments",          choices=NULL, multiple=TRUE),
        hr(),
        h6("Significance thresholds"),
        sliderInput("gsc_pval","p-value  ≤", 0,1,0.05,step=0.01),
        sliderInput("gsc_lfc", "|log2FC| ≥",0,5,0.25,step=0.25),
        actionButton("gsc_search","Search",class="btn-primary w-100 mt-2"),
        hr(),
        downloadButton("dl_gsc","Download results",
                       class="btn-sm btn-outline-primary")
      ),
      navset_tab(
        nav_panel("Dot Plot",      plotlyOutput("gsc_dotplot", height="560px")),
        nav_panel("Summary Table", DTOutput("gsc_table"))
      )
    )
  )
)

# ===========================================================================
# SERVER
# ===========================================================================
server <- function(input, output, session) {

  # Tracks which omics types the user has chosen after a search
  gs_selected_types <- reactiveVal(BULK_TYPES)

  # Render the current selection as a compact tag list with a "Change" button
  output$gs_omics_ui <- renderUI({
    sel <- gs_selected_types()
    tagList(
      tags$small(style="color:#555",
        if (length(sel)==0) "None selected"
        else paste(sel, collapse=", ")),
      tags$br(),
      actionButton("gs_omics_btn","Change data types…",
                   class="btn-sm btn-outline-secondary mt-1 w-100")
    )
  })

  # Show data-type modal when button clicked
  observeEvent(input$gs_omics_btn, {
    avail <- gs_selected_types()
    showModal(modalDialog(
      title = "Select data types to display",
      checkboxGroupInput("gs_omics_modal", NULL,
        choices  = BULK_TYPES,
        selected = avail),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("gs_omics_apply","Apply",class="btn-primary")
      )
    ))
  })

  # Apply modal selection
  observeEvent(input$gs_omics_apply, {
    sel <- input$gs_omics_modal
    gs_selected_types(if (length(sel)==0) BULK_TYPES else sel)
    removeModal()
  })

  # After a new gene search, auto-detect available types for those genes
  # and restrict display to those (user can open modal to change)
  observeEvent(input$gs_search, {
    req(nzchar(trimws(input$gs_genes)))
    genes <- parse_genes(input$gs_genes)
    df    <- FLAT_BULK
    sel_m <- input$gs_model[nzchar(input$gs_model)]
    sel_t <- input$gs_treat[nzchar(input$gs_treat)]
    sel_h <- input$gs_time[ nzchar(input$gs_time)]
    sel_p <- input$gs_pkgs[ nzchar(input$gs_pkgs)]
    if (length(sel_m)>0) df <- df[df$model_abbrev    %in% sel_m, ]
    if (length(sel_t)>0) df <- df[df$treatment_clean %in% sel_t, ]
    if (length(sel_h)>0) df <- df[df$time_h_clean    %in% sel_h, ]
    if (length(sel_p)>0) df <- df[df$pkg_id          %in% sel_p, ]
    matched <- match_genes(df, genes)
    avail   <- intersect(BULK_TYPES, unique(matched$omics_type))
    if (length(avail)==0) avail <- BULK_TYPES
    gs_selected_types(avail)
  }, ignoreInit=TRUE, priority=10)   # runs before gs_all


  # ── Overview ────────────────────────────────────────────────
  output$meta_table <- renderDT({
    datatable(META_RAW, filter="top", rownames=FALSE,
              options=list(pageLength=25, scrollX=TRUE))
  })

  # ── Bulk Gene Search: cascading filters ─────────────────────
  # Drive entirely from FLAT_BULK (actual loaded data), not metadata table.
  # Uses model_abbrev / treatment_clean / time_h_clean normalised columns.

  all_models <- reactive({
    sort(unique(na.omit(FLAT_BULK$model_abbrev)))
  })
  observe({
    updateSelectizeInput(session, "gs_model",
      choices = c("(all models)" = "", all_models()), server = TRUE)
  })

  # Treatments filtered by selected models
  avail_treats <- reactive({
    df <- FLAT_BULK
    sel_m <- input$gs_model[nzchar(input$gs_model)]
    if (length(sel_m) > 0) df <- df[df$model_abbrev %in% sel_m, , drop = FALSE]
    sort(unique(na.omit(df$treatment_clean)))
  })
  observe({
    updateSelectizeInput(session, "gs_treat",
      choices = c("(all treatments)" = "", avail_treats()), server = TRUE)
  })

  # Times filtered by model + treatment
  avail_times <- reactive({
    df <- FLAT_BULK
    sel_m <- input$gs_model[nzchar(input$gs_model)]
    sel_t <- input$gs_treat[nzchar(input$gs_treat)]
    if (length(sel_m) > 0) df <- df[df$model_abbrev    %in% sel_m, , drop = FALSE]
    if (length(sel_t) > 0) df <- df[df$treatment_clean %in% sel_t, , drop = FALSE]
    sort(unique(na.omit(df$time_h_clean)))
  })
  observe({
    updateSelectizeInput(session, "gs_time",
      choices = c("(all durations)" = "", avail_times()), server = TRUE)
  })

  # Packages filtered by model + treatment + time — show only "Pck001" label
  avail_pkgs <- reactive({
    df <- FLAT_BULK
    sel_m <- input$gs_model[nzchar(input$gs_model)]
    sel_t <- input$gs_treat[nzchar(input$gs_treat)]
    sel_h <- input$gs_time[ nzchar(input$gs_time) ]
    if (length(sel_m) > 0) df <- df[df$model_abbrev    %in% sel_m, , drop = FALSE]
    if (length(sel_t) > 0) df <- df[df$treatment_clean %in% sel_t, , drop = FALSE]
    if (length(sel_h) > 0) df <- df[df$time_h_clean    %in% sel_h, , drop = FALSE]
    sort(unique(na.omit(df$pkg_id)))
  })
  observe({
    p <- avail_pkgs()
    updateSelectizeInput(session, "gs_pkgs",
      choices = c("(all packages)" = "", setNames(p, p)), server = TRUE)
  })

  # ── Bulk Gene Search: core reactive ──────────────────────────
  parse_genes <- function(raw) {
    g <- trimws(unlist(strsplit(raw,"[,\n]+")))
    g[nzchar(g)]
  }

  gs_all <- eventReactive(input$gs_search, {
    raw <- input$gs_genes; req(nzchar(trimws(raw)))
    genes <- parse_genes(raw); req(length(genes)>0)

    df <- FLAT_BULK

    # Apply cascading filters using normalised columns
    sel_m <- input$gs_model[nzchar(input$gs_model)]
    sel_t <- input$gs_treat[nzchar(input$gs_treat)]
    sel_h <- input$gs_time[ nzchar(input$gs_time) ]
    sel_p <- input$gs_pkgs[ nzchar(input$gs_pkgs) ]
    if (length(sel_m) > 0) df <- df[df$model_abbrev    %in% sel_m, ]
    if (length(sel_t) > 0) df <- df[df$treatment_clean %in% sel_t, ]
    if (length(sel_h) > 0) df <- df[df$time_h_clean    %in% sel_h, ]
    if (length(sel_p) > 0) df <- df[df$pkg_id          %in% sel_p, ]

    matched <- match_genes(df, genes)
    if (nrow(matched)==0) return(list(bulk=NULL, ptm=NULL, genes=genes))

    matched <- add_sig(matched, input$gs_pval, input$gs_pval, input$gs_lfc)
    # Factor levels: queried genes first, then any PTM variants
    all_lev <- unique(c(genes, as.character(matched$gene_name)))
    matched$gene_name <- factor(matched$gene_name, levels=rev(all_lev))

    # Split bulk (non-PTM) vs PTM
    bulk_df <- matched[matched$omics_type %in% gs_selected_types(),]
    # KEY FIX: drop PTM-named factor levels so y-axis only shows genes with data
    bulk_df$gene_name <- droplevels(bulk_df$gene_name)

    list(
      bulk  = bulk_df,
      ptm   = matched[matched$omics_type %in% PTM_TYPES,],
      genes = genes
    )
  })

  output$gs_dotplot <- renderPlotly({
    res <- gs_all(); req(!is.null(res))
    df  <- res$bulk
    if (is.null(df)||nrow(df)==0)
      return(ggplotly(ggplot()+
               labs(title="No bulk data found for these genes / filters")+
               theme_minimal()))
    p <- make_bulk_dotplot(df, xlabel_fmt=input$gs_xlabel)
    ggplotly(p, tooltip="text") %>%
      layout(legend=list(orientation="v", x=1.02, y=1))
  })

  output$gs_bar <- renderPlotly({
    res <- gs_all(); req(!is.null(res))
    df  <- res$bulk
    if (is.null(df)||nrow(df)==0)
      return(ggplotly(ggplot()+labs(title="No results")+theme_minimal()))
    p <- make_bulk_bar(df)
    ggplotly(p, tooltip=c("text","x"))
  })

  output$gs_table <- renderDT({
    res <- gs_all()
    df  <- if (!is.null(res)) as.data.frame(res$bulk) else NULL
    if (is.null(df)||nrow(df)==0) return(empty_dt("No results — enter genes and click Search"))
    show <- c("gene_name","omics_type","treatment","model","log2FC",
              "pvalue","padj","sig_status","comparison","pkg_id","replicate_averaged")
    show <- intersect(show, names(df))
    datatable(df[,show], filter="top", rownames=FALSE,
              options=list(pageLength=30, scrollX=TRUE)) %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })

  output$gs_ptm_table <- renderDT({
    res <- gs_all()
    df  <- if (!is.null(res)) as.data.frame(res$ptm) else NULL
    if (is.null(df)||nrow(df)==0)
      return(empty_dt("No PTM data found for these genes"))
    show <- c("gene_name","omics_type","treatment","model","log2FC",
              "pvalue","padj","sig_status","comparison","pkg_id")
    show <- intersect(show, names(df))
    datatable(df[,show], filter="top", rownames=FALSE,
              caption="gene_name format: PROTEIN_SITE (e.g. STAT1_T701)",
              options=list(pageLength=20, scrollX=TRUE)) %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })

  output$dl_gs <- downloadHandler(
    filename=function() paste0("bulk_search_",
      gsub("[, \n]+","_",trimws(input$gs_genes)),".csv"),
    content=function(f) {
      res <- gs_all()
      if (!is.null(res)) {
        all_df <- rbind(
          if (!is.null(res$bulk)&&nrow(as.data.frame(res$bulk))>0)
            as.data.frame(res$bulk) else NULL,
          if (!is.null(res$ptm) &&nrow(as.data.frame(res$ptm))>0)
            as.data.frame(res$ptm)  else NULL
        )
        if (!is.null(all_df)) write.csv(all_df, f, row.names=FALSE)
      }
    }
  )

  # ── Single Cell: filters ─────────────────────────────────────
  observe({
    models <- sort(unique(na.omit(FLAT_SC$model)))
    updateSelectizeInput(session,"gsc_model",
      choices=c("(all models)"="", models), server=TRUE)
  })

  observe({
    df <- FLAT_SC
    sel_m <- input$gsc_model[nzchar(input$gsc_model)]
    if (length(sel_m)>0) df <- df[df$model %in% sel_m,]
    cts   <- sort(unique(na.omit(gsub("_markers","",df$cell_type,ignore.case=TRUE))))
    comps <- sort(unique(na.omit(df$comparison)))
    updateSelectizeInput(session,"gsc_ct",   choices=cts,   selected=cts,   server=TRUE)
    updateSelectizeInput(session,"gsc_comp", choices=comps, selected=comps, server=TRUE)
  })

  # ── Single Cell: core reactive ────────────────────────────────
  gsc_all <- eventReactive(input$gsc_search, {
    raw <- input$gsc_genes; req(nzchar(trimws(raw)))
    genes <- parse_genes(raw); req(length(genes)>0)

    df <- FLAT_SC
    sel_m <- input$gsc_model[nzchar(input$gsc_model)]
    if (length(sel_m)>0) df <- df[df$model %in% sel_m,]
    if (length(input$gsc_ct)>0)
      df <- df[gsub("_markers","",df$cell_type,ignore.case=TRUE) %in% input$gsc_ct,]
    if (length(input$gsc_comp)>0)
      df <- df[df$comparison %in% input$gsc_comp,]

    matched <- match_genes(df, genes)
    if (nrow(matched)==0) return(list(data=NULL, genes=genes))

    matched <- add_sig(matched, input$gsc_pval, input$gsc_pval, input$gsc_lfc)
    all_lev <- unique(c(genes, as.character(matched$gene_name)))
    matched$gene_name <- factor(matched$gene_name, levels=rev(all_lev))
    matched$gene_name <- droplevels(matched$gene_name)

    list(data=matched, genes=genes)
  })

  output$gsc_dotplot <- renderPlotly({
    res <- gsc_all(); req(!is.null(res))
    df  <- res$data
    if (is.null(df)||nrow(df)==0)
      return(ggplotly(ggplot()+labs(title="No scRNA results found")+theme_minimal()))
    p <- make_sc_dotplot(df)
    ggplotly(p, tooltip="text") %>%
      layout(legend=list(orientation="v"))
  })

  output$gsc_table <- renderDT({
    res <- gsc_all()
    df  <- if (!is.null(res)) as.data.frame(res$data) else NULL
    if (is.null(df)||nrow(df)==0)
      return(empty_dt("No results — enter genes and click Search"))
    show <- c("gene_name","cell_type","comparison","treatment","model",
              "log2FC","pvalue","padj","sig_status","pkg_id")
    show <- intersect(show, names(df))
    datatable(df[,show], filter="top", rownames=FALSE,
              options=list(pageLength=30, scrollX=TRUE)) %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })

  output$dl_gsc <- downloadHandler(
    filename=function() paste0("sc_search_",
      gsub("[, \n]+","_",trimws(input$gsc_genes)),".csv"),
    content=function(f) {
      res <- gsc_all()
      if (!is.null(res$data))
        write.csv(as.data.frame(res$data), f, row.names=FALSE)
    }
  )
}

shinyApp(ui, server)
