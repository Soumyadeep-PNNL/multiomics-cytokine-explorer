# ============================================================
# Multiomics In Vitro Cytokine Explorer  — v4
# Tabs: Overview | Bulk Omics | scRNA / snRNA | Gene Search | Multi-Dataset
# ============================================================

# ── Working directory auto-fix ──────────────────────────────
local({
  sd <- tryCatch(
    normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)),
    error = function(e) ""
  )
  if (!nzchar(sd)) {
    sd <- tryCatch(
      normalizePath(dirname(sys.frame(1)$ofile)),
      error = function(e) ""
    )
  }
  if (nzchar(sd) && file.exists(file.path(sd, "data_loader.R"))) setwd(sd)
})

library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(scales)
library(grid)

source("data_loader.R")

# ── Load data (cache-aware) ──────────────────────────────────
CACHE_FILE <- "data_cache.rds"
if (file.exists(CACHE_FILE)) {
  cat("Loading from cache...\n")
  DATA <- readRDS(CACHE_FILE)
} else {
  cat("Parsing all packages (first run — may take 1-2 min)...\n")
  DATA <- load_all_data()
  saveRDS(DATA, CACHE_FILE, compress = "xz")
  cat("Cache saved.\n")
}

FLAT      <- DATA$flat
FLAT_BULK <- DATA$flat_bulk
FLAT_SC   <- DATA$flat_sc
PKGS      <- names(DATA$packages)
META      <- DATA$metadata

# Pkg → display label
pkg_label <- function(id, meta) {
  row <- meta[grepl(id, meta[[1]], fixed = TRUE), , drop = FALSE]
  if (nrow(row) == 0) return(id)
  # Try to find model and treatment columns
  cols <- tolower(names(row))
  m_col <- grep("model|cell", cols)[1]; t_col <- grep("treat|cyto", cols)[1]
  parts <- c(id,
             if (!is.na(m_col)) as.character(row[[m_col]]),
             if (!is.na(t_col)) as.character(row[[t_col]]))
  parts <- parts[!is.na(parts) & nzchar(parts)]
  paste(parts, collapse = " — ")
}
PKG_CHOICES <- setNames(PKGS, vapply(PKGS, pkg_label, character(1), meta = META))

# Omics type palette
OMICS_COLORS <- c(
  "RNA-seq"      = "#2196F3",
  "scRNA-seq"    = "#9C27B0",
  "ATAC-seq"     = "#4CAF50",
  "ChIP-seq"     = "#FF9800",
  "Proteomics"   = "#F44336",
  "Lipidomics"   = "#00BCD4",
  "Metabolomics" = "#795548",
  "UMI-4C"       = "#607D8B"
)

MAX_PLOT_ROWS <- 20000

cap_for_plot <- function(df, max_ns = MAX_PLOT_ROWS) {
  if (is.null(df) || nrow(df) == 0) return(df)
  sig <- df[df$sig_status != "NS", ]
  ns  <- df[df$sig_status == "NS",  ]
  if (nrow(ns) > max_ns) ns <- ns[sample(nrow(ns), max_ns), ]
  rbind(sig, ns)
}

add_sig <- function(df, pval_thr, padj_thr, lfc_thr) {
  if (is.null(df) || nrow(df) == 0) return(df)
  # Use padj if available, fall back to pvalue
  use_adj <- !is.na(df$padj)
  p_eff   <- ifelse(use_adj, df$padj, df$pvalue)
  sig     <- !is.na(p_eff) & p_eff < ifelse(use_adj, padj_thr, pval_thr) &
             !is.na(df$log2FC) & abs(df$log2FC) >= lfc_thr
  df$sig_status <- dplyr::case_when(
    !sig           ~ "NS",
    df$log2FC > 0  ~ "Up",
    TRUE           ~ "Down"
  )
  df
}

# ── Volcano helper ───────────────────────────────────────────
make_volcano <- function(df, title = "", pval_thr = 0.05, lfc_thr = 1,
                         top_n = 20) {
  if (is.null(df) || nrow(df) == 0)
    return(ggplot() + labs(title = "No data") + theme_minimal())
  df$neglog10p <- -log10(pmax(df$pvalue, 1e-300))
  df$color <- df$sig_status
  top <- df %>% filter(sig_status != "NS") %>%
    arrange(desc(neglog10p)) %>% slice_head(n = top_n)
  ggplot(df, aes(x = log2FC, y = neglog10p,
                 color = color, text = gene_name)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(pval_thr), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = top, aes(label = gene_name),
                    size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c("NS" = "grey70", "Up" = "#E53935", "Down" = "#1E88E5"),
                       name = "Status") +
    labs(title = title, x = "log2(Fold Change)", y = "-log10(p-value)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "right")
}

# ── MA plot helper ───────────────────────────────────────────
make_ma <- function(df, title = "", pval_thr = 0.05, lfc_thr = 1) {
  if (is.null(df) || nrow(df) == 0)
    return(ggplot() + labs(title = "No data") + theme_minimal())
  # MA: Y = log2FC, X = mean expression proxy (-log10 pvalue)
  df$xval <- -log10(pmax(df$pvalue, 1e-300))
  top <- df %>% filter(sig_status != "NS") %>%
    arrange(desc(abs(log2FC))) %>% slice_head(n = 15)
  ggplot(df, aes(x = xval, y = log2FC, color = sig_status, text = gene_name)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = top, aes(label = gene_name),
                    size = 3, max.overlaps = 15, show.legend = FALSE) +
    scale_color_manual(values = c("NS"="#cccccc","Up"="#E53935","Down"="#1E88E5")) +
    labs(title = title, x = "-log10(p-value)", y = "log2(Fold Change)") +
    theme_bw(base_size = 13)
}

# ── Heatmap helper ───────────────────────────────────────────
make_heatmap <- function(df, top_n = 50, title = "") {
  sig <- df %>% filter(sig_status != "NS") %>%
    arrange(desc(abs(log2FC))) %>% slice_head(n = top_n)
  if (nrow(sig) == 0) sig <- df %>% arrange(desc(abs(log2FC))) %>% slice_head(n = top_n)
  if (nrow(sig) == 0) return(NULL)
  mat <- as.matrix(sig[, "log2FC", drop = FALSE])
  rownames(mat) <- sig$gene_name
  colnames(mat) <- "log2FC"
  breaks <- seq(-max(abs(mat), na.rm = TRUE), max(abs(mat), na.rm = TRUE), length.out = 101)
  pheatmap(mat, color = colorRampPalette(c("#1E88E5","white","#E53935"))(100),
           breaks = breaks, cluster_cols = FALSE, fontsize_row = 8,
           main = title, silent = TRUE)
}

# ── Gene dot plot (search results) ──────────────────────────
make_gene_dotplot <- function(df, genes) {
  if (is.null(df) || nrow(df) == 0)
    return(ggplot() + labs(title="No results") + theme_minimal())
  df$neglog10p <- pmin(-log10(pmax(df$pvalue, 1e-300)), 20)
  df$gene_name <- factor(df$gene_name, levels = rev(genes))
  df$label     <- paste0(df$pkg_id, "\n", df$layer_name)

  ggplot(df, aes(x = label, y = gene_name,
                 size  = neglog10p,
                 color = log2FC,
                 shape = sig_status)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(low="#1E88E5", mid="white", high="#E53935",
                          midpoint = 0, name = "log2FC") +
    scale_size_continuous(range = c(2, 8), name = "-log10(p)") +
    scale_shape_manual(values = c("NS"=1,"Up"=16,"Down"=16), name="Status") +
    labs(x = NULL, y = NULL, title = "Gene × Dataset dot plot") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          panel.grid.major = element_line(color = "grey92"))
}

# ── Multi-dataset heatmap ────────────────────────────────────
make_multi_heatmap <- function(df, top_n = 40) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  sig <- df %>% filter(sig_status != "NS")
  if (nrow(sig) == 0) sig <- df
  top_genes <- sig %>%
    group_by(gene_name) %>%
    summarise(max_lfc = max(abs(log2FC), na.rm = TRUE)) %>%
    arrange(desc(max_lfc)) %>% slice_head(n = top_n) %>% pull(gene_name)
  mat_df <- df %>%
    filter(gene_name %in% top_genes) %>%
    group_by(gene_name, pkg_id) %>%
    summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = pkg_id, values_from = log2FC, values_fill = 0)
  mat <- as.matrix(mat_df[, -1])
  rownames(mat) <- mat_df$gene_name
  if (nrow(mat) < 2 || ncol(mat) < 1) return(NULL)
  lim <- max(abs(mat), na.rm = TRUE)
  breaks <- seq(-lim, lim, length.out = 101)
  pheatmap(mat, color = colorRampPalette(c("#1E88E5","white","#E53935"))(100),
           breaks = breaks, fontsize_row = 8, fontsize_col = 9,
           main = "Top variable genes across datasets", silent = TRUE)
}

# ===========================================================================
# UI
# ===========================================================================
ui <- page_navbar(
  title = "Multiomics Cytokine Explorer",
  theme = bs_theme(bootswatch = "flatly", base_font = font_google("Inter")),
  bg    = "#2c3e50",

  # ── TAB 1: Overview ─────────────────────────────────────────────────────
  nav_panel("Overview",
    layout_sidebar(
      sidebar = sidebar(width = 280,
        h5("Filter packages"),
        selectizeInput("ov_model", "Cell/Tissue model", choices=NULL, multiple=TRUE),
        selectizeInput("ov_treat", "Treatment",         choices=NULL, multiple=TRUE),
        selectizeInput("ov_omics", "Omics type",        choices=NULL, multiple=TRUE),
        checkboxInput("ov_sc", "Show scRNA/snRNA only", value=FALSE)
      ),
      card(
        card_header("Package metadata — 25 data packages"),
        DTOutput("meta_table")
      )
    )
  ),

  # ── TAB 2: Bulk Omics ───────────────────────────────────────────────────
  nav_panel("Bulk Omics",
    layout_sidebar(
      sidebar = sidebar(width = 290,
        h5("Select data"),
        selectizeInput("s_model", "Model", choices=c("All"=""), multiple=FALSE),
        selectizeInput("s_treat", "Treatment", choices=c("All"=""), multiple=FALSE),
        selectizeInput("s_pkg",   "Package", choices=PKG_CHOICES, multiple=FALSE),
        selectizeInput("s_layer", "Omics layer", choices=NULL, multiple=FALSE),
        selectizeInput("s_comp",  "Comparison", choices=NULL, multiple=FALSE),
        hr(),
        h5("Significance thresholds"),
        sliderInput("s_pval", "p-value", 0, 1, 0.05, step=0.01),
        sliderInput("s_padj", "adj p-value", 0, 1, 0.05, step=0.01),
        sliderInput("s_lfc",  "|log2FC|", 0, 5, 1, step=0.25),
        hr(),
        checkboxInput("s_sig_only", "Show significant only", FALSE),
        checkboxInput("s_hide_genomic", "Hide genomic coordinates", TRUE),
        downloadButton("dl_bulk", "Download table", class="btn-sm btn-outline-primary")
      ),
      navset_tab(
        nav_panel("Volcano",  plotlyOutput("bulk_volcano", height="520px")),
        nav_panel("MA Plot",  plotlyOutput("bulk_ma",      height="520px")),
        nav_panel("Heatmap",  plotOutput( "bulk_heatmap",  height="560px")),
        nav_panel("Data Table", DTOutput("bulk_table"))
      )
    )
  ),

  # ── TAB 3: scRNA / snRNA ────────────────────────────────────────────────
  nav_panel("scRNA / snRNA",
    layout_sidebar(
      sidebar = sidebar(width = 290,
        h5("Select data"),
        selectizeInput("sc_pkg",   "Package", choices=NULL, multiple=FALSE),
        selectizeInput("sc_ct",    "Cell type", choices=NULL, multiple=TRUE),
        selectizeInput("sc_comp",  "Comparison (treatment)", choices=NULL, multiple=TRUE),
        hr(),
        h5("Significance thresholds"),
        sliderInput("sc_pval", "p-value",    0, 1, 0.05, step=0.01),
        sliderInput("sc_padj", "adj p-value",0, 1, 0.05, step=0.01),
        sliderInput("sc_lfc",  "|log2FC|",   0, 5, 0.5,  step=0.25),
        downloadButton("dl_sc", "Download table", class="btn-sm btn-outline-primary")
      ),
      navset_tab(
        nav_panel("Volcano",    plotlyOutput("sc_volcano", height="520px")),
        nav_panel("MA Plot",    plotlyOutput("sc_ma",      height="520px")),
        nav_panel("Heatmap",    plotOutput( "sc_heatmap",  height="560px")),
        nav_panel("Data Table", DTOutput("sc_table"))
      )
    )
  ),

  # ── TAB 4: Gene Search ──────────────────────────────────────────────────
  nav_panel("Gene Search",
    layout_sidebar(
      sidebar = sidebar(width = 290,
        h5("Search genes / proteins"),
        textAreaInput("gs_genes", "Gene names (comma or newline separated)",
                      placeholder="STAT1, IRF1, MX1, IFIT1", rows=5),
        checkboxGroupInput("gs_omics", "Include omics types",
          choices  = c("RNA-seq","Proteomics","Lipidomics","Metabolomics","scRNA-seq","ATAC-seq","ChIP-seq","UMI-4C"),
          selected = c("RNA-seq","Proteomics","Lipidomics","Metabolomics","scRNA-seq")),
        sliderInput("gs_pval", "p-value threshold", 0, 1, 0.05, step=0.01),
        sliderInput("gs_lfc",  "|log2FC| min",      0, 5, 0.5,  step=0.25),
        actionButton("gs_search", "Search", class="btn-primary w-100 mt-2"),
        hr(),
        downloadButton("dl_gs", "Download results", class="btn-sm btn-outline-primary")
      ),
      navset_tab(
        nav_panel("Dot Plot",   plotlyOutput("gs_dotplot",  height="600px")),
        nav_panel("Summary Table", DTOutput("gs_table")),
        nav_panel("Bar Chart",  plotlyOutput("gs_bar",      height="480px"))
      )
    )
  ),

  # ── TAB 5: Multi-Dataset ────────────────────────────────────────────────
  nav_panel("Multi-Dataset",
    layout_sidebar(
      sidebar = sidebar(width = 290,
        h5("Compare packages"),
        selectizeInput("m_model", "Filter by model",     choices=c("All"=""), multiple=FALSE),
        selectizeInput("m_treat", "Filter by treatment", choices=c("All"=""), multiple=FALSE),
        selectizeInput("m_omics", "Omics type",          choices=NULL),
        selectizeInput("m_pkgs",  "Packages (≥2)",       choices=PKG_CHOICES, multiple=TRUE),
        hr(),
        sliderInput("m_pval", "p-value",    0, 1, 0.05, step=0.01),
        sliderInput("m_padj", "adj p-value",0, 1, 0.05, step=0.01),
        sliderInput("m_lfc",  "|log2FC|",   0, 5, 1,    step=0.25)
      ),
      navset_tab(
        nav_panel("Heatmap",  plotOutput("m_heatmap", height="600px")),
        nav_panel("Summary",  DTOutput("m_table"))
      )
    )
  )
)

# ===========================================================================
# SERVER
# ===========================================================================
server <- function(input, output, session) {

  # ── Debounced slider inputs ──────────────────────────────────────────────
  d_sp  <- debounce(reactive(input$s_pval),  400)
  d_sa  <- debounce(reactive(input$s_padj),  400)
  d_sl  <- debounce(reactive(input$s_lfc),   400)
  d_scp <- debounce(reactive(input$sc_pval), 400)
  d_sca <- debounce(reactive(input$sc_padj), 400)
  d_scl <- debounce(reactive(input$sc_lfc),  400)
  d_mp  <- debounce(reactive(input$m_pval),  400)
  d_ma  <- debounce(reactive(input$m_padj),  400)
  d_ml  <- debounce(reactive(input$m_lfc),   400)

  # ── Overview tab ────────────────────────────────────────────────────────
  observe({
    models  <- sort(unique(na.omit(FLAT$model)))
    treats  <- sort(unique(na.omit(FLAT$treatment)))
    omics   <- sort(unique(na.omit(FLAT$omics_type)))
    updateSelectizeInput(session, "ov_model", choices=c("All"="", models), server=TRUE)
    updateSelectizeInput(session, "ov_treat", choices=c("All"="", treats), server=TRUE)
    updateSelectizeInput(session, "ov_omics", choices=c("All"="", omics),  server=TRUE)
  })

  output$meta_table <- renderDT({
    df <- META
    datatable(df, filter="top", rownames=FALSE,
              options=list(pageLength=25, scrollX=TRUE))
  })

  # ── Bulk Omics: cascading filters ───────────────────────────────────────
  observe({
    models <- sort(unique(na.omit(FLAT_BULK$model)))
    updateSelectizeInput(session, "s_model", choices=c("All"="", models), server=TRUE)
  })

  avail_pkgs_bulk <- reactive({
    df <- FLAT_BULK
    if (nzchar(input$s_model)) df <- df[df$model == input$s_model, ]
    if (nzchar(input$s_treat)) df <- df[df$treatment == input$s_treat, ]
    sort(unique(df$pkg_id))
  })

  observe({
    pkgs   <- avail_pkgs_bulk()
    treats <- sort(unique(na.omit(FLAT_BULK$treatment[FLAT_BULK$pkg_id %in% pkgs])))
    choices_t <- c("All"="", treats)
    updateSelectizeInput(session, "s_treat", choices=choices_t, server=TRUE)
    choices_p <- PKG_CHOICES[names(PKG_CHOICES) %in% pkgs | unname(PKG_CHOICES) %in% pkgs]
    # Rebuild from available pkgs
    valid <- PKG_CHOICES[PKG_CHOICES %in% pkgs]
    updateSelectizeInput(session, "s_pkg", choices=valid, server=TRUE)
  })

  observe({
    req(input$s_pkg)
    layers <- sort(unique(FLAT_BULK$layer_name[FLAT_BULK$pkg_id == input$s_pkg]))
    updateSelectizeInput(session, "s_layer", choices=layers, server=TRUE)
  })

  observe({
    req(input$s_pkg, input$s_layer)
    comps <- sort(unique(FLAT_BULK$comparison[
      FLAT_BULK$pkg_id == input$s_pkg & FLAT_BULK$layer_name == input$s_layer]))
    updateSelectizeInput(session, "s_comp",
      choices=c("All comparisons"="", comps), server=TRUE)
  })

  bulk_data <- reactive({
    req(input$s_pkg)
    df <- FLAT_BULK[FLAT_BULK$pkg_id == input$s_pkg, ]
    if (nzchar(input$s_layer)) df <- df[df$layer_name == input$s_layer, ]
    if (nzchar(input$s_comp))  df <- df[df$comparison  == input$s_comp, ]
    if (input$s_hide_genomic)  df <- df[!df$is_genomic, ]
    df <- add_sig(df, d_sp(), d_sa(), d_sl())
    if (input$s_sig_only) df <- df[df$sig_status != "NS", ]
    df
  })

  bulk_plot_data <- reactive({ cap_for_plot(bulk_data()) })

  output$bulk_volcano <- renderPlotly({
    df <- bulk_plot_data(); req(nrow(df) > 0)
    p  <- make_volcano(df, title=paste(input$s_pkg, input$s_layer),
                       pval_thr=d_sp(), lfc_thr=d_sl())
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="v"))
  })

  output$bulk_ma <- renderPlotly({
    df <- bulk_plot_data(); req(nrow(df) > 0)
    p  <- make_ma(df, title=paste(input$s_pkg, input$s_layer),
                  pval_thr=d_sp(), lfc_thr=d_sl())
    ggplotly(p, tooltip="text")
  })

  output$bulk_heatmap <- renderPlot({
    df <- bulk_data(); req(nrow(df) > 0)
    ph <- make_heatmap(df, title=paste(input$s_pkg, input$s_layer))
    if (!is.null(ph)) print(ph)
  })

  output$bulk_table <- renderDT({
    df <- bulk_data()
    show_cols <- c("gene_name","log2FC","pvalue","padj","comparison",
                   "sig_status","omics_type","replicate_averaged")
    show_cols <- intersect(show_cols, names(df))
    datatable(df[, show_cols], filter="top", rownames=FALSE,
              options=list(pageLength=20, scrollX=TRUE)) %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })

  output$dl_bulk <- downloadHandler(
    filename = function() paste0(input$s_pkg, "_", input$s_layer, ".csv"),
    content  = function(f) write.csv(bulk_data(), f, row.names=FALSE)
  )

  # ── scRNA tab ────────────────────────────────────────────────────────────
  observe({
    pkgs <- sort(unique(FLAT_SC$pkg_id))
    updateSelectizeInput(session, "sc_pkg", choices=pkgs, server=TRUE)
  })

  observe({
    req(input$sc_pkg)
    cts <- sort(unique(na.omit(FLAT_SC$cell_type[FLAT_SC$pkg_id == input$sc_pkg])))
    updateSelectizeInput(session, "sc_ct", choices=cts, selected=cts, server=TRUE)
  })

  observe({
    req(input$sc_pkg)
    comps <- sort(unique(FLAT_SC$comparison[FLAT_SC$pkg_id == input$sc_pkg]))
    updateSelectizeInput(session, "sc_comp", choices=comps, selected=comps, server=TRUE)
  })

  sc_data <- reactive({
    req(input$sc_pkg)
    df <- FLAT_SC[FLAT_SC$pkg_id == input$sc_pkg, ]
    if (length(input$sc_ct)   > 0) df <- df[df$cell_type %in% input$sc_ct, ]
    if (length(input$sc_comp) > 0) df <- df[df$comparison %in% input$sc_comp, ]
    df <- add_sig(df, d_scp(), d_sca(), d_scl())
    df
  })

  output$sc_volcano <- renderPlotly({
    df <- cap_for_plot(sc_data()); req(nrow(df) > 0)
    p  <- make_volcano(df, title=paste(input$sc_pkg, "scRNA"),
                       pval_thr=d_scp(), lfc_thr=d_scl())
    ggplotly(p, tooltip="text")
  })

  output$sc_ma <- renderPlotly({
    df <- cap_for_plot(sc_data()); req(nrow(df) > 0)
    p  <- make_ma(df, title=paste(input$sc_pkg, "scRNA"),
                  pval_thr=d_scp(), lfc_thr=d_scl())
    ggplotly(p, tooltip="text")
  })

  output$sc_heatmap <- renderPlot({
    df <- sc_data(); req(nrow(df) > 0)
    ph <- make_heatmap(df, title=paste(input$sc_pkg, "scRNA"))
    if (!is.null(ph)) print(ph)
  })

  output$sc_table <- renderDT({
    df <- sc_data()
    show_cols <- c("gene_name","log2FC","pvalue","padj","cell_type",
                   "comparison","sig_status")
    show_cols <- intersect(show_cols, names(df))
    datatable(df[, show_cols], filter="top", rownames=FALSE,
              options=list(pageLength=20, scrollX=TRUE)) %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })

  output$dl_sc <- downloadHandler(
    filename=function() paste0(input$sc_pkg,"_scRNA.csv"),
    content =function(f) write.csv(sc_data(), f, row.names=FALSE)
  )

  # ── Gene Search tab ──────────────────────────────────────────────────────
  gs_results <- eventReactive(input$gs_search, {
    raw <- input$gs_genes
    req(nzchar(trimws(raw)))
    genes <- trimws(unlist(strsplit(raw, "[,\n]+")))
    genes <- genes[nzchar(genes)]
    req(length(genes) > 0)

    df <- FLAT
    if (length(input$gs_omics) > 0)
      df <- df[df$omics_type %in% input$gs_omics, ]
    df <- df[!df$is_genomic, ]

    # Case-insensitive gene match
    pattern <- paste0("^(", paste(genes, collapse="|"), ")$")
    matched <- df[grepl(pattern, df$gene_name, ignore.case=TRUE), ]
    if (nrow(matched) == 0) return(NULL)

    matched <- add_sig(matched, input$gs_pval, input$gs_pval, input$gs_lfc)
    matched$gene_name <- factor(matched$gene_name,
                                levels=genes[genes %in% matched$gene_name])
    matched
  })

  output$gs_dotplot <- renderPlotly({
    df <- gs_results(); req(!is.null(df) && nrow(df)>0)
    genes <- levels(df$gene_name)
    p <- make_gene_dotplot(df, genes)
    ggplotly(p, tooltip=c("text","size","color"))
  })

  output$gs_table <- renderDT({
    df <- gs_results()
    if (is.null(df)) return(datatable(data.frame(message="No results found")))
    show_cols <- c("gene_name","pkg_id","layer_name","omics_type",
                   "comparison","log2FC","pvalue","padj","sig_status","replicate_averaged")
    show_cols <- intersect(show_cols, names(df))
    datatable(as.data.frame(df)[, show_cols], filter="top", rownames=FALSE,
              options=list(pageLength=30, scrollX=TRUE)) %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })

  output$gs_bar <- renderPlotly({
    df <- gs_results(); req(!is.null(df) && nrow(df)>0)
    df2 <- as.data.frame(df) %>%
      group_by(gene_name, pkg_id, omics_type) %>%
      summarise(log2FC=mean(log2FC, na.rm=TRUE),
                pvalue=min(pvalue, na.rm=TRUE), .groups="drop")
    p <- ggplot(df2, aes(x=reorder(paste0(pkg_id,"\n",omics_type), log2FC),
                          y=log2FC, fill=omics_type, text=gene_name)) +
      geom_col() +
      geom_hline(yintercept=0) +
      facet_wrap(~gene_name, scales="free_y") +
      coord_flip() +
      scale_fill_manual(values=OMICS_COLORS) +
      labs(x=NULL, y="log2FC") +
      theme_bw(base_size=11)
    ggplotly(p, tooltip=c("text","y"))
  })

  output$dl_gs <- downloadHandler(
    filename=function() paste0("gene_search_",
      gsub("[, \n]+","_",trimws(input$gs_genes)),".csv"),
    content=function(f) {
      df <- gs_results()
      if (!is.null(df)) write.csv(as.data.frame(df), f, row.names=FALSE)
    }
  )

  # ── Multi-Dataset tab ─────────────────────────────────────────────────────
  observe({
    omics_types <- sort(unique(FLAT_BULK$omics_type))
    updateSelectizeInput(session, "m_omics", choices=omics_types, server=TRUE)
    models  <- sort(unique(na.omit(FLAT_BULK$model)))
    treats  <- sort(unique(na.omit(FLAT_BULK$treatment)))
    updateSelectizeInput(session, "m_model", choices=c("All"="", models), server=TRUE)
    updateSelectizeInput(session, "m_treat", choices=c("All"="", treats), server=TRUE)
  })

  m_data <- reactive({
    req(length(input$m_pkgs) >= 1)
    df <- FLAT_BULK[FLAT_BULK$pkg_id %in% input$m_pkgs, ]
    if (nzchar(input$m_omics)) df <- df[df$omics_type == input$m_omics, ]
    df <- add_sig(df, d_mp(), d_ma(), d_ml())
    df
  })

  output$m_heatmap <- renderPlot({
    df <- m_data(); req(nrow(df) > 0)
    ph <- make_multi_heatmap(df)
    if (!is.null(ph)) print(ph)
  })

  output$m_table <- renderDT({
    df <- m_data()
    show_cols <- c("gene_name","pkg_id","layer_name","omics_type",
                   "log2FC","pvalue","padj","comparison","sig_status")
    show_cols <- intersect(show_cols, names(df))
    datatable(df[, show_cols], filter="top", rownames=FALSE,
              options=list(pageLength=25, scrollX=TRUE)) %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })
}

shinyApp(ui, server)
