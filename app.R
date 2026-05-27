# ============================================================
# Multiomics In Vitro Cytokine Explorer — v5
# Tabs: Overview | Bulk Gene Search | Single Cell
# ============================================================
# GIT REMINDER — run after every session:
#   cd "/Users/sark224/Library/CloudStorage/OneDrive-PNNL/Documents/Projects/In-vitroCT treated exp-Shinnyapp"
#   git add App/app.R App/data_loader.R
#   git commit -m "describe changes here"
#   git push
# ============================================================

local({
  sd <- tryCatch(
    normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)),
    error = function(e) "")
  if (!nzchar(sd))
    sd <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile)), error=function(e) "")
  if (nzchar(sd) && file.exists(file.path(sd,"data_loader.R"))) setwd(sd)
})

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(dplyr); library(tidyr)
  library(ggplot2); library(plotly); library(DT)
  library(ggrepel); library(RColorBrewer); library(viridis); library(scales); library(grid)
})

source("data_loader.R")

# ── Load data (cache-aware) ──────────────────────────────────
CACHE_FILE <- "data_cache.rds"
if (file.exists(CACHE_FILE)) {
  cat("Loading from cache...\n")
  DATA <- readRDS(CACHE_FILE)
} else {
  cat("Parsing all packages (first run — may take 1–2 min)...\n")
  DATA <- load_all_data()
  saveRDS(DATA, CACHE_FILE, compress="xz")
  cat("Cache saved.\n")
}

FLAT      <- DATA$flat
FLAT_BULK <- DATA$flat_bulk
FLAT_SC   <- DATA$flat_sc
META      <- DATA$metadata

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

MAX_PLOT_ROWS <- 20000

cap_for_plot <- function(df, max_ns=MAX_PLOT_ROWS) {
  if (is.null(df)||nrow(df)==0) return(df)
  sig <- df[df$sig_status!="NS",]; ns <- df[df$sig_status=="NS",]
  if (nrow(ns)>max_ns) ns <- ns[sample(nrow(ns),max_ns),]
  rbind(sig,ns)
}

add_sig <- function(df, pval_thr, padj_thr, lfc_thr) {
  if (is.null(df)||nrow(df)==0) return(df)
  use_adj <- !is.na(df$padj)
  p_eff   <- ifelse(use_adj, df$padj, df$pvalue)
  sig     <- !is.na(p_eff) & p_eff < ifelse(use_adj,padj_thr,pval_thr) &
             !is.na(df$log2FC) & abs(df$log2FC) >= lfc_thr
  df$sig_status <- dplyr::case_when(!sig~"NS", df$log2FC>0~"Up", TRUE~"Down")
  df
}

# ── Short dataset label: "Treatment [#N]" ─────────────────────
ds_label <- function(treatment, pkg_id) {
  n <- sub("^Pck0*","",pkg_id)
  paste0(treatment,"\n[#",n,"]")
}

# ── Bulk dot plot — faceted by omics type, PTM excluded ───────
make_bulk_dotplot <- function(df) {
  if (is.null(df)||nrow(df)==0)
    return(ggplot()+labs(title="No results — enter gene names and click Search")+theme_minimal())

  df$neglog10p <- pmin(-log10(pmax(df$pvalue,1e-300)),20)
  df$ds_label  <- ds_label(df$treatment, df$pkg_id)

  # Sort: within each facet, order x by treatment then pkg
  df <- df %>% arrange(omics_type, treatment, pkg_id)
  df$ds_label  <- factor(df$ds_label, levels=unique(df$ds_label))
  # Preserve gene order from caller
  if (!is.factor(df$gene_name)) df$gene_name <- factor(df$gene_name)

  df$tip <- paste0("<b>",df$gene_name,"</b>\n",
                   df$pkg_id," | ",df$omics_type,"\n",
                   "Treatment: ",df$treatment,"\n",
                   "log2FC: ",round(df$log2FC,2),"\n",
                   "p: ",signif(df$pvalue,3),
                   ifelse(!is.na(df$padj), paste0("  padj: ",signif(df$padj,3)),""))

  ggplot(df, aes(x=ds_label, y=gene_name,
                 size=neglog10p, color=log2FC, text=tip)) +
    geom_point(alpha=0.85) +
    facet_grid(.~omics_type, scales="free_x", space="free_x") +
    scale_color_gradient2(low="#1565C0", mid="white", high="#C62828",
                          midpoint=0, name="log2FC", na.value="grey70") +
    scale_size_continuous(range=c(2,9), name="-log10(p)") +
    labs(x=NULL, y=NULL) +
    theme_bw(base_size=12) +
    theme(axis.text.x    = element_text(angle=50,hjust=1,size=8),
          axis.text.y    = element_text(size=10),
          strip.text     = element_text(face="bold",size=10,color="#2c3e50"),
          strip.background = element_rect(fill="#ecf0f1",color="#bdc3c7"),
          panel.grid.major = element_line(color="grey92"),
          panel.spacing    = unit(0.5,"lines"))
}

# ── Bulk bar chart — grouped by omics type ────────────────────
make_bulk_bar <- function(df) {
  if (is.null(df)||nrow(df)==0)
    return(ggplot()+labs(title="No results")+theme_minimal())

  df2 <- as.data.frame(df) %>%
    mutate(ds=ds_label(treatment,pkg_id)) %>%
    group_by(gene_name, ds, omics_type) %>%
    summarise(log2FC=mean(log2FC,na.rm=TRUE),.groups="drop") %>%
    arrange(omics_type, log2FC)

  df2$ds_om <- factor(paste0("[",df2$omics_type,"] ",df2$ds),
                       levels=unique(paste0("[",df2$omics_type,"] ",df2$ds)))

  ggplot(df2, aes(x=log2FC, y=ds_om, fill=omics_type, text=gene_name)) +
    geom_col() +
    geom_vline(xintercept=0, color="black", linewidth=0.4) +
    facet_wrap(~gene_name, scales="free_x") +
    scale_fill_manual(values=OMICS_COLORS, na.value="grey50") +
    labs(x="log2FC", y=NULL, fill="Omics type") +
    theme_bw(base_size=11) +
    theme(strip.text=element_text(face="bold"),
          axis.text.y=element_text(size=8))
}

# ── Single-cell dot plot — genes × cell type, faceted by treatment ──
make_sc_dotplot <- function(df) {
  if (is.null(df)||nrow(df)==0)
    return(ggplot()+labs(title="No results — enter gene names and click Search")+theme_minimal())

  df$neglog10p <- pmin(-log10(pmax(df$pvalue,1e-300)),20)
  df$ct_label  <- gsub("_markers","",df$cell_type,ignore.case=TRUE)
  if (!is.factor(df$gene_name)) df$gene_name <- factor(df$gene_name)

  df$tip <- paste0("<b>",df$gene_name,"</b>\n",
                   "Cell: ",df$ct_label,"\n",
                   "Treatment: ",df$comparison,"\n",
                   "log2FC: ",round(df$log2FC,2),"\n",
                   "p: ",signif(df$pvalue,3),
                   ifelse(!is.na(df$padj),paste0("  padj: ",signif(df$padj,3)),""))

  ggplot(df, aes(x=ct_label, y=gene_name,
                 size=neglog10p, color=log2FC, text=tip)) +
    geom_point(alpha=0.85) +
    facet_wrap(~comparison, scales="free_x") +
    scale_color_gradient2(low="#1565C0", mid="white", high="#C62828",
                          midpoint=0, name="log2FC", na.value="grey70") +
    scale_size_continuous(range=c(2,9), name="-log10(p)") +
    labs(x="Cell type", y=NULL) +
    theme_bw(base_size=12) +
    theme(axis.text.x   = element_text(angle=50,hjust=1,size=9),
          axis.text.y   = element_text(size=10),
          strip.text    = element_text(face="bold",size=10,color="#2c3e50"),
          strip.background = element_rect(fill="#ecf0f1",color="#bdc3c7"),
          panel.grid.major = element_line(color="grey92"),
          panel.spacing    = unit(0.5,"lines"))
}

# ── Gene match helper ─────────────────────────────────────────
match_genes <- function(df, genes) {
  gene_pat  <- paste(vapply(genes,
    function(g) gsub("([.\\^$*+?|(){}\\[\\]])","\\\\\\1",g),""),collapse="|")
  pat_exact <- paste0("(?i)^(",gene_pat,")$")
  pat_ptm   <- paste0("(?i)^(",gene_pat,")_")
  df[grepl(pat_exact,df$gene_name,perl=TRUE)|
     grepl(pat_ptm,  df$gene_name,perl=TRUE),]
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
      color: #2c3e50 !important; background-color: #ffffff !important;
      border-color: #dee2e6 #dee2e6 #fff !important;
    }
    .nav-tabs .nav-link       { color: #495057 !important; }
    .nav-tabs .nav-link:hover { color: #2c3e50 !important; }
    .sidebar { overflow-y: auto; }
    .ptm-section { background:#fff8e1; border-left:4px solid #FFC107;
                   padding:12px 16px; border-radius:4px; margin-top:18px; }
    .ptm-section h5 { color:#795548; margin-bottom:8px; font-weight:700; }
  "))),

  # ── TAB 1: Overview ─────────────────────────────────────────────────────
  nav_panel("Overview",
    layout_sidebar(
      sidebar = sidebar(width=280,
        h5("Filter packages"),
        selectizeInput("ov_model","Cell / Tissue model",choices=NULL,multiple=TRUE),
        selectizeInput("ov_treat","Treatment",          choices=NULL,multiple=TRUE),
        selectizeInput("ov_omics","Omics type",         choices=NULL,multiple=TRUE),
        checkboxInput("ov_sc","Show scRNA / snRNA only",value=FALSE),
        hr(),
        h6("Excluded / partial packages"),
        tags$small(
          tags$b("Not shown in plots:"),tags$br(),
          tags$b("Pck012, Pck017")," — alternative splicing (ΔPSI, no FC)",tags$br(),
          tags$b("Pck014")," — phospho-intensity matrix (no gene IDs)",tags$br(),
          tags$b("Pck015")," — top-down proteomics catalog (no stats)",tags$br(),
          tags$b("Pck021 Metabolomics")," — raw intensities only (no FC)",tags$br(),
          tags$br(),
          tags$b("Replicate-averaged (exploratory):"),tags$br(),
          "Pck011, Pck018–021 Proteomics — no per-gene p-value; treat as exploratory."
        )
      ),
      card(card_header("Package metadata — 25 data packages"), DTOutput("meta_table"))
    )
  ),

  # ── TAB 2: Bulk Gene Search ──────────────────────────────────────────────
  nav_panel("Bulk Gene Search",
    layout_sidebar(
      sidebar = sidebar(width=290,
        h5("Search genes / proteins"),
        textAreaInput("gs_genes","Gene names (comma or newline separated)",
                      placeholder="STAT1, IRF1, MX1, IFIT1", rows=5),
        hr(),
        h6("Filters"),
        selectizeInput("gs_model","Experimental model",choices=NULL,multiple=TRUE),
        checkboxGroupInput("gs_omics","Omics types to include",
          choices  = BULK_TYPES,
          selected = BULK_TYPES),
        hr(),
        h6("Significance thresholds"),
        sliderInput("gs_pval","p-value  ≤", 0,1,0.05,step=0.01),
        sliderInput("gs_lfc", "|log2FC| ≥",0,5,0.5, step=0.25),
        actionButton("gs_search","Search",class="btn-primary w-100 mt-2"),
        hr(),
        downloadButton("dl_gs","Download all results",class="btn-sm btn-outline-primary")
      ),

      tagList(
        navset_tab(
          nav_panel("Dot Plot",      plotlyOutput("gs_dotplot",  height="560px")),
          nav_panel("Bar Chart",     plotlyOutput("gs_bar",      height="520px")),
          nav_panel("Summary Table", DTOutput("gs_table"))
        ),
        # PTM section — shown below the main plot
        div(class="ptm-section",
          h5(icon("flask"), " PTM Results — Phosphoproteomics / Redox / Acetylomics"),
          p(tags$small(style="color:#888",
            "PTM data is excluded from the dot plot above to avoid mixing site-level ",
            "data with protein/gene-level data. Rows here correspond to individual ",
            "modification sites (gene_name = PROTEIN_SITE, e.g. STAT1_T701).")),
          DTOutput("gs_ptm_table")
        )
      )
    )
  ),

  # ── TAB 3: Single Cell ──────────────────────────────────────────────────
  nav_panel("Single Cell",
    layout_sidebar(
      sidebar = sidebar(width=290,
        h5("Search genes in scRNA / snRNA"),
        textAreaInput("gsc_genes","Gene names (comma or newline separated)",
                      placeholder="STAT1, IRF1, IFIT1", rows=5),
        hr(),
        h6("Filters"),
        selectizeInput("gsc_model","Experimental model",choices=NULL,multiple=TRUE),
        selectizeInput("gsc_ct",  "Cell types",          choices=NULL,multiple=TRUE),
        selectizeInput("gsc_comp","Treatments / comparisons",choices=NULL,multiple=TRUE),
        hr(),
        h6("Significance thresholds"),
        sliderInput("gsc_pval","p-value  ≤", 0,1,0.05,step=0.01),
        sliderInput("gsc_lfc", "|log2FC| ≥",0,5,0.25,step=0.25),
        actionButton("gsc_search","Search",class="btn-primary w-100 mt-2"),
        hr(),
        downloadButton("dl_gsc","Download results",class="btn-sm btn-outline-primary")
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

  # ── Overview ──────────────────────────────────────────────────────────────
  observe({
    updateSelectizeInput(session,"ov_model",choices=c("All"="",sort(unique(na.omit(FLAT$model)))),server=TRUE)
    updateSelectizeInput(session,"ov_treat",choices=c("All"="",sort(unique(na.omit(FLAT$treatment)))),server=TRUE)
    updateSelectizeInput(session,"ov_omics",choices=c("All"="",sort(unique(na.omit(FLAT$omics_type)))),server=TRUE)
  })

  output$meta_table <- renderDT({
    datatable(META, filter="top", rownames=FALSE,
              options=list(pageLength=25, scrollX=TRUE))
  })

  # ── Bulk Gene Search: populate model choices ──────────────────────────────
  observe({
    models <- sort(unique(na.omit(FLAT_BULK$model)))
    updateSelectizeInput(session,"gs_model",
      choices=c("(all models)"="", models), server=TRUE)
  })

  # ── Bulk Gene Search: core reactive ───────────────────────────────────────
  parse_genes <- function(raw) {
    g <- trimws(unlist(strsplit(raw,"[,\n]+")))
    g[nzchar(g)]
  }

  gs_all <- eventReactive(input$gs_search, {
    raw <- input$gs_genes; req(nzchar(trimws(raw)))
    genes <- parse_genes(raw); req(length(genes)>0)

    df <- FLAT[!FLAT$is_genomic,]
    if (length(input$gs_model)>0 && any(nzchar(input$gs_model)))
      df <- df[df$model %in% input$gs_model[nzchar(input$gs_model)],]

    matched <- match_genes(df, genes)
    if (nrow(matched)==0) return(list(bulk=NULL, ptm=NULL, genes=genes))

    matched <- add_sig(matched, input$gs_pval, input$gs_pval, input$gs_lfc)

    # Compute factor levels: queried genes first, then PTM site variants
    all_levels <- unique(c(genes, as.character(matched$gene_name)))
    matched$gene_name <- factor(matched$gene_name, levels=rev(all_levels))

    list(
      bulk  = matched[matched$omics_type %in% input$gs_omics,],
      ptm   = matched[matched$omics_type %in% PTM_TYPES,],
      genes = genes
    )
  })

  # Dot plot — bulk only (no PTM)
  output$gs_dotplot <- renderPlotly({
    res <- gs_all(); req(!is.null(res))
    df <- res$bulk
    if (is.null(df)||nrow(df)==0)
      return(ggplotly(ggplot()+labs(title="No bulk results found")+theme_minimal()))
    p <- make_bulk_dotplot(df)
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="v"))
  })

  # Bar chart — bulk only
  output$gs_bar <- renderPlotly({
    res <- gs_all(); req(!is.null(res))
    df <- res$bulk
    if (is.null(df)||nrow(df)==0)
      return(ggplotly(ggplot()+labs(title="No results")+theme_minimal()))
    p <- make_bulk_bar(df)
    ggplotly(p, tooltip=c("text","x"))
  })

  # Summary table — bulk
  output$gs_table <- renderDT({
    res <- gs_all()
    df  <- if (!is.null(res)) res$bulk else NULL
    if (is.null(df)||nrow(df)==0)
      return(datatable(data.frame(message="No results — enter genes and click Search")))
    show <- c("gene_name","omics_type","treatment","model","log2FC","pvalue","padj",
              "sig_status","comparison","pkg_id","layer_name","replicate_averaged")
    show <- intersect(show, names(df))
    datatable(as.data.frame(df)[,show], filter="top", rownames=FALSE,
              options=list(pageLength=30, scrollX=TRUE)) %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })

  # PTM table — only PTM types
  output$gs_ptm_table <- renderDT({
    res <- gs_all()
    df  <- if (!is.null(res)) res$ptm else NULL
    if (is.null(df)||nrow(df)==0)
      return(datatable(data.frame(message="No PTM data found for these genes")))
    show <- c("gene_name","omics_type","treatment","model","log2FC","pvalue","padj",
              "sig_status","comparison","pkg_id","layer_name")
    show <- intersect(show, names(df))
    datatable(as.data.frame(df)[,show], filter="top", rownames=FALSE,
              options=list(pageLength=20, scrollX=TRUE),
              caption="gene_name format: PROTEIN_SITE (e.g. STAT1_T701). Search by protein name to retrieve all its sites.") %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })

  output$dl_gs <- downloadHandler(
    filename=function() paste0("gene_search_",gsub("[, \n]+","_",trimws(input$gs_genes)),".csv"),
    content =function(f) {
      res <- gs_all()
      if (!is.null(res)) {
        all_df <- rbind(
          if (!is.null(res$bulk)) as.data.frame(res$bulk) else NULL,
          if (!is.null(res$ptm))  as.data.frame(res$ptm)  else NULL
        )
        write.csv(all_df, f, row.names=FALSE)
      }
    }
  )

  # ── Single Cell: populate filters ─────────────────────────────────────────
  observe({
    models <- sort(unique(na.omit(FLAT_SC$model)))
    updateSelectizeInput(session,"gsc_model",
      choices=c("(all models)"="", models), server=TRUE)
  })

  # Cell types and comparisons depend on model selection
  observe({
    df <- FLAT_SC
    if (length(input$gsc_model)>0 && any(nzchar(input$gsc_model)))
      df <- df[df$model %in% input$gsc_model[nzchar(input$gsc_model)],]
    cts   <- sort(unique(na.omit(gsub("_markers","",df$cell_type,ignore.case=TRUE))))
    comps <- sort(unique(na.omit(df$comparison)))
    updateSelectizeInput(session,"gsc_ct",   choices=cts,   selected=cts,   server=TRUE)
    updateSelectizeInput(session,"gsc_comp", choices=comps, selected=comps, server=TRUE)
  })

  # ── Single Cell: core reactive ─────────────────────────────────────────────
  gsc_all <- eventReactive(input$gsc_search, {
    raw <- input$gsc_genes; req(nzchar(trimws(raw)))
    genes <- parse_genes(raw); req(length(genes)>0)

    df <- FLAT_SC
    if (length(input$gsc_model)>0 && any(nzchar(input$gsc_model)))
      df <- df[df$model %in% input$gsc_model[nzchar(input$gsc_model)],]
    if (length(input$gsc_ct)>0)
      df <- df[gsub("_markers","",df$cell_type,ignore.case=TRUE) %in% input$gsc_ct,]
    if (length(input$gsc_comp)>0)
      df <- df[df$comparison %in% input$gsc_comp,]

    matched <- match_genes(df, genes)
    if (nrow(matched)==0) return(list(data=NULL, genes=genes))

    matched <- add_sig(matched, input$gsc_pval, input$gsc_pval, input$gsc_lfc)
    all_levels <- unique(c(genes, as.character(matched$gene_name)))
    matched$gene_name <- factor(matched$gene_name, levels=rev(all_levels))

    list(data=matched, genes=genes)
  })

  output$gsc_dotplot <- renderPlotly({
    res <- gsc_all(); req(!is.null(res))
    df  <- res$data
    if (is.null(df)||nrow(df)==0)
      return(ggplotly(ggplot()+labs(title="No scRNA results found")+theme_minimal()))
    p <- make_sc_dotplot(df)
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="v"))
  })

  output$gsc_table <- renderDT({
    res <- gsc_all()
    df  <- if (!is.null(res)) res$data else NULL
    if (is.null(df)||nrow(df)==0)
      return(datatable(data.frame(message="No results — enter genes and click Search")))
    show <- c("gene_name","cell_type","comparison","treatment","model",
              "log2FC","pvalue","padj","sig_status","pkg_id")
    show <- intersect(show, names(df))
    datatable(as.data.frame(df)[,show], filter="top", rownames=FALSE,
              options=list(pageLength=30, scrollX=TRUE)) %>%
      formatRound(c("log2FC","pvalue","padj"), digits=4)
  })

  output$dl_gsc <- downloadHandler(
    filename=function() paste0("sc_search_",gsub("[, \n]+","_",trimws(input$gsc_genes)),".csv"),
    content =function(f) {
      res <- gsc_all()
      if (!is.null(res$data)) write.csv(as.data.frame(res$data), f, row.names=FALSE)
    }
  )
}

shinyApp(ui, server)
