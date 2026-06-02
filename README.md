# Multiomics Cytokine Explorer

This is a dataset resource from previously published multi-omics studies of pancreatic islets and β-cell models exposed to pro-inflammatory cytokines. The omics datasets are aggregated and quality-controlled into unified “data packages” for sharing via PNNL DataHub. The resource currently comprises 70 datasets from 23 studies organized into 25 data packages (Pck001–Pck025) spanning data from primary human islets, human β-cell lines (EndoC-βH1, iPSC-derived islet-like cells), mouse islets, and mouse β-cell lines (MIN6, β-TC-3, β-TC-6).

---

## Overview

This interactive R Shiny application provides a unified interface to explore how pro-inflammatory cytokines (IFNα, IFNγ, IL-1β, TNFα) remodel the molecular landscape of pancreatic β-cells across multiple omics layers and in vitro model systems. It enables cross-package gene and protein queries without requiring programming expertise, supporting both exploratory data analysis and hypothesis generation relevant to type 1 diabetes and β-cell dysfunction.

---

## Data Packages

| Package | Model | Treatment | Time | Omics |
|---------|-------|-----------|------|-------|
| Pck001 | Human pancreatic islets (HPI) | IFNγ + IL-1β | 48h | RNA-seq, ATAC-seq, UMI-4C, ChIP-seq |
| Pck002 | EndoC-βH1 | IFNγ + IL-1β | 48h | RNA-seq, ATAC-seq, Proteomics, Lipidomics, Metabolomics |
| Pck003 | HPI / EndoC-βH1 | IFNα | 2h, 8h, 24h | RNA-seq, ATAC-seq, Proteomics, Lipidomics, Metabolomics |
| Pck004 | HPI | BFA, TG, IL-1β+IFNγ+TNFα, IFNγ, IL-1β, TNFα | 48h | scRNA-seq |
| Pck005 | βTC-6 | IFNγ + IL-1β | 4h, 8h, 24h | Proteomics, Redox Proteomics, Phosphoproteomics, Acetylomics |
| Pck006 | HPI | IFNγ + IL-1β | 24h | Proteomics, Lipidomics, Metabolomics |
| Pck007 | EndoC-βH1 | IFNα, IFNγ, IL-1β, TNFα | 24h | RNA-seq |
| Pck008 | HPI | IFNγ + IL-1β | 24h | RNA-seq, Proteomics |
| Pck009 | iPSC-derived islets | IFNα | 24h | RNA-seq, scRNA-seq |
| Pck010 | EndoC-βH1 | IFNα | 8h, 24h | RNA-seq |
| Pck011 | HPI | IFNγ + IL-1β | 48h | RNA-seq |
| Pck014 | HPI | IFNγ + IL-1β | 24h | Phosphoproteomics |
| Pck015 | HPI | IFNγ + IL-1β | 24h | Proteomics (top-down) |
| Pck016 | βTC-3 | IL-1β + IFNγ + TNFα | 24h | RNA-seq |
| Pck018 | MIN6 | IL-1β + IFNγ + TNFα | 24h | Proteomics |
| Pck019 | MIN6 | IL-1β + IFNγ + TNFα | 24h | Proteomics |
| Pck020 | MIN6 | IL-1β + IFNγ + TNFα | 24h | Proteomics |
| Pck021 | MIN6 | IL-1β + IFNγ + TNFα | 8h, 24h | Proteomics, Metabolomics |
| Pck022 | Mouse pancreatic islets (MPI) | IFNγ, IL-1β, IFNγ+IL-1β | 6h | scRNA-seq |
| Pck023 | MPI | IFNγ + IL-1β | 18h | scRNA-seq |
| Pck024 | HPI | IFNγ, IL-1β, IFNγ+IL-1β | 6h, 18h | scRNA-seq |
| Pck025 | HPI | IFNγ + IL-1β | 18h | RNA-seq (Discovery + Validation), Proteomics |

> Packages Pck012 and Pck017 contain alternative splicing data (ΔPSI) and are shown in the Overview table but excluded from gene search plots.

---

## Application Features

### Bulk Gene Search
- Enter any number of gene or protein names (comma or newline separated)
- Cascading filters for **in vitro model**, **cytokine treatment**, **duration**, and **data package**
- Treatment filter shows cytokines only (IFNα, IFNγ, IL-1β, TNFα and combinations); non-cytokine additives are stripped automatically
- Duration filter shows clean time points (2h, 4h, 6h, 8h, 18h, 24h, 48h) derived per-row from the data
- After search, a **data type popup** shows which omics types returned results — select only those you want to display
- **Dot plot** with hierarchical x-axis ordering (treatment → duration → cell model), significance-coded point size (−log10 p-value) and colour (log2FC)
- **Four x-axis label formats** selectable from sidebar (Treatment·Time·Model, Treatment·Model·Time, Model·Treatment·Time, compact Treatment·Time)
- **Bar chart** of mean log2FC per gene per dataset
- **Summary table** with download

### Single Cell
- Search genes across scRNA-seq / snRNA-seq data from Pck004, Pck009, Pck022, Pck023, Pck024
- Filter by model, cell type, and treatment comparison
- Dot plot coloured by log2FC, sized by −log10(p-value), faceted by cell type

### Overview
- Searchable metadata table listing all 25 data packages with model, treatment, timepoint, omics types, replicates, PMID, and repository accession

---

## Repository Structure

```
App/
├── app.R               # Shiny UI + server
├── data_loader.R       # Per-package Excel parsers (Pck001–Pck025)
├── build_cache.R       # One-time cache builder (run before first launch)
├── install_packages.R  # R package installer
├── HOW_TO_RUN.txt      # Quick-start instructions
└── data_cache.rds      # Pre-built data cache (auto-generated, not tracked in git)

Files/
├── Table 1-Metadata.xlsx        # Master metadata table
└── Pck_Analysis/
    ├── Pck001.xlsx
    ├── Pck002.xlsx
    └── ... (Pck001–Pck025)
```

---

## Quick Start

**1. Install R dependencies** (first time only):
```r
source("App/install_packages.R")
```

**2. Run the app** from RStudio — open `App/app.R` and click **Run App**, or from the R console:
```r
shiny::runApp("App", launch.browser = TRUE)
```

On first launch the app parses all 25 Excel files and saves a cache (`data_cache.rds`). This takes approximately 1–2 minutes. Subsequent launches load instantly from cache.

---

## Technical Notes

### Data Normalization
- **Treatments** are normalized to canonical cytokine names using regex matching; non-cytokine additives (DFMO, palmitate, NMMA, ML351, BFA, TG) are stripped or excluded
- **Time points** are extracted per-row: the `comparison` column is scanned for a `\bNh\b` pattern first; the metadata `time_h` field is used as fallback for single-value entries
- **Model names** are abbreviated: Human pancreatic islets → HPI, EndoC-βH1 → EndoC, βTC-6 → βTC-6, Mouse pancreatic islets → MPI, iPSC-derived → iPSC, all MIN6 variants → MIN6

### Statistical Methods
- **Fisher's combined probability test** (χ² = −2 Σ ln pᵢ, df = 2k) used for Pck011 to merge p-values across 5 patient samples
- **One-sample t-test** (H₀: μ = 0) used for Pck020 replicate log2FC values
- **Welch's t-test** used for Pck021 Metabolomics (CT_Eth vs NCT_Eth on log₂-transformed peak areas)
- Packages without per-gene statistics report replicate-averaged log2FC and are flagged as exploratory

### PTM Data
Phosphoproteomics and redox proteomics entries use a `GENE_POSITION` naming convention (e.g., `STAT1_701`) to distinguish modification sites. These appear in a dedicated PTM table below the main dot plot.

---

## Citation

If you use this resource, please cite the original publications associated with each data package (PMIDs listed in the Overview tab).

---

## Contact

Soumyadeep Sarkar — Pacific Northwest National Laboratory (PNNL)
