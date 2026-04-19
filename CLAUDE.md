# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Computational immunology research project profiling the Non-Small Cell Lung Cancer (NSCLC) tumor microenvironment by integrating bulk and single-cell RNA-seq data from three cohorts (LungPredict, Vanderbilt, TCGA). Published in *Frontiers in Immunology* (2024).

## Environment Setup

All code is R-based. Scripts live in `Scripts/` and use `renv` for reproducibility (R 4.3.1).

```r
# One-time setup — run inside Scripts/ in R console
install.packages('renv')
renv::activate()
renv::restore()

# Deactivate when done
renv::deactivate()
```

Open `Scripts/LungPredict.Rproj` in RStudio to work within the project.

## Running Analyses

Analyses are RMarkdown notebooks (`.Rmd`) run sequentially in RStudio. Intended execution order:

1. `Data_preprocessing.Rmd` — load counts, filter protein-coding genes, convert ENSEMBL→SYMBOL, QC/PCA
2. `LP_all_analysis.Rmd` — LungPredict cohort, all stages (I–IV)
3. `LP_early_analysis.Rmd` — LungPredict early stage (I–II), most comprehensive analysis
4. `Vanderbilt_analysis.Rmd` — Vanderbilt early-stage adenocarcinoma (bulk RNA-seq)
5. `Vanderbilt_scRNAseq_analysis.Rmd` — Vanderbilt single-cell RNA-seq (15 patients)
6. `Vanderbilt_survival_analysis.Rmd` — survival/prognostic analysis
7. `TCGA_data_preprocessing.Rmd` — download and preprocess TCGA data
8. `TCGA_analysis.Rmd` — TCGA validation

Results are written to `Results/<Cohort>/` as CSVs and PDFs. All analyses use `set.seed(123)`.

## Code Architecture

### `src/` — shared helper functions sourced by each notebook

- **`environment_set.R`**: Loads 50+ libraries and defines 20+ reusable functions covering normalization, PCA, transcription factor (TF) activity inference (decoupleR/viper), pathway/enrichment analysis (GSVA, ReactomePA, fgsea), WGCNA module identification, survival analysis (coxph, KM), and Boruta feature selection.
- **`cell_deconvolution.R`**: Preprocessing for immune cell deconvolution — cell type splitting, correlation filtering, subgrouping logic used by multiple deconvolution methods (Quantiseq, MCPcounter, XCell, EpiDISH, DeconRNASeq).
- **`immunoscores.R`**: Computes 10+ published immune signature scores (CYT, Roh_IS, Davoli_IS, IFNγ, Ayers_expIS, T-cell inflamed, MSI, TLS, chemokines, ICB genes).

### Key dependencies

| Category | Packages |
|----------|----------|
| Differential expression | DESeq2 |
| Deconvolution | EpiDISH, DeconRNASeq, Quantiseq, MCPcounter, XCell |
| TF activity | decoupleR, viper, OmnipathR |
| Pathway analysis | GSVA, ReactomePA, clusterProfiler, fgsea, msigdbr |
| Co-expression | WGCNA |
| Survival | survival, rms, survminer |
| Visualization | ggplot2, ComplexHeatmap, pheatmap, ggstatsplot, PCAtools |
| Feature selection | Boruta |
