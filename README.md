# Transcriptomics Analysis (limma-based RNA-seq DE pipeline)

This repository contains an R pipeline for bulk RNA-seq analysis using the limma package.  
Focus: Differential gene expression associated with survival in breast cancer, with confounder adjustment.

## Features

- Preprocessing and matching of expression, survival, and clinical datasets
- Quality control (log transform, zero variance filtering)
- PCA and batch/higher-level covariate assessment
- Differential expression (simple and confounder-adjusted; DE gene export)
- Volcano plots, heatmaps, Venn diagrams, PCA visualizations
- Clean code structure and extensive comments

## Files

- `transcriptomics_analysis.R`: Main analysis script
- `README.md`: Project documentation

## Requirements

- R and RStudio
- Packages: `limma`, `pheatmap`, `ggplot2`, `VennDiagram`  
  (Install with: `install.packages(c("limma","pheatmap","ggplot2","VennDiagram"))`)

## Data

**No raw data is uploaded (privacy/policy).**  
- The script expects three input files:  
  - `exp`: expression matrix (genes as rows, samples as columns)
  - `survival`: survival/phenotype table (sample IDs, status)
  - `breast`: clinical annotation table (sample covariates)
- All data files should be in plain text, as described in script comments.

## Usage

1. Install the required R packages.
2. Place your compatible data files in the working directory named as above.
3. Open and run `transcriptomics_analysis.R` in R or RStudio.
4. Review output: DE gene tables, volcano plots, heatmaps, PCA and correlation/ANOVA summaries.

*For questions or to request a sample (fake) data structure, contact the author.*

## Contact

Yasser Abdelrahim — yasserabdelrahim6@gmail.com

