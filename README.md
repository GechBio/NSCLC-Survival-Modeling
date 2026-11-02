# Reproducibility Package for: Comparative Survival Modeling of Clinical, Methylation, and Expression Predictors in NSCLC

## Overview
This repository contains the complete R analysis pipeline used to generate all results, figures, and tables for the manuscript. The pipeline evaluates six survival modeling strategies across expression and methylation modalities in TCGA NSCLC data.

## System Requirements
* R version 4.5.2 (or compatible)
* RStudio (recommended)

## Environment Setup
This project uses `renv` to ensure exact package versions. To restore the environment:
1. Open R in this directory.
2. Install renv if you haven't already: `install.packages("renv")`
3. Restore the locked environment: `renv::restore()`

## How to Run the Pipeline
To execute the entire analysis from start to finish:
```R
source("run_all.R")
```
This will sequentially execute Scripts 00 through 08, generating all tables (`.csv`) and figures (`.pdf`) in the working directory.

## Data Sources
Raw RNA-seq and 450K methylation data were obtained from the TCGA GDC Portal (https://portal.gdc.cancer.gov/).
Due to size constraints, raw data files are not included in this archive.
Script `00_base_setup.R` assumes the raw data are present in the specified local directories.

## Contact
Corresponding Author: Yuan-Chii Gladys Lee, PhD (ycgl@tmu.edu.tw)
