options("instal.lock" = FALSE)

### 1. Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install()

### 2. Install required packages
install.packages("pacman")
library(pacman)
pacman::p_load(Seurat, Matrix, ggplot2, dplyr, SingleCellExperiment, scran, scuttle, celldex, SingleR)
