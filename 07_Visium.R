install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork", "Matrix", "remotes", "scatterpie"))
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)

library(Seurat)
library(spacexr)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(1234)   # 재현성을 위해