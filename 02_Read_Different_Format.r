# install.packages("rlang")
# packageVersion("rlang")

library(Seurat)

getwd()

### Read single cell matrix

# This section contents scripts to read single cell matrix from different formats and to convert it into Seurat object for further analysis. 

# 1. 10X Cell-Ranger .HDF5 format - Cell ranger outputs the feature barcode matrix in .hdf5 format

data <- Read10X_h5("single_cell/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
seurat.obj <- CreateSeuratObject(data)

rds_obj <- readRDS('single_cell/sc_embryo.rds')
seurat.obj1 <- CreateSeuratObject(rds_obj)

mtx_obj <- ReadMtx(mtx = "single_cell/raw_feature_bc_matrix/matrix.mtx.gz",
        features = "single_cell/raw_feature_bc_matrix/features.tsv.gz",
        cells = "single_cell/raw_feature_bc_matrix/barcodes.tsv.gz")
seurat_mtx <- CreateSeuratObject(counts = mtx_obj)