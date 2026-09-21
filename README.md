# Single-Cell RNA-seq Analysis Pipeline (Seurat-based)

> **Note**: This repository contains practice code written while following an external single-cell RNA-seq analysis course/workshop. It is shared here as a personal learning record, not as original methodology. Course-specific materials (slides, etc.) are not included.

A step-by-step Seurat (v5) pipeline for single-cell RNA-seq analysis, from raw data loading to cell type annotation, subclustering, and visualization. Developed and tested on a mouse pancreas dataset (10x Genomics, Cell Ranger output).

## Pipeline Overview

```
01_Required_Packages.R
        │
02_Read_Different_Format.R
        │
03_Pre-processing_and_Clustering.R
        │
04_Cell_Type_Annotation.R
        │
        ├── 05_Customized_sc_RNA_data_analysis.R  (subclustering branch)
        │
06_Visualization_basedOn_Seurat_Pipeline.R
```

## Scripts

### `01_Required_Packages.R`
Installs and loads all required R/Bioconductor packages via `pacman`:
`Seurat`, `Matrix`, `ggplot2`, `dplyr`, `SingleCellExperiment`, `scran`, `scuttle`, `celldex`, `SingleR`.

### `02_Read_Different_Format.R`
Demonstrates loading single-cell count matrices into a Seurat object from three common formats:
- 10x Cell Ranger `.h5` (`Read10X_h5`)
- Pre-saved `.rds` matrix
- MTX/TSV triplet (`ReadMtx`: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)

### `03_Pre-processing_and_Clustering.R`
Core QC and clustering workflow:
1. Load 10x `.h5` data, create Seurat object
2. Compute mitochondrial percentage (`^mt-` pattern; mouse gene naming)
3. QC visualization (violin/density plots) before and after filtering
4. Filter cells (`nFeature_RNA` 250–2500, `percent.mt` < 5, `nCount_RNA` 500–20000)
5. Normalize (`LogNormalize`), find variable features (top 2000, `vst`)
6. Scale data, run PCA, elbow plot for dimensionality
7. `FindNeighbors` (dims 1:35), `FindClusters` (resolution 0.6)
8. `RunUMAP` and `RunTSNE` for visualization
9. Saves → `single_cell/Seurat_Clustering.RDS`

### `04_Cell_Type_Annotation.R`
Three complementary annotation strategies:
1. **Reference-based (TabulaMuris, SmartSeq2/EH1618)** via `SingleR`, using `cell_ontology_class` labels
2. **Reference-based (MouseRNAseq, `celldex`)** with broad + fine label prediction, combined via `combineCommonResults`
3. **DEG-based manual annotation**: `FindAllMarkers` for top marker genes per cluster; ambiguous clusters (e.g., cluster 16, which SingleR could not confidently assign) are resolved manually by inspecting canonical marker gene expression (e.g., `Calca`) via `VlnPlot`/`FeaturePlot`
4. Final cluster identities are assigned with `RenameIdents` and stored in metadata (`celltypes` column)
5. Saves → `single_cell/Session 2/Clusterwithcompletelabels.RDS`

### `05_Customized_sc_RNA_data_analysis.R`
Focused subclustering analysis on fibroblast-like clusters (identified as clusters 10 and 17 from the annotated object) to resolve cancer-associated fibroblast (CAF) subtypes:
- Subsets and re-runs the standard pipeline (normalize → HVG → scale → PCA → cluster → UMAP) on the fibroblast subset
- Screens for three CAF subtypes using curated marker gene panels:
  - **i-CAF** (inflammatory): `Vim, Fap, Col3a1, Il6, Cxcl2, ...`
  - **my-CAF** (myofibroblastic): `Dcn, Postn, Tpm1/2, Thbs2, Acta2, ...`
  - **ap-CAF** (antigen-presenting): `H2-Ab1, H2-Aa, Cd74, Cd83, Krt8/18`
- Visualizes marker expression (violin, feature, dot plots, heatmap via `dittoSeq`)
- Subsets cells by marker gene expression (e.g., `Acta2 > 0`), merges subsets, and runs `FindAllMarkers` (with `JoinLayers` for Seurat v5 layer compatibility) to compare marker-defined groups
- Saves → `single_cell/Session 3/Session 3/Fibroblast_Clustering.RDS`

### `06_Visualization_basedOn_Seurat_Pipeline.R`
A reference collection of Seurat's built-in visualization functions applied to the annotated and subclustered objects:
- `RidgePlot`, `VlnPlot`, `FeaturePlot` (including quantile-based contrast scaling and two-gene co-expression blending), `DotPlot`, `DoHeatmap`, `dittoHeatmap`, `DimPlot`/`TSNEPlot`
- Custom plot theming (`ggmin::theme_powerpoint()`, `DarkTheme()`)
- Example QC metric scatter plot colored by mitochondrial percentage

## Requirements

- R (Seurat v5 — note: `JoinLayers()` is required after `merge()` before running `FindAllMarkers()`)
- Packages listed in `01_Required_Packages.R`
- `dittoSeq` (Bioconductor) for heatmap visualization
- `ggmin` (GitHub: `sjessa/ggmin`) for optional plot theming

## Data

Raw input data (10x Cell Ranger output: `filtered_feature_bc_matrix.h5`, `raw_feature_bc_matrix/`) is **not included** in this repository. Scripts assume a working directory containing a `single_cell/` folder with the expected input files.

## Notes

- All file paths are relative to a `single_cell/` working directory and contain session-specific subfolders (`Session 2/`, `Session 3/Session 3/`) — update paths to match your own directory structure before running.
- Mitochondrial gene pattern (`^mt-`) and marker gene panels are specific to **mouse** data; adjust for human data (`^MT-`) accordingly.
