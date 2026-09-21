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

## Suggested Next Steps (추천 후속 분석)

기존 파이프라인(QC → clustering → annotation → CAF subclustering)에서 확장해볼 만한 세 가지 분석 방향을 정리했습니다. 우선순위 순서로 배치했습니다.

### 1. 싱글셀 → 공간전사체 Deconvolution 연결

**목적**
지금까지 만든 싱글셀 cell type reference를 이용해서, 공개된 공간전사체(spatial transcriptomics) 데이터의 각 spot에 어떤 세포 타입이 얼마나 섞여 있는지 추정하는 분석입니다. 공간전사체 데이터(예: 10x Visium)는 spot 하나에 여러 세포가 섞여 들어가는 경우가 많아서(single-cell 해상도가 아님), 싱글셀 데이터를 reference로 활용해 각 spot의 세포 구성비를 통계적으로 추정하는 **deconvolution** 과정이 필요합니다.

**방법 / 도구**
1. 공개 Visium 데이터셋 다운로드 (예: 10x Genomics의 유방암 조직 데모 데이터)
2. 본 파이프라인의 `Clusterwithcompletelabels.RDS`(cell type 라벨이 붙은 최종 Seurat object)를 reference로 사용
3. `cell2location`(Python, Bayesian 기반) 또는 `RCTD`(R, `spacexr` 패키지)로 deconvolution 수행
4. 각 spot에 대해 세포 타입 구성비를 추정하고, 공간 좌표 위에 파이 차트 또는 색상 그라데이션으로 시각화

**기대 효과**
- 종양미세환경(TME) 내에서 특정 세포 타입(예: CAF, 면역세포)이 공간적으로 어디에 몰려 있는지 확인 가능
- 싱글셀 분석과 공간 데이터 분석을 연결하는 실전 경험 확보

**English summary**
Use the annotated single-cell reference from this pipeline to deconvolve public Visium spatial transcriptomics data (via `cell2location` or `RCTD`), estimating the cell-type composition of each spot and visualizing spatial distribution of specific populations (e.g., CAFs, immune cells) within the tumor microenvironment.

---

### 2. Cell-Cell Communication 분석

**목적**
어노테이션이 완료된 Seurat object를 대상으로, 세포 타입 간에 어떤 리간드-리셉터(ligand-receptor) 상호작용이 일어나는지 추론하는 분석입니다. 특히 CAF 아형(i-CAF/my-CAF/ap-CAF)별로 종양세포·면역세포와 다른 신호 경로를 쓰는지 비교하면, 단순 marker gene 나열을 넘어선 기능적 해석이 가능합니다.

**방법 / 도구**
1. `CellChat` 또는 `CellPhoneDB` 패키지를 annotation이 완료된 `seurat.obj`(및 필요시 fibroblast subcluster 결과)에 적용
2. 세포 타입 간 상호작용 강도를 계산 (`computeCommunProb` 등)
3. CAF 아형별로 어떤 신호 경로(signaling pathway)가 유의하게 활성화되는지 비교
4. 결과를 circle plot, heatmap, bubble plot 등으로 시각화

**기대 효과**
- 단순 세포 구성 파악을 넘어, 세포 간 상호작용이라는 기능적 레이어 추가
- 종양미세환경 연구에서 자주 요구되는 분석 유형이라 실전 활용도가 높음

**English summary**
Apply `CellChat` or `CellPhoneDB` to the annotated Seurat object to infer ligand-receptor interactions between cell types, comparing signaling pathways used by different CAF subtypes (i-CAF/my-CAF/ap-CAF) with tumor and immune cells — moving beyond marker gene lists to a functional interpretation of the tumor microenvironment.
