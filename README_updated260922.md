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

> **Note**: Unlike the practice code above (`01`–`06`), the analyses proposed in this section are written and designed by the repository author as original extensions of the pipeline, not adapted from the course/workshop.

기존 파이프라인(QC → clustering → annotation → CAF subclustering)에서 확장해볼 만한 두 가지 분석 방향을 정리했습니다. 우선순위 순서로 배치했습니다.

### 1. 싱글셀 → 공간전사체 Deconvolution 연결

> **Status (2026-09-22): ✅ Implemented.** RCTD 기반 파이프라인을 GSE279507 마우스 췌장암 전암병변 데이터에 적용했습니다. 상세 단계, 결과, 해석은 아래 **[Update: Visium Deconvolution Pipeline (RCTD)](#update-2026-09-22-visium-deconvolution-pipeline-rctd)** 섹션을 참고하세요.

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

---

## Update (2026-09-22): Visium Deconvolution Pipeline (RCTD)

> Source: *싱글셀 → Visium Deconvolution 단계별 가이드 (RCTD)*, 26p internal working doc (2026-09-21). This section condenses that guide into README form and records what was actually run, the decisions made along the way, and the resulting findings.

### Overview

싱글셀 reference(`Clusterwithcompletelabels.RDS`)와 Visium spot 발현 행렬을 [RCTD](https://github.com/dmcable/spacexr)(`spacexr`, R)에 넣어, spot마다 세포 타입 구성비를 추정하는 파이프라인입니다. GPU가 필요 없고, reference가 이미 Seurat 객체(R)이므로 R 기반 RCTD를 기본 도구로 선택했습니다(Python 기반 `cell2location`은 대안으로 Step 10에 정리).

```
Step 1  Visium 다운로드·로드
Step 2  Reference 점검·변환         ─┐
Step 3  SpatialRNA + RCTD 실행       │
Step 4  구성비 추출                  ├─→ Step 5  공간 시각화 → Step 6  검증·해석
Step 5  공간 시각화
Step 6  검증·해석
Step 7  5개 샘플 전체로 확장
Step 8  샘플별 리포트 그림 일괄 생성
Step 9  5개 샘플 비교·결론
Step 10 Acinar 거리 × 전체 세포 타입 + 통계 검정 + cell2location 대안
```

**시작 전 체크리스트**
1. **종이 같아야 함** — 마우스 reference ↔ 마우스 Visium (다르면 ortholog 변환 필요, 정확도 저하)
2. **조직이 같아야 함** — 췌장암 reference라면 췌장암(PDAC) Visium (다르면 reference에 없는 세포 타입이 엉뚱한 타입으로 배정됨)
3. **라벨이 충분히 세분화되어 있어야 함** — 세포 타입당 최소 25개 이상 (`CELL_MIN_INSTANCE = 25`)

### Data used: GSE279507

실제 사용 데이터는 **[GSE279507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279507)** — *p48Cre/+;LSL-KrasG12D/+* 마우스(6개월령) 췌장 5개 샘플에서 10x Visium(CytAssist)과 scRNA-seq를 함께 생산한 연구로, 췌장암 **전암병변**(acinar-to-ductal metaplasia, ADM → mucinous tubular complex, MTC) 진행 모델입니다. 종·조직 모두 본 pipeline의 reference와 일치합니다.

| 파일 | 내용 | 비고 |
|---|---|---|
| `GSE279507_RAW.tar` | Visium 원본 (mtx/tsv/json/png/tiff), 샘플별 미분리·flat 묶음 | 샘플별 `spatial/` 폴더 구조로 재구성 필요 |
| `GSE279507_adata_whole_scRNA.h5ad` | 전체 scRNA-seq (Python AnnData) | R 사용 시 변환 단계 필요 |
| `GSE279507_adata_acinar_pseudotime.h5ad` | acinar 세포만 모은 pseudotime 부분집합 | reference용 아님 (whole scRNA 사용) |

5개 샘플: `GSM8573641` – `GSM8573645`. GEO summary에 "Integration with the Mouse Cell Atlas (MCA) 2.0 single-cell dataset was performed"라는 문구가 있으나 통합 비율·처리 방식은 미공개(관련 논문 미출판)이므로, 본 파이프라인은 **자체 reference(`Clusterwithcompletelabels.RDS`)를 주 reference로, `adata_whole_scRNA.h5ad`를 비교용 보조 reference로** 사용하는 것으로 결정했습니다.

### Step 0 — R 환경 준비

`Seurat`, `spacexr`(GitHub: `dmcable/spacexr`), `Matrix`, `ggplot2`, `dplyr`, `patchwork`, `scatterpie` 설치. `set.seed(1234)`로 재현성 고정. Visium 한 장당 8–16GB 메모리면 충분.

### Step 1 — Visium 다운로드·로드

`GSE279507_RAW.tar`는 표준 10x 폴더 구조가 아니라 `GSM<ID>_<시료번호>_<표준파일명>.gz` 형식으로 평평하게 묶여 있어, 샘플별 `spatial/` 하위 폴더로 재구성하는 전처리가 필요합니다. deconvolution에는 샘플당 아래 7개 파일만 있으면 됩니다:

```
barcodes.tsv.gz · features.tsv.gz · matrix.mtx.gz · scalefactors_json.json.gz ·
tissue_hires_image.png.gz · tissue_lowres_image.png.gz · tissue_positions_list.csv.gz
```

`Read10X()`는 `barcodes.tsv.gz`/`matrix.mtx.gz`를 압축된 상태로 찾고, `Read10X_Image()`는 이미지·JSON·CSV를 압축 해제된 상태로만 읽으므로 두 종류를 다르게 처리해야 합니다(압축 유지 vs `gunzip`). h5 파일이 없으므로 `Load10X_Spatial()` 대신 `Read10X()` + `Read10X_Image()`로 직접 Seurat 객체를 조립합니다.

**기본 QC**: 이 데이터는 probe 기반 CytAssist(FFPE) 패널이라 **미토콘드리아 유전자가 패널에 아예 없습니다** (`mt_genes` 벡터가 빔). 따라서 `percent.mt` QC는 생략하고 `nCount_Spatial`/`nFeature_Spatial`만으로 QC를 진행하며, spot당 UMI가 거의 0인 spot(`nCount_Spatial > 100`)만 제거합니다. (probe 패널은 전체 유전체가 아니라 수천 개 수준의 큐레이션된 유전자만 담고 있어 Step 2 reference와 공통 유전자가 일부만 겹치지만, RCTD는 자동으로 공통 유전자만 사용하므로 실행 자체엔 문제없음.)

### Step 2 — Reference 점검과 RCTD Reference 생성

- Reference: `Clusterwithcompletelabels.RDS`, 라벨 컬럼 `celltypes`.
- **라벨 세트(14개, `04_Cell_Type_Annotation.R` 기준)**: `Acinar, CD4_T, NK, RegT, NK_T, B, Myleoid, Granulocyte, Fibroblast, Macrophage, Mastcell, Endothelial, T, Neutrophils`
- 정상 조직 atlas(TabulaMuris/MouseRNAseq) 기준 라벨이라 **"Tumor"/"Malignant" 라벨은 존재하지 않음** → 두 가지 명칭 보정 필요:
  - `CAF` → `Fibroblast` (라벨명 통일)
  - `Tumor`에 대응하는 라벨이 없음 → GSE279507이 전암병변(ADM→PanIN) 모델이라, 이 파이프라인에서는 **`prop_Acinar`를 `prop_Tumor`의 현실적 대안(proxy)**으로 사용 (단, "생물학적 종양"이 아니라 "acinar 기원 상피가 몰려있는 자리"를 보는 것임을 명시)
- 세포 타입당 25개 미만인 타입은 합치거나 제외 (`CELL_MIN_INSTANCE = 25`)
- 라벨에 `/`나 공백이 있으면 RCTD 오류 → `gsub("[/ +]", "_", ...)`로 정리
- RCTD는 **정규화되지 않은 raw count**를 요구 (`counts`, not `data`/`scale.data`)
- `Reference(counts, cell_types, nUMI)` 객체 생성 후 `results/RCTD_reference.rds`로 저장 (재사용용)

### Step 3 — SpatialRNA 객체 생성 + RCTD 실행

Visium을 RCTD 형식(`SpatialRNA` = 좌표 + count 행렬 + spot별 UMI)으로 변환 → `create.RCTD(puck, reference, max_cores = 4)` → `run.RCTD(myRCTD, doublet_mode = "full")`.

| `doublet_mode` | 의미 | 사용 여부 |
|---|---|---|
| `full` | spot당 임의 개수 타입 허용, 타입별 가중치 추정 | **사용** — 가중치 행렬이 바로 나와 시각화가 쉬움 |
| `multi` | spot당 최대 4개 타입으로 제한 | 잡음 줄이고 싶을 때 |
| `doublet` | spot당 1–2개 타입 | Slide-seq류 고해상도용, Visium엔 부적합 |

실행 시간: spot 2–4천 개 × 타입 10–20개 기준 수 분–수십 분.

### Step 4 — 결과 추출과 Seurat에 붙이기

`myRCTD@results$weights` → `normalize_weights()`로 행 합을 1로 맞춰 구성비(`prop_*` 컬럼)로 변환 → `AddMetaData()`로 Seurat 객체에 부착 → `results/spot_celltype_proportions.csv`, `results/visium_with_deconv.rds`로 저장.

검증 포인트: (1) 타입별 평균 비율이 상식과 부합하는지(이 데이터는 acinar/상피 계열과 fibroblast/면역세포가 우세한 패턴이 자연스러움), (2) 한 타입이 90% 이상을 차지하거나 어떤 타입이 0에 가까우면 reference 라벨 문제나 종 불일치를 의심.

### Step 5 — 공간 시각화

세 가지 플롯:
1. **관심 타입 색 그라데이션** — `SpatialFeaturePlot(features = c("prop_Fibroblast", "prop_Acinar", "prop_Macrophage", "prop_Myleoid"))`
2. **spot별 우세 타입 지도** — 각 spot에서 비율이 가장 높은 타입을 `dominant` 컬럼으로 뽑아 `SpatialDimPlot`
3. **파이 차트** — `scatterpie::geom_scatterpie()`로 spot마다 전체 구성비를 한 번에 표시 (spot이 많으면 상위 4–6개 타입만 남기고 나머지는 `Other`로 묶는 것을 권장; 조직 이미지 위에 얹지 않으므로 5-1이 논문 figure로는 더 적합)

### Step 6 — 결과 검증과 해석

**6-1. 마커 유전자로 교차확인.** 세포 타입별 잘 알려진 마우스 마커와 RCTD 비율의 spot 단위 Spearman 상관을 확인 (Fibroblast: `Col1a1, Dcn, Pdgfra` / Acinar: `Cpa1, Prss2, Cela1` / T세포: `Cd3` / 대식세포: `Adgre1, Cd68` / 혈관내피: `Pecam1`). 상관이 뚜렷하게 양수(대략 0.4 이상)면 정상 작동.

**6-2. 세포 타입 간 co-localization.** 구성비 행렬의 타입 간 상관행렬 → heatmap. 양의 상관 = 같은 spot에 자주 공존, 음의 상관 = 서로 다른 영역에 분리.

**6-3. Acinar 거리에 따른 Fibroblast/면역세포 분포.** `prop_Acinar > 0.5`인 spot을 "acinar-dominant"로 정의하고, 나머지 spot에서 최근접 acinar-dominant spot까지의 거리(`RANN::nn2`)를 계산 → 거리 대비 `prop_Fibroblast` 산점도 + LOESS. (Tumor 라벨이 없어 Acinar로 대체한 분석이며, 좌표를 실제 μm 단위로 바꾸려면 `scalefactors_json.json`의 `spot_diameter_fullres`를 사용.)

### Step 7–8 — 5개 샘플 전체로 확장·자동화

Step 2의 reference는 재사용하고, Step 1(QC)·Step 3(RCTD)·Step 4(비율 추출)만 `run_rctd_pipeline()` 함수로 묶어 5개 샘플(`GSM8573641`–`GSM8573645`)에 반복 적용. `make_sample_plots()`로 샘플마다 QC violin + 공간 시각화 3종 + 마커 검증 + 상관 히트맵 + Acinar 거리 플롯까지 총 7개 그림을 일괄 생성 (`results/` 폴더에 샘플당 7개, 총 35개 PNG).

### Step 9 — 5개 샘플 결과 종합 비교

각 샘플을 `marker_cor`(마커-비율 내적 타당성), `acinar_dist_cor`/`acinar_dist_slope`(Acinar 거리-Fibroblast 비율 상관/기울기)로 요약 → `results/all_samples_summary.csv`.

| sample | n_spots | marker_cor | acinar_dist_cor | acinar_dist_slope |
|---|---|---|---|---|
| GSM8573641 | 2,961 | 0.740 | −0.029 | −5.06e-05 |
| GSM8573642 | 2,672 | 0.570 | −0.106 | −1.16e-05 |
| GSM8573643 | 3,306 | 0.656 | +0.132 | +1.06e-05 |
| GSM8573644 | 2,419 | 0.579 | −0.024 | −7.07e-06 |
| GSM8573645 | 3,274 | 0.684 | +0.406 | +4.97e-05 |

**해석 기준**: `marker_cor`가 5개 샘플 모두 양수·일정 수준 이상이면 RCTD 추정치의 내적 타당성이 확보된 것(본 데이터는 0.57–0.74로 양호). `acinar_dist_cor`가 5개 샘플 모두 같은 방향이면 "Acinar에서 멀어질수록 Fibroblast/면역세포가 몰린다"는 생물학적으로 일관된 패턴으로 해석 가능하지만, 이 데이터에서는 **방향이 샘플마다 갈리고(음 3개, 양 2개) 절댓값도 대체로 작아, Fibroblast 하나만으로는 일관된 결론을 내리기 어렵다는 것이 Step 9의 1차 결론**이었습니다. 5개를 하나의 평균으로 억지로 묶지 않고, 개체 간 이질성 자체를 결과로 보는 방향으로 Step 10에서 확장했습니다.

### Step 10 — Acinar 거리 × 전체 세포 타입 + 통계 검정

Fibroblast 하나만 보던 Step 9를 reference의 **모든 세포 타입(Acinar 제외 13종)**으로 확장하고, 각 (샘플 × 타입) 상관에 `cor.test()`로 p-value를 붙여 5개 샘플 × 13개 타입 다중검정을 BH(FDR) 보정했습니다.

**⚠️ 통계적 주의사항**: 이 p-value는 spot들을 서로 독립이라 가정하지만, 실제로는 이웃 spot끼리 조성이 비슷한 **공간적 자기상관**(spatial autocorrelation)이 있어 유효 표본크기가 과대평가 → p-value가 지나치게 유의하게 나올 수 있습니다. 여기서는 "어느 타입이 상대적으로 더 강하고 일관된 신호를 보이는가"를 랭킹하는 용도로만 사용했고, 엄밀한 검정을 위해서는 공간 블록 permutation이나 Moran's I 기반 보정이 추가로 필요합니다.

**결과 요약 (신뢰도순)**

| 신뢰도 | 세포 타입 | 패턴 |
|---|---|---|
| 가장 강하고 일관됨 | **Myeloid** | 5/5 샘플 유의, 모두 양의 상관(+0.36~+0.66) — Acinar에서 멀어질수록 증가 |
| 가장 강하고 일관됨 | **Macrophage** | 5/5 샘플 유의, 모두 음의 상관(−0.15~−0.56) — Acinar 근처에 몰림 |
| 대체로 일관(예외 1개) | NK_T(+), Neutrophils(+), Granulocyte(+) | 4~5개 샘플 동일 방향 |
| 대체로 일관(예외 1개) | CD4_T(−), RegT(−), Endothelial(+) | 4개 샘플 동일 방향, **GSM8573642만 반대** (절편 위치·개체 특성 차이로 별도 점검 권장) |
| 유의하지만 신뢰도 낮음 | NK, Mastcell, T, B | 방향이 샘플마다 뒤섞임, 유의성이 spot 수(공간적 자기상관)에 의한 것일 가능성 |
| **Fibroblast** | **5개 중 3개만 유의, 부호도 갈림**(642 −0.11 / 643 +0.13 / 645 +0.41) | "desmoplasia(섬유화)가 병변 주변에 형성된다"는 가설이 이 데이터에서는 가장 약하고 비일관적인 축으로 확인됨 — 전암병변 단계라 문헌(대부분 침습암 PDAC/KPC 기준)만큼 뚜렷하지 않을 수 있음 |

**결론**: 이 파이프라인에서 재현성 있게 잡히는 신호는 **골수계 세포의 대비(Macrophage vs Myeloid)**와 GSM8573642를 제외한 4개 샘플의 T세포/혈관 패턴이며, Fibroblast를 포함한 나머지는 샘플 간 일관성이 낮아 추가 검증(예: 병리학적 병변 면적 데이터와의 교차 확인) 없이는 결론을 내리기 어렵습니다.

**타입별 문헌 근거 대조** (주의: 대부분 사람 침습성 췌장암(PDAC)이나 마우스 KPC 모델 대상 논문이며, GSE279507은 침습암이 아닌 전암병변(ADM/PanIN) 단계이므로 문헌보다 약하거나 아직 나타나지 않은 형태로 관찰될 수 있음):

| 세포 타입 | 문헌 기반 예상 패턴 | 근거 |
|---|---|---|
| Fibroblast (CAF) | Acinar/PanIN 병변 주변에 desmoplastic stroma 형성 (PanIN·IPMN 단계부터) | Masugi, *Cancers* 2022 |
| Macrophage / Myeloid | 대식세포 분비 사이토카인(RANTES, TNF)이 NF-κB 통해 ADM을 직접 유도, 억제성 골수계 세포가 병변 초기부터 우세 | Liou et al., *J Cell Biol* 2013; Clark et al., *Cancer Res* 2007 |
| T세포 (CD4_T/NK_T/RegT/T) | FOXP3+ Treg는 PanIN→침습암 진행에 따라 증가, 세포독성 CD8+ T세포는 감소 | Hiraoka et al., *Clin Cancer Res* 2006 |
| Endothelial (혈관내피) | 종양 조직은 미세혈관 밀도(MVD)가 높지만 분화도와 단순 비례하진 않음 | Bărău et al., *Virchows Arch* 2013 |

### 자주 겪는 문제

| 증상 | 원인 | 해결 |
|---|---|---|
| `Reference` 생성 시 라벨 관련 오류 | 라벨에 `/` 또는 공백, 혹은 factor의 names가 없음 | Step 2-2처럼 `gsub` 하고 `names(cell_types) <- colnames(ref)` |
| `create.RCTD`에서 공통 유전자가 거의 없다는 경고 | 종 불일치, 유전자 ID(Ensembl vs symbol) 불일치 | `head(rownames())`로 두 객체의 형식을 맞춤 |
| 한 세포 타입이 모든 spot을 덮음 | reference에 없는 세포가 있거나 라벨이 너무 성기게 나뉘어 있음 | 라벨을 합치거나, 빠진 세포 타입 추가 |
| `Load10X_Spatial`이 좌표 파일을 못 찾음 | GEO에서 낱개로 받은 파일명이 표준과 다름 | `spatial/` 폴더 안에 표준 파일명으로 재구성 (Step 1) |
| 실행이 너무 오래 걸림 | reference 세포가 너무 많음 | 타입당 수백~수천 개로 다운샘플링, `max_cores` 증가 |

분석 품질을 올리려면 reference가 대상 Visium 샘플과 같은 종류의 조직이어야 가장 좋습니다(서로 다른 시료·실험실 유래 싱글셀은 배치 효과를 만듦). RCTD는 이를 어느 정도 보정하지만 완전히 없애지는 못하므로, 발표 그림은 되도록 같은 모델(예: KPC 자료)로 맞추는 것을 권장합니다.

### RCTD vs cell2location — 대안 검토

RCTD 결과의 검증용으로 두 번째 방법을 돌리거나, 파이프라인 전체를 Python으로 통일하고 싶을 때는 `cell2location`을 씁니다: (1) Seurat reference를 `SeuratDisk::SaveH5Seurat()` + `Convert()`로 `.h5ad` 변환, (2) `cell2location.models.RegressionModel`로 세포 타입별 발현 시그니처 학습(250+ epoch, GPU 권장), (3) `cell2location.models.Cell2location`에 Visium AnnData + 시그니처를 넣고 조직당 예상 세포 수(`N_cells_per_location`)와 스팟 변동성(`detection_alpha`) 설정해 학습, (4) 결과는 세포 타입별 **예상 세포 수(절대량)**로 나오며 행별로 정규화하면 비율이 됨.

| 항목 | RCTD (R) | cell2location (Python) |
|---|---|---|
| 설치·실행 난이도 | 낮음 | 높음 (환경 설정, 하이퍼파라미터) |
| GPU | 불필요 | 권장 |
| 출력 | 세포 타입 비율 | 타입별 예상 세포 수(절대량) |
| 대표적 장점 | 빠르고 Seurat과 연결이 쉬움 | 데이터가 많을 때 안정적이고 미세 영역까지 잘 잡는다고 보고됨 |

처음에는 RCTD로 한 바퀴 끝까지 돌리고, 논문이나 발표에서 결과의 강건함을 보여줘야 할 때 cell2location을 검증용 두 번째 방법으로 추가하는 순서를 권장합니다.

### Output files (`results/`)

```
RCTD_reference.rds                              # Step 2, reference 객체 (재사용 가능)
myRCTD.rds                                       # Step 3, 단일 샘플 RCTD 원본 결과
spot_celltype_proportions.csv                    # Step 4, spot × 세포타입 비율
visium_with_deconv.rds                           # Step 4, 비율이 붙은 Seurat 객체
spatial_proportions.png / spatial_pie.png / ...  # Step 5, 공간 시각화
celltype_correlation_heatmap.png                 # Step 6-2
<GSM ID>_*.png  (× 5 샘플 × 7종)                  # Step 8, 샘플별 리포트 그림 (35개)
all_samples_summary.csv                          # Step 9, 샘플별 marker_cor/acinar_dist_cor 요약
all_samples_composition.png / all_samples_acinar_cor.png / all_samples_acinar_distance_facet.png
all_samples_celltype_distance.csv                # Step 10, 타입 × 샘플 상관·기울기
all_samples_celltype_distance_heatmap.png
all_samples_celltype_distance_pvalues.csv        # Step 10, BH 보정 p-value 포함
```
