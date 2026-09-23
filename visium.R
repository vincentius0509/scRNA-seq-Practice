getwd()

install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork", "Matrix", "remotes", "scatterpie"))
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)

options(download.file.method = "libcurl", timeout = 600)
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)

library(Seurat)
library(spacexr)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(1234)   # 재현성을 위해

untar("single_cell/visium_deconv/GSE279507_RAW.tar", list = TRUE, tar = "internal")

dir.create("single_cell/visium_deconv/data", recursive = TRUE, showWarnings = FALSE)
dir.create("single_cell/visium_deconv/results", recursive = TRUE, showWarnings = FALSE)

setwd("single_cell/visium_deconv")
getwd()

## Step 1. Visium 데이터 다운로드와 로드

# 진행 상황 (2026-09-22). 실제로 사용 중인 데이터는 GSE279507입니다 — p48Cre/+-LSL-KrasG12D/+ 마우스(6개월령) 취장 5개 샘플에서 10x Visium과 scRNA-seq를 함께 생산한 연구
# GSE279507_RAW.tar: Visium 원본(mtx, tsv, json, png/jpg/tiff)이 샘플별로 섞이지 않고 하나의 tar에 묶임
# GSE279507_adata_whole_scRNA.h5ad: 전체 scRNA-seq, Python AnnData(h5ad) 형식
# GSE279507_adata_acinar_pseudotime.h5ad: acinar 세포만 모은 pseudotime 부분집합

# barcodes.tsv.gz · features.tsv.gz · matrix.mtx.gz · scalefactors_json.json.gz · tissue_hires_image.png.gz · tissue_lowres_image.png.gz · tissue_positions_list.csv.gz
# 압축 해제와 샘플별 폴더 재구성. 접두어를 때고 spatial/ 하위 폴더로 나누는 작업입니다. 이미지와 JSON은 Read10X_Image()가 압축 푸린 파일만 읽을 수 있어서, 발현 행렬을 포함해 모두 압축을 해제

tar_file    <- "GSE279507_RAW.tar"
extract_dir <- "GSE279507_extracted"
sample_root <- "visium_samples"

dir.create(extract_dir, showWarnings = FALSE)
untar(tar_file, exdir = extract_dir, tar = "internal")

files <- list.files(extract_dir)
gsm_ids <- unique(sub("^(GSM[0-9]+)_.*", "\\1", files))
gsm_ids
#> [1] "GSM8573641" "GSM8573642" "GSM8573643" "GSM8573644" "GSM8573645"

# Read10X()는 CellRanger v3+ 형식을 인식할 때 barcodes.tsv.gz / matrix.mtx.gz를
# 압축된 상태 그대로 찾습니다 — 압축을 풀면 오히려 못 읽습니다.
matrix_trio <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
# 반대로 Read10X_Image()는 압축 파일을 못 읽으므로 이쪽만 푸니다.
spatial_gz  <- c("scalefactors_json.json.gz", "tissue_hires_image.png.gz",
                 "tissue_lowres_image.png.gz", "tissue_positions_list.csv.gz")
needed <- c(matrix_trio, spatial_gz)

# 외부 패키지 없이 base R로 gz 해제 (spatial 파일에만 쓸)
gunzip_file <- function(src, dst) {
  con_in  <- gzfile(src, "rb"); on.exit(close(con_in), add = TRUE)
  con_out <- file(dst, "wb");   on.exit(close(con_out), add = TRUE)
  repeat {
    buf <- readBin(con_in, what = "raw", n = 5e7)
    if (length(buf) == 0) break
    writeBin(buf, con_out)
  }
}

dir.create(sample_root, showWarnings = FALSE)

for (gsm in gsm_ids) {
  s_dir  <- file.path(sample_root, gsm)
  sp_dir <- file.path(s_dir, "spatial")
  dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)

  gsm_files <- files[startsWith(files, paste0(gsm, "_"))]
  gsm_files <- gsm_files[sub("^GSM[0-9]+_[^_]+_", "", gsm_files) %in% needed]

  for (f in gsm_files) {
    suffix <- sub("^GSM[0-9]+_[^_]+_", "", f)   # GSM<ID>_<시료번호>_ 제거 (.gz 포함)

    if (suffix %in% matrix_trio) {
      # 압축 유지한 채 이름만 바끈서 복사
      file.copy(file.path(extract_dir, f), file.path(s_dir, suffix), overwrite = TRUE)
    } else {
      # spatial 파일은 압축 해제
      base <- sub("\\.gz$", "", suffix)
      gunzip_file(file.path(extract_dir, f), file.path(sp_dir, base))
    }
  }
  message(gsm, " -> ", s_dir)
}

list.files(sample_root, recursive = TRUE)

# Seurat으로 로드하기. h5 파일이 없으므로 Load10X_Spatial() 대신 Read10X() + Read10X_Image()로 직접 조립합니다.
sample_dirs <- list.dirs(sample_root, recursive = FALSE)
names(sample_dirs) <- basename(sample_dirs)

load_visium_sample <- function(d) {
  counts <- Read10X(data.dir = d)                         # barcodes/features/matrix
  obj <- CreateSeuratObject(counts = counts, assay = "Spatial", project = basename(d))

  image <- Read10X_Image(image.dir = file.path(d, "spatial"))
  image <- image[Cells(obj)]      # 바코드 집합을 세포 객체에 맞춤
  DefaultAssay(image) <- "Spatial"
  obj[["slice1"]] <- image
  obj
}

visium_list <- lapply(sample_dirs, load_visium_sample)

# 처음에는 샘플 하나로 진행 (5개 모두 돌리려면 이 이름을 for/lapply로 감싸면 됨)
visium <- visium_list[[1]]
visium
#> An object of class Seurat
#> ... features across ... samples

dim(visium)
head(rownames(visium))   # 마우스라 Cd8a처럼 첫 글자만 대문자면 정상
View(rownames(visium))
# 기본 QC. 조직 밖의 빈 spot은 이미 제외된 상태이므로, 사이즈와 미토콘드리아 비율만 가볼게 봅니다.
# 먼저 실제 미토콘드리아 유전자가 패널에 있는지 확인
mt_genes <- grep("mt-", rownames(visium), ignore.case = TRUE, value = TRUE)
mt_genes
length(mt_genes) 

has_mt <- length(mt_genes) > 0
if (has_mt) {
  mt_pat <- if (any(grepl("^mt-", rownames(visium)))) "^mt-" else "^MT-"
  visium[["percent.mt"]] <- PercentageFeatureSet(visium, pattern = mt_pat)
} else {
  message("이 패널에는 미토콘드리아 유전자가 없습니다 (probe 기반 CytAssist 데이터에서 흔함). percent.mt 컬럼은 만들지 않고 nCount/nFeature만으로 QC합니다.")
}
# probe 기반(CytAssist FFPE) 이라 패널에 미토콘드리아 유전자가 아예 포함되지 않았을

# percent.mt가 없으면 VlnPlot에서도 빼야 함 (전부 NA면 VlnPlot이 에러남)
qc_features <- c("nCount_Spatial", "nFeature_Spatial")
if (has_mt) qc_features <- c(qc_features, "percent.mt")

VlnPlot(visium, features = qc_features, pt.size = 0.1)
SpatialFeaturePlot(visium, features = "nCount_Spatial")

hist(visium$nCount_Spatial)
sum(visium$nCount_Spatial <= 100)  # 10

visium <- subset(visium, subset = nCount_Spatial > 100)
nrow(visium)  # 19465


## Step 2. Reference 점검과 RCTD Reference 만들기
# Clusterwithcompletelabels.RDS(single cell 수행 결과)를 주 reference
# adata_whole_scRNA.h5ad는 (후속 h5ad → Seurat 변환 후) 비교용 보조 reference

# 2-1. Reference 불러오기와 라벨 확인
getwd()
ref <- readRDS("../../single_cell/Session 2/Clusterwithcompletelabels.RDS")

colnames(ref@meta.data)
table(Idents(ref))

label_col <- "celltypes"            # 04_Cell_Type_Annotation.R의 AddMetaData(col.name = "celltypes")와 일치
table(ref@meta.data[[label_col]])

# 종 확인: 마우스면 Cd8a처럼 첫 글자만 대문자 (TabulaMuris/MouseRNAseq 기반이라 마우스가 맞음)
head(rownames(ref))
table(ref$celltypes)

# 2-2. 라벨 정리
ref$ct <- as.character(ref@meta.data[[label_col]])
ref$ct <- gsub("[/ +]", "_", ref$ct)             # CD4+T, NK+T의 '+'도 함께 치환

keep <- names(which(table(ref$ct) >= 25))         # 25개 이상만 유지
ref  <- subset(ref, subset = ct %in% keep)
table(ref$ct)

# 2-3. RCTD Reference 객체 생성
counts <- GetAssayData(ref, assay = "RNA", layer = "counts")   # Seurat 4면 slot = "counts"
cell_types <- factor(ref$ct)
names(cell_types) <- colnames(ref)
nUMI <- colSums(counts)

reference <- Reference(counts, cell_types, nUMI)

# 저장해 두면 다음에 재사용 가능
saveRDS(reference, "results/RCTD_reference.rds")

## Step 3. SpatialRNA 객체 만들고 RCTD 실행 (GPU 사용이 어려워 cell2location 말고 RCTD를 사용)
# Visium을 RCTD가 이해하는 SpatialRNA 형식(좌표 + count 행렬 + spot별 UMI)으로 바꿉니다. 좌표의 행 이름과 count 행렬의 열 이름이 같은 barcode여야 합니다.
# raw count 행렬, 세포별 타입 라벨, 세포별 총 UMI. Reference() 함수를 부르는 건 이 세 가지만 Seurat 객체에서 뽑아내서 RCTD가 이해하는 형식으로 "번역"해주는 작업
# RCTD에게 "각 세포 타입이 순수하게 있을 때 유전자 발현이 어떤 모양인지"를 가르쳐주는 교과서 역할
# T세포는 Cd3e를 많이 쓴다는 식의 유전자별 기댓값)을 추정합니다. 그다음 Visium의 각 spot을 볼 때, "이 spot에서 관찰된 발현이 Acinar 프로파일 x%, Fibroblast 프로파일 y%, T세포 프로파일 z%를 섞으면 가장 잘 재현되는가"를 통계적으로 풀어서 구성비를 역산
# visium spot data의 발현을 설명해보려함
counts_sp <- GetAssayData(visium, assay = "Spatial", layer = "counts")

coords <- GetTissueCoordinates(visium)[, 1:2]      # 앞 2열 = x, y (Seurat 4는 imagerow, imagecol)
colnames(coords) <- c("x", "y")
coords <- coords[colnames(counts_sp), ]            # 순서 맞추기
stopifnot(identical(rownames(coords), colnames(counts_sp)))

nUMI_sp <- colSums(counts_sp)
puck <- SpatialRNA(coords, counts_sp, nUMI_sp)

# RCTD를 생성하고 실행
myRCTD <- create.RCTD(puck, reference, max_cores = 4)   # Windows는 max_cores = 1
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

# full: spot에 임의의 세포 타입 수를 허용하고 타입별 가중치를 추정 (시각화가 쉬움)
# multi: spot당 최대 4개 타입으로 제한
# doublet: spot당 1~2개 타입

saveRDS(myRCTD, "results/myRCTD.rds")

## Step 4. 결과 추출과 Seurat에 붙이기
# full 모드의 결과는 spot × 세포타입 가중치 행렬입니다. 행별로 합을 1로 정규화하면 구성비(비율)가 됩니다.
weights      <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)          # 행 합 = 1

prop <- as.data.frame(as.matrix(norm_weights))
colnames(prop) <- paste0("prop_", make.names(colnames(prop)))

# 검사: 모든 spot의 합이 1이고, 타입별 평균 비율이 상식적인지
summary(rowSums(prop))
round(sort(colMeans(prop), decreasing = TRUE), 3)
#  prop_Fibroblast     prop_Myleoid      prop_Acinar           prop_B 
#            0.543            0.190            0.096            0.052 
# prop_Neutrophils        prop_NK_T  prop_Macrophage prop_Endothelial 
#            0.028            0.025            0.014            0.013 
#          prop_NK           prop_T        prop_RegT       prop_CD4_T 
#            0.009            0.009            0.007            0.007 
#    prop_Mastcell prop_Granulocyte 
#            0.004            0.004 

# Seurat 객체에 연결 (UMI 조건으로 RCTD가 제외한 spot은 NA)
visium <- AddMetaData(visium, metadata = prop)
sum(is.na(visium$prop_Fibroblast))                  # 이 reference는 CAF 대신 Fibroblast 라벨

write.csv(prop, "results/spot_celltype_proportions.csv")
saveRDS(visium, "results/visium_with_deconv.rds")

## Step 5. 공간 시각화
# 색 그라데이션으로 보는 지도, 여러 타입을 바로 비교하는 지도, 그리고 spot마다 전체 구성을 보여주는 파이 차트
# 5-1.색 그라데이션으로 (CAF, 면역세포 등)
types_of_interest <- c("prop_Fibroblast", "prop_Acinar", "prop_Macrophage", "prop_Myleoid")  # 이 reference의 실제 라벨

p1 <- SpatialFeaturePlot(
  visium, features = types_of_interest,
  pt.size.factor = 1.6, alpha = c(0.1, 1), ncol = 2
)
ggsave("results/spatial_proportions.png", p1, width = 10, height = 9, dpi = 300)

# 5-2. spot별 우세 세포타입 지도
prop_cols <- grep("^prop_", colnames(visium@meta.data), value = TRUE)
m <- as.matrix(visium@meta.data[, prop_cols]); m[is.na(m)] <- 0
visium$dominant <- sub("^prop_", "", prop_cols[max.col(m, ties.method = "first")])

p2 <- SpatialDimPlot(visium, group.by = "dominant", pt.size.factor = 1.6)
ggsave("results/spatial_dominant_celltype.png", p2, width = 8, height = 7, dpi = 300)

# 5-3. 파이 차트 (spot마다 전체 구성비를 한 번에)
# install.packages("RANN")   # 최초 1회
library(RANN)

common <- intersect(rownames(coords), rownames(prop))
pie_df <- cbind(coords[common, ], prop[common, ])
type_cols <- colnames(prop)

# 파이 반지름 = 이웃 spot 간 최소 거리의 절반 정도
nn <- RANN::nn2(pie_df[, c("x", "y")], k = 2)$nn.dists[, 2]
pie_df$r <- 0.45 * median(nn)

p3 <- ggplot() +
  scatterpie::geom_scatterpie(data = pie_df, aes(x = x, y = y, r = r), cols = type_cols,
                              pie_scale = 1, color = NA) +
  scale_y_reverse() + coord_fixed() + theme_void() +
  labs(fill = "Cell type")

ggsave("results/spatial_pie.png", p3, width = 9, height = 8, dpi = 300)

# Spot 위치 
# Visium 슬라이드는 10x Genomics가 공장에서 미리 찍어놓은 고정된 격자입니다. 
# 슬라이드 한 장에 약 5,000개의 spot이 육각형(벌집) 패턴으로 정확히 55 µm 지름, spot 간 100 µm 간격으로 인쇄돼 있고, 이 좌표는 Visium 버전마다 항상 동일합니다. 
# 조직을 그 위에 올려놓으면 조직이 이 고정 격자의 일부만 덮게 되는데, Space Ranger가 조직 이미지와 슬라이드 모서리의 정렬 마커(fiducial)를 대조
# "어느 고정 spot이 실제 조직 아래에 있었는지"를 계산해 tissue_positions 파일에 기록

# 14개 타입 전부에 대해 "이 타입이 몇 %일 것이다"라는 비율을 동시에 내놓고, 그 14개 숫자를 다 더하면 1이 되는 식입니다. 
# 그래서 "이 spot은 Fibroblast다"가 아니라 "이 spot은 Fibroblast 62%, Myleoid 15%, Acinar 8%, ... 식으로 섞여있다"는 방식으로 계산
# 수천 개 유전자의 발현 패턴 전체가 어느 조합에서 제일 잘 맞아떨어지는지를 통계적으로 계산

## Step 6. 결과 검증과 해석
# 두 가지를 확인합니다. 추정한 비율 지도가 해당 세포의 마커 유전자 발현 지도와 겹치는지, 그리고 조직 구조가 상식과 맞는지
# 6-1. 마커 유전자로 교차확인 (각 세포별 주요 마커 유전자)
# Fibroblast (섬유아세포, 이 reference의 CAF 역할): Col1a1, Dcn, Pdgfra
# Acinar (이 모델의 전암병변 상피가 가장 가까이 잡힐 타입): Cpa1, Prss2, Cela1
# T세포 (CD4_T/RegT/NK_T/T 공통): Cd3e
# 대식세포: Adgre1, Cd68
# 혈관내피: Pecam1

visium <- NormalizeData(visium, assay = "Spatial")     # 또는 SCTransform(visium, assay = "Spatial")

p_marker <- SpatialFeaturePlot(visium, features = c("Col1a1", "Dcn", "prop_Fibroblast"),
                               pt.size.factor = 1.6, ncol = 3)

# 수치로도 확인: 비율과 마커 발현의 spot 단위 상관
expr <- FetchData(visium, vars = c("Col1a1", "Dcn", "prop_Fibroblast"))
cor(expr, use = "complete.obs", method = "spearman")
#                    Col1a1       Dcn prop_Fibroblast
# Col1a1          1.0000000 0.6361623       0.6890379
# Dcn             0.6361623 1.0000000       0.6234062
# prop_Fibroblast 0.6890379 0.6234062       1.0000000

# 6-2. 세포 타입 간 co-localization
# spot 단위로 비율을 서로 비교하면 어떤 세포끼리 같은 자리에 모이는지
prop_mat <- as.matrix(prop)
cor_mat  <- cor(prop_mat, method = "spearman")

heatmap(cor_mat, symm = TRUE, margins = c(12, 12))   # 또는 pheatmap::pheatmap(cor_mat)

png("results/celltype_correlation_heatmap.png", width = 8, height = 7, units = "in", res = 300)
heatmap(cor_mat, symm = TRUE, margins = c(12, 12))
dev.off()

# 6-3. Acinar 영역으로부터의 거리에 따른 Fibroblast/면역세포 분포
# 데이터 자체가 acinar 세포에서 시작된 전암병변 모델이므로, Acinar 영역에서 멀어질수록 Fibroblast/면역세포가 늘어나는지 확인
# 따로 Tumor 세포에 대한 정보가 없어서 Acinar 영역을 기준으로 거리 계산
meta <- visium@meta.data
meta$x <- coords[rownames(meta), "x"]
meta$y <- coords[rownames(meta), "y"]
meta <- meta[!is.na(meta$prop_Acinar), ]

# 이 reference에는 "Tumor" 라벨이 없으므로, acinar 비율이 높은 spot을
# 전암병변 상피가 몰려있는 자리의 proxy로 쓴다
acinar_spots <- meta[meta$prop_Acinar > 0.5, c("x", "y")]

meta$dist_to_acinar <- RANN::nn2(acinar_spots, meta[, c("x", "y")], k = 1)$nn.dists[, 1]

p<- ggplot(meta, aes(dist_to_acinar, prop_Fibroblast)) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_smooth(method = "loess") +
  labs(x = "Distance to nearest acinar-dominant spot (px)", y = "Fibroblast proportion") +
  theme_classic()
ggsave("results/fibroblast_distance_to_acinar.png", p, width = 7, height = 5, dpi = 300)

## Step 7. 5개 샘플 전체로 확장하기
# Step 2의 reference는 그대로 재사용하고, Step 1(QC)·Step 3(RCTD)·Step 4(비율 추출)만 샘플별로 반복하는 함수로 묶습니다.

run_rctd_pipeline <- function(visium_sample, reference, min_count = 100) {
  v <- visium_sample

  # Step 1: QC 필터 (거의 빈 spot만 제거)
  v <- subset(v, subset = nCount_Spatial > min_count)

  # Step 3: SpatialRNA puck 생성
  counts_sp <- GetAssayData(v, assay = "Spatial", layer = "counts")
  coords <- GetTissueCoordinates(v)[, 1:2]
  colnames(coords) <- c("x", "y")
  coords <- coords[colnames(counts_sp), ]
  nUMI_sp <- colSums(counts_sp)
  puck <- SpatialRNA(coords, counts_sp, nUMI_sp)

  # Step 3: RCTD 실행
  myRCTD <- create.RCTD(puck, reference, max_cores = 4)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

  # Step 4: 비율 추출
  weights      <- myRCTD@results$weights
  norm_weights <- normalize_weights(weights)
  prop <- as.data.frame(as.matrix(norm_weights))
  colnames(prop) <- paste0("prop_", make.names(colnames(prop)))

  v <- AddMetaData(v, metadata = prop)

  list(visium = v, prop = prop, coords = coords)   # myRCTD 자체는 용량이 커서 리스트에는 안 담고 따로 저장
}

# 5개 샘플 전체에 대해
dir.create("results", showWarnings = FALSE)

# Step 1~6에서 이미 계산한 GSM8573641 결과를 다른 샘플과 같은 형식으로 저장
saveRDS(list(visium = visium, prop = prop, coords = coords),
        "results/GSM8573641_deconv.rds")

results_list <- list(GSM8573641 = prop)   # 다시 계산하지 않고 그대로 재사용

samples_to_run <- setdiff(names(visium_list), "GSM8573641")   # 이미 돌린 샘플 제외

for (nm in samples_to_run) {
  message("===== ", nm, " 처리 중 =====")
  res <- run_rctd_pipeline(visium_list[[nm]], reference)
  saveRDS(res, paste0("results/", nm, "_deconv.rds"))
  results_list[[nm]] <- res$prop
}

sample_means <- sapply(results_list, colMeans, na.rm = TRUE)
round(sample_means, 3)
#                  GSM8573641 GSM8573642 GSM8573643 GSM8573644 GSM8573645
# prop_Acinar           0.096      0.022      0.077      0.087      0.229
# prop_B                0.052      0.059      0.006      0.008      0.018
# prop_CD4_T            0.007      0.004      0.004      0.009      0.005
# prop_Endothelial      0.013      0.022      0.021      0.013      0.012
# prop_Fibroblast       0.543      0.613      0.625      0.620      0.500
# prop_Granulocyte      0.004      0.006      0.006      0.007      0.007
# prop_Macrophage       0.014      0.011      0.004      0.009      0.019
# prop_Mastcell         0.004      0.007      0.012      0.009      0.006
# prop_Myleoid          0.190      0.191      0.143      0.137      0.108
# prop_Neutrophils      0.028      0.019      0.027      0.027      0.027
# prop_NK               0.009      0.006      0.011      0.013      0.014
# prop_NK_T             0.025      0.013      0.031      0.034      0.034
# prop_RegT             0.007      0.007      0.010      0.008      0.009
# prop_T                0.009      0.020      0.023      0.021      0.014

## Step 8. 샘플별로 Step 1·5·6 그림 한 번에 뽑기
library(RANN)
library(scatterpie)

make_sample_plots <- function(sample_name, res) {
  visium <- res$visium
  prop   <- res$prop
  coords <- res$coords

  # Step 1. QC violin
  mt_genes <- grep("mt-", rownames(visium), ignore.case = TRUE, value = TRUE)
  qc_features <- c("nCount_Spatial", "nFeature_Spatial")
  if (length(mt_genes) > 0) {
    mt_pat <- if (any(grepl("^mt-", rownames(visium)))) "^mt-" else "^MT-"
    visium[["percent.mt"]] <- PercentageFeatureSet(visium, pattern = mt_pat)
    qc_features <- c(qc_features, "percent.mt")
  }
  p_qc <- VlnPlot(visium, features = qc_features, pt.size = 0.1)
  ggsave(paste0("results/", sample_name, "_qc_violin.png"), p_qc, width = 10, height = 4, dpi = 300)

  # Step 5-1. 관심 타입 비율
  types_of_interest <- c("prop_Fibroblast", "prop_Acinar", "prop_Macrophage", "prop_Myleoid")
  types_of_interest <- intersect(types_of_interest, colnames(visium@meta.data))
  p1 <- SpatialFeaturePlot(visium, features = types_of_interest,
                            pt.size.factor = 1.6, alpha = c(0.1, 1), ncol = 2)
  ggsave(paste0("results/", sample_name, "_spatial_proportions.png"), p1, width = 10, height = 9, dpi = 300)

  # Step 5-2. spot별 우세 타입
  prop_cols <- grep("^prop_", colnames(visium@meta.data), value = TRUE)
  m <- as.matrix(visium@meta.data[, prop_cols]); m[is.na(m)] <- 0
  visium$dominant <- sub("^prop_", "", prop_cols[max.col(m, ties.method = "first")])
  p2 <- SpatialDimPlot(visium, group.by = "dominant", pt.size.factor = 1.6)
  ggsave(paste0("results/", sample_name, "_dominant_type.png"), p2, width = 7, height = 6, dpi = 300)

  # Step 5-3. 파이 차트
  common <- intersect(rownames(coords), rownames(prop))
  pie_df <- cbind(coords[common, ], prop[common, ])
  type_cols <- colnames(prop)
  nn <- RANN::nn2(pie_df[, c("x", "y")], k = 2)$nn.dists[, 2]
  pie_df$r <- 0.45 * median(nn)
  p3 <- ggplot() +
    scatterpie::geom_scatterpie(data = pie_df, aes(x = x, y = y, r = r), cols = type_cols,
                                 pie_scale = 1, color = NA) +
    scale_y_reverse() + coord_fixed() + theme_void() +
    labs(fill = "Cell type")
  ggsave(paste0("results/", sample_name, "_spatial_pie.png"), p3, width = 9, height = 8, dpi = 300)

  # Step 6-1. 마커 유전자 vs Fibroblast 비율
  visium <- NormalizeData(visium, assay = "Spatial")
  p_marker <- SpatialFeaturePlot(visium, features = c("Col1a1", "Dcn", "prop_Fibroblast"),
                                  pt.size.factor = 1.6, ncol = 3)
  ggsave(paste0("results/", sample_name, "_marker_overlay.png"), p_marker, width = 12, height = 4, dpi = 300)

  # Step 6-2. 세포 타입 간 상관 히트맵
  prop_mat <- as.matrix(prop)
  cor_mat  <- cor(prop_mat, method = "spearman")
  png(paste0("results/", sample_name, "_celltype_correlation_heatmap.png"),
      width = 8, height = 7, units = "in", res = 300)
  heatmap(cor_mat, symm = TRUE, margins = c(12, 12))
  dev.off()

  # Step 6-3. Acinar 영역으로부터의 거리 vs Fibroblast 비율
  meta <- visium@meta.data
  meta$x <- coords[rownames(meta), "x"]
  meta$y <- coords[rownames(meta), "y"]
  meta <- meta[!is.na(meta$prop_Acinar), ]
  acinar_spots <- meta[meta$prop_Acinar > 0.5, c("x", "y")]

  if (nrow(acinar_spots) >= 1) {
    meta$dist_to_acinar <- RANN::nn2(acinar_spots, meta[, c("x", "y")], k = 1)$nn.dists[, 1]
    p_dist <- ggplot(meta, aes(dist_to_acinar, prop_Fibroblast)) +
      geom_point(alpha = 0.2, size = 0.6) +
      geom_smooth(method = "loess") +
      labs(x = "Distance to nearest acinar-dominant spot (px)", y = "Fibroblast proportion") +
      theme_classic()
    ggsave(paste0("results/", sample_name, "_acinar_distance.png"), p_dist, width = 6, height = 5, dpi = 300)
  } else {
    message(sample_name, ": prop_Acinar > 0.5인 spot이 없어 6-3 그림을 건너뜁니다.")
  }

  invisible(NULL)
}

for (nm in names(visium_list)) {
  res <- readRDS(paste0("results/", nm, "_deconv.rds"))
  message("===== ", nm, " 그림 생성 중 =====")
  make_sample_plots(nm, res)
}

## Step 9. 5개 샘플 결과 한눈에 비교
# 각 샘플을 숫자 마지막(marker_cor, acinar_dist_cor 등)로 요약해서 표 하나와 비교 그림 세 개로 모아보는 코드
# (1) RCTD가 매긴 Fibroblast 비율이 실제 마커 유전자 발현과 맞는가(내적 타당성)
# (2) Acinar 영역에서 멀어질수록 Fibroblast 비율이 증가하는 경향이 5개 샘플 모두에서 재현
# https://www.mdpi.com/2072-6694/14/13/3293 => desmoplastic stroma는 침습암뿐 아니라 PanIN/IPMN 전암병변 단계부터 이미 형성됨

library(RANN)

samples <- sub("_deconv\\.rds$", "", list.files("results", pattern = "_deconv\\.rds$"))

summary_rows <- list()
comp_rows    <- list()
dist_rows    <- list()

for (nm in samples) {
  res    <- readRDS(paste0("results/", nm, "_deconv.rds"))
  visium <- res$visium
  coords <- res$coords
  prop   <- res$prop

  # --- (a) Acinar 거리 - Fibroblast 비율 관계 (Step 6-3과 동일한 계산) ---
  meta <- visium@meta.data
  meta$x <- coords[rownames(meta), "x"]
  meta$y <- coords[rownames(meta), "y"]
  meta <- meta[!is.na(meta$prop_Acinar) & !is.na(meta$prop_Fibroblast), ]
  acinar_spots <- meta[meta$prop_Acinar > 0.5, c("x", "y")]

  if (nrow(acinar_spots) == 0) {
    dist_cor <- NA_real_; dist_slope <- NA_real_
  } else {
    meta$dist_to_acinar <- RANN::nn2(acinar_spots, meta[, c("x", "y")], k = 1)$nn.dists[, 1]
    dist_cor   <- cor(meta$dist_to_acinar, meta$prop_Fibroblast, method = "spearman")
    dist_slope <- unname(coef(lm(prop_Fibroblast ~ dist_to_acinar, data = meta))[2])
    dist_rows[[nm]] <- data.frame(sample = nm, dist_to_acinar = meta$dist_to_acinar,
                                   prop_Fibroblast = meta$prop_Fibroblast)
  }

  # --- (b) 마커 발현 vs RCTD 비율 상관 (Step 6-1과 동일한 검증 논리) ---
  v_norm  <- NormalizeData(visium, assay = "Spatial", verbose = FALSE)
  present <- intersect(c("Col1a1", "Dcn"), rownames(v_norm))
  marker_cor <- if (length(present) > 0) {
    expr <- FetchData(v_norm, vars = present)
    cor(rowMeans(expr), visium$prop_Fibroblast, use = "complete.obs", method = "spearman")
  } else NA_real_

  summary_rows[[nm]] <- data.frame(
    sample = nm, n_spots = nrow(prop),
    marker_cor = marker_cor,
    acinar_dist_cor = dist_cor, acinar_dist_slope = dist_slope
  )

  # --- (c) 샘플별 평균 세포 구성 (전체 비교용) ---
  comp_rows[[nm]] <- data.frame(
    sample = nm, celltype = sub("^prop_", "", colnames(prop)),
    proportion = colMeans(prop, na.rm = TRUE)
  )
}

summary_df <- do.call(rbind, summary_rows); rownames(summary_df) <- NULL
comp_long  <- do.call(rbind, comp_rows);    rownames(comp_long)  <- NULL
dist_long  <- do.call(rbind, dist_rows);    rownames(dist_long)  <- NULL

summary_df
write.csv(summary_df, "results/all_samples_summary.csv", row.names = FALSE)

# 1) 샘플별 평균 세포 구성 비교
p_comp <- ggplot(comp_long, aes(x = sample, y = proportion, fill = celltype)) +
  geom_col(position = "stack") +
  labs(x = NULL, y = "평균 비율", fill = "Cell type", title = "샘플별 평균 세포 구성 비교") +
  theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("results/all_samples_composition.png", p_comp, width = 8, height = 6, dpi = 300)

# 2) 샘플별 Acinar 거리-Fibroblast 상관의 일관성
p_cor <- ggplot(summary_df, aes(x = sample, y = acinar_dist_cor)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Spearman r (거리 vs Fibroblast 비율)", title = "샘플별 Acinar 거리-Fibroblast 상관") +
  theme_classic()
ggsave("results/all_samples_acinar_cor.png", p_cor, width = 7, height = 5, dpi = 300)

# 3) 거리-비율 관계를 5개 샘플 한 화면에서 비교
p_facet <- ggplot(dist_long, aes(dist_to_acinar, prop_Fibroblast)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = "loess", color = "firebrick") +
  facet_wrap(~ sample, scales = "free_x") +
  labs(x = "Distance to nearest acinar-dominant spot (px)", y = "Fibroblast proportion",
       title = "Acinar 거리에 따른 Fibroblast 비율 — 5개 샘플 비교") +
  theme_classic()
ggsave("results/all_samples_acinar_distance_facet.png", p_facet, width = 12, height = 7, dpi = 300)

cat("=== 요약 ===\n")
cat("마커(Col1a1/Dcn)-RCTD Fibroblast 비율 상관 범위:",
    round(range(summary_df$marker_cor, na.rm = TRUE), 2), "\n")
    # 마커(Col1a1/Dcn)-RCTD Fibroblast 비율 상관 범위: 0.57 0.74 
cat("Acinar 거리-Fibroblast 비율 상관 범위:",
    round(range(summary_df$acinar_dist_cor, na.rm = TRUE), 2), "\n")
    # Acinar 거리-Fibroblast 비율 상관 범위: -0.11 0.41 

# Fibroblast: 병변 진행/섬유화 확산 정도가 개체마다 달라서 공간적 관계의 강도 자체가 달라진다

### Step 10. Acinar 거리에 따른 다른 세포 타입 비교
type_dist_rows <- list()

for (nm in samples) {
  res    <- readRDS(paste0("results/", nm, "_deconv.rds"))
  visium <- res$visium
  coords <- res$coords

  meta <- visium@meta.data
  meta$x <- coords[rownames(meta), "x"]
  meta$y <- coords[rownames(meta), "y"]
  meta <- meta[!is.na(meta$prop_Acinar), ]
  acinar_spots <- meta[meta$prop_Acinar > 0.5, c("x", "y")]
  if (nrow(acinar_spots) == 0) next
  meta$dist_to_acinar <- RANN::nn2(acinar_spots, meta[, c("x", "y")], k = 1)$nn.dists[, 1]

  other_types <- setdiff(grep("^prop_", colnames(meta), value = TRUE), "prop_Acinar")

  for (ct in other_types) {
    v    <- meta[[ct]]
    keep <- !is.na(v)
    if (sum(keep) < 10) next   # spot 수가 너무 적으면 상관계수가 불안정해지므로 건너뜁니다
    cor_val   <- suppressWarnings(cor(meta$dist_to_acinar[keep], v[keep], method = "spearman"))
    slope_val <- tryCatch(unname(coef(lm(v[keep] ~ meta$dist_to_acinar[keep]))[2]), error = function(e) NA_real_)
    type_dist_rows[[paste(nm, ct)]] <- data.frame(
      sample = nm, celltype = sub("^prop_", "", ct),
      cor = cor_val, slope = slope_val
    )
  }
}

type_dist_df <- do.call(rbind, type_dist_rows); rownames(type_dist_df) <- NULL
write.csv(type_dist_df, "results/all_samples_celltype_distance.csv", row.names = FALSE)

# 타입 × 샘플 상관 히트맵 — 한 행(타입)이 5개 샘플에서 같은 색이면 패턴이 일관된 것
 p_heat <- ggplot(type_dist_df, aes(x = sample, y = celltype, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0,
                        limits = c(-1, 1), name = "Spearman r") +
  labs(x = NULL, y = NULL,
       title = "Acinar 거리와 세포 타입 비율의 상관 — 타입 × 샘플") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("results/all_samples_celltype_distance_heatmap.png", p_heat, width = 9, height = 7, dpi = 300)

# 타입별 평균 상관 + 양수 샘플 개수(일관성 확인용)
type_summary <- do.call(rbind, lapply(split(type_dist_df, type_dist_df$celltype), function(d) {
  data.frame(celltype = unique(d$celltype),
             mean_cor = mean(d$cor, na.rm = TRUE),
             n_positive = sum(d$cor > 0, na.rm = TRUE),
             n_samples = sum(!is.na(d$cor)))
}))
rownames(type_summary) <- NULL
type_summary <- type_summary[order(-type_summary$mean_cor), ]
type_summary

# 통계적 유의성 검증 (p-value + 다중검정 보정)
sig_rows <- list()

for (nm in samples) {
  res    <- readRDS(paste0("results/", nm, "_deconv.rds"))
  visium <- res$visium
  coords <- res$coords

  meta <- visium@meta.data
  meta$x <- coords[rownames(meta), "x"]
  meta$y <- coords[rownames(meta), "y"]
  meta <- meta[!is.na(meta$prop_Acinar), ]
  acinar_spots <- meta[meta$prop_Acinar > 0.5, c("x", "y")]
  if (nrow(acinar_spots) == 0) next
  meta$dist_to_acinar <- RANN::nn2(acinar_spots, meta[, c("x", "y")], k = 1)$nn.dists[, 1]

  other_types <- setdiff(grep("^prop_", colnames(meta), value = TRUE), "prop_Acinar")

  for (ct in other_types) {
    v    <- meta[[ct]]
    keep <- !is.na(v)
    if (sum(keep) < 10) next
    test <- suppressWarnings(cor.test(meta$dist_to_acinar[keep], v[keep], method = "spearman"))
    sig_rows[[paste(nm, ct)]] <- data.frame(
      sample = nm, celltype = sub("^prop_", "", ct),
      cor = unname(test$estimate), p_value = test$p.value, n_spots = sum(keep)
    )
  }
}

sig_df <- do.call(rbind, sig_rows); rownames(sig_df) <- NULL

# 5개 샘플 × 여러 타입 = 다중검정 → BH(FDR) 보정
sig_df$p_adj       <- p.adjust(sig_df$p_value, method = "BH")
sig_df$significant <- sig_df$p_adj < 0.05
head(sig_df)
write.csv(sig_df, "results/all_samples_celltype_distance_pvalues.csv", row.names = FALSE)

# 타입별로 "유의 + 방향 일관"인 샘플 수를 집계
sig_summary <- do.call(rbind, lapply(split(sig_df, sig_df$celltype), function(d) {
  data.frame(
    celltype       = unique(d$celltype),
    mean_cor       = mean(d$cor, na.rm = TRUE),
    n_sig          = sum(d$significant, na.rm = TRUE),
    n_sig_positive = sum(d$significant & d$cor > 0, na.rm = TRUE),
    n_sig_negative = sum(d$significant & d$cor < 0, na.rm = TRUE),
    n_samples      = nrow(d)
  )
}))
rownames(sig_summary) <- NULL
sig_summary[order(-sig_summary$n_sig), ]

# 가장 강하고 완전히 일관된 신호 — 골수계 내에서 정반대 방향

# Myleoid: 5개 샘플 모두 유의, 모두 양의 상관(+0.36~+0.66) — Acinar에서 멀어질수록 뚜렷하게 증가
# Macrophage: 5개 샘플 모두 유의, 모두 음의 상관(-0.15~-0.56) — Acinar 근처에 뚜렷하게 몰림
# 같은 "골수계"라도 두 라벨이 정확히 반대로 갈린다는 게 이번 분석에서 가장 명확한 결과입니다. Acinar 주변엔 조직 상주성 대식세포가, 병변부로 갈수록 침윤성 Myleoid가 우세