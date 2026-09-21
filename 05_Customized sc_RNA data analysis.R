getwd()

input.file="single_cell/Seurat_Clustering.RDS"
seurat.obj = readRDS(input.file)
DimPlot(seurat.obj)
DimPlot(seurat.obj, label.size = 4, label = TRUE) + NoLegend()

try_subset = subset(x=seurat.obj, idents = c(10,17))

try_subset = NormalizeData(try_subset,normalization.method = "LogNormalize", scale.factor = 10000)

try_subset = FindVariableFeatures(try_subset)

try_subset <- ScaleData(try_subset)

try_subset = RunPCA(try_subset, features = VariableFeatures(object = try_subset))
ElbowPlot(try_subset, ndims = 50)

try_subset = FindNeighbors(try_subset, dims = 1:30, force.recalc = T)
try_subset = FindClusters(object = try_subset, resolution = 1.0)

try_subset = RunUMAP(object = try_subset, dims = 1:30)
DimPlot(object = try_subset, reduction = "umap", label = T, pt.size = 1.0, label.size = 5)

saveRDS(try_subset, file = paste0("single_cell/Session 3/Session 3/Fibroblast_Clustering.RDS"))

# i_CAF specific markers (Pancreas) - "Vim", "Fap", "Col3a1", "Il6", "Cxcl2", "C3","Pdgfra", "Cfb", "Cfh", "Cxcl12", "Cxcl1", "Cxcl10", "Ccl2", "Il1r1"
# My_CAF specific markers (Pancreas) - "Dcn","Postn","Tpm1","Tpm2","Thbs2","Thy1", "Mmp11","Bgn","Col8a1","Col15a1","Igfbp7", "Acta2"
# Ap_CAF specific markers (Pancreas) - "H2-Ab1", "H2-Aa", "Cd74", "Cd83","Krt8","Krt18"

Allmarkers = FindAllMarkers(try_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Allmarkers, row.names = TRUE, "single_cell/Session 3/Session 3/Allmarkers.csv", quote=TRUE)

###### Pulling out CAF marker data from all marker file (based of differential expression patterns)####
# **Prepare a dataframe with respective markers to pull**

Gene=c("COL3A1","Serpine2","Cxcl14", "Crlf1", "Igfbp3", "Acta2", "Tagln", "Thy1", "Col8a1",
       "Cthrc1", "Sfrp1", "Tnc,Sparcl1", "Col15a1", "Col12a1", "Tgfb1", 
       "Col1a1", "Sdc1", "Cilp", "H19", "Thbs2", "CALD1", "MMP11", "HOPX", "BGN",
       "MYL9", "IGFBP7", "TPM2", "INHBA", "COL10A1", "TPM1", "POSTN", "GRP", "CST1")

# **Create a dtaframe**
df = data.frame(Gene)
df
My_caf_markers =Allmarkers[Allmarkers$gene %in% unique(df$Gene), ] 
My_caf_markers
write.csv(My_caf_markers, row.names = TRUE, col.names = TRUE, "single_cell/Session 3/Session 3/My_caf_markers.csv")

# **But this method is difficult**
## we can import the marker list from csv file double quote it, however dQuote (need to transpose the genes "t") and then collapse with ', '
# **For ap_CAF**

df1 = read.csv("single_cell/Session 3/Session 3/AP_CAF.csv")
paste(dQuote(t(df1)), collapse=', ')
Gene = c(df1)
Gene

Allmarkers = read.csv("single_cell/Session 3/Session 3/Allmarkers.csv")
head(Allmarkers)
My_I_markers = Allmarkers[Allmarkers$gene %in% unique(df1$Gene), ]
My_I_markers
write.csv(My_I_markers, row.names = TRUE, "single_cell/Session 3/Session 3/AP_CAFpull.csv")

# **For i_CAF**
df1 = read.csv("single_cell/Session 3/Session 3/I_CAF.csv")
paste(dQuote(t(df1)), collapse=', ')
Gene = c(df1)
Gene

Allmarkers = read.csv("single_cell/Session 3/Session 3/Allmarkers.csv")
head(Allmarkers)
My_I_markers = Allmarkers[Allmarkers$gene %in% unique(df1$Gene), ]
My_I_markers
write.csv(My_I_markers, row.names = TRUE, "single_cell/Session 3/Session 3/I_CAFpull.csv")

# **For my_CAF**
df1 = read.csv("single_cell/Session 3/Session 3/my_CAF.csv")
paste(dQuote(t(df1)), collapse=', ')
Gene = c(df1)
Gene

Allmarkers = read.csv("single_cell/Session 3/Session 3/Allmarkers.csv")
head(Allmarkers)
My_I_markers = Allmarkers[Allmarkers$gene %in% unique(df1$Gene), ]
My_I_markers
write.csv(My_I_markers, row.names = TRUE, "single_cell/Session 3/Session 3/my_CAFpull.csv")

## Visulaization for Ap-CAF specific markers 
VlnPlot(try_subset, features = c("H2-Ab1", "Krt8", "Krt18", "H2-Aa","Cd74","Cd83"))
FeaturePlot(try_subset, features = c("H2-Ab1", "Krt8", "Krt18", "H2-Aa","Cd74","Cd83"))

## Visulaization for i-CAF specific markers 
VlnPlot(try_subset, features = c("Cxcl2", "Col3a1", "Il6", "Cfh", "Pdgfra", "Cxcl12"))
FeaturePlot(try_subset, features = c("Cxcl2", "Col3a1", "Il6", "Cfh", "Pdgfra", "Cxcl12"))

## Visulaization for my-CAF specific markers 
VlnPlot(try_subset, features = c("Bgn","Dcn", "Igfbp7","Postn"))
FeaturePlot(try_subset, features = c("Bgn","Dcn", "Igfbp7","Postn"))

##Dotplot
cell.markers = c("H2-Ab1", "H2-Aa", "Cd74", "Cd83","Krt8","Krt18","Dcn","Postn","Tpm1","Tpm2","Thbs2","Thy1",
          "Mmp11","Bgn","Col8a1","Col15a1","Igfbp7", "Acta2","Vim", "Fap", "Col3a1", "Il6", 
          "Cxcl2", "C3","Pdgfra", "Cfb", "Cfh", "Cxcl12", "Cxcl1",
          "Cxcl10", "Ccl2", "Il1r1")
DotPlot(try_subset, features = cell.markers, cols = c("blue", "red"),scale = FALSE) + RotatedAxis()

## Heatmap
BiocManager::install("dittoSeq")
library(dittoSeq)
cell.markers = c("H2-Ab1", "H2-Aa", "Cd74", "Cd83","Krt8","Krt18","Dcn","Postn","Tpm1","Tpm2","Thbs2","Thy1",
          "Mmp11","Bgn","Col8a1","Col15a1","Igfbp7", "Acta2","Vim", "Fap", "Col3a1", "Il6", 
          "Cxcl2", "C3","Pdgfra", "Cfb", "Cfh", "Cxcl12", "Cxcl1",
          "Cxcl10", "Ccl2", "Il1r1")
dittoHeatmap (try_subset, cell.markers, annot.by = "seurat_clusters")

##Blending two genes to visualize expression
FeaturePlot(try_subset, features = c("Col3a1", "Cd74"), combine = FALSE, blend = TRUE, cols = c("grey","green", "red"))

##### Subsetting cells based on genes
# Acta = subset(seurat.obj, subset= Acta2 >0)
# Fap = subset(seurat.obj, subset= Fap >0)
# Il6 = subset(seurat.obj, subset= Il6 >0)

Acta = subset(try_subset, subset= Acta2 >0)
Fap = subset(try_subset, subset= Fap >0)
Il6 = subset(try_subset, subset= Il6 >0)

#### Adding new coloumn in meta data 
Acta@meta.data$Gene=paste0(Acta@meta.data$Gene, 'Acta')
Fap@meta.data$Gene=paste0(Fap@meta.data$Gene, 'Fap')
Il6@meta.data$Gene=paste0(Il6@meta.data$Gene, 'Il6')

## Merging three objects
Try4 = merge(x=Acta, y=list(Fap, Il6))
Try4 = JoinLayers(Try4) 
Idents(Try4) = "Gene"
head(Try4)

## Identifying DEG
ACTA_FAP_Il6=FindAllMarkers(Try4, logfc.threshold = 0.25, only.pos = TRUE)
head(ACTA_FAP_Il6)
write.csv(ACTA_FAP_Il6, row.names = TRUE, "single_cell/Session 3/Session 3/ACTA_FAP_Il6.csv")

Try = merge(x = Acta, y = Fap)
Try = JoinLayers(Try)
Idents(Try) = "Gene"
ACTA_FAP=FindAllMarkers(Try, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(ACTA_FAP, row.names = TRUE, "single_cell/Session 3/Session 3/ACTA_FAP.csv")

Try2 = merge(x=Acta, y=Il6)
Try2 = JoinLayers(Try2)
Idents(Try2) = "Gene"
ACTA_Il6=FindAllMarkers(Try2, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(ACTA_Il6, row.names = TRUE, "single_cell/Session 3/Session 3/ACTA_Il6.csv")

Try3 = merge(x=Fap, y=Il6)
Try3 = JoinLayers(Try3)
Idents(Try3) = "Gene"
Fap_Il6=FindAllMarkers(Try3, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(Fap_Il6, row.names = TRUE, "single_cell/Session 3/Session 3/Fap_Il6.csv")
