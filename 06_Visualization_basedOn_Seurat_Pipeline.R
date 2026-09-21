## I) Inbuilt functions (RidgePlot,VlnPlot, FeaturePlot, DotPlot,DoHeatmap, DimPlot) 
# **First lets read files (.RDS) saved from our earlier analysis (cell-type annotation and sub-setting fibroblast cells)**

getwd()

input.file = "single_cell/Session 2/Clusterwithcompletelabels.RDS"
seurat.obj = readRDS(input.file)
DimPlot(seurat.obj, label.size = 4, label = TRUE) + NoLegend()

input.file = "single_cell/Session 3/Session 3/Fibroblast_Clustering.RDS"
seurat.obj1 = readRDS(input.file)
DimPlot(seurat.obj1, label.size = 4, label = TRUE) + NoLegend()

# **Create a data frame to save few of the marker genes to use for later stage**
Ap_CAF= c("H2-Ab1", "Krt8", "Krt18", "H2-Aa","Cd74","Cd83")

# 1. Ridge Plot - Shows a graphical distribution of gene expression across either clusters or labelled cell types.You can visualize expression of gene at different levels and in which cluster or cell types it is highly expressing.
RidgePlot(seurat.obj1, features = Ap_CAF)

# 2. Another method is using violin plot which is similar to a box plot, but instead of showing the median and interquartile range of the data, it shows the entire distribution of the gene expression.
VlnPlot(seurat.obj1, features = Ap_CAF)

# 3. Feature plot - It shows pattern of gene expression or abundance  across different cells 
FeaturePlot(seurat.obj1, features = Ap_CAF)
FeaturePlot(seurat.obj1, features = Ap_CAF, pt.size = 2.0)

# **Calculate feature-specific contrast levels based on quantiles of non-zero expression**
FeaturePlot(seurat.obj1, features = Ap_CAF, min.cutoff = "q10", max.cutoff = "q90")

# **Visualize co-expression of two features simultaneously**
FeaturePlot(seurat.obj1, features = c("Col3a1", "Cd74"), blend = TRUE)
FeaturePlot(seurat.obj1, features = c("Col3a1", "Cd74"), combine = FALSE, blend = TRUE, cols = c("grey", "green", "red"))

# 4. Dot plot 
DotPlot(seurat.obj1, features = Ap_CAF, cols = c("blue", "red"),scale = FALSE) + RotatedAxis()

# 5. Heatmap
DoHeatmap(seurat.obj1, features = Ap_CAF ,size=3)

# **To visualize how variable genes are expressing across the clusters**
DoHeatmap(seurat.obj1, features = VariableFeatures(seurat.obj1)[1:20], cells = 1:200)

# **For heatmap prefer visualizing using dittoSeq package**
library(dittoSeq)
dittoHeatmap (seurat.obj, Ap_CAF, annot.by = "seurat_clusters")

# 6. Dimension reduction plot
DimPlot(seurat.obj, label.size = 4, label = TRUE) + NoLegend()

seurat.obj = RunTSNE(seurat.obj, reduction.use = "pca", dims.use = 1:35, do.fast = T)
p2 = TSNEPlot(seurat.obj)
p2

## Applying themes to plots
baseplot= DimPlot(seurat.obj, label.size = 4, label = TRUE) + NoLegend()

baseplot + labs(title = "Clustering analysis with annotation")

remotes::install_github('sjessa/ggmin')
baseplot + ggmin::theme_powerpoint()

baseplot + DarkTheme()

library(dplyr)
qc.metrics = seurat.obj@meta.data

qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Example of plotting QC metrics") +
  geom_hline(yintercept = 750) +
  geom_hline(yintercept = 2000) 
