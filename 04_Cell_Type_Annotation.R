options("install.lock"=FALSE)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

install.packages("pacman")
library(pacman)
pacman::p_load(Seurat, Matrix, ggplot2, dplyr, SingleCellExperiment, scran, scuttle, celldex, SingleR, TabulaMurisData)

################################################################################################
### Unsupervised Reference based annotation using TabulaMurius reference ## 
################################################################################################



## Step I -  Prepare the reference data 
# install.packages("dbplyr", repos = "https://cloud.r-project.org", type = "binary")

library(ExperimentHub)####Experiment hub provides a central location to store curated data.
eh <- ExperimentHub() #### Create an ExperimentHub object
query(eh, "TabulaMurisData")

# **Two types data based on different technique EH1617 10x (droplet) and EH1618 SmartSeq2 (on FACS-sorted cells)**
eh[['EH1618']] ### The individual data set can be accessed using title number by this script
Complete_ref <- eh[['EH1618']] #### I have selected EH618 because it contains cells for Pancreas (based on the reference)

Complete_ref
colnames(colData(Complete_ref))
colnames(rowData(Complete_ref))
View(as.data.frame(Complete_ref@colData))
head(colData(Complete_ref))
head((rowData(Complete_ref)))
write.csv(Complete_ref$cell_ontology_class, "single_cell/Session 2/check.csv")
write.csv(Complete_ref$tissue, "single_cell/Session 2/check.csv")

# Eliminate NA
Complete_ref <- Complete_ref[,!is.na(Complete_ref$cell_ontology_class)]

library(scuttle)
Complete_ref <- logNormCounts(Complete_ref)

## Step II - Prepare the query data set based on Seurat pipeline 

input.file="single_cell/Seurat_Clustering.RDS"
seurat.obj = readRDS(input.file)
DimPlot(object = seurat.obj, reduction = "umap", label = T, pt.size = 0.5, label.size = 5)

## Step III - Annotation using Single R tool
sce.obj = as.SingleCellExperiment(seurat.obj)

results <- SingleR(test = sce.obj, ref = Complete_ref, labels = Complete_ref$cell_ontology_class)

Summary = table(results$labels)
Summary

write.csv(Summary, row.names = TRUE, "single_cell/Session 2/Summaryexample.csv")

## Step IV - Prepare the output table to check the cells assigned to the clusters
Clustertocell=table(results$labels, seurat.obj$seurat_clusters)
Clustertocell
write.csv(Clustertocell, "single_cell/Session 2/Clustertocell.csv", row.names = TRUE)

################################################################################################
### Unsupervised Reference based annotation using Mouse-RNAseq reference ## 
################################################################################################

mimd.se = celldex::MouseRNAseqData()
sceP  = as.SingleCellExperiment(seurat.obj)

# **Predict using the Broad labels**
pred.mimd <- SingleR(test = sceP, 
                     ref = mimd.se, 
                     labels = mimd.se$label.main)

# # Predict using the Fine-grained labels
pred.mimd.fine <- SingleR(test = sceP, 
                          ref = mimd.se, 
                          labels = mimd.se$label.fine)

# Combine the labels for improved accuracy
pred.comb <- combineCommonResults(
  list(
    "Broad" = pred.mimd,
    "Fine" = pred.mimd.fine
  )
)
head(pred.comb)

# **Add the predicted labels and compare to clusters**
seurat.obj$predicted_id <- pred.comb$pruned.labels
Clustertocell_Mouse_RNA_seq=table(seurat.obj$predicted_id, seurat.obj$seurat_clusters)
Clustertocell_Mouse_RNA_seq
write.csv(Clustertocell_Mouse_RNA_seq, "single_cell/Session 2/Clustertocell_Mouse_RNA_seq.csv", row.names = TRUE)

################################################################################################
### DEG based Semi-supervised manual annotation ## 
################################################################################################

All.markers = FindAllMarkers(seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

topDEGs <- All.markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC) # n is the number of top marker genes
topDEGs [1:20, ]

# -**p_val** is the p-value from the Wilcoxon rank sum test
# -**avg_logFC** is the log fold change between the cluster and the aggregate of others
# -**pct.1** is the percentage of cells in the cluster with detectable expression of the gene
# -**pct.2** is the percentage of cells in the aggregate of other clusters that express it
# -**p_val_adj** is the adjusted p-value
# -**cluster** is the cluster being tested
# -**gene** is the gene

write.csv(All.markers, row.names = TRUE, "single_cell/Session 2/Allmarkers_MAST2.csv", quote=TRUE)
write.csv(topDEGs, row.names = TRUE, col.names = TRUE, "single_cell/Session 2/TOP10_MAST1.csv")

Clustertocell <- table(results$labels, seurat.obj$seurat_clusters)
Clustertocell

# 다른 클러스터에 비해서 16이 비교적 균등한 라벨로 이루어짐
################################################################################################
### Annotation of an unknown cluster (#16) ## 
################################################################################################

Cluster16_CalcaVln <- VlnPlot(seurat.obj, features = "Calca")
Cluster16_CalcaFP <- FeaturePlot(seurat.obj, features = "Calca")

# cluster 16에서 모여있음
ggsave("single_cell/Session 2/Cluster16_CalcaVln.png", plot = Cluster16_CalcaVln, width = 8, height = 6, dpi = 300)
ggsave("single_cell/Session 2/Cluster16_CalcaFP.png", plot = Cluster16_CalcaFP, width = 8, height = 6, dpi = 300)


################################################################################################
### Adding final annotation to the clusters ## 
################################################################################################

# **Add the predicted labels- Prepare final IDs**
cluster.celltype = c("Acinar", "CD4+T", "NK", "RegT", "NK+T", 
                     "NK+T", "B", "Myleoid", "Acinar", "Granulocyte", "Fibroblast",
                     "Acinar", "Macrophage", "Acinar", "Mastcell", "Macrophage",
                     "Endothelial", "Fibroblast", "T", "B","Neutrophils")
cluster.celltype
names(cluster.celltype) = levels(seurat.obj)
cluster.celltype

seurat.obj = RenameIdents(seurat.obj, 
                          cluster.celltype)
DimPlot(seurat.obj)
DimPlot(seurat.obj, label.size = 4, label = TRUE) + NoLegend() 

DimPlot_Label <- DimPlot(seurat.obj, label.size = 4, label = TRUE) + NoLegend()
ggsave("single_cell/Session 2/DimPlot_Label.png", plot = DimPlot_Label, width = 8, height = 6, dpi = 300)

# **Add the final IDs to the metadata**
seurat.obj = AddMetaData(object = seurat.obj, 
                         metadata = as.character(Idents(seurat.obj)), 
                         col.name = "celltypes") 

saveRDS(object = seurat.obj, file = paste0("single_cell/Session 2/Clusterwithcompletelabels.RDS")) 
