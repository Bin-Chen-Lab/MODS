setwd("/mnt/ufs18/home-047/ramashan/data/MODS_scRNA/")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(Signac)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
#load("final_analysis.RData")
samples= dir(path = ".", pattern="output")

CT_sample= samples[c(1:5)]
MODS_sample=samples[c(6:10)]


#samples= c("MODS001_output","MODS002_output", "CT001_output", "CT003_output" )

#sam="CT003_out"

CT.list=list()
for (sam in CT_sample){
  Sample1_data=Read10X(data.dir = paste(sam,"/outs/filtered_feature_bc_matrix/", sep = ""))
  sample_type=strsplit(sam, "_")[[1]][1]
  Sample1 <- CreateSeuratObject(counts = Sample1_data, project = sample_type, min.cells = 3, min.features = 100)
  #mito.genes <- grep(pattern = "^MT-", x = rownames(Sample1@assays[["RNA"]]), value = TRUE)
  Sample1[["percent.mt"]] <- PercentageFeatureSet(Sample1, pattern = "^MT-")
  
  # Visualize QC metrics as a violin plot
  VlnPlot(Sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste(sample_type,"_feature_RNA_MT.pdf", sep = ""), width = 6, height = 4)
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  
  plot1 <- FeatureScatter(Sample1, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Sample1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  
  ggsave(paste(sample_type,"_correlation_feature_RNA_MT.pdf", sep = ""), width = 6, height = 4)
  ### Filter the cells 
  Sample1 <- subset(Sample1, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 20)
  
  plot1 <- FeatureScatter(Sample1, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Sample1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  ggsave(paste(sample_type, "_filtered_correlation_feature_RNA_MT.pdf", sep = ""), width = 6, height = 4)
  
  ### Normalization of data
  #Sample1 <- NormalizeData(Sample1)
  
  Sample1@meta.data$sample=sample_type
  CT.list[[sample_type]]=Sample1
  
}



####  For MODS samples

MODS.list=list()
for (sam in MODS_sample){
  Sample1_data=Read10X(data.dir = paste(sam,"/outs//filtered_feature_bc_matrix/", sep = ""))
  sample_type=strsplit(sam, "_")[[1]][1]
  Sample1 <- CreateSeuratObject(counts = Sample1_data, project = sample_type, min.cells = 3, min.features = 100)
  #mito.genes <- grep(pattern = "^MT-", x = rownames(Sample1@assays[["RNA"]]), value = TRUE)
  Sample1[["percent.mt"]] <- PercentageFeatureSet(Sample1, pattern = "^MT-")
  
  # Visualize QC metrics as a violin plot
  VlnPlot(Sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste(sample_type,"_feature_RNA_MT.pdf", sep = ""), width = 6, height = 4)
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  
  plot1 <- FeatureScatter(Sample1, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Sample1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  
  ggsave(paste(sample_type,"_correlation_feature_RNA_MT.pdf", sep = ""), width = 6, height = 4)
  ### Filter the cells 
  Sample1 <- subset(Sample1, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 20)
  
  plot1 <- FeatureScatter(Sample1, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Sample1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  ggsave(paste(sample_type, "_filtered_correlation_feature_RNA_MT.pdf", sep = ""), width = 6, height = 4)
  
  ### Normalization of data
  #Sample1 <- NormalizeData(Sample1)
  
  Sample1@meta.data$sample=sample_type
  MODS.list[[sample_type]]=Sample1
  
}





### Integration of both
MODSCT.list <- c(CT001, CT002, CT003, CT004, CT005, MODS001, MODS002, MODS003, MODS004, MODS005)
names(MODSCT.list)=c("CT001", "CT002", "CT003", "CT004", "CT005", "MODS001", "MODS002", "MODS003", "MODS004", "MODS005")
#saveRDS(MODSCT.list, file = "MODSCT.list.RDS")

MODSCT.list <- lapply(X = MODSCT.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

#Next, select features for downstream integration, and run PCA on each object in the list, 
#which is required for running the alternative reciprocal PCA workflow.

features <- SelectIntegrationFeatures(object.list = MODSCT.list)
MODSCT.list <- lapply(X = MODSCT.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})



anchors <- FindIntegrationAnchors(object.list = MODSCT.list, reference = c(1, 2), reduction = "rpca",
                                  dims = 1:50)
CT.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)



CT.integrated <- ScaleData(CT.integrated, verbose = FALSE)
CT.integrated <- RunPCA(CT.integrated, verbose = FALSE)

for (i in 1:nrow(CT.integrated@meta.data)){
  if (grepl("MODS",CT.integrated@meta.data$orig.ident[i])){
    CT.integrated@meta.data$sample.type[i]="MODS"
  }else{
    CT.integrated@meta.data$sample.type[i]="CT"
  }
  
}

#saveRDS(CT.integrated, file = "CT_MODS.RDS")
#CT.integrated= readRDS("CT_MODS.RDS")
CT.integrated <- FindNeighbors(CT.integrated, dims = 1:50)
CT.integrated <- FindClusters(CT.integrated, resolution = 0.8)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
CT.integrated <- RunUMAP(CT.integrated, dims = 1:50)
saveRDS(CT.integrated, file = "CT_MODS_after_UMAP.RDS")

DimPlot(CT.integrated, group.by = c("sample.type", "seurat_clusters"),reduction = "umap", ncol = 2)
# note that you can set `label = TRUE` or use the LabelClusters function to help label

DimPlot(CT.integrated, reduction = "umap", group.by = "sample")
ggsave("MODS_CT_UMAP_samplewise.pdf", width = 6, height = 5)

Idents(CT.integrated)= CT.integrated$seurat_clusters
DimPlot(CT.integrated, reduction = "umap", label = T)+NoLegend()
ggsave("MODS_CT_UMAP_cluster.pdf", width = 5, height = 5)



Idents(CT.integrated)=CT.integrated$sample.type
DimPlot(CT.integrated, group.by = "sample.type",reduction = "umap")+NoLegend()
ggsave("sample_type_Umap_plot_CT_MODS.pdf", width = 6, height = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
CT.markers <- FindAllMarkers(CT.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

### assign the cell names to each cluster
new.cluster.ids= c("Neutrophil", "CD14.monocyte", "cDC1", "CD14+CD16.monocyte", "Helper.T.cells", "cDC2", "M1.Macrophages", "M2.Macrophages", "CD16.monocyte", "Mature.B.cell", 
           "Memory.B", "GZMK.NK", "T.reg", "Memory.T.cells", "CD8.T.cells", "GZMK.NK", "immature.Neutrophils", "Naive.B", "Neutrophil", "Neutrophil", "Neutrophil", "pDC+high.TCF4",
           "GZMK.NK", "Mature.B.cell", "Gamma.delta.T.cells", "Platelet", "pDC+high.PLD4", "Neutrophil", "Plasma.cell")


names(new.cluster.ids) = levels(CT.integrated)
CT.integrated = RenameIdents(CT.integrated, new.cluster.ids)
CT.integrated@meta.data$Cell.type=Idents(CT.integrated)

save(CT.integrated, "CT_MODS_umap_cell_annotation.RData")



### Count fraction of cells in each sample
#Idents(CT.integrated)=CT.integrated$monaco.main
Cell.number.each.type=as.data.frame.matrix(table(Idents(CT.integrated), CT.integrated$orig.ident))
cell_count_fraction = apply(Cell.number.each.type, 2, function(x){x/sum(x)})
write.csv(cell_count_fraction, "cellcount.fraction.csv")
cell_count_fraction= t(cell_count_fraction)
cell_count_fraction= data.frame(cell_count_fraction)
cell_count_fraction$sample= c(rep("CT", 5), rep("MODS",5))

source("../../plot_boxplot.R")

for (c in colnames(cell_count_fraction)[-(length(cell_count_fraction))]){
  p=plot_boxplot(df =cell_count_fraction, x = "sample", y = c, ylab = "cell.fraction", test.method = "t.test",title = c)
  print(p)
  ggsave(paste(c, "plot.pdf", sep = "_"), height = 4, width = 4)
  
}





