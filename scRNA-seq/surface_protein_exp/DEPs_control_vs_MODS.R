library(dplyr)
library(Seurat)
library(SPIDER)

load('/home/ubuntu/single_cell/MODS/CT_MODS_umap_cell_annotation.RData')

RNA_MODS <- CT.integrated[, filter(CT.integrated@meta.data, sample.type == 'MODS') %>% rownames()]
RNA_control <- CT.integrated[, filter(CT.integrated@meta.data, sample.type == 'CT') %>% rownames()] 

RNA_MODS[["study"]] = RNA_MODS@meta.data$orig.ident #36865 cells
RNA_control[["study"]] = RNA_control@meta.data$orig.ident #49974 cells

ADT_MODS_seen <- read.csv('/home/ubuntu/single_cell/MODS/SPIDER_MODS/all_seen_proteins_predicted.csv', stringsAsFactors = F, row.names = 1, check.names = F)
ADT_MODS_unseen <- read.csv('/home/ubuntu/single_cell/MODS/SPIDER_MODS/all_unseen_proteins_predicted.csv', stringsAsFactors = F, row.names = 1, check.names = F) 
ADT_control_seen <- read.csv('/home/ubuntu/single_cell/MODS/SPIDER_control/all_seen_proteins_predicted.csv', stringsAsFactors = F, row.names = 1, check.names = F)
ADT_control_unseen <- read.csv('/home/ubuntu/single_cell/MODS/SPIDER_control/all_unseen_proteins_predicted.csv', stringsAsFactors = F, row.names = 1, check.names = F)

highly_confident_seen <- read.csv('/home/ubuntu/single_cell/SPIDER_python/SPIDER_python/SPIDER_weight/cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', stringsAsFactors = F)
highly_confident_seen = highly_confident_seen$X[highly_confident_seen$pearson > 0.6] #163 proteins
highly_confident_unseen_MODS <- read.csv('/home/ubuntu/single_cell/MODS/SPIDER_MODS/confidence_score_all_unseen_proteins.csv', stringsAsFactors = F, row.names = 1, check.names = F)
highly_confident_unseen_MODS = rownames(highly_confident_unseen_MODS)[highly_confident_unseen_MODS$max_inferred_coef > 0.85]
highly_confident_unseen_MODS = highly_confident_unseen_MODS[!is.na(highly_confident_unseen_MODS)] #454 proteins
highly_confident_unseen_control <- read.csv('/home/ubuntu/single_cell/MODS/SPIDER_control/confidence_score_all_unseen_proteins.csv', stringsAsFactors = F, row.names = 1, check.names = F)
highly_confident_unseen_control = rownames(highly_confident_unseen_control)[highly_confident_unseen_control$max_inferred_coef > 0.85]
highly_confident_unseen_control = highly_confident_unseen_control[!is.na(highly_confident_unseen_control)] #407 proteins

use_highly_confident_unseen_control_MODS <- intersect(highly_confident_unseen_control, highly_confident_unseen_MODS) #219 proteins
use_highly_confident_unseen_control_MODS = setdiff(use_highly_confident_unseen_control_MODS, highly_confident_seen) #209 proteins
use_RNA_MODS_control <- merge(RNA_MODS, RNA_control)

ADT_MODS_seen = t(ADT_MODS_seen)
ADT_control_seen = t(ADT_control_seen)
ADT_MODS_control_seen <- cbind(ADT_MODS_seen[highly_confident_seen, ], ADT_control_seen[highly_confident_seen, ])[, colnames(use_RNA_MODS_control)] #163 86839

ADT_MODS_unseen = t(ADT_MODS_unseen)
ADT_control_unseen = t(ADT_control_unseen)
ADT_MODS_control_unseen <- cbind(ADT_MODS_unseen[use_highly_confident_unseen_control_MODS, ], ADT_control_unseen[use_highly_confident_unseen_control_MODS, ])[, colnames(use_RNA_MODS_control)] #209 86839

ADT_MODS_control_seen_unseen <- rbind(ADT_MODS_control_seen, ADT_MODS_control_unseen)
use_RNA_MODS_control[["ADT"]] <- CreateAssayObject(counts = ADT_MODS_control_seen_unseen)
use_RNA_MODS_control@assays$ADT@data = ADT_MODS_control_seen_unseen
use_RNA_MODS_control <- ScaleData(use_RNA_MODS_control, assay = "ADT")

celltype_MODS <- unique(as.character(RNA_MODS@meta.data$Cell.type))
celltype_control <- unique(as.character(RNA_control@meta.data$Cell.type))
celltype = intersect(celltype_MODS, celltype_control)

#DEPs between MODS and control, for each cell type:
for (c in celltype) {
  
  tmp <- use_RNA_MODS_control[, colnames(use_RNA_MODS_control)[use_RNA_MODS_control@meta.data$Cell.type == c]]
  
  if((length(colnames(tmp)[grep('MODS', tmp@meta.data$study)])>=3) & (length(colnames(tmp)[grep('CT', tmp@meta.data$study)])>=3)){
  adt_markers <- FindMarkers(tmp, ident.1 = colnames(tmp)[grep('MODS', tmp@meta.data$study)], ident.2 = colnames(tmp)[grep('CT', tmp@meta.data$study)], assay = "ADT")
  }
  
  assign(paste0(c, '_DEPs'), adt_markers)

}

#DEPs for each cell type:
tmp2 <- FindAllMarkers(use_RNA_MODS_control, assay = 'ADT', only.pos = TRUE, verbose = FALSE)
celltype_DEPs <- tmp2
write.csv(celltype_DEPs, '/home/ubuntu/single_cell/MODS/celltype_DEPs.csv')

save.image('/home/ubuntu/single_cell/MODS/DEPs.RData')

#-----------------------------------------------------------------------------------------------------------
#Visualize by volcano plot:
library(ggrepel)

#
#volcano plots for cell type markers:
for (i in unique(as.character(celltype_DEPs$cluster))){
  
dz_signature = celltype_DEPs[celltype_DEPs$cluster == i, ]
dz_signature$diffexpressed='NO'
dz_signature[dz_signature$p_val_adj<0.05 & dz_signature$avg_log2FC>0.5,'diffexpressed']='UP'
dz_signature[dz_signature$p_val_adj<0.05 & dz_signature$avg_log2FC<(-0.5),'diffexpressed']='DOWN'
table(dz_signature$diffexpressed)
#label DE genes
dz_signature$delabel <- NA
dz_signature$delabel[dz_signature$diffexpressed != "NO" & dz_signature$avg_log2FC ] <- as.character(rownames(dz_signature)[dz_signature$diffexpressed != "NO"])

#fancy volcano plot
pdf(paste0('/home/ubuntu/single_cell/MODS/celltype_', i, '_DEPs_high_confidence_volcano_plot.pdf'))
#pdf('/Users/ruoqiaochen/Downloads/Control_DEPs_gd1T_gd2T_high_confidence_volcano_plot.pdf')
p=ggplot(data=dz_signature, aes(x=avg_log2FC , y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
print(p)
dev.off()

}

#
#volcano plots for disease markers:
for (i in celltype){
  
  dz_signature = get(paste0(i, '_DEPs'))
  write.csv(dz_signature, paste0('/home/ubuntu/single_cell/MODS/MODS_DEPS_', i, '.csv'))
  dz_signature$diffexpressed='NO'
  dz_signature[dz_signature$p_val_adj<0.05 & dz_signature$avg_log2FC>0.5,'diffexpressed']='UP'
  dz_signature[dz_signature$p_val_adj<0.05 & dz_signature$avg_log2FC<(-0.5),'diffexpressed']='DOWN'
  table(dz_signature$diffexpressed)
  #label DE genes
  dz_signature$delabel <- NA
  dz_signature$delabel[dz_signature$diffexpressed != "NO" & dz_signature$avg_log2FC ] <- as.character(rownames(dz_signature)[dz_signature$diffexpressed != "NO"])
  
  #fancy volcano plot
  pdf(paste0('/home/ubuntu/single_cell/MODS/MODS_', i, '_DEPs_high_confidence_volcano_plot.pdf'))
  #pdf('/Users/ruoqiaochen/Downloads/Control_DEPs_gd1T_gd2T_high_confidence_volcano_plot.pdf')
  p=ggplot(data=dz_signature, aes(x=avg_log2FC , y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.5, 0.5), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  print(p)
  dev.off()
  
}
