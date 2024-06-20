### Subtypes analysis
#names(CT.integrated@meta.data)[9]="Cell.type"
cells=c("CD14.monocyte", "M2.Macrophages", "CD14+CD16.monocyte", "M1.Macrophages", "CD16.monocyte") #### Assign the group of cells 
Monocytes= subset(CT.integrated, idents=cells) 
Monocytes <- RunUMAP(Monocytes, dims = 1:50)
#CT001 <- FindNeighbors(CT001,  dims = 1:30)
Monocytes <- FindNeighbors(Monocytes, dims = 1:50)
Monocytes <- FindClusters(Monocytes, resolution = 0.1)
Idents(Monocytes)=Monocytes$Cell.type
DimPlot(Monocytes, reduction = "umap", pt.size = 0.5, label = T, label.size = 5, repel = T, split.by = "sample.type")+NoLegend()
ggsave("Monocytes_cells_umap.pdf", width = 10, height = 5)

#Idents(CT.integrated)= CT.integrated$assign.ident
Monocytes_marker=cell_type_marker[cell_type_marker$cluster %in% cells,]


Monocytes_marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

## plot heatmaps
Idents(Monocytes)= Monocytes$Cell.type
data_for_Hm= subset(Monocytes, downsample= 200)
DoHeatmap(data_for_Hm, features = top10$gene, size = 4) + NoLegend()+theme(text = element_text(size = 6))
ggsave("Monocytes_top10_HM.pdf", width = 7, height = 8)


Idents(Monocytes)= Monocytes$Phase
DimPlot(Monocytes, reduction = "umap", pt.size = 0.5)
ggsave("Monocytes_cell_cycle_umap.pdf", width = 6, height = 5)

#### Find the DE genes in each celltype and their role in pathways

immature_Neutrophils= subset(CT.integrated, ident="immature.Neutrophils")
Idents(immature_Neutrophils)=immature_Neutrophils@meta.data$sample.type
DE_genes=FindMarkers(immature_Neutrophils, ident.1 = "CT", ident.2 = "MODS", logfc.threshold = 0.5, min.pct = 0.5)
DE_gene_sig= DE_genes[DE_genes$p_val_adj < 0.05,]
DE_gene_sig$gene=rownames(DE_gene_sig)

DE_gene_sig %>%
  top_n(n = 12, wt = avg_log2FC) -> top10

#DefaultAssay(immature_Neutrophils)="RNA"
VlnPlot(immature_Neutrophils, features = top10$gene, split.by = "sample.type", pt.size = 0.01)
ggsave("immature_Neutrophils_DE_MODS_CT_vln.pdf", width = 12, height = 9)
# we want the log2 fold change 
original_gene_list <- DE_gene_sig$avg_log2FC

# name the vector
names(original_gene_list) <- DE_gene_sig$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
library(clusterProfiler)
library(org.Hs.eg.db)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
# #gse <- setReadable(gse, OrgDb = org.Hs.eg.db)
# dotplot(gse, showCategory=10, font.size=8)
# #dotplot(gse, showCategory=10, font.size=8, split=".sign") + facet_grid(.~.sign)
# 
# #emapplot(gse, showCategory = 10)
# ridgeplot(gse, showCategory=10) + labs(x = "enrichment distribution")

library(ggplot2)


up_bp=gse@result[gse@result$enrichmentScore > 0,]
down_bp=gse@result[gse@result$enrichmentScore < 0,]

up_bp_sign = up_bp[up_bp$p.adjust <= 0.05,]

if(nrow(up_bp_sign) >1){
  # up_bp_sign=up_bp_sign %>% separate(Overlap, c("A", "B"))
  # up_bp_sign$GeneRatio= as.numeric(up_bp_sign$A)/as.numeric(up_bp_sign$B)
  up_bp_sign$log10_adjp = -log10(up_bp_sign$p.adjust)
  up_bp_sign$reg="up"
  #write.csv(up_bp_sign, "GSE156063_Male_CT_vs_SARS2_Up_uniq.csv")
  up_bp_sign= up_bp_sign[c(1:5),]
}else{
  print("Nothing in up")
}
# down_enriched <- enrichr(GSE156063_MF_venn_b60_gene$Male_CT_vs_SARS2_Down, dbs)
# # printEnrich(down_enriched, paste(GSE,"down_female_CT_vs_SARS_GO.txt",sep="_") , sep = "\t", columns = c(1:9))
# down_bp <- down_enriched[["GO_Biological_Process_2018"]]
down_bp_sign = down_bp[down_bp$p.adjust <= 0.05,]

if(nrow(down_bp_sign) >1){
  # down_bp_sign=down_bp_sign %>% separate(Overlap, c("A", "B"))
  # down_bp_sign$GeneRatio= as.numeric(down_bp_sign$A)/as.numeric(down_bp_sign$B)
  down_bp_sign$log10_adjp <- -log10(down_bp_sign$p.adjust)
  down_bp_sign$reg="down"
  #write.csv(down_bp_sign, "GSE156063_Male_CT_vs_SARS2_Down_uniq.csv")
  down_bp_sign = down_bp_sign[c(1:5),]}
total_go = rbind(up_bp_sign,down_bp_sign)
total_go = total_go[complete.cases(total_go),]
total_go$reg=factor(total_go$reg, levels = c("up", "down"))

#group_color= c("up"="red","down"="blue")

if(length(unique(total_go$reg)) == 1 && "down" %in% unique(total_go$reg)){
  ggplot(data = total_go, aes(x = Description, y = log10_adjp)) +
    geom_bar(stat = "identity", fill = "cyan3", width = 0.4) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 14, face = "bold", lineheight = 0.9),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
          axis.ticks = element_line(linewidth = 0.5)) +
    labs(title = "DE genes in MODS vs CT")
} else if(length(unique(total_go$reg)) == 1 && "up" %in% unique(total_go$reg)){
  ggplot(data = total_go, aes(x = Description, y = log10_adjp)) +
    geom_bar(stat = "identity", fill = "coral1", width = 0.4) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 14, face = "bold", lineheight = 0.9),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
          axis.ticks = element_line(linewidth = 0.5)) +
    labs(title = "DE genes in MODS vs CT")
} else {
  ggplot(data = total_go, aes(x = Description, y = log10_adjp, fill = reg)) +
    geom_bar(stat = "identity", width = 0.4) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 14, face = "bold", lineheight = 0.9),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
          axis.ticks = element_line(linewidth = 0.5)) +
    labs(title = "DE genes in MODS vs CT") +
    scale_color_manual(values = c("cyan3", "coral1"))
}


ggsave("immature_Neutrophils_DE_MODS_CTGO.pdf", width = 8, height = 3)



library(enrichplot)
barplot(gse, showCategory=20) 

