### Bulk RNA-seq gene expression
bulk_MODS_CT_up= read.csv("MODS_CT_base_DGE.csv")
data_for_Hm= subset(CT.integrated, downsample= 200)
DoHeatmap(data_for_Hm, features = bulk_MODS_CT_up$gene, size = 4) + NoLegend()+theme(text = element_text(size = 12, face = "bold"))
ggsave("MODS_CT_bulk_Up_genes.pdf", width = 7, height = 5)


#### Cell cycle analysis

CT.integrated <- CellCycleScoring(
  object = CT.integrated,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

#### Visualize the cell cycle status
DimPlot(CT.integrated, group.by = "Phase",label = F, split.by = "sample.type")
ggsave("MODS_CT_image_by_cell_cycle.pdf", width = 12, height = 5)

### Comparison of cell cycle between control and MODS
cells= unique(CT.integrated$Cell.type)
test_sign = list(c("CT", "MODS"))
library(ggpubr)
vp_case1 <- function(gene_signature, file_name, test_sign, y_max, cell_name){
  plot_case1 <- function(signature){
    VlnPlot(CT.integrated, features = signature, idents = cell_name,
            pt.size = 0, 
            group.by = "sample.type", 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, "_cell_cycle.pdf")
  ggsave(file_name, width = 8, height = 4)
}



gene_sig <- c("S.Score", "G2M.Score")
comparisons <- list(c("CT", "MODS"))
for (cell in cells){
  vp_case1(gene_signature = gene_sig, file_name = cell, cell_name =cell ,test_sign = comparisons, y_max = 1)
}



