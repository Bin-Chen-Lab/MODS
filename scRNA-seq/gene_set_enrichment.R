### Cell status and other enrichment analysis
immune_cell_marker= read.csv("Immune_cells_markers.csv")
immune_features= unique(immune_cell_marker$cells)
immune_panel= list()
for (i in 1:length(immune_features)){
  immune_panel[["gene_set"]][[i]]=immune_cell_marker$gene[immune_cell_marker$cells==immune_features[i]]
  immune_panel[["gene_names"]][[i]]=immune_features[i]
}

immune_features= immune_panel$gene_set
names(immune_features)= immune_panel$gene_names

tmp.data <- as.matrix(GetAssayData(CT.integrated,
                                   slot = "data"))
#tmp.data <- tmp.data[VariableFeatures(se.merged), ]
library(sigatlas)
t.scores <- t(sigatlas(expr = tmp.data, cellmarker = immune_features, tissue = "Other", organism = "Other"))
#names(t.scores)="malignancy"
#t.scores= data.frame(t.scores)
CT.integrated@meta.data=cbind(CT.integrated@meta.data, t.scores)


FeaturePlot(CT.integrated, features = c("Cytotoxic.cell", "Resting.cell", "Activating.cell", "T.regulatory", "Exhaustive.cells"), split.by = "sample.type", cols = c("gray", "blue"), ncol = 2)
ggsave("immunecell_feature.pdf",width = 8, height = 20)
#### Fraction of cells in each samples

vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(CT.integrated, features = signature,
            pt.size = 0, 
            group.by = "sample.type", 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, "_cell_status.pdf")
  ggsave(file_name, width = 8, height = 7)
}



gene_sig <- c("Cytotoxic.cell", "Resting.cell", "Activating.cell", "T.regulatory", "Exhaustive.cells")
comparisons <- list(c("CT", "MODS"))

library(ggpubr)
vp_case1(gene_signature = gene_sig, file_name = "immune_response" ,test_sign = comparisons, y_max = 0.6)

### specific genes expression

gene_sig= immune_features$T.regulatory  ### change the name 

vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(CT.integrated, features = signature,
            pt.size = 0, 
            group.by = "sample.type",
            y.max = y_max# add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, "gene_exp.pdf")
  ggsave(file_name, width = 8, height = 8)
}


vp_case1(gene_signature = gene_sig, file_name = "Regulatory" ,test_sign = comparisons, y_max = 6)


for (i in 1:length(immune_features)){
  DotPlot(CT.integrated, features = immune_features$Cytotoxic.cell, split.by = "sample.type") + RotatedAxis()
  ggsave("Dotplot_Cytotoxic_cell_marker_exp.pdf", width = 8, height = 10)
}

