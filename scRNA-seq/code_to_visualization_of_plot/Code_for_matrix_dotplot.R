
### Code to plot the DotPlot similar to Seurat DaoPlot

library(ggplot2)

#library(ggplot2)

set.seed(42)

liver_cell_atlas_cells= read_xls("../atlas_data/liver_cell_atlas_cell_fraction.xls", sheet = 1)
liver_cell_atlas_cells= data.frame(liver_cell_atlas_cells)
rownames(liver_cell_atlas_cells)= liver_cell_atlas_cells[,1]
liver_cell_atlas_cells=liver_cell_atlas_cells[,-1]
liver_cell_atlas_cells_log=log(liver_cell_atlas_cells+1)


library(reshape2)
liver_cell_atlas_cells_log$sample=rownames(liver_cell_atlas_cells_log)
#d= melt(liver_cell_atlas_cells_log, id.vars = "sample")


#liver_cell_atlas_cells_log$sample= rownames(liver_cell_atlas_cells_log)

liver_cell_atlas_cells_log$sample= factor(liver_cell_atlas_cells_log$sample, 
                               levels = c("GSE115469_Healthy", "GSE136103_Healthy",  "GSE162616_Healthy", "GSE234241_Healthy" , 
                                          "GSE124395_Control" , "GSE140228_Control" ,"GSE149614_Control", "GSE162616_Control", "GSE169446_Control","GSE182159_Control",
                                          "GSE202642_Control",  "GSE242889_Control", "GSE236382_ALD", "GSE217235_NAFLD","GSE169446_NASH" ,"GSE217235_NASH", "GSE234241_CHB" ,
                                          "GSE136103_Cirrhosis", "GSE125449_HCC", "GSE140228_HCC", "GSE149614_HCC", "GSE151530_HCC", "GSE162616_HCC", "GSE202642_HCC", "GSE182159_HBV_HCC",
                                          "GSE242889_MVI-HCC" ,  "GSE151530_IHC", "GSE125449_ICC"))

d= melt(liver_cell_atlas_cells_log)

ggplot(d, aes(variable, forcats::fct_rev(sample), fill = value, size = value)) +
  geom_point(shape = 21, stroke = 0) +
  #geom_hline(yintercept = seq(.5, 4.5, 1), size = .2) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(2, 15)) +
  scale_fill_gradient(low = "blue", high = "orange", breaks = c(0,  0.50,  1.00), labels = c("0", "0.50", "1.00"), limits = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 6, angle = 45)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "right", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "Percentage.of.cells", fill = "Percentage:", x = NULL, y = NULL)


ggsave("~/Desktop/Liver_sc_Atlas/Cell_fraction_dotplot.pdf", width = 15, height = 7)
