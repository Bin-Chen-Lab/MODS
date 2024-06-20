### Cell cycle related work
## For Monocytes
#Monocytes$Cell.type= Monocytes@active.ident
Monocytes$Sample_celltype= paste(Monocytes$sample.type, Monocytes$Cell.type, sep = "_")
Cell.cycle=as.data.frame.matrix(table(Monocytes$Sample_celltype, Monocytes$Phase))
Cell.cycle= data.frame(t(Cell.cycle))
Cell_cycle_fraction = data.frame(apply(Cell.cycle, 2, function(x){x/sum(x)}))
#VlnPlot(Monocytes, features = c("phase"), split.by = "sample.type")
#Cell_cycle_fraction= data.frame(Cell_cycle_fraction)
Cell_cycle_fraction$phase=rownames(Cell_cycle_fraction) 

library(reshape2)
library(ggpattern)


d=melt(Cell_cycle_fraction, id.vars = "phase")

d$group=ifelse(grepl("CT_", d$variable), "CT", "MODS")
# d$variable= factor(d$variable, levels = c("CT_CD8.T.cells", "MODS_CD8.T.cells", "CT_Gamma.delta.T.cells", "MODS_Gamma.delta.T.cells", "CT_GZMK.NK",
#                                           "MODS_GZMK.NK",  "CT_Helper.T.cells","MODS_Helper.T.cells","CT_Memory.T.cells", "MODS_Memory.T.cells",
#                                           "CT_T.reg", "MODS_T.reg"))
unique(d$variable)
#d$variable= factor(d$variable, levels = c("CT_cDC1", "MODS_cDC1", "CT_cDC2", "MODS_cDC2","CT_pDC.high.PLD4", "MODS_pDC.high.PLD4",
#                                          "CT_pDC.high.TCF4", "MODS_pDC.high.TCF4")) 

d$variable= factor(d$variable, levels = c("CT_CD14.monocyte", "MODS_CD14.monocyte", "CT_CD14.CD16.monocyte", "MODS_CD14.CD16.monocyte", "CT_CD16.monocyte",
                                          "MODS_CD16.monocyte", "CT_M1.Macrophages", "MODS_M1.Macrophages", "CT_M2.Macrophages", "MODS_M2.Macrophages"))


ggplot(data = d, aes(x=variable, y=value, fill=phase))+geom_bar(stat = "identity", width = 0.8, aes(color=group))+
    theme_bw()+theme(axis.text = element_text( size=12, face="bold"),
                                          axis.text.x = element_text( size=12, face="bold", angle = 45, hjust = 1),
                                          axis.title.x = element_text(size=12, face="bold"),
                                         axis.title.y = element_blank(),
                                          legend.position = "bottom",
                                          plot.title = element_text(size=14, face="bold",lineheight = 0.9),
                                         panel.border = element_rect(fill=NA, colour = "black", linewidth=1),
                                          axis.ticks = element_line(linewidth = 0.5))+
    labs(title = "Monocytes")+scale_color_manual(values = c("black", "red1"))+
    scale_fill_manual(values = c("burlywood4", "ivory3", "darkseagreen2"))
ggsave("Monocytes_cellcycle_phase.pdf", width = 6, height = 5)

Idents(Monocytes)= Monocytes$Phase
DimPlot(Monocytes, reduction = "umap")
ggsave("Monocytes_cell_cycle_umap.pdf", width = 4.5, height = 4)












set.seed(35)
df <- data.frame(Class = factor(rep(c(1,2),times = 80), labels = c("Math","Science")),
                 StudyTime = factor(sort(sample(1:4, 16, prob = c(0.25,0.3,0.3,0.15), replace = TRUE)),labels = c("<5","5-10","10-20",">20")),
                 Nerd = factor(sapply(rep(c(0.1,0.3,0.5,0.8),c(30,50,50,30)), function(x)sample(c("Nerd","NotNerd"),size = 1, prob = c(x,1-x))),levels = c("NotNerd","Nerd")))

ggplot(data = df, aes(x = Class, fill = StudyTime, pattern = Nerd)) +
  geom_bar_pattern(position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(4)) +
  scale_pattern_manual(values = c(Nerd = "stripe", NotNerd = "none")) +
  labs(x = "Class", y = "Number of Students", pattern = "Nerd?") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))

