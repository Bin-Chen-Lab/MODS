###Heatmaps for Key genes in all the patients at baseline
#ECMO_MODS_gene_ids_base : dataframe contains FPKM value of key genes
#rows genes
#columns sample ids
library(dplyr)
library(gplots)
library(RColorBrewer)
mod_ECMO_MODS_gene_ids_base <- ECMO_MODS_gene_ids_base+1
key_base_validation <- as.data.frame(apply(t(mod_ECMO_MODS_gene_ids_base), 2, scale), row.names = row.names(t(mod_ECMO_MODS_gene_ids_base)))
heatmap.2(as.matrix(t(key_base_validation)),symm=FALSE,trace="none", dendrogram = "row",Rowv=T, Colv=F, col=rev(brewer.pal(11,"RdBu")),scale="none", na.rm = T, cexRow = 1, cexCol = 1, margins=c(6,8))



####Boxplot with dots
#Key gene : dataframe contains FPKM of key genes with ECMO or MODS labele 
#rows samples ids
#columns genes 
library(ggplot2)
library (ggthemes)
ggplot(key_gene, aes(Samples, ENSG00000278637)) + 
  geom_boxplot(outlier.shape=NA, color=c("blue","#999999","#56B4E9","#999999","#56B4E9","#999999","#56B4E9"))+#avoid plotting outliers twice
  ylim (0,9)+geom_jitter(position=position_jitter(width=.1, height=0))+scale_x_discrete(limits=c("CT","ECMO", "MODS", "ECMO_72h", "MODS_72h","ECMO_8d", "MODS_8d"))+
  theme_few()+theme(axis.text = element_text(size=12, face="bold"))+border(color = "black", size = 1.2)
#calculate the p-value using the t.test for each gene
t.test(gene~ Samples, data = key_gene)