###Clustering of patients based on key genes identified by logistic regression
library(ggplot2)
library(dplyr)
library(plotly)
library(ggthemes)
#base_FPKM.csv: dataframe contains FPKM values for all the patients (CT, ECMO and MODS) at baseline  
#columns: samples name
#Rows: gene
#ECMO_MODS_key_genes: dataframe containes list of genes strongly associated with ECMO####
load("~/Desktop/ECMO_Vs_MODS_data/data/ECMO_MODS_key_genes.RData")
load("~/Desktop/ECMO_Vs_MODS_data/data/base_FPKM_count.RData")
ECMO_MODS_key_genes <- ECMO_MODS_key_genes %>% left_join(base_FPKM_count, by=c('gene_Ids'='ensembl'))
row.names(ECMO_MODS_key_genes) <- ECMO_MODS_key_genes[,1]
ECMO_MODS_key_genes[,1] <- NULL
pca_value_base <-prcomp(t(ECMO_MODS_key_genes), scale. = T)
scores_base = as.data.frame(pca_value_base$x)

#create dataframe contains sample ids list
group_samples_base <- as.data.frame(t(colnames(base_FPKM_count)))
group_samples_base <- t(group_samples_base)
colnames(group_samples_base) <- "Ids"
group_samples_base <- as.data.frame(group_samples_base)
group_samples_base <- as.data.frame(group_samples_base[-28,])

#Add sample class and combined the it in PCA values file####
base_line_class <- as.data.frame(c(rep("CT", times=4), rep("ECMO", times=6), rep("MODS", times=17)))
colnames(base_line_class) <- "class"
group_samples_base <-cbind(group_samples_base,base_line_class)
colnames(group_samples_base) <- c("sample", "class")
score_base_ids <- cbind(scores_base, group_samples_base)

#Create and save the cluster image
png(filename = "baseline_patients_cluster.png")
cluster <- ggplot(data = score_base_ids, aes(x=PC1, y=PC2, label=rownames(score_base_ids)), alpha = 0.9, size = 25)+
  geom_point(color=c(rep("mediumseagreen", 4), rep("sienna1", 6), rep("steelblue", 17)), size=3)+theme_few()+theme(axis.text = element_text(size=10, face="bold"))
plot(cluster)
dev.off()

#Create and save the clustering image in 3d

if (!require("processx")) install.packages("processx")
plotly_base <- plot_ly(score_base_ids, x=~PC1, y=~PC2, z=~PC3, color =factor(score_base_ids$class))
orca(plotly_base, "Baseline_3d_clusterplot.png")

# write the PCA score in a file 
write.table(score_base_ids, file = "PCA_score_for_baseline_clustering.txt", row.names = T)


####Clustering of total patients belong to different time points (0h, 72h and 8days) using the key genes###
#ECMO_MODS_key_genes: dataframe containes list of genes strongly associated with ECMO
#combined_FPKM.csv :dataframe containes the FPKM value for all the genes in all the patients at different time points
#columns: sample ids
#Rows: gene
load("~/Desktop/ECMO_Vs_MODS_data/data/combined_FPKM_count.RData")
load("~/Desktop/ECMO_Vs_MODS_data/data/class.Rdata")
combined_FPKM_count$ensembl <- row.names(combined_FPKM_count)
Key_gene_ids_combined <- ECMO_MODS_key_genes %>% left_join(combined_FPKM_count, by=c('gene_Ids'='ensembl'))
row.names(Key_gene_ids_combined) <- Key_gene_ids_combined[,1]
Key_gene_ids_combined[,1] <- NULL
pca_value_combined <-prcomp(t(Key_gene_ids_combined), scale. = T)
scores_combined = as.data.frame(pca_value_combined$x)

#class: dataframe with two colums (sample ids and sample categories)
score_combined_ids <- cbind(scores_combined, class)
png(filename = "comined_patients_cluster.png")
combined_cluster <- ggplot(data = score_combined_ids, aes(x=PC1, y=PC2, label=rownames(score_combined_ids)), alpha = 0.9, size = 25)+
  geom_point(aes(color=Group), size=2)+
  theme_few()+theme(axis.text = element_text(size=10, face="bold"))
dev.off()
#Create 3d clutering plot for all the samples
if (!require("processx")) install.packages("processx")
plotly_combined <- plot_ly(score_combined_ids, x=~PC1, y=~PC2, z=~PC3, color =factor(score_combined_ids$Group))
orca(plotly_combined, "Combined_3d_clusterplot.png")
