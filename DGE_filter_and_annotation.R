# Differential gene expression (DGE) analysis analysis analysis and annotation
# EdgeR pipeline was used for DGE analysis provided in "https://github.com/Bin-Chen-Lab/octad_desktop/blob/master/Code/deprecated/Diff_Exp.R"
# Further expression of genes was filtered to identify the Differentially expressed genes
#column: Gene
# rows: sample ids
load("~/Desktop/ECMO_Vs_MODS_data/ECMO_MODS_base_DE_genes.RData")
load("~/Desktop/ECMO_Vs_MODS_data/ECMO_MODS_72h_DE_genes.RData")
load("~/Desktop/ECMO_Vs_MODS_data/ECMO_MODS_8d_DE_genes.RData")
load("~/Desktop/ECMO_Vs_MODS_data/MODS_CT_base_DE_genes.RData")
load("~/Desktop/ECMO_Vs_MODS_data/MODS_CT_72h_DE_genes.RData")
load("~/Desktop/ECMO_Vs_MODS_data/MODS_CT_8d_DE_genes.RData")
load("~/Desktop/ECMO_Vs_MODS_data/gene_description.RData")
######
#First analysis
####
ECMO_MODS_base_DGE <- filter(ECMO_MODS_base_DE_genes.csv, log2FoldChange >= 1| log2FoldChange <= -1, padj <= 0.25)
colnames(ECMO_MODS_base_DGE)[1] <- "Ids" 
ECMO_MODS_72h_DGE <- filter(ECMO_MODS_72h_DE_genes.csv, log2FoldChange >= 1| log2FoldChange <= -1, padj <= 0.25)
colnames(ECMO_MODS_72h_DGE)[1] <- "Ids"
ECMO_MODS_8d_DGE <- filter(ECMO_MODS_8d_DE_genes.csv, log2FoldChange >= 1| log2FoldChange <= -1, padj <= 0.25)
colnames(ECMO_MODS_8d_DGE)[1] <- "Ids"
MODS_CT_base_DGE <- filter(MODS_CT_base_DE_genes.csv, log2FoldChange >= 1| log2FoldChange <= -1, padj <= 0.25)
colnames(MODS_CT_base_DGE)[1] <- "Ids"
MODS_CT_72h_DGE <- filter(MODS_CT_72h_DE_genes.csv, log2FoldChange >= 1| log2FoldChange <= -1, padj <= 0.25)
colnames(MODS_CT_72h_DGE)[1] <- "Ids"
MODS_CT_8d_DGE <- filter(MODS_CT_8d_DE_genes.csv, log2FoldChange >= 1| log2FoldChange <= -1, padj <= 0.25)
colnames(MODS_CT_8d_DGE)[1] <- "Ids"

# Add gene name and annotation in the DGE
# "merged_ensembl_geneinfo.csv" file contains information of all genes of human
library(dplyr)
ECMO_MODS_base_DGE <- ECMO_MODS_base_DGE %>% left_join(gene_description %>% select(gene, description, ensembl), by=c('Ids'='ensembl'))
ECMO_MODS_72h_DGE <- ECMO_MODS_72h_DGE %>% left_join(gene_description %>% select(gene, description, ensembl), by=c('Ids'='ensembl'))
ECMO_MODS_8d_DGE <- ECMO_MODS_8d_DGE %>% left_join(gene_description %>% select(gene, description, ensembl), by=c('Ids'='ensembl'))
MODS_CT_base_DGE <- MODS_CT_base_DGE %>% left_join(gene_description %>% select(gene, description, ensembl), by=c('Ids'='ensembl'))
MODS_CT_72h_DGE <- MODS_CT_72h_DGE %>% left_join(gene_description %>% select(gene, description, ensembl), by=c('Ids'='ensembl'))
MODS_CT_8d_DGE <- MODS_CT_8d_DGE %>% left_join(gene_description %>% select(gene, description, ensembl), by=c('Ids'='ensembl'))

# Go enrichment analysis 
library(org.Hs.eg.db)
library(clusterProfiler)
png("~/Desktop/ECMO_Vs_MODS_data/GO_ECMO_MODS_base.png")
GO_ECMO_MODS_base <- enrichGO(gene = as.character(ECMO_MODS_base_DGE$Ids), OrgDb= org.Hs.eg.db, keyType= 'ENSEMBL', ont= "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
GO_ECMO_MODS_base1 <- setReadable(GO_ECMO_MODS_base, OrgDb = org.Hs.eg.db)
dotplot(GO_ECMO_MODS_base1, showCategory=30, font.size=14)
dev.off()
png("~/Desktop/ECMO_Vs_MODS_data/GO_MODS_CT_base_DGE.png")
GO_MODS_CT_base_DGE <- enrichGO(gene = as.character(MODS_CT_base_DGE$Ids), OrgDb= org.Hs.eg.db, keyType= 'ENSEMBL', ont= "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
GO_MODS_CT_base_DGE1 <- setReadable(GO_MODS_CT_base_DGE, OrgDb = org.Hs.eg.db)
dotplot(GO_MODS_CT_base_DGE1, showCategory=30, font.size=14)
dev.off()