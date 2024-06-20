
## Cell cell coomunications
#CT.integrated= readRDS("CT_MODS_after_UMAP.RDS")
#setwd("data/MODS_scRNA/GSE171524/")
library(CellChat)
library(patchwork)
library(data.table)
library(Matrix)
library(Seurat)
library(ComplexHeatmap)
CT.integrated= readRDS("CT_MODS_manual_annotation.RDS")
cell_names= data.frame("cells"=Idents(CT.integrated))
CT.integrated@meta.data$Cell.type=cell_names$cells
set.seed(1234)
Control <- subset(x = CT.integrated, subset = sample.type == "CT")
MODS <- subset(x = CT.integrated, subset = sample.type == "MODS")

MODS_exp= MODS@assays[["RNA"]]@counts
#MODS_exp=fread("MODS_Gene_Count_per_Cell.tsv")
#MODS_exp= data.frame(MODS_exp)
# rownames(MODS_exp)= MODS_exp[,1]
# MODS_exp= MODS_exp[,-1]
#MODS_Norm_exp=NormalizeData(MODS_exp)
#data=Matrix(MODS_Norm_exp, sparse = T)
MODS_meta= MODS@meta.data
names(MODS_meta)[9]= "cells"
MODS_meta$cells <- droplevels(MODS_meta$cells)
MODS_meta$cells=factor(MODS_meta$cells)
# MODS_meta= data.frame(MODS_meta)
# rownames(MODS_meta)=MODS_meta[,1]

#load(url("https://ndownloader.figshare.com/files/25950872"))
#data.HS = data_humanSkin$data

#data.input = data_humanSkin$data # normalized data matrix
#meta = data_humanSkin$meta # a dataframe with rownames containing cell mata data
#for MODS data
#cell.use = meta$cells[meta$case == "MODS"] # extract the cell names from disease data

# Prepare input data for CelChat analysis
#data.input = data[, cell.use]
#MODS_meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
#MODS_meta$cells = droplevels(MODS_meta$cells, exclude = setdiff(levels(MODS_meta$cells),unique(MODS_meta$cells)))

MODS_cellchat <- createCellChat(object = MODS_exp, meta = MODS_meta, group.by = "cells")

### Add cell meta information to cell chat
#MODS_cellchat <- addMeta(MODS_cellchat, meta = MODS_meta)
MODS_cellchat <- setIdent(MODS_cellchat, ident.use = "cells") # set "labels" as default cell identity
levels(MODS_cellchat@idents) # show factor levels of the cell labels


CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
MODS_cellchat@DB <- CellChatDB.use
MODS_cellchat <- subsetData(MODS_cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = 4)
MODS_cellchat <- identifyOverExpressedGenes(MODS_cellchat)
MODS_cellchat <- identifyOverExpressedInteractions(MODS_cellchat)
# project gene expression data onto PPI network (optional)
MODS_cellchat <- projectData(MODS_cellchat, PPI.human)


MODS_cellchat <- computeCommunProb(MODS_cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
MODS_cellchat <- filterCommunication(MODS_cellchat, min.cells = 10)
#Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(MODS_cellchat)

#Infer the cell-cell communication at a signaling pathway level
MODS_cellchat <- computeCommunProbPathway(MODS_cellchat)

#Calculate the aggregated cell-cell communication network
MODS_cellchat <- aggregateNet(MODS_cellchat)
groupSize <- as.numeric(table(MODS_cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("MODS_intereaction.pdf")
netVisual_circle(MODS_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(MODS_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

#edge weights between differet networks
mat <- MODS_cellchat@net$weight
par(mfrow = c(1,4), xpd=TRUE)
pdf("MODS_indiv_cells_intereaction.pdf")
for (i in 1:nrow(mat)) {
  mat2 <- matrix(1, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()

MODS_cellchat <- netAnalysis_computeCentrality(MODS_cellchat, slot.name = "netP")

save(MODS_cellchat, file = "MODS_cellchat.RData")

#CT_cellchat <- netAnalysis_computeCentrality(MODS_cellchat, slot.name = "netP")







### For control data

CT_exp=Control@assays[["RNA"]]@counts
#CT_exp= data.frame(CT_exp)
# rownames(CT_exp)= CT_exp[,1]
# CT_exp= CT_exp[,-1]
options(future.globals.maxSize = 8000 * 1024^2)
#CT_Norm_exp=NormalizeData(CT_exp)
#data=Matrix(MODS_Norm_exp, sparse = T)
library(dplyr)
#CT_norm_final=CT_Norm_exp[,which(colnames(CT_exp) %in% CT_meta$cells)]
CT_meta= Control@meta.data
names(CT_meta)[9]= "cells"
CT_meta= data.frame(CT_meta)
CT_meta$cells= factor(CT_meta$cells)
# rownames(CT_meta)=CT_meta[,1]

#CT_cell.use = meta$cells[meta$case == "CT"] # extract the cell names from disease data

# Prepare input data for CelChat analysis
#CT_data.input = data[, CT_cell.use]
#CT_meta = meta[CT_cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels

CT_cellchat <- createCellChat(object = CT_exp, meta = CT_meta, group.by = "cells")

### Add cell meta information to cell chat
#CT_cellchat <- addMeta(CT_cellchat, meta = CT_meta)
CT_cellchat <- setIdent(CT_cellchat, ident.use = "cells") # set "labels" as default cell identity
levels(CT_cellchat@idents)

CT_CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CT_CellChatDB.use <- subsetDB(CT_CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
CT_cellchat@DB <- CT_CellChatDB.use
CT_cellchat <- subsetData(CT_cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = 4)
CT_cellchat <- identifyOverExpressedGenes(CT_cellchat)
CT_cellchat <- identifyOverExpressedInteractions(CT_cellchat)
# project gene expression data onto PPI network (optional)
CT_cellchat <- projectData(CT_cellchat, PPI.human)


CT_cellchat <- computeCommunProb(CT_cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
CT_cellchat <- filterCommunication(CT_cellchat, min.cells = 10)
#Extract the inferred cellular communication network as a data frame
CT_df.net <- subsetCommunication(CT_cellchat)

#Infer the cell-cell communication at a signaling pathway level
CT_cellchat <- computeCommunProbPathway(CT_cellchat)

#Calculate the aggregated cell-cell communication network
CT_cellchat <- aggregateNet(CT_cellchat)
CT_groupSize <- as.numeric(table(CT_cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
pdf("CT_intereaction.pdf")
source("netVisual_circle.R")
netVisual_circle(CT_cellchat@net$count, vertex.weight = CT_groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(CT_cellchat@net$weight, vertex.weight = CT_groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

#edge weights between differet networks
mat <- CT_cellchat@net$weight
par(mfrow = c(1,3), xpd=TRUE)
pdf("CT_indiv_cells_intereaction.pdf")
for (i in 1:nrow(mat)) {
  mat2 <- matrix(1, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()

CT_cellchat <- netAnalysis_computeCentrality(CT_cellchat, slot.name = "netP")

##### individual group cell interaction
# Define the cell labels to lift up
group.new = levels(CT_cellchat@idents)
MODS_cellchat <- liftCellChat(MODS_cellchat, group.new)




### Compare control vs MODS interaction

object.list <- list(CT = CT_cellchat, MODS = MODS_cellchat)


cellchat_merge <- mergeCellChat(object.list, add.names = names(object.list))
####Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

ggsave("CT_MODS_intereaction_level.pdf")

#save(CT_cellchat, cellchat, file = "male_cellchat_CT_tumor.RData")

#Differential number of interactions or interaction strength among different cell populations
#where red colored edges represent increased and blue colored represents the decreased signaling in the second dataset compared to the first one.
par(mfrow = c(1,2), xpd=TRUE)
pdf("CT_MODS_interactions.pdf")
source("netVisual_diffInteraction.R")
netVisual_diffInteraction(cellchat_merge, weight.scale = T)
netVisual_diffInteraction(cellchat_merge, weight.scale = T, measure = "weight")
dev.off()

gg1 <- netVisual_heatmap(cellchat_merge)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_merge, measure = "weight")
#> Do heatmap based on a merged object
pdf("CT_MODS_intereaction_strength.pdf", width = 10, height = 6)
gg1 + gg2
dev.off()

##### Get maximum number of intereactions
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("CT_MODS_maximum_interections.pdf")
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


####Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
#net=netAnalysis_computeCentrality(object = object.list, slot.name = "netP")
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("CT_MODS_agrregated_cell_cell_intereaction.pdf")
patchwork::wrap_plots(plots = gg)
dev.off()


#### find the important cells interections
gg1 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = "GZMK.NK", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = "CD8.T", signaling.exclude = "MIF")
#gg3 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = "Airway epithelial cells", signaling.exclude = c("MIF"))
#gg4 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = "Cycling NK/T cells", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
pdf("CT_MODS_NK_CD8T_cell.pdf",width = 7, height = 4)
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


pdf("CT_MODS_plasma_cell.pdf",width = 4, height = 4)
patchwork::wrap_plots(plots = list(gg2))
dev.off()

##Identify signaling groups based on their functional similarity
cellchat_merge <- computeNetSimilarityPairwise(cellchat_merge, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat_merge <- netEmbedding(cellchat_merge, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat_merge <- netClustering(cellchat_merge, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("CT_MODS_interaction_cluster.pdf")
netVisual_embeddingPairwise(cellchat_merge, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
dev.off()

save.image("CT_MODS_analysis.Rdata")
load("CT_MODS_analysis.Rdata")

#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat_merge, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merge, mode = "comparison", stacked = F, do.stat = TRUE)
pdf("CT_MODS_signalling_pathways.pdf", width = 7, height = 4)
gg1 + gg2
dev.off()

##Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
pdf("CT_MODS_signalling_pathways_out_HM.pdf", width = 7, height = 4) 
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

#### signalling pathways going in to the cells

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
pdf("CT_MODS_signalling_pathways_in_HM.pdf", width = 7, height = 4) 
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


##### overall siganlling pathways in cells
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
pdf("CT_MODS_signalling_pathways_all_HM.pdf", width = 7, height = 4)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

netVisual_bubble(cellchat_merge, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MODS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat_merge <- identifyOverExpressedGenes(cellchat_merge, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_merge, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat_merge, net = net, datasets = "MODS",ligand.logFC = 0.1, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat_merge, net = net, datasets = "CT",ligand.logFC = -0.1, receptor.logFC = -0.1)

write.csv(net.up, "MODS_up_genes_signaling.csv")
write.csv(net.down, "MODS_down_genes_signaling.csv")
#Since the signaling genes in the net.up and net.down might be complex with multi-subunits, 
#we can do further deconvolution to obtain the individual signaling genes.


gene.up <- extractGeneSubsetFromPair(net.up, cellchat_merge)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_merge)

write.csv(gene.up, "MODS_genes_up_in_signal.csv")
write.csv(gene.down, "MODS_genes_down_in_signal.csv")

load("CT_MODS_analysis.Rdata")
##Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
# Chord diagram
par(mfrow = c(1,1), xpd=TRUE)
pdf("MODS_CT_Up_sign_pathways.pdf", width = 16, height = 8)
netVisual_chord_gene(object.list[[2]], slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 1.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()
#> Note: The first link end is drawn out of sector 'MIF'.
pdf("MODS_CT_Down_sign_pathways.pdf", width = 16, height = 8)
netVisual_chord_gene(object.list[[1]], slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 1.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()