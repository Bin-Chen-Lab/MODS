###Compute geometic mean and AUC representing the risk score###
#ECMO_MODS_gene_ids_base: dataframe contains the information of Key genes expression in all the ECMO and MODS patients at baseline
#Columns:samples
#Rows: genes
load("~/Desktop/ECMO_Vs_MODS_data/data/plos_key_gene_exp.RData")
load("~/Desktop/ECMO_Vs_MODS_data/data/ECMO_MODS_gene_ids_base.RData")
ECMO_MODS_gene_ids_base[nrow(ECMO_MODS_gene_ids_base)+ 1,] <- apply(ECMO_MODS_gene_ids_base, 2, function(x){prod(x)^(1/length(x))})
geomean_ECMO_MODS <- ECMO_MODS_gene_ids_base[14,]
geomean_ECMO_MODS[,1:4] <- NULL
rownames(geomean_ECMO_MODS) <- "geomean"
geomean_ECMO_MODS1 <- as.data.frame(t(geomean_ECMO_MODS))
write.csv(geomean_ECMO_MODS1, file="geomean_ECMO_MODS.csv", row.names = T)
geomean_ECMO_MODS1$Survival <- c(rep(1, times=6), rep(0, times=17))

#compute geometic mean for Plos Medicine data
#plos_key_gene_exp: dataframe contains key gene expression in MODS and no-MODS patients at base line
#Columns: genes
#Rows: sample ids
plos_key_gene_exp$geomean <- apply(plos_key_gene_exp, 1, function(x){prod(x)^(1/length(x))})
plos_key_geomean <- as.data.frame(plos_key_gene_exp[,-c(1:10)], row.names = row.names(plos_key_gene_exp))
colnames(plos_key_geomean) <- "prediction"
plos_key_geomean$Survival <- c(rep(1, times=11), rep(0, times=15))

#Compute AUC for ECMO and MODS data
library(pROC)
ECMO_MODS_roc_obj <- roc(geomean_ECMO_MODS1$Survival, geomean_ECMO_MODS1$geomean)
ECMO_MODS_auc <- auc(ECMO_MODS_roc_obj)
ECMO_MODS_roc_df <- data.frame(
  TPR=rev(ECMO_MODS_roc_obj$sensitivities), 
  FPR=rev(1 - ECMO_MODS_roc_obj$specificities))

#Compute AUC for validation data
roc_obj <- roc(plos_key_geomean$Survival, plos_key_geomean$prediction)
auc <- auc(roc_obj)
roc_df <- data.frame(
  TPR=rev(roc_obj$sensitivities), 
  FPR=rev(1 - roc_obj$specificities))

#generate AUC Plot for our data and validation data 
png(filename = "~/Desktop/ECMO_Vs_MODS_data/AUC_Plot.png")
plot(roc_df$FPR, y=roc_df$TPR, type = "l", col = "blue", xlab = "FPR", ylab = "TPR",font=2, main = "ECMO_MODS")
lines(ECMO_MODS_roc_df$FPR, y=ECMO_MODS_roc_df$TPR, col = "red", xlab = "FPR", ylab = "TPR",font=2, main = "ECMO_MODS")
axis(1, seq(0.0,1.0,0.1),lwd = 3, cex.axis=1, font = 2)
axis(2, seq(0.0,1.0,0.1),lwd = 3, cex.axis=1, font = 2)
abline(a=0.0, b=1.0)
text(x=0.5,y=0.3, label= "ECMO Vs MODS, AUC=0.99", cex=0.7, adj = 0.1, font=2)
text(x=0.4,y=0.3, label="--", cex=1, adj=0, col="red")
text(x=0.5,y=0.25, label= "Validation data, AUC=0.80", cex=0.7, adj = 0.1, font=2)
text(x=0.4,y=0.25, label="--", cex=1, adj=0, col="blue")
box(lwd = 3) 
dev.off()