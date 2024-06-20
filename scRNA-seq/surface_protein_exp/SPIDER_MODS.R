library(dplyr)
library(Seurat)
library(SPIDER)

load('/home/ubuntu/single_cell/MODS/CT_MODS_umap_cell_annotation.RData')

RNA_MODS <- CT.integrated[, filter(CT.integrated@meta.data, sample.type == 'MODS') %>% rownames()]
RNA_control <- CT.integrated[, filter(CT.integrated@meta.data, sample.type == 'CT') %>% rownames()]

RNA_MODS[["study"]] = RNA_MODS@meta.data$orig.ident
RNA_control[["study"]] = RNA_control@meta.data$orig.ident

SPIDER_predict ( seurat_data = RNA_MODS,
                 tissue = 'blood',
                 disease = 'MODS',
                 SPIDER_model_file_path = '/home/ubuntu/single_cell/SPIDER/SPIDER_python/SPIDER_weight/', 
                 use_cell_type = 'SingleR',
                 query_cell_type = NULL,
                 protein = 'All',
                 use_pretrain = 'T', #Using pretrained SPIDER
                 save_path = '/home/ubuntu/single_cell/MODS/SPIDER_MODS/', 
                 use_python_path = '/home/ubuntu/chenlab_deeplearning/chenlab_deeplearning_V2/anaconda3/envs/SPIDER/bin/python', 
                 scarches_path = '/home/ubuntu/single_cell/SPIDER/scarches-0.4.0/',
                 all_trainable_proteins_gene_names = NULL, 
                 file_A_B_C_matching_table = NULL ) 



