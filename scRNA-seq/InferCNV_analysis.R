setwd('/mnt/home/koiralas/infercnv/run_tests/march_7_mods/2k_hvg_subset')

library(Seurat)
raw_counts_2k_hvg <- FindVariableFeatures(CT.integrated, selection.method = "vst", nfeatures = 2000)
#load('./raw_counts_2k_hvg.RData')

library(infercnv)
Control <- subset(x = CT.integrated, subset = sample.type == "CT")
#MODS <- subset(x = CT.integrated, subset = sample.type == "MODS")

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_counts_2k_hvg,
                                    annotations_file='./filtered_ann_ct_mods_cell_type.txt',
                                    delim='\t',
                                    gene_order_file='./gencode.txt',
                                    ref_group_names= Control$Cell.type)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir='./ct_m_cell_type_results/', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)
