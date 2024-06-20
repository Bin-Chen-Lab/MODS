# MODS
### MODS RNA-seq data analysis
1. Author: ramashan@msu.edu; ramashankar12@gmail.com

2. Steps to reproduce the analysis and figures:


3. The expression value of each gene in all the samples was used for differential gene expression analysis using the EdgeR package workflow provided in "https://github.com/Bin-Chen-Lab/octad_desktop/blob/master/Code/deprecated/Diff_Exp.R"
4. Create a folder and named it DGE_analysis
5. Create a Rstudio project file under DGE_analysis folder named script_code
6. Download all the code in the script_code folder
7. All the relevant files for analysis has been provided in data.zip. 
8. A absolute path "~/Desktop/ECMO_Vs_MODS_data" was used in the code. This can be modified to DGE_analysis
9. Run the script as follows:
	a. DGE_annotation_filter.R
	b. OR_calculation_and_forestplot.R
	c. PCA_plot_for_clustering_patients.R
	d. Heatmaps_and_boxplot.R
	e. Risk_score_and_AUC_plot.R


### MODS scRNA-seq data analysis
1. Download the Cell Ranger processed data from GSE269751.
2. Import all the HDF5 files into R.
3. Use the `scRNA-seq_processing.R` script to process all files and integrate them into one Seurat object with cell annotation. It also provides the number of differentially expressed genes in MODS for each cell type.
4. The `Bulk_gene_Expression_cell_cycle_analysis.R` script visualizes upregulated genes from bulk RNA-seq data and estimates the cell cycle phase for all cell types.
5. The `Cell_type_DE_genes_pathways.R` script performs differential expression (DE) analysis for each cell type and provides bar plots for upregulated and downregulated genes.
6. The `Cell_subtype_cycle_analysis.R` script generates cell type-specific cell cycle phase distributions and creates a figure showing the percentage of cells in G1, G2M, and S phases.
7. The `Cell_cell_communication.R` script performs cell-cell interaction analysis and generates all the relevant figures.
8. The `InferCNV_analysis.R` script generates the InferCNV data.
9. The Surface Protein directory contains three scripts to estimate cell surface protein abundance in control and MODS samples and compute their differential expression in MODS.
10. The `Code_to_visualization_of_plot` directory contains five scripts for visualizing results. The `plot_boxplot`, `plot_violinplot`, and `plot_dotplot` scripts compare two groups, while `Boxplot_multilevel.R` can compare more than two groups. Additionally, `Code_for_matrix_dotplot.R` can plot matrix dot plots, such as correlation plots.

# Citations:
### For RNA-seq data
Shankar R, Leimanis ML, Newbury PA, Liu K, Xing J, Nedveck D, Kort EJ, Prokop JW, Zhou G, Bachmann AS, Chen B, Rajasekaran S. Gene expression signatures identify paediatric patients with multiple organ dysfunction who require advanced life support in the intensive care unit. EBioMedicine. 2020 Dec;62:103122. doi: 10.1016/j.ebiom.2020.103122. Epub 2020 Nov 25. PMID: 33248372; PMCID: PMC7704404.

