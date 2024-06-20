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
1. Download the Cell ranger process data from GSE269751.
2. Import all the hd5 files into R
3. Process the data using the scRNA-seq_processing.R. It will process all the files and integrate it into one seurat object with the cell annotation. It also provide the number of differentially expressed genes in MODS for each cell type. 
4. Bulk_gene_Expression_cell_cycle_analysis.R help to visualized he upregulated genes from our bulk RNA-seq data and estimate the cell cycle phase for all the cell types.
5. Cell_type_DE_genes_pathways.R helps to perform the DE analysis for each cell type and provide the barplot for upregulated and downregulated genes.
6. Cell_subtype_cycle_analysis.R code will generate the cell type specific cell cycle phase distribution and generate a figure for percentage of cells with G1, G2M and S phase.
7. Cell_cell_communication.R will help to perform the cell cell interactions analysis and generate all the figures.
8. InferCNV_analysis.R generate the InferCNV data.
9. Surface protein directory has three code to estimate the cell surface protein aboundance in control and MODS and compute their differential expression in MODS.
10. Code_to_visualization_of_plot directory has five codes for visualization of results. The plot_boxplot, plot_voilinplot, and plot_dotplot can be used to compare two groups, while Boxplot_multilevel.R can be used to compare more than two groups. Addtionally, Code_for_matrix_dotplot.R can be used to plot the matrix dotplot such as correlation plot. 
