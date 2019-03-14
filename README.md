# MODS
MODS data analysis
Author: ramashankar@hc.msu.edu ramashankar12@gmail.com

Steps to reproduce the analysis and figures

1. The expression value of each gene in all the samples was used for differential gene expression analysis using the 
   EdgeR package workflow provided in "https://github.com/Bin-Chen-Lab/octad_desktop/blob/master/Code/deprecated/Diff_Exp.R"
2. Create a folder and named it DGE_analysis
3. Create a Rstudio project file under DGE_analysis folder named script_code
4. Download all the code in the script_code folder
5. Download all the files provided in another subfolder under DGE_analysis named as data
6. A absolute path "~/Desktop/ECMO_Vs_MODS_data" was used in the code. This can be modified to DGE_analysis
7. Run the script as follows
	a) DGE_annotation_filter.R
	b) OR_calculation_and_forestplot.R
	c) PCA_plot_for_clustering_patients.R
	d) Heatmaps_and_boxplot.R
	e) Risk_score_and_AUC_plot.R
