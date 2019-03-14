# MODS
MODS data analysis
Author: ramashankar@hc.msu.edu ramashankar12@gmail.com

Steps to reproduce the analysis and figures:


The expression value of each gene in all the samples was used for differential gene expression analysis using the EdgeR package workflow provided in "https://github.com/Bin-Chen-Lab/octad_desktop/blob/master/Code/deprecated/Diff_Exp.R"
Create a folder and named it DGE_analysis
Create a Rstudio project file under DGE_analysis folder named script_code
Download all the code in the script_code folder
All the relevant files for analysis has been provided in data.zip. 
A absolute path "~/Desktop/ECMO_Vs_MODS_data" was used in the code. This can be modified to DGE_analysis
Run the script as follows:
	a) DGE_annotation_filter.R
	b) OR_calculation_and_forestplot.R
	c) PCA_plot_for_clustering_patients.R
	d) Heatmaps_and_boxplot.R
	e) Risk_score_and_AUC_plot.R
