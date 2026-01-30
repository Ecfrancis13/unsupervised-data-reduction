Placental Insulin-mTOR and Stress–Inflammatory Signaling Patterns are Associated with Childhood Adiposity

This repository contains the R code and sample data used for the analysis presented in " Placental Insulin-mTOR and Stress–Inflammatory Signaling Patterns are Associated with Childhood Adiposity". 
Prerequisites

To replicate this analysis, please ensure you have R installed. You can install the required packages using the `requirements.txt` file or by running the following command in R:

Example of how to install key dependencies
install.packages(c("tidyverse", "pheatmap", "ggpubr"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")

Data Input 
A sample dataset for deriving latent variables that captures patterns is provided in the /data folder: 
•	sample_data.csv: A numeric matrix of [e.g., simulated protein levels] with samples as rows and features as columns. 

Analysis Workflow 
Please execute the scripts in the following order: 
1.	0_dataPrep.R: Cleans the raw data, handles missing, and performs scaling.
	    Input: sample_data.csv
	    Output: bdf1_scalelg_rr.csv
2.	1_ConsensusClust.R: Runs ConsensusClusterPlus to determine optimal 𝑘
and cluster assignments, as well as visualizations
3.	1_PCA.R: Runs principal component analysis, determines the optimal number of PCs and includes visualizations
4.	1_WGCNA: Runs WCNA, extracts modules, and includes visualizations
5.	2_MergeData: Combines the latent patterns derived with Consensus clustering, WGCNA, and PCA with the sample meta data for association testing. It also examines patterns of missingness in meta data. 
6.	3_Associations: Runs LASSO regression on the three sets of patterns (consensus cluster, WGCNA, and PCA). It then runs generalized linear models to test associations with the primary outcome fat mass, and secondary outcomes metabolic health and BMI percentile. 

•	/scripts: Contains all R scripts listed above.
•	/data: Contains sample input files.
•	requirements.txt: List of all package versions used. 

