# Biomarkers of NSCLC project
This repository contains the code files and data used for a gene expression analysis project based on the GEO dataset. The project involves processing gene expression data and performing various analyses. Below, we describe the main components of this project.

## Project Overview
In this project, we analyze gene expression data from the GEO dataset. The goal is to gain insights into gene expression patterns in different tumor types. We perform data preprocessing, create relevant CSV files, and visualize gene expression patterns.

### Code Files
* main.py- The main.py script serves as the entry point for the project. It performs the following tasks:

  1) Data Processing: Reads and processes the CSV files from the GEO dataset and cleans the data leaving only gene expressions with P-value of less than 0.05

  2) CSV File Creation:
     
      a. **combined_output.csv**: This CSV file holds information about gene expression, tumor profiles, mean expression, standard deviation, and z-scores for specific genes and samples.

      b. **gene_per_patient.csv**: In this CSV file each row represents a unique patient and each column a specific gene. The values in this matrix are binary,
     indicating the presence (1) or absence (0) of each gene's expression in the patient's sample and the last column holds the tumor type of the patient.

* plot.py- The plot.py script is designed to create box plots that visualize gene expression patterns for different tumor types based on the data in the 'combined_output.csv' file.
* tsne.py- The tsne.py script is responsible for performing dimensionality reduction and clustering on gene expression data and generating visualizations to explore the relationships between genes and tumor types.
* final_project.rmd- The final_project.rmd script is a ascript in R that includes data analysis and visualization. It performs the following:

  1) PCA analysis with gene colors and tumor labels
  2) PCA analysis with gene colors and tumor labels while omitting NA values
  3) Correlation heatmap of gene expressions that exports to PDF file
  4) PCA analysis of genes per patient CSV file
  5) PCA analysis of genes per patient CSV file while using alternative appraoch
  6) UMAP analysis of genes per patient CSV file

## Data Requirements

Before running the code in this project, you will need to obtain specific data files from the GEO dataset available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135304.
Due to their size, these data files are not included in this repository.

Please download the following files and convert them into CSV files:

1) **GSE135304_raw.txt.gz-** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135304
2) **GSE135304_series_matrix.txt.gz-** https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135304/matrix/
3) **GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz-** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558
