# Predicting gene distribution in ammonia-oxidizing archaea using phylogenetic signals- code

This repository contains the scripts and data re-run the analyses of the manuscript [ _Predicting gene distribution of Ammonia oxidising archaea using phylogenetic signals_](https://doi.org/10.1093/ismeco/ycaf087). 

## Requirements
The code is written in RMarkdown and is intended to be run in RStudio. Before running the code, make sure you have the necessary R packages installed (see each .Rmd file for package requirements).

## Contents
The repository contains the following files and folders:

- `AOA_gene_predictions.Rproj`: This is the R project file that contains the structure of the repository. Opening this file in RStudio automatically sets the correct working directory.
- `01_elastic_net_predictions.Rmd`: This file contains the code to run the phylogenetic signal analyses.
- `02_AOA_case_study.Rmd`: This file contains the code to run the phylogenetic signal analyses with a different approach.
- `functions/`, `data/`, `taxonomy_genomes/`, etc: Folders containing the necessary functions and data files used in the analyses.

## Running the Code
To run the code, follow these steps:

1. Download or clone the repository.
2. Unzip it (if downloaded as a ZIP).
3. Open AOA_gene_predictions.Rproj in RStudio. This makes sure the working directory is correctly set.
4. Open and run the RMarkdown files (`01_elastic_net_predictions.Rmd` and `02_AOA_case_study.Rmd`) without modifying the folder structure. All required data will be read in as expected.

