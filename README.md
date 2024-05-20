# Parallelism in ecotype differenciation between two hybrid zones of _Littorina fabalis_


This directory contains all the scripts allowing to replicate analysis showing parallelism in ecotype differenciation between two ecotypes of _Littorina fabalis_ in two locations: Sweden and France. As the Swedish population was already studied in a [previous paper](https://academic.oup.com/evlett/advance-article/doi/10.1093/evlett/qrae014/7656805), this population is taken as a reference and always presented first.
This directory is composed of 4 folders.


## 1. [Snakemake_processins_of_short-read_data](https://github.com/PAJOT-Basile/L_fabalis/tree/main/Snakemake_processing_of_short-read_data#snakemake-processing-of-short-read-sequencing-data-from-littorina-snails)


This folder is a duplicate of the repository [Snakemake_processing_of_short-read_data](https://github.com/PAJOT-Basile/Snakemake_processing_of_short-read_data). 
It contains all the necessary scripts to run the analysis of raw data for shrot read sequencing on a high-performance computing (HPC) cluster. It allows to use at its best the ressources of the cluster to clean raw data, map it on a reference genome and extract the parts we want to use for a genomics analysis.
This folder contains all the necessary instructions to run it on any short-read sequencing data (supposing you have a reference genome to map it on).

## 2. [Phenotypic_analysis](https://github.com/PAJOT-Basile/L_fabalis/tree/main/Phenotypic_analysis#phenotypic_analysis)

This folder contains scripts allowing to replicate the analysis of the phenotypic variations along the transect for both studied locations.

## 3. [General_scripts](https://github.com/PAJOT-Basile/L_fabalis/blob/main/General_scripts/README.md#general_scripts)

This folder contains two scripts containing all the functions that we created to analyse, calculate and represent genomic data in R, in coordination with the [tidyverse](https://www.tidyverse.org/) and the adegenet [(Thibaut Jombart et al, 2008)](https://pubmed.ncbi.nlm.nih.gov/18397895/) packages.

## 4. [Report](https://github.com/PAJOT-Basile/L_fabalis/tree/main/Report#report)

This folder contains all the code required to reproduce the figures and tables used in the report.
