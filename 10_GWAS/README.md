# 10_GWAS


## What is in the directory?

This directory contains the scripts used to run the Genome Wide Association studies (GWAS) on size and colour. The scripts are ordered in the order in which they should be run.

`01_Relatedness.r`: this first script contains the code to compute the relatedness between individuals. These calculations are run separately for the two populations. These computations are used in the GWAS analyses.

`02_Run_GWAS.r`: this script contains the code to run the GWAS for the traits (size and colour).

`03_Boxplots_candidate_SNPs.r`: this script contains code to plot the distribution of the traits (size and colour) depending on the genotypes of individuals for QTL SNPs.

## How to use these scripts?

To run the analyses, open [RStudio](https://www.posit.co) and execute the scripts line by line.


## What output to expect?

`01_Relatedness.r`: this first script produces a table for each population containing the relatedness between individuals. The script also produces a graph that is not saved automatically.

`02_Run_GWAS.r`: this script produces a table containing the results of the GWAS for the two traits. The code that is used to represent these results is in the [Figures folder](../Figures/).

`03_Boxplots_candidate_SNPs.r`: this script produces two graphs (one for each trait) that are saved automatically.


## Dependencies


```
- vcftools = 0.1.16

- R             # v.3.3.1 at least
  - tidyverse
  - adegenet
  - vcfR
  - readxl
  - statgenGWAS
  - ggforce
  - ggh4x
  - ggpubr
  - reshape2
  - patchwork
```
