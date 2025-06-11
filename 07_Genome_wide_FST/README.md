# 07_Genome_wide_FST


## What is in the directory?

This directory contains the script used to compute the genome wide FST

## How to use these scripts?

To run the script, open it using [RStudio](https://posit.co/) and execute it line by line. 


## What output to expect?

This script will a table in which the pairwise Fst is calculated between the ecotypes of each population and a second one is produced containing the p-values of these Fst values.

## Dependencies


```
  - vcftools = 0.1.16

  - R             # v.3.3.1 at least
    - pacman
    - tidyverse
    - adegenet
    - vcfR
    - StAMPP
```
