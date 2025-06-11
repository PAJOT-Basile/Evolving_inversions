# 11_Stats


## What is in the directory?

This directory contains the scripts used to do some stats on the genetic data and the PCA output. There is no particular order in which to execute the scripts.

`Hobs_triangle.r`: this script contains the code to represent the mean heterozygosity on the genome for each individual to have an idea of the population structure (are there hybrids? Backcrosses? ...).

`Correlation_delta_freqs_Fr_Sw.r`: this script contains the code to represent the correlation in differentiation between locations (Sweden vs France) and see how they are linked to the pca outputs.

`Hobs_inversions.r`: this script contains code to compute and start to represent the heterozygosity in the inversions using the genetic data.

## How to use these scripts?

To run the analyses, open [RStudio](https://www.posit.co) and execute the scripts line by line.


## What output to expect?

`Hobs_triangle.r`: this first script produces a graph that is not saved automatically.

`Correlation_delta_freqs_Fr_Sw.r`: this script produces a graph that is not saved automatically.

`Hobs_inversions.r`: this script produces two graphs that are saved automatically.


## Dependencies


```
- R             # v.3.3.1 at least
  - tidyverse
  - adegenet
  - vcfR
  - readxl
  - ggforce
  - ggh4x
  - see
```
