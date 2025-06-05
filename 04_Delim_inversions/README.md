# 04_Delim_inversions


## What is in the directory?

This directory contains scripts to detect chromosomal inversions. These scripts use as input the HMM outputs that were produced in the previous steps. These scripts were not automated so their manual execution is required using R. The scripts are ordered from 1 to 4, giving the order in which to use them.

## How to use these scripts?

To run these scripts, open them using [RStudio](https://posit.co/) and execute them line by line. 

The first script `01_Analysis_HMM.r` allows to analyse the output of the HMM and define regions of the genome that are highly differentiated (islands of differentiation). Once this is done, it writes a csv file that contains the information of the islands of high differentiation (chromosome, start position, end position, length) to use in other scripts (for example figures).
#### Warning !! The positions of the inversions were also determined by hand so starting line 220, the delimiations of other inversions may not be the same.

The second script `02_Local_PCA_inversions.r` runs a local PCA on the delimitations of the islands of high differentiation identified earlier. If we find three clusters, the island of differentiation is kept for the next step. Again, the selection of the islands of differentiation are done by visual ispection of the local PCA results. The selection process might be different in other cases.

The third script `03_Stats_inversions.r` computes basic stats (heterozygosity and mean linkage disequilibrium calculated on the inversion) on the identified islands of differentiation kept in the previous script. The linkage disequilibrium part requires to run scripts in the [next part](../05_Linkage_disequilibrium/), so you can look ahead and come back to this script later. This step produces a table that keeps only islands of differentiation that have a higher heterozygosity in the heterokaryotype than in the homokaryotypes. We call these resulting islands of differentiation inversions.

Finally, the fourth script `04_Figure_local_pca_inside_differetiation_islands.r` reruns the local PCA analyses on the two populations just in the inversions that were kept after all these steps to make a nice figure.


## What output to expect?

These scripts will create tsv files that contain the coordinates of the islands of diffeerntiation that are kept after each step. The temporary files can be deleted once all the analyses are run, but it is not automatic.

## Dependencies


```
- R             # v.3.3.1 at least
    - pacman
    - tidyverse
    - vcfR
    - adegenet
    - readxl
```
