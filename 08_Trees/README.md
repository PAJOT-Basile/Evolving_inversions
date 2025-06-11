# 08_Trees


## What is in the directory?

This directory contains the scripts used to make some trees on the inversions (or other parts of the genome).

## How to use these scripts?

To run the scripts, use the `launcher.sh` script to run the analysis in the background:
```
sbatch launcher.sh
```


## What output to expect?

This script will produce some png figures of the trees and a table in which the divergence between ecotypes is calculated.

## Dependencies


```
  - vcftools = 0.1.16
  - bcftools = 1.9

  - R             # v.3.3.1 at least
    - tidyverse
    - vcfR
    - adegenet
    - ggpubr
    - phangorn
    - reshape2
    - ggtree
    - MetBrewer
```
