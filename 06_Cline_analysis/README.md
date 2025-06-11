# 06_Cline_analysis


## What is in the directory?

This directory contains scripts used to run cline analyses on the genetic data. There are two parts to this analysis. First, do the cline fitting and second, analyse these fits. The cline fitting is done using the snakemake and the analysis is done with the `Cline_analysis.r` script. Thus, the Rscript takes the output of the snakemake as input.

## How to use these scripts?
### The Snakemake
To use these scripts, there are several steps. 
1. First, modify the `configuration_file.yaml` file so the paths correspond to some paths you can use on your machine. Full paths are prefered, but relative paths can also be used.
Modify the variables to optimise the use of Snakemake.

2. Then, type:
```
sbatch launcher.sh -f configuration_file.yaml -s Cline_analysis.snk
```

If you are lost, you can display an assistance message with: 
```
sbatch launcher.sh -h
```
### The Rscript
To run the last script (`Cline_analysis.r`), open it using [RStudio](https://posit.co/) and execute it line by line. 


## What output to expect?

This script will create new directories in the output path that was given in the `configuration_file.yaml`. The files that appear in these directories are the final output. The script will also create directories in the given `tmp_path`, but these files will not be conserved .They will dissapear before the end of the Snakemake as they are defined as temporary files.

The Rscript will not produce any file. It will simply print the results in the Console.

## Dependencies


```
- python = 3.12

- snakemake     # v.7.25.0 with a conda environment

- R             # v.4.4.1
  - pacman
  - tidyverse
  - adegenet
  - vcfR
  - readxl
  - StAMPP
  - dartR
  - argparse
```
