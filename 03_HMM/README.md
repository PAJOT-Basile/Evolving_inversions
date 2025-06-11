# 03_ HMM


## What is in the directory?

This directory contains the scripts needed to run a Hidden Markow Model (HMM) to determine islands of high differentiation. The code to run the HMM was adapted from [Hofer et al, 2012](https://github.com/marqueda/HMM-detection-of-genomic-islands/tree/master), [Soria-Carrasco et al, 2014](https://pubmed.ncbi.nlm.nih.gov/24833390/) and [Marques et al, 2016](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13774). In addition, we chose to implement this code into a snakemake to run the analyses on multiple populations at one on a cluster. 

There is also a configuration file (`configuration_file.yaml`) that is used to prepare the environment for the Snakemake to run **as well as variables necessary for the Snakemkae to work**.


This directory also contains a directory called `Scripts_snk` that contains custom scripts to prepare the conda ennvironments and the configuration files for the Snakemake to run on its own. It also contains Rscripts used to filter VCFs, and plot stats on the VCF files produced after each filtration step. The last script contains custom python functions that are called in the Snakemake to produce the workflow.


The snakefile requires a list of individuals that are in each population. The indications for these files are the same as the ones given previously for the `Pop_map` directory (in the [previous folder](https://github.com/PAJOT-Basile/Evolving_inversions/tree/main/02_Filter_VCF_File)): Each file has to contain one name of sample per line. In addition:
1. The file name will be used to name the output directory containing the results for the population.
2. The path to the directory in which the population files are located is to be modified in the configuation file.

## How to use these scripts?

To use these scripts, there are several steps. 
1. First, modify the `configuration_file.yaml` file so the paths correspond to some paths you can use on your machine. Full paths are prefered, but relative paths can also be used.
Modify the variables to optimise the use of Snakemake.

2. Then, type:
```
sbatch launcher.sh -f configuration_file.yaml -s HMM_prep_and_execution.snk
```

If you are lost, you can display an assistance message with: 
```
sbatch launcher.sh -h
```

## What output to expect?

This script will create new directories in the output path that was given in the `configuration_file.yaml`. The files that appear in these directories are the final output. The script will also create directories in the given `tmp_path`, but these files will not be conserved .They will dissapear before the end of the Snakemake as they are defined as temporary files.

## Dependencies


```
- vcftools = 0.1.16
  
- python3:
    - os

- snakemake     # v.7.25.0 with a conda environment

- R             # v.3.3.1 at least
    - pacman
    - tidyverse
    - adegenet
    - pegas
    - vcfR
    - readxl
    - HiddenMarkov
    - foreach
    - doParallel
    - argparse
```
