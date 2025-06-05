# 05_Linkage_disequilibrium


## What is in the directory?

This directory contains scripts to compute the linkage disequilibrium (LD) along each chromosome of each population. Some scripts are used to run the LD calculations (statistical LD) in parallel using snakemake. The script `LD_indivs.r` is used to calculate a hand LD as such: LD = percentage of individuals that have the same karyotype on several inversions. For example, if all individuals that have the most common karyotype in the exposed part on chromosome 3 also have the most common inversion karyotype in the exposed part on chromosome 11, we say that they have a LD of 1. 

The rest of the scripts are used in a snakemake to run the LD computations in parallel for all populations and all chromosomes. This snakemake thus requires that a `Chromosome` directory be specified with a file containing the names of all the chromosomes and a `Pop_map` directory with the names of the individuals in the different populations. For this directory, see section [02](../02_Filter_VCF_File/) for indications on how to fill the `Pop_map` directory.

## How to use these scripts?

To use these scripts, there are several steps. 
1. First, modify the `configuration_file.yaml` file so the paths correspond to some paths you can use on your machine. Full paths are prefered, but relative paths can also be used.
Modify the variables to optimise the use of Snakemake.

2. Then, type:
```
sbatch launcher.sh -f configuration_file.yaml -s Comuting_LD_by_pop.snk
```

If you are lost, you can display an assistance message with: 
```
sbatch launcher.sh -h
```

To run the last script (`LD_indivs.r`), open it using [RStudio](https://posit.co/) and execute it line by line. 


## What output to expect?

This script will create new directories in the output path that was given in the `configuration_file.yaml`. The files that appear in these directories are the final output. The script will also create directories in the given `tmp_path`, but these files will not be conserved .They will dissapear before the end of the Snakemake as they are defined as temporary files.

## Dependencies


```
  - bcftools = 1.16
  - vcftools = 0.1.16

  - R             # v.3.3.1 at least
    - pacman
    - tidyverse
    - argparse
```
