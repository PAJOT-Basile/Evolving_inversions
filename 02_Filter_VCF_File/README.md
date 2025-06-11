# 02_Filter_VCF_File

## What is in the directory?

This directory contains the scripts needed to filter the VCF file obtained in the previous step. It contains a launcher that is used to start running the Snakemake [(Mölder et al, 2021)](https://doi.org/10.12688/f1000research.29032.2), the snakefile that is called `Filter_VCF.snk`. This file contains a workflow that is executed for each of the 142 VCF files obtained in the separate 142 genomic regions produced in the previous step. 
There is also a configuration file (`configuration_file.yaml`) that is used to prepare the environment for the Snakemake to run **as well as variables necessary for the Snakemkae to work**.


This directory also contains a directory called `Scripts_snk` that contains custom scripts to prepare the conda ennvironments and the configuration files for the Snakemake to run on its own. It also contains Rscripts used to filter VCFs, and plot stats on the VCF files produced after each filtration step. The last script contains custom python functions that are called in the Snakemake to produce the workflow.


The last directory called `Pop_map` is required to run the snakemake and must contain a text file in which there is a list of all the indiviudals to keep from the raw vcf files. If all the individuals are required, add the names of all individuals in this file. There must be one name per line.

## How to use these scripts?

To use these scripts, there are several steps. 
1. First, modify the `configuration_file.yaml` file so the paths correspond to some paths you can use on your machine. Full paths are prefered, but relative paths can also be used.
Modify the variables to optimise the use of Snakemake.

2. Then, type:
```
sbatch launcher.sh -f configuration_file.yaml -s Filter_VCF.snk
```

If you are lost, you can display an assistance message with: 
```
sbatch launcher.sh -h
```

## What output to expect?

This script will create new directories in the output path that was given in the `configuration_file.yaml`. The files that appear in these directories are the final output. The script will also create directories in the given `tmp_path`, but these files will not be conserved .They will dissapear before the end of the Snakemake as they are defined as temporary files.

## Dependencies


```
- bcftools = 1.16
- vcftools = 0.1.16

- python3:
    - math
    - pandas
    - os
    - numpy

- snakemake     # v.8.9.0 with a conda environment

- R             # v.3.3.1 at least
    - tidyverse
    - adegenet
    - pegas
    - vcfR
    - argparse
    - ggpubr
    - pacman

```
