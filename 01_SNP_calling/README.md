# 01_SNP_calling

## What is in the directory?

This directory contains the scripts needed to do the SNP calling. It contains a launcher that is used to start running the Snakemake [(Mölder et al, 2021)](https://doi.org/10.12688/f1000research.29032.2), the snakefile that is called `Processing_of_short-read_data.snk`. This file contains the workflow for each of the samples. 
There is also a configuration file that is used to prepare the environment for the Snakemake to run **as well as variables necessary for the Snakemkae to work**.
Finally, there is a directory calle `Scripts_snk` that contains custom scripts to prepare the conda ennvironments and the configuration files for the Snakemake to run on its own. It also contains an Rscript used to plot stats on the VCF file once the SNP calling is done and a script containing custom python functions that are called in the Snakemake.

## How to use these scripts?

To use these scripts, there are several steps. 
1. First, modify the `configuration_file.yaml` file so the paths correspond to some paths you can use on your machine. Full paths are prefered, but relative paths can also be used.
Modify the variables to optimise the use of Snakemake.

2. Then, type:
```
sbatch launcher.sh -f configuration_file.yaml -s Processing_of_short-read_data.snk
```

If you are lost, you can display an assistance message with: 
```
sbatch launcher.sh -h
```

## What output to expect?

This script will create new directories in the output path that was given in the `configuration_file.yaml`. The files that appear in these directories are the final output. The script will also create directories in the given `tmp_path`, but these files will not be conserved .They will dissapear before the end of the Snakemake as they are defined as temporary files.

## Dependencies

This uses Snakemake v.7.25.0 [(Mölder et al, 2021)](https://doi.org/10.12688/f1000research.29032.2) and relies on a conda environment.
Additionally, R v 3.3.1 and a few R packages are needed: `tidyverse, ggpubr, argparse`. Their installation should be automated, but the conda environment creation can sometimes get blocked if the packages are not installed in the right version of R. Should this happen, load the packages by hand before executing the workflow.
Finally, the `math` package is necessary in python3.
