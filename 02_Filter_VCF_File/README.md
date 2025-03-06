# 02.1_Filter_VCF_Files

This folder contains all the standalone scripts and parameters to run the snakemake on several samples to filter a vcf file on several parameters and do some stats on the filtration steps.


The folder contains one sub_folder, two scripts and one example of a configuration file.



## [launcher.sh](./launcher.sh) (script)

This script, as its name indicates it is used as a launcher to start running the snakemake. It is used to call the scripts needed before running the snakemake, preparing the environment for it. It takes as input the configuration file (`configuration_file.yaml`) and creates the environment necessary to run the snakemake.
The second argument is the name of the Snkafile allowing to run the analysis.


## [Filter_VCF.snk](./Filter_VCF.snk) (script)

This script is the workflow. It filters a raw VCF file and returns the final filtered vcf file as well as all the vcf files after the filtration steps. It also returns some stats (Allele frequencies, sequencing depth per site, missing data, sequencing quality and the strand bias per site, and number of SNPs counted per chromosomal region).

## [Scripts_snk](./Scripts_snk/) (directory)

This folder contains the scripts and functions that are used as accessories in the launcher and the snakemake during their execution.
There are three files in this folder containing the scripts for these functions.

### [Configuration.sh](./Scripts_snk/Configuration.sh) (script)

This script allows to parse the configuration file so as to prepare the run of the snakemake.

### [snakemake_functions.py](./Scripts_snk/snakemake_functions.py) (script)

This script contains custom functions that allow to run the snakemake and parallelise it at best.

### [Graph_quanlity.r](./Scripts_snk/Graph_quality.r) (script)

This is an R script that allows to plot the quality of the VCF file after a filtration process. 

### [Table_maker_SP.r](./Scripts_snk/Table_maker_SP.r) (script)

This is an R script that makes a list of the positions to keep in the VCF file according to a strand bias threshold value indicated by the user in the configuration file.

### [Filter_Hobs.r](./Scripts_snk/Filter_Hobs.r) (script)

This is an R script that makes a list of the positions to keep to remove paralogs. To do this, we filter the VCF file on the observed heterozygosity. The used Hobs threshold is 55%, but, this value can be changed by the user in the configuration file.


## [configuration_file.yaml](./configuration_file.yaml) (configuration)

This file contains all the configurations you need to change to adapt the workflow to your data. __Please go through this file and change the paths and parameters to match your directories/clusters parameters!__
This file is structured in a particular way. It is separated into two parts with the use of:
```
---------------------------------------------------------
-----------------------  Profile  -----------------------
---------------------------------------------------------
name: Profile
```
These flags (with the "----") are used to parse this file so please keep them. In addition, the `name` tag is reserved to name the configuration files that are used by the snakemake. The two ones here are the only necessary ones. Of course, if you want to add some configuration files in the snakemake, you are welcome to add some flags in the `configuration_file.yaml` with a new `name` tag.
__There are 5 paths to change in the "Profile" part of the configuration file for you to run the snakemake.__
## How to run the script ?

1. Change the paths in the places that are indicated before in the `README.md` file.
1. Finally, you can place yourself in this terminal and type:
```
sbatch launcher.sh -f configuration_file.yaml -s Filter_VCF.snk
```
