#! /bin/bash

#SBATCH --mem=50GB
#SBATCH --cpus-per-task=15
#SBATCH --partition=long

# Here, we load the r module
module load r/4.3.1
# And execute the HMM_log10FST+1_3norm.R script
Rscript ../HMM/HMM_log10FST+1_3norm.R \
    /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_whole_genome/Data/Whole_genome_sweden_freqs.csv \
    15
