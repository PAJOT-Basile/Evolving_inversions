#! /bin/sh
#SBATCH --mem=100G

module load r/4.4.1
Rscript Trace_trees_inversions.R
module unload r/4.4.1
