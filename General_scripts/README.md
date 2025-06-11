# General scripts


## What is in the directory?

This directory contains scripts with custom functions that are used in nearly every analyse. The script `Cline_functions.R` was written and published in [Le Moan *et al*., (2024)](https://doi.org/10.1093/evlett/qrae014). The functions in this script were slighly modified to be used in coordination with the other functions and scripts.

The functions in the `Functions_optimise_plot_clines.r` script were designed to be useful in different cases:
1. General useful functions that can be used in different cases in R
2. Functions that use the output format of other functions or packages to transform/filter them easily
3. Functions that use the adegenet format to compute stats (delta_freqs, *f*4, get the genotypes, ...)
4. Represent genomic output (manhattan plots). The functions are used to make a manhattan plot in different cases and represent threshold lines on these graphs.
5. Run analyses (optimise the best model to represent allele variations along transects, run local PCAs, subset genetic data to run phylogenies in parts of the genome)

## How to use these scripts?

These scripts only contain functions so they will not work on their own. They are called in the different scripts that require them so the functions can be imported.