# Parallelism in ecotype differenciation between two hybrid zones of _Littorina fabalis_


This directory contains all the scripts allowing to replicate analyses to study parallel differenciation between two pairs of ecotypes of _Littorina fabalis_ in two locations: Sweden and France. As the Swedish population was already studied in a previous paper [(Le Moan et al, 2024)](https://academic.oup.com/evlett/advance-article/doi/10.1093/evlett/qrae014/7656805), this population is taken as a reference and always presented first.
This directory is composed of 14 folders.


## 1. [01_SNP_calling](https://github.com/PAJOT-Basile/L_fabalis/tree/main/01_SNP_calling)


This folder contains the scrips used to run the SNP calling from the raw fastq files. It uses Snakemake [(Mölder et al, 2021)](https://f1000research.com/articles/10-33/v1) to parallelise the process for each of the indiviudals selected in the study. The SNP calling is done on 142 non-overlapping windows along the genome that are run in parallel to go faster.


## 2. [02_Filter_VCF_File]([https://github.com/PAJOT-Basile/L_fabalis/tree/main/Phenotypic_analysis#phenotypic_analysis](https://github.com/PAJOT-Basile/L_fabalis/tree/main/02_Filter_VCF_File))

This folder contains the scripts and Pop map file to run the filtering of the raw VCF file that was obtained in step 01_SNP_calling. It uses the VCF file cut in 142 windows along the genome as input and filters each of these windows on the selected parameters. This step also uses Snakemake [(Mölder et al, 2021)](https://f1000research.com/articles/10-33/v1) to go faster. The outputs are several folders with a VCF and a corresponding stat file for each filtration step. Some of these filtration steps are defined as temporary and are suppressed automatically in the snakemake. 

## 3. [03_HMM](https://github.com/PAJOT-Basile/L_fabalis/tree/main/03_HMM)

This folder contains the scripts to run the HMM from [Hofer et al, 2012](https://github.com/marqueda/HMM-detection-of-genomic-islands/tree/master), [Soria-Carrasco et al, 2014](https://pubmed.ncbi.nlm.nih.gov/24833390/) and [Marques et al, 2016](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13774). We run this using Snakemake [(Mölder et al, 2021)](https://f1000research.com/articles/10-33/v1) to go faster and run the same analysis on the two populations separately. The script using the HMM was slightly modified to return only two levels of differenciation rather than the three levels that were used in the orginal paper.

## 4. [04_Delim_inversions](https://github.com/PAJOT-Basile/L_fabalis/tree/main/04_Delim_inversions)

This folder contains the scripts that take the output of the HMM to delimitate large chromosomal inversions (> than 1e6 bp). The first step is to use the HMM output to define large islands of differenciation (> 10e5 bp) between ecotypes separately in the two populations. Then, a local PCA is calculated on the defined island of differenciation. If we observe three clusters that are what is expected for a polymorphic chromosomal inversion in a population, we validate that the island of differentiation is an inversion. 
Finally, we calculate the heterozygosity of the SNPs within the identified inversions. If the previously identified inversions do not have a heterozygosity that is greater in the heterozygote genotype of the inversion than in the homozygote genotypes, we do not consider it to contribute to the differentiation between ecotypes.

## 5. [05_Linkage_disequilibrium](https://github.com/PAJOT-Basile/L_fabalis/tree/main/05_Linkage_disequilibrium)

This folder contains the scripts to compute the mean LD between windows of 1e5 bp. It uses Snakemake [(Mölder et al, 2021)](https://f1000research.com/articles/10-33/v1) to parallelise on all the chromosomes at once to go faster and on the different populations to run them separately.
It also contains the script to calculate the percentage of individuals that share the same rearrangements for the different inversions. For example, if an individual is homozygote for the colinear version of the inversion, we look if it is homozygote for the colinear version on the other inversions. We run this process for each individual on all the inversions. This measure allows us to approximate linkage disequilibrium (LD) between inversions on the whole genome. This measure allows to regroup inversions located on the same chromosome as one inversion if they have a LD of 1.

## 6. [06_Cline_analysis](https://github.com/PAJOT-Basile/L_fabalis/tree/main/06_Cline_analysis)

This folder contains the scripts to run the cline analyses on the highly differentiated SNPs between transects. It uses Snakemake [(Mölder et al, 2021)](https://f1000research.com/articles/10-33/v1) to parallelise the cline analyses on batches of 1000 SNPs to go faster. The cline functions that were used here were published in [Westram et al, 2018]( https://doi.org/10.1002/evl3.74), but were slightly modified to be more easily usable in my code and explanations were added to describe the function's functionning.

## 7. [07_Genome_wide_FST](https://github.com/PAJOT-Basile/L_fabalis/tree/main/07_Genome_wide_FST)

This folder contains a script to calculate the genome wide Fst between ecotypes.

## 8. [08_Trees](https://github.com/PAJOT-Basile/L_fabalis/tree/main/08_Trees)

This folder contains the scripts used to compute phylogenetic trees inside the identified inversions. A supplementary file is provided to give an example of what the phylogenetic trees look like outside of the chromosomal inversions.

## 9. [09_Cline_Rstan](https://github.com/PAJOT-Basile/L_fabalis/tree/main/09_Cline_Rstan)

This folder contains the scripts used to calculate the Fis along the transect on the genotypes of the inversions. This analysis uses Rstan and the scripts were first published in [Le Moan et al, 2024](https://doi.org/10.1093/evlett/qrae014). The Rscript was modified to be used with our data, but the stan files were not modified.

## 10. [10_GWAS](https://github.com/PAJOT-Basile/L_fabalis/tree/main/10_GWAS)

This folder contains the scripts to run the genome wide association study (GWAS) on shell size and shell color. In these analyses, the populations were treated separately. It also contains the scripts used to get the boxplots of the QTL SNPs. It also contains a script to compute the relatedness between individuals along each transect.

## 11. [11_Stats](https://github.com/PAJOT-Basile/L_fabalis/tree/main/11_Stats)

This folder contains the scripts to compute the heterozygosity of the inversion genotypes for each individual and scripts to compute the correlation between the differences in allele frequency in France and in Sweden. Finally, it contains the script to compute the mean heterozygosity calculated on all SNPs in individual.

## 11. [Figures](https://github.com/PAJOT-Basile/L_fabalis/tree/main/Figures)

This folder contains the scripts to create the figures used in the paper. The first figure is missing as it was done using Google maps and Inkscape.

## 12. [General_scripts](https://github.com/PAJOT-Basile/L_fabalis/tree/main/General_scripts)

This folder contains two scripts containing useful functions created to analyse, compute and represent genomic data in this study. Other functions were also created to be used everywhere. All the functions were created to work in coordination with the [tidyverse](https://www.tidyverse.org/) and the adegenet [(Thibaut Jombart et al, 2008)](https://pubmed.ncbi.nlm.nih.gov/18397895/) packages.
At the beginnnig of aach function, a documentation was written.

## 13. [Phenotypic_analysis](https://github.com/PAJOT-Basile/L_fabalis/tree/main/Phenotypic_analysis)

This folder contains the scripts to run the phenotypic analysis of shell size and shell color variations along the transects in Sweden and in France. It is used to compute the tables of the best model and the cline parameters described in the paper.

