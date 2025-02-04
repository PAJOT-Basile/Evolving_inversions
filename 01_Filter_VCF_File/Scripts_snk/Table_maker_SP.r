# Install libraries if needed and load them
libraries <- c("tidyverse", "argparse")
if (!require("pacman")) install.packages("pacman")
for (lib in libraries){
  pacman::p_load(lib, character.only = TRUE)
}

# Take inputs from the snakemake program
parser <- ArgumentParser(description = "This program is used to create a table of sites to keep on the filtration on Strand Bias")

# Add the arguments that are used
parser$add_argument("--input", "-i", help= "The input file of this script is the name of the stats file that is calculated from the previous step")
parser$add_argument("--output", "-o", help= "The name that you want to give to the output file (the path)")
parser$add_argument("--cutoff", "-c", help= "Value to use as a filtering threshold on the Strand Bias (SP) parameter", type= "double")

xargs <- parser$parse_args()

read.table(xargs$input, col.names = c("CHROM", "POS", "SP")) %>%
    filter(SP < xargs$cutoff) %>%
    select(-SP) %>%
    write.table(xargs$output, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
