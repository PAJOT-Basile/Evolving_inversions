# Install libraries if needed and load them
libraries <- c("tidyverse", "adegenet", "pegas", "vcfR", "argparse")
if (!require("pacman")) install.packages("pacman")
for (lib in libraries){
  pacman::p_load(lib, character.only = TRUE)
}

############################ Parse and use arguments ##########################
# Take inputs from the snakemake program
parser <- ArgumentParser(description = "This program is used to filter out the Hobs")

# Add the arguments that are used
parser$add_argument("--input", "-i", help = "The input folder (where to find the folders containing all the stats outputs)")
parser$add_argument("--output", "-o", help = "The path and name to give to the output file. A png extension will be added")
parser$add_argument("--Hobs", "-H", help = "The value of heterozygosity to use as theshold")
parser$add_argument("--tmp", "-t", help = "The temporary folder")

xargs <- parser$parse_args()

input_path <- xargs$input
output_path <- xargs$output
Hobs_thres <- xargs$Hobs
tmp <- xargs$tmp

############################ Define variables ##########################
# VCFtools program
vcftools <- "/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools"
# Get the name of the output directory
outdir__ <- sub("[^/]+$", "", output)
# get the name of the region that is to be used
region_name <- (output %> %>%
    str_split_fixed(., "/", 9))[, 9] %>%
    str_remove_all("VCF_File_|.vcf.gz")
############################ Import vcf file ##########################
data <- read.vcfR(input) %>% 
  vcfR2genind()

# We use the summary function to calculate the Hobs of the data we imported
data_info <- summary(data)
Hobs <- data_info$Hobs

############################ Filter on the Observed heterozygosity ##########################
# We make a list of all the positions that have a higher Hobs than 55% that we use to filter the filtered and thinned vcf
cbind(
  (names(Hobs[which(Hobs < Hobs_thres)]) %>% str_split_fixed(., "_", 4))[, 1:3],
  (names(Hobs[which(Hobs < Hobs_thres)]) %>% str_split_fixed(., "_", 4))[, 4]
  ) %>% 
  as.data.frame() %>% 
  rename(Chrom = V1,
         Super = V2,
         Super_frag = V3,
         Position_chrom = V4) %>% 
  mutate(Where = paste(Chrom, Super, Super_frag, sep="_")) %>% 
  select(Where, Position_chrom) %>% 
  write.table(paste0(outdir__, "Filter_vcf_Hobs", region_name, ".tsv"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")



system2(vcftools,
        args = paste0("--gzvcf ",
                      vcf_file,
                      " --positions", outdir__, "Filter_vcf_Hobs.tsv",
                      " --recode",
                      " --stdout",
                      " --temp ", tmp,
                      " | gzip -c > "output))
