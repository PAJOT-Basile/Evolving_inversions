################## Libraries  ##################
libraries <- c("tidyverse", "vcfR", "adegenet", "StAMPP")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(characters = libraries, character.only = TRUE)
rm(libraries)


################## Useful function  ##################

################## Import vcf  ##################
# data <- read.vcfR("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz") %>% 
#   vcfR2genind()

# data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 7))[, 3] %>% str_remove_all("n") %>% as.factor

# ################## Compute genome wide for all the population  ##################
# # First, import the VCF file as a genlight object
# data2 <- read.vcfR("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz") %>% 
#   vcfR2genlight()

# data2@pop <- (data2@ind.names %>% str_split_fixed(., "_", 7))[, 3] %>% str_remove_all("n") %>% as.factor
# Fst <- stamppFst(data2, nboots = 10000, percent = 95, nclusters = 3)

# Fst$Fsts %>% 
#   write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Population_wide_Fst.tsv",
#               sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
# Fst$Pvalues %>% 
#   write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Population_wide_Pvalues.tsv",
#               sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

################## Re-sample a vcf file with the reference individuals  ##################
# system2("/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools",
#         args = paste0("--gzvcf /shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz",
#                       " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/France_exposed.txt",
#                       " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/France_sheltered.txt",
#                       " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/Sweden_exposed.txt",
#                       " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/Sweden_sheltered.txt",
#                       " --recode",
#                       " --stdout | gzip -c > /shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Reference_indivs.vcf.gz"))

# First, import the VCF file as a genlight object
data3 <- read.vcfR("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Reference_indivs.vcf.gz") %>% 
  vcfR2genlight()

ref_indivs <- read.table("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/France_exposed.txt", header = FALSE,
                         sep = "\t", col.names = "Sample_Name") %>% 
  mutate(Exposition = "Fr_Exp") %>% 
  rbind(read.table("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/France_sheltered.txt", header = FALSE,
                   sep = "\t", col.names = "Sample_Name") %>% 
          mutate(Exposition = "Fr_Shl")) %>% 
  rbind(read.table("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/Sweden_exposed.txt", header = FALSE,
                   sep = "\t", col.names = "Sample_Name") %>% 
          mutate(Exposition = "Sw_Exp")) %>% 
  rbind(read.table("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/Sweden_sheltered.txt", header = FALSE,
                   sep = "\t", col.names = "Sample_Name") %>% 
          mutate(Exposition = "Sw_Shl"))

data3@pop <- data3@ind.names %>%
  as.data.frame %>% 
  rename(Sample_Name = ".") %>% 
  left_join(ref_indivs, by = "Sample_Name") %>% 
  pull(Exposition) %>% 
  as.factor


Fst_ecotypes <- stamppFst(data3, nboots = 10000, percent = 95, nclusters = 3)

Fst_ecotypes$Fsts %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Ecotype_Fst.tsv",
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

Fst_ecotypes$Pvalues %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Ecotype_Pvalues.tsv",
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
