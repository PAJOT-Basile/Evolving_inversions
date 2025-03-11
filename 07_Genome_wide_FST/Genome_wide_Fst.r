################## Libraries  ##################
libraries <- c("tidyverse", "vcfR", "adegenet", "StAMPP")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(characters = libraries, character.only = TRUE)
rm(libraries)

################## Sample a vcf file with the reference individuals  ##################
 system2("/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools",
         args = paste0("--gzvcf /shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz",
                       " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/France_exposed.txt",
                       " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/France_sheltered.txt",
                       " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/Sweden_exposed.txt",
                       " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/Sweden_sheltered.txt",
                       " --recode",
                       " --stdout | gzip -c > /shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Reference_indivs.vcf.gz"))

# First, import the VCF file as a genlight object
data <- read.vcfR("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Reference_indivs.vcf.gz") %>% 
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

data@pop <- data@ind.names %>%
  as.data.frame %>% 
  rename(Sample_Name = ".") %>% 
  left_join(ref_indivs, by = "Sample_Name") %>% 
  pull(Exposition) %>% 
  as.factor

############### Compute Fst between reference individuals ###############
Fst_ecotypes <- stamppFst(data, nboots = 10000, percent = 95, nclusters = 3)

Fst_ecotypes$Fsts %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Ecotype_Fst.tsv",
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

Fst_ecotypes$Pvalues %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Ecotype_Pvalues.tsv",
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
