# Import libraries
require("anyLib")
anyLib(c("tidyverse"))


################## Run vcftools to calculate the frequency of each allele for the exposed parts of the transects  ##################
vcf_file <- "/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/filtering_vcf_files/Final_outputs/Fully_filtered.vcf.gz"
vcftools <- "/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools"

# First, for the french exposed population
system2(vcftools,
        args = paste0("--gzvcf ",
                      vcf_file,
                      " --keep /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_preparation/French_pop_exposed.txt",
                      " --freq",
                      " --out /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_whole_genome/French_expos_freqs",
                      " --maf 0.05"))

# Then, for the Swedish exposed population
system2(vcftools,
        args = paste0("--gzvcf ",
                      vcf_file,
                      " --keep /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_preparation/Swedish_pop_exposed.txt",
                      " --freq",
                      " --out /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_whole_genome/Swedish_expos_freqs",
                      " --maf 0.05"))

################## Run vcftools to calculate the frequency of each allele for the sheltered parts of the transects  ##################
# First for the French sheltered population
system2(vcftools,
        args = paste0("--gzvcf ",
                      vcf_file,
                      " --keep /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_preparation/French_pop_sheltered.txt",
                      " --freq",
                      " --out /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_whole_genome/French_shelt_freqs",
                      " --maf 0.05"))

# Then for the Swedish sheltered population
system2(vcftools,
        args = paste0("--gzvcf ",
                      vcf_file,
                      " --keep /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_preparation/Swedish_pop_sheltered.txt",
                      " --freq",
                      " --out /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_whole_genome/Swedish_shelt_freqs",
                      " --maf 0.05"))


################## Import created files  ##################
# Now that we calculated the allelic frequencies for each allele using VCFtools, we have to calculate the differences in allelic
# frequencies for each SNP to run the HMM. Therefore, we import each table we just created
fr_exp_freqs <- read_csv("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_whole_genome/French_expos_freqs.frq",
                     col_names = TRUE) %>% 
  separate("CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}",
           c("Chromosome", "Position", "Num_alleles", "N_indivs", "First_allele_freq", "Second_allele_freq"), "\t") %>% 
  separate(First_allele_freq, c("First_allele", "First_allele_freq"), sep = ":") %>% 
  separate(Second_allele_freq, c("Second_allele", "Second_allele_freq"), sep = ":") %>% 
  unite("Position", Chromosome, Position) %>% 
  rename(First_all_freq_expos = First_allele_freq,
         Second_all_freq_expos = Second_allele_freq) %>% 
  mutate(First_all_freq_expos = First_all_freq_expos %>% as.numeric,
         Second_all_freq_expos = Second_all_freq_expos %>% as.numeric) %>% 
  select(-c(Num_alleles, N_indivs))

fr_shelt_freqs <- read_csv("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_whole_genome/French_shelt_freqs.frq",
                         col_names = TRUE) %>% 
  separate("CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}",
           c("Chromosome", "Position", "Num_alleles", "N_indivs", "First_allele_freq", "Second_allele_freq"), "\t") %>% 
  separate(First_allele_freq, c("First_allele", "First_allele_freq"), sep = ":") %>% 
  separate(Second_allele_freq, c("Second_allele", "Second_allele_freq"), sep = ":") %>% 
  unite("Position", Chromosome, Position) %>% 
  rename(First_all_freq_shelt = First_allele_freq,
         Second_all_freq_shelt = Second_allele_freq) %>% 
  mutate(First_all_freq_shelt = First_all_freq_shelt %>% as.numeric,
         Second_all_freq_shelt = Second_all_freq_shelt %>% as.numeric) %>% 
  select(-c(Num_alleles, N_indivs))

sw_exp_freqs <- read_csv("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_whole_genome/Swedish_expos_freqs.frq",
                         col_names = TRUE) %>% 
  separate("CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}",
           c("Chromosome", "Position", "Num_alleles", "N_indivs", "First_allele_freq", "Second_allele_freq"), "\t") %>% 
  separate(First_allele_freq, c("First_allele", "First_allele_freq"), sep = ":") %>% 
  separate(Second_allele_freq, c("Second_allele", "Second_allele_freq"), sep = ":") %>% 
  unite("Position", Chromosome, Position) %>% 
  rename(First_all_freq_expos = First_allele_freq,
         Second_all_freq_expos = Second_allele_freq) %>% 
  mutate(First_all_freq_expos = First_all_freq_expos %>% as.numeric,
         Second_all_freq_expos = Second_all_freq_expos %>% as.numeric) %>% 
  select(-c(Num_alleles, N_indivs))

sw_shelt_freqs <- read_csv("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/HMM_whole_genome/Swedish_shelt_freqs.frq",
                           col_names = TRUE) %>% 
  separate("CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}",
           c("Chromosome", "Position", "Num_alleles", "N_indivs", "First_allele_freq", "Second_allele_freq"), "\t") %>% 
  separate(First_allele_freq, c("First_allele", "First_allele_freq"), sep = ":") %>% 
  separate(Second_allele_freq, c("Second_allele", "Second_allele_freq"), sep = ":") %>% 
  unite("Position", Chromosome, Position) %>% 
  rename(First_all_freq_shelt = First_allele_freq,
         Second_all_freq_shelt = Second_allele_freq) %>% 
  mutate(First_all_freq_shelt = First_all_freq_shelt %>% as.numeric,
         Second_all_freq_shelt = Second_all_freq_shelt %>% as.numeric) %>% 
  select(-c(Num_alleles, N_indivs))

################## Compute allelic frequencies  ##################
# And then, we calculate the allelic frequenciesfor all the SNPs
Total_pops_freqs <- fr_exp_freqs %>% 
  inner_join(fr_shelt_freqs) %>% 
  mutate(Population = "France") %>% 
  rbind(sw_exp_freqs %>%
          inner_join(sw_shelt_freqs) %>% 
          mutate(Population = "Sweden")) %>% 
  rename(First_all_freq_shelt = Fisrt_all_freq_shelt,
         First_all_freq_expos = Fisrt_all_freq_expos) %>% 
  mutate(Delta_freq_allele_1 = First_all_freq_expos - First_all_freq_shelt,
         Delta_freq_allele_2 = Second_all_freq_expos - Second_all_freq_shelt)


################## Export to text files  ##################
# We then separate the two populations and export the calculated delta frequencies into 
# specific tables for each location.
Total_pops_freqs  %>% 
  filter(Population == "Sweden") %>% 
  select(Position, Delta_freq_allele_1) %>% 
  mutate(x = Delta_freq_allele_1 %>% abs) %>% 
  select(-Delta_freq_allele_1) %>% 
  column_to_rownames("Position") %>% 
  write.table("./HMM/Data/Whole_genome_sweden_freqs.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

Total_pops_freqs  %>% 
  filter(Population == "France") %>% 
  select(Position, Delta_freq_allele_1) %>% 
  mutate(x = Delta_freq_allele_1 %>% abs) %>% 
  select(-Delta_freq_allele_1) %>% 
  column_to_rownames("Position") %>% 
  write.table("./HMM/Data/Whole_genome_france_freqs.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
