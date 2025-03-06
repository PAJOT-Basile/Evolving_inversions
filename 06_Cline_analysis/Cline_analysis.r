# Import libraries
libraries <- c("tidyverse", "adegenet", "vcfR", "readxl", "StAMPP", "dartR")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)


################## Useful functions  ##################
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

count_nb_unique_positions <- function(df){
  df %>% 
    pull(Position) %>% 
    unique %>% 
    length %>% 
    return()
}

################## Import SNPs with deltafreq > 30 %  ##################
High_delta_freqs <- read.table("/shared/home/bpajot/fabalis/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Priors/High_delta_freqs.tsv",
                               sep = "\t", header = TRUE)

################### Import the cline parameters  ##################
Params_clines_high_delta_freqs <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Cline_parameters.tsv",
                                             sep = "\t", header = TRUE) %>% 
  unique %>% 
  select_clinal_SNPs()

################### Import the inversion delimitations  ##################
Delim_inversions <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
                               sep = "\t", header = TRUE) 

################### Nb SNPs highly diff (deltafreq > 0.3)  ##################
Highly_diff_SNPs <- High_delta_freqs %>% 
  count_nb_unique_positions

# See how many are located within chromosomal inversions
Highly_diff_SNPs_invs <- High_delta_freqs %>%
  mutate(pos = Position) %>%
  transform_position_ade2tidy() %>%
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>%
  unite(Chromosome, Chromosome, SUPER, SUPER_frag, sep = "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>%
  drop_na %>%
  filter(Position > Start & Position < End) %>% 
  count_nb_unique_positions


################### Clines in at least one pop  ##################
# All SNPs making at least one cline
Clines_at_least_one_pop <- Params_clines_high_delta_freqs %>% 
  count_nb_unique_positions

# How many located in the inversions
Clines_at_least_one_pop_invs <- Params_clines_high_delta_freqs %>% 
  mutate(pos = Position) %>%
  transform_position_ade2tidy() %>%
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>%
  unite(Chromosome, Chromosome, SUPER, SUPER_frag, sep = "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>%
  drop_na %>%
  filter(Position > Start & Position < End) %>% 
  count_nb_unique_positions
  

################### Clines both pops  ##################
# Clinal SNPs in France
clines_fr <- Params_clines_high_delta_freqs %>%
  filter(Population == "France",
         Best_model == "Clinal",
         abs(Delta_AIC_second_best_model) >= 4)

# Clinal SNPs in Sweden
clines_sw <- Params_clines_high_delta_freqs %>%
  filter(Population == "Sweden",
         Best_model == "Clinal",
         abs(Delta_AIC_second_best_model) >= 4)

# SNPs with clines in France and Sweden
clines_fr_and_sw <- clines_fr %>%
  filter(Position %in% clines_sw$Position) %>%
  left_join(High_delta_freqs %>%
              select(Position, F4_stat), by = "Position")

Clines_both_locations <- clines_fr_and_sw %>% 
  count_nb_unique_positions

# How many in the inversions ?
Clines_both_locations_invs <- clines_fr_and_sw %>%
  mutate(pos = Position) %>%
  transform_position_ade2tidy() %>%
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>%
  unite(Chromosome, Chromosome, SUPER, SUPER_frag, sep = "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>%
  drop_na %>%
  filter(Position > Start & Position < End) %>% 
  count_nb_unique_positions

################### Parallel clines  ##################
# SNPs that make parallel clines
para_clines <- clines_fr_and_sw %>%
  filter(F4_stat >= 0.09)

Parallel_clines <- para_clines %>% 
  count_nb_unique_positions

# How many in inversions ?
Parallel_clines_invs <- para_clines %>%
  mutate(pos = Position) %>%
  transform_position_ade2tidy() %>%
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>%
  unite(Chromosome, Chromosome, SUPER, SUPER_frag, sep = "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>%
  drop_na %>%
  filter(Position > Start & Position < End) %>% 
  count_nb_unique_positions

################### Anti-parallel clines  ##################
# SNPs that make anti parallel clines
antipara_clines <- clines_fr_and_sw %>%
  filter(F4_stat <= -0.09)

Anti_parallel_clines <- antipara_clines %>% 
  count_nb_unique_positions

# How many in inversions
Anti_parallel_clines_invs <- antipara_clines %>%
  mutate(pos = Position) %>%
  transform_position_ade2tidy() %>%
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>%
  unite(Chromosome, Chromosome, SUPER, SUPER_frag, sep = "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>%
  drop_na %>%
  filter(Position > Start & Position < End) %>% 
  count_nb_unique_positions

################### Remaining clines in both populations  ##################
# Remaining
remaining <- clines_fr_and_sw %>%
  filter(Position %!in% para_clines$Position & Position %!in% antipara_clines$Position)

Remaining_clines <- remaining %>% 
  count_nb_unique_positions

# How many in the inversions
Remaining_clines_invs <- remaining %>%
  mutate(pos = Position) %>%
  transform_position_ade2tidy() %>%
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>%
  unite(Chromosome, Chromosome, SUPER, SUPER_frag, sep = "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>%
  drop_na %>%
  filter(Position > Start & Position < End) %>% 
  count_nb_unique_positions
################### Total private clines  ##################
# SNPs that do not make clines in two populations
all_private <- Params_clines_high_delta_freqs %>%
  filter(Position %!in% clines_fr_and_sw$Position)


Private_clines_tot <- all_private %>% 
  count_nb_unique_positions

# How many in the inversions
Private_clines_tot_invs <- all_private %>%
  mutate(pos = Position) %>%
  transform_position_ade2tidy() %>%
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>%
  unite(Chromosome, Chromosome, SUPER, SUPER_frag, sep = "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>%
  drop_na %>%
  filter(Position > Start & Position < End) %>% 
  count_nb_unique_positions
################### Clines private France  ##################
# SNPs that make clines only in France
private_fr <- Params_clines_high_delta_freqs %>%
  filter(Population == "France",
         Best_model == "Clinal",
         abs(Delta_AIC_second_best_model) >= 4,
         Position %!in% clines_sw$Position)

Private_clines_france <- private_fr %>% 
  count_nb_unique_positions

# How many in the inversions
Private_clines_france_invs <- private_fr %>%
  mutate(pos = Position) %>%
  transform_position_ade2tidy() %>%
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>%
  unite(Chromosome, Chromosome, SUPER, SUPER_frag, sep = "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>%
  drop_na %>%
  filter(Position > Start & Position < End) %>% 
  count_nb_unique_positions

################### Clines private Sweden  ##################
# SNPs that make clines only in France
private_sw <- Params_clines_high_delta_freqs %>%
  filter(Population == "Sweden",
         Best_model == "Clinal",
         abs(Delta_AIC_second_best_model) >= 4,
         Position %!in% clines_fr$Position)


Private_clines_sweden <- private_sw %>% 
  count_nb_unique_positions

# How many in inversions
Private_clines_sweden_invs <- private_sw %>%
  mutate(pos = Position) %>%
  transform_position_ade2tidy() %>%
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>%
  unite(Chromosome, Chromosome, SUPER, SUPER_frag, sep = "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>%
  drop_na %>%
  filter(Position > Start & Position < End) %>% 
  count_nb_unique_positions


###################  Summarise everything  ##################
Summary_stats <- data.frame(
  "Stat" = c("Delta_freq > 0.3", "Clines at least one pop", "Clines both pops",
             "Parallel", "Anti-parallel", "Remaining", "Total private",
             "Private_Sweden", "Private_France"),
  "Nb_SNPs" = c(Highly_diff_SNPs, Clines_at_least_one_pop, Clines_both_locations,
                Parallel_clines, Anti_parallel_clines, Remaining_clines,
                Private_clines_tot, Private_clines_sweden, Private_clines_france),
  "In_inversions" = c(Highly_diff_SNPs_invs, Clines_at_least_one_pop_invs,
                      Clines_both_locations_invs, Parallel_clines_invs,
                      Anti_parallel_clines_invs, Remaining_clines_invs,
                      Private_clines_tot_invs, Private_clines_sweden_invs,
                      Private_clines_france_invs)
)

Summary_stats %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Summary_stats.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE) 

###################  Save names of parallel and anti-parallel SNPs  ##################
Para_cline_SNPs <- para_clines %>% 
  select(Position) %>% 
  unique %>% 
  transform_position_ade2tidy() %>% 
  unite(Position, c(Chromosome, Position), sep = "_") %>% 
  pull(Position)

Antipara_cline_SNPs <- antipara_clines %>% 
  select(Position) %>% 
  unique %>% 
  transform_position_ade2tidy() %>% 
  unite(Position, c(Chromosome, Position), sep = "_") %>% 
  pull(Position)

###################  Calculate mean FST between ecotypes on parallel SNPs  ##################
# Import the vcf file
data <- read.vcfR("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/07_Genome_wide_FST/Reference_indivs.vcf.gz") %>% 
  vcfR2genind()

# Import the reference individuals
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

# Attribute a population to each individual
data@pop <- data@tab %>%
  rownames %>% 
  as.data.frame %>% 
  rename(Sample_Name = ".") %>% 
  left_join(ref_indivs, by = "Sample_Name") %>% 
  pull(Exposition) %>% 
  as.factor

# Subset the data to keep the parallel clines
data_para <- data[loc = Para_cline_SNPs] %>% 
  gi2gl()

Fst_parallel <- stamppFst(data_para, nboots = 10000, percent = 95, nclusters = 3)

Fst_parallel$Fsts %>% 
  as.data.frame %>% 
  rownames_to_column("Ecotype_1") %>% 
  pivot_longer(cols = !Ecotype_1,
               names_to = "Ecotype_2",
               values_to = "Fst") %>% 
  drop_na %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Parallel_SNPs_Fst.tsv")

Fst_parallel$Pvalues %>% 
  as.data.frame %>% 
  rownames_to_column("Ecotype_1") %>% 
  pivot_longer(cols = !Ecotype_1,
               names_to = "Ecotype_2",
               values_to = "P_value") %>% 
  drop_na %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Parallel_SNPs_P_value.tsv")

###################  Calculate mean FST between ecotypes on anti-parallel SNPs  ##################
data_antipara <- data[loc = Antipara_cline_SNPs] %>% 
  gi2gl()

Fst_antiparallel <- stamppFst(data_antipara, nboots = 10000, percent = 95, nclusters = 3)

Fst_antiparallel$Fsts %>% 
  as.data.frame %>% 
  rownames_to_column("Ecotype_1") %>% 
  pivot_longer(cols = !Ecotype_1,
               names_to = "Ecotype_2",
               values_to = "Fst") %>% 
  drop_na %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Anti_parallel_SNPs_Fst.tsv")

Fst_antiparallel$Pvalues %>% 
  as.data.frame %>% 
  rownames_to_column("Ecotype_1") %>% 
  pivot_longer(cols = !Ecotype_1,
               names_to = "Ecotype_2",
               values_to = "P_value") %>% 
  drop_na %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Anti_parallel_SNPs_P_value.tsv")

###################  Calculate mean FST between ecotypes on private_france  ##################
Private_fr <- private_fr %>% 
  select(Position) %>% 
  unique %>% 
  transform_position_ade2tidy() %>% 
  unite(Position, c(Chromosome, Position), sep = "_") %>% 
  pull(Position)

data_private_fr <- data[loc = Private_fr] %>% 
  gi2gl()

Fst_private_fr <- stamppFst(data_private_fr, nboots = 10000, percent = 95, nclusters = 3)

Fst_private_fr$Fsts %>% 
  as.data.frame %>% 
  rownames_to_column("Ecotype_1") %>% 
  pivot_longer(cols = !Ecotype_1,
               names_to = "Ecotype_2",
               values_to = "Fst") %>% 
  drop_na %>% 
  filter(grepl("Fr", Ecotype_1) & grepl("Fr", Ecotype_2)) %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Private_France_SNPs.tsv")

Fst_private_fr$Pvalues %>% 
  as.data.frame %>% 
  rownames_to_column("Ecotype_1") %>% 
  pivot_longer(cols = !Ecotype_1,
               names_to = "Ecotype_2",
               values_to = "P_value") %>% 
  drop_na %>% 
  filter(grepl("Fr", Ecotype_1) & grepl("Fr", Ecotype_2)) %>%
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Private_France_SNPs_P_value.tsv")


###################  Calculate mean FST between ecotypes on private_sweden  ##################
Private_sw <- private_sw %>% 
  select(Position) %>% 
  unique %>% 
  transform_position_ade2tidy() %>% 
  unite(Position, c(Chromosome, Position), sep = "_") %>% 
  pull(Position)

data_private_sw <- data[loc = Private_sw] %>% 
  gi2gl()

Fst_private_sw <- stamppFst(data_private_sw, nboots = 10000, percent = 95, nclusters = 3)

Fst_private_sw$Fsts %>% 
  as.data.frame %>% 
  rownames_to_column("Ecotype_1") %>% 
  pivot_longer(cols = !Ecotype_1,
               names_to = "Ecotype_2",
               values_to = "Fst") %>% 
  drop_na %>% 
  filter(grepl("Sw", Ecotype_1) & grepl("Sw", Ecotype_2)) %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Private_Sweden_SNPs_Fst.tsv")

Fst_private_sw$Pvalues %>% 
  as.data.frame %>% 
  rownames_to_column("Ecotype_1") %>% 
  pivot_longer(cols = !Ecotype_1,
               names_to = "Ecotype_2",
               values_to = "P_value") %>% 
  drop_na %>% 
  filter(grepl("Sw", Ecotype_1) & grepl("Sw", Ecotype_2)) %>%
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Private_Sweden_SNPs_P_value.tsv")

