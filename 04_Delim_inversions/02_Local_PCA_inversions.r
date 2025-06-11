############################ Libraries ##########################
libraries <- c("tidyverse", "adegenet", "vcfR")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(libraries, character.only = TRUE)
rm(libraries)

############################ Useful functions ##########################
source("../General_scripts/Functions_optimise_plot_clines.r")

my_theme <- theme_bw() +
  theme(text = element_text(size = 30))

# Import the candidate regions for the inversions
## Sweden
Sweden_candidate_regions <- read_tsv("../../Output/Sweden_France_parallelism/04_Inversions/Candidate_inversions.tsv",
                                     col_names = TRUE) %>% 
  mutate(Chromosome = Chromosome %>%
           factor(levels = read_tsv("../../Output/Sweden_France_parallelism/04_Inversions/Candidate_inversions.tsv",
                                    col_names = TRUE) %>% 
                    select(Chromosome) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                    as.vector %>% unlist %>% unname)) %>% 
  filter(Population == "Sweden")
## France
France_candidate_regions <- read_tsv("../../Output/Sweden_France_parallelism/04_Inversions/Candidate_inversions.tsv",
                                     col_names = TRUE) %>% 
  mutate(Chromosome = Chromosome %>%
           factor(levels = read_tsv("../../Output/Sweden_France_parallelism/04_Inversions/Candidate_inversions.tsv",
                                    col_names = TRUE) %>% 
                    select(Chromosome) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                    as.vector %>% unlist %>% unname)) %>% 
  filter(Population == "France") %>% 
  mutate(Inversion = case_when(
    Inversion == "Inv_2.2" ~ "Inv_2.1",
    Inversion == "Inv_2.3" ~ "Inv_2.2",
    TRUE ~ Inversion))
# Import the vcf file
data <- read.vcfR("../../Output/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz") %>% 
  vcfR2genind()

data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor


############### Run local PCA and clusterisation on candidate inversion localisations ###############
## Sweden
local_pca_inversions_sweden <- Sweden_candidate_regions %>% 
  run_loca_pca_inversions(population = "LOKn")
pcas_candidate_inversions_sweden <- local_pca_inversions_sweden$PC_scores
# Represent all the potential inversions for PCA
pcas_candidate_inversions_sweden %>% 
  mutate(Habitat = str_split_fixed(Sample_Name, "_", 6)[, 5],
         Habitat = ifelse(Habitat == "EXPOS", "Exposed", ifelse(Habitat == "SHELT", "Sheltered", "Transition")) %>% 
           factor(levels = c("Exposed", "Transition", "Sheltered")),
         Inversion = Inversion %>% 
           factor(levels = pcas_candidate_inversions_sweden %>% 
                    select(Inversion) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist),
         
         Group = Group %>% factor(levels = c("1", "2", "3"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis2, colour = Habitat), size = 3, alpha = 0.7) +
  facet_wrap(vars(Inversion), scales = "free", ncol = 4) +
  scale_color_manual(values = c("Exposed" = "orange2", "Transition" = "deeppink", "Sheltered" = "dodgerblue2")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis1, y = Axis2, colour = Group), lwd = 1.2) +
  my_theme

list_inversions_sweden <- c("Isl_diff_1.3", "Isl_diff_3.1", "Isl_diff_3.2", "Isl_diff_4.1",
                            "Isl_diff_4.2", "Isl_diff_4.3", "Isl_diff_6.1",
                            "Isl_diff_6.2", "Isl_diff_6.3", "Isl_diff_7.1",
                            "Isl_diff_7.4", "Isl_diff_8.4", "Isl_diff_11.1",
                            "Isl_diff_13.2", "Isl_diff_14.1", "Isl_diff_14.2",
                            "Isl_diff_15.1", "Isl_diff_16.1", "Isl_diff_16.2")


## France
local_pca_inversions_france <- France_candidate_regions %>% 
  run_loca_pca_inversions(population = "LAMn")
pcas_candidate_inversions_france <- local_pca_inversions_france$PC_scores
# Represent all the potential inversions for PCA
pcas_candidate_inversions_france %>% 
  filter(!is.na(Inversion)) %>% 
  mutate(Habitat = str_split_fixed(Sample_Name, "_", 6)[, 5],
         Habitat = ifelse(Habitat == "EXPOS", "Exposed", ifelse(Habitat == "SHELT", "Sheltered", "Transition")) %>% 
           factor(levels = c("Exposed", "Transition", "Sheltered")),
         Inversion = Inversion %>% 
           factor(levels = pcas_candidate_inversions_france %>% 
                    select(Inversion) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist),
         
         Group = Group %>% factor(levels = c("1", "2", "3"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis2, colour = Habitat), size = 3, alpha = 0.7) +
  facet_wrap(vars(Inversion), scales = "free", ncol = 4) +
  scale_color_manual(values = c("Exposed" = "orange2", "Transition" = "deeppink", "Sheltered" = "dodgerblue2")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis1, y = Axis2, colour = Group), lwd = 1.2) +
  my_theme

list_inversions_france <- c("Isl_diff_1.2", "Isl_diff_2.3", "Isl_diff_3.1",
                            "Isl_diff_3.2", "Isl_diff_4.1", "Isl_diff_4.2",
                            "Isl_diff_4.3", "Isl_diff_6.1", "Isl_diff_6.2",
                            "Isl_diff_6.3", "Isl_diff_7.1", "Isl_diff_7.3",
                            "Isl_diff_8.3", "Isl_diff_11.1", "Isl_diff_13.1",
                            "Isl_diff_14.1", "Isl_diff_14.2", "Isl_diff_16.2",
                            "Isl_diff_16.3")


############### Regroup inversions that have three clusters ###############
correspondance_table <- Sweden_candidate_regions %>% 
  filter(Inversion %in% list_inversions_sweden) %>% 
  mutate(Population = "Sweden",
         Isl_diff = Inversion,
         Inversion = case_when(Inversion == "Isl_diff_1.3" ~ "Inv_1.1",
                               Inversion == "Isl_diff_7.4" ~ "Inv_7.2",
                               Inversion == "Isl_diff_8.4" ~ "Inv_8.1",
                               Inversion == "Isl_diff_13.2" ~ "Inv_13.1",
                               TRUE ~ Inversion) %>%
           str_replace_all("Isl_diff", "Inv")) %>%
  rbind(
    France_candidate_regions %>% 
      filter(Inversion %in% list_inversions_france) %>% 
      mutate(Population = "France",
             Isl_diff = Inversion,
             Inversion = case_when(
               Inversion == "Isl_diff_1.2" ~ "Inv_1.1",
               Inversion == "Isl_diff_2.3" ~ "Inv_2.1",
               Inversion == "Isl_diff_7.3" ~ "Inv_7.2",
               Inversion == "Isl_diff_8.3" ~ "Inv_8.1",
               Inversion == "Isl_diff_16.2" ~ "Inv_16.1",
               Inversion == "Isl_diff_16.3" ~ "Inv_16.2",
               TRUE ~ Inversion
             ) %>% str_replace_all("Isl_diff", "Inv"))) %>% 
  select(Chromosome, Isl_diff, Inversion, Population)


Confirmed_inversion <- Sweden_candidate_regions %>% 
  mutate(Population = "Sweden") %>% 
  rbind(France_candidate_regions %>% 
          mutate(Population = "France")) %>% 
  rename(Isl_diff = Inversion) %>% 
  inner_join(correspondance_table, by = c("Chromosome", "Isl_diff", "Population")) %>% 
  select(-Isl_diff) %>% 
  relocate(Inversion, .after = Chromosome)

# Write the output
Confirmed_inversion %>% 
  write.table("../../Output/Sweden_France_parallelism/04_Inversions/Candidate_inversion_post_pca.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

############### Save output ###############
# Save the local pca output
pcas_candidate_inversions_france %>% 
  rename(Isl_diff = Inversion) %>%
  mutate(Population = "France") %>% 
  rbind(pcas_candidate_inversions_sweden %>% 
          rename(Isl_diff = Inversion) %>% 
          mutate(Population = "Sweden")) %>% 
  inner_join(correspondance_table, by = c("Chromosome", "Isl_diff", "Population"), relationship = "many-to-one") %>% 
  select(-Isl_diff) %>% 
  relocate(Inversion, .after = Chromosome) %>% 
  write.table("../../Output/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

