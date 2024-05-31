# Import libraries
#require("anyLib")
anyLib(c("tidyverse", "adegenet", "vcfR", "readxl"))

################## Useful functions  ##################
source("../General_scripts/Functions_optimisation_visualisation.r")

################## Import the vcf file  ##################
data <- read.vcfR("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/filtering_vcf_files/Final_outputs/Fully_filtered_thinned_Hobs.vcf.gz") %>% 
  vcfR2genind()

# We add the poputations to the vcf object
data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor
data@other$exposition <- (data@tab %>% rownames %>% str_split_fixed(., "_", 6))[, 5] %>% as.factor
data@other$exposition[which(data@other$exposition == "TRANSI")] <- "TRANS" %>% as.factor
data@other$exposition <- data@other$exposition %>% droplevels()

############################ PCA on whole genome ##########################

# Scale the genome to get rid of missing data
X <- scaleGen(data, NA.method="mean", scale=FALSE, center=TRUE)

# Run the PCA
pca <- dudi.pca(X, scale=TRUE, nf=5, scannf=FALSE)

rm(X)
################## Import the metadata  ##################
metadata <- read_excel(path = "../Data/data_Fabalis_resequencing_Basile.xlsx",
                       sheet = 1,
                       col_names = TRUE,
                       trim_ws = TRUE) %>%
  
  # Convert to the correct formats
  mutate(Species = as.factor(Species),
         ID_number = as.factor(ID_number),
         Population = as.factor(Population),
         Transect = as.factor(TRANSECT),
         Id = as.factor(ID),
         Shell_colour = factor(Shell.colour %>% str_to_title, levels = c("Black", "Black/Square", "Brown", "Brown/Square", "Dark", "Yellow", "Yellow/Brown", "Yellow/Square", "Grey", "White", "Banded", NA)),
         LCmeanDist = as.numeric(LCmeanDist),
         Mreads = as.numeric(Mreads),
         Gbp = as.numeric(Gbp),
         Q30 = as.numeric(Q30),
         x = as.numeric(x),
         y = as.numeric(y),
         Length = as.numeric(length),
         Bi_size = as.factor(biSIZE %>% str_to_title),
         Habitat = ifelse(Habitat %in% c("EXPOS"), "Exposed", Habitat),
         Habitat = ifelse(Habitat %in% c("HARB", "SHELT"), "Sheltered", Habitat),
         Habitat = ifelse(Habitat %in% c("TRANS", "TRANSI"), "Transition", Habitat),
         Habitat = as.factor(Habitat)
  ) %>%
  
  # Select only the necessary columns for the analysis
  select(-c(length, biSIZE, Shell.colour, ID, TRANSECT)) %>% 
  
  # select only the data we need (the one on fabalis just in the transects from LAM and LOKn)
  filter(Species == "FAB",
         Population != "BRE",
         Transect == "n") %>% 
  
  # Modify the population column to get only the name of the country
  mutate(Population = ifelse(Population == "LOK", "Sweden", "France") %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  # Drop unused levels
  droplevels %>% 
  # Change the color information to have it into two simple colors: yellow or brown
  mutate(Shell_color_naive = (Shell_colour %>% str_split_fixed(., "/", 2))[, 1],
         Shell_color_morphology = (Shell_colour %>% str_split_fixed(., "/", 2))[, 2],
         Shell_color_naive = ifelse(Shell_color_naive %in% c("Yellow", "White", "Grey"), "Yellow", Shell_color_naive),
         Shell_color_morphology = ifelse(Shell_color_naive == "Banded", "Banded", Shell_color_morphology),
         Shell_color_naive = ifelse(Shell_color_naive %in% c("Black", "Brown", "Dark"), "Brown", Shell_color_naive),
         Shell_color_naive = ifelse(Shell_color_naive == "Banded", "Brown", Shell_color_naive),
         Shell_color_naive = Shell_color_naive %>% factor(levels = c("Yellow", "Brown")),
         Shell_color_morphology = ifelse(! Shell_color_morphology %in% c("Banded", "Square"), "Uniform", "Banded"))

################## Calculate the delta frequencies along the genome  ##################
Delta_freqs_whole_genome <- get_delta_freqs_and_F4(genetic_data = data,
                                                   SNP_subset = NULL,
                                                   nb_extreme_indivs = 30,
                                                   nb_indivs_to_keep = 20,
                                                   meta_data = metadata) %>% 
  select(-c(starts_with("p_"))) %>% 
  mutate(across(contains("Delta_"), .fns = abs, .names = "{.col}_abs")) %>% 
  # Calculate the mean delta freq
  mutate(Mean_delta_freq = rowMeans(select(., ends_with("_abs")), na.rm = TRUE)) %>% 
  select_good_SNPs(Delta_freq_Sweden)

################## Filter SNPs for which delta freqs are bigger than 30%  ##################
# First, we isolate all the delta freqs that are greater than 0.3
High_delta_freqs <- Delta_freqs_whole_genome %>% 
  filter(abs(Delta_freq_France) > 0.3 | abs(Delta_freq_Sweden) > 0.3) %>% 
  mutate(pos = Position) %>% 
  separate(pos, c("SNP_name", "allele"), "\\.") %>% 
  select(-allele)

Allele_freqs_high_delta_freqs <- get_allelic_frequencies(
  genetic_data = data,
  SNP_subset = High_delta_freqs,
  Extreme_values = pca$li,
  var = "Axis2",
  meta_data = metadata
)

################## Fit the cline models  ##################
# Prepare the priors for the filtered SNPs
Priors_high_delta_freqs <- Allele_freqs_high_delta_freqs %>% 
  filter(Position %in% High_delta_freqs$Position) %>% 
  mutate(Centre_prior = ifelse(Population == "France", 200, 70),
         Width_prior = ifelse(Population == "France", 50, 40),
         Centre_max = ifelse(Population == "France", 300, 150),
         Centre_min = ifelse(Population == "France", 100, 0),
         Width_max = ifelse(Population == "France", 700, 360),
         Width_min = 2)

Params_clines_high_delta_freqs <- optimise_clines(Priors = Priors_high_delta_freqs,
                                                  logarithm = TRUE,
                                                  batch_size = 100,
                                                  genetic_data = data,
                                                  SNP_subset = High_delta_freqs,
                                                  meta_data = metadata) %>% 
    select_clinal_SNPs()


# All SNPs making at least one cline
print(paste0("We have ", (Params_clines_high_delta_freqs %>% nrow) / 2, " SNPs making a cline in at least one population"))


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
  left_join(Delta_freqs_whole_genome %>% select(Position, F4_stat), by = "Position")

print(paste0("We have ", 
             (clines_fr_and_sw %>% nrow),
             " SNPs making clines in both populations"))

# SNPs that make parallel clines
para_clines <- clines_fr_and_sw %>% 
  filter(F4_stat >= 0.09)

print(paste0("We have ",
             (para_clines %>% nrow),
             " making parallel clines in both populations"))

# SNPs that make anti parallel clines
antipara_clines <- clines_fr_and_sw %>% 
  filter(F4_stat <= -0.09)

print(paste0("We have ",
             (antipara_clines %>% nrow),
             " SNPs making anti parallel clines"))

# Remaining
remaining <- clines_fr_and_sw %>% 
  filter(Position %!in% para_clines$Position & Position %!in% antipara_clines$Position)

print(paste0("We finally have ",
             nrow(remaining),
             " SNPs that form parallel clines but that are not highly differenciated"))



# SNPs that do not make clines in two populations
all_private <- Params_clines_high_delta_freqs %>% 
  filter(Position %!in% clines_fr_and_sw$Position)

print(paste0("We have ",
             (all_private %>% nrow) / 2,
             " SNPs making clines in exactly one location"))

# SNPs that make clines only in France
private_fr <- Params_clines_high_delta_freqs %>% 
  filter(Population == "France",
         Best_model == "Clinal",
         abs(Delta_AIC_second_best_model) >= 4,
         Position %!in% clines_sw$Position) 

print(paste0("We have ",
             (private_fr %>% nrow),
             " SNPs making clines only in France"))

# SNPs that make clines only in France
private_sw <- Params_clines_high_delta_freqs %>% 
  filter(Population == "Sweden",
         Best_model == "Clinal",
         abs(Delta_AIC_second_best_model) >= 4,
         Position %!in% clines_fr$Position) 

print(paste0("We have ",
             (private_sw %>% nrow),
             " SNPs making clines only in Sweden"))

# Find if private snps are in inversions
## France
private_fr_in_invs <- private_fr %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>% 
  drop_na %>% 
  filter(Position > Pos_min_inv & Position < Pos_max_inv)

print(paste0("We have ",
             private_fr_in_invs %>% nrow,
             " SNPs making clines only in France located in the inversions"))
## Sweden
private_sw_in_invs <- private_sw %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>% 
  drop_na %>% 
  filter(Position > Pos_min_inv & Position < Pos_max_inv)

print(paste0("We have ",
             private_sw_in_invs %>% nrow,
             " SNPs making clines only in Sweden located in the inversions"))

# Find if clinal snps in both populations are in inversions
clines_in_inversions <- clines_fr_and_sw %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>% 
  drop_na %>% 
  filter(Position > Pos_min_inv & Position < Pos_max_inv)

print(paste0("We have ",
             clines_in_inversions %>% nrow,
             " SNPs making clines in both populations located in the inversions"))

# Find if parallel clinal snps are in inversions
para_in_inversions <- para_clines %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>% 
  drop_na %>% 
  filter(Position > Pos_min_inv & Position < Pos_max_inv)

print(paste0("We have ",
             para_in_inversions %>% nrow,
             " SNPs making parallel clines located in the inversions"))

# Find if anti-parallel clinal snps are in inversions
antipara_in_inversions <- antipara_clines %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  separate(pos, c("Chromosome", "SUPER", "SUPER_frag", "pos"), "_") %>% 
  left_join(Delim_inversions, by = "Chromosome", relationship = "many-to-many") %>% 
  drop_na %>% 
  filter(Position > Pos_min_inv & Position < Pos_max_inv)

print(paste0("We have ",
             antipara_in_inversions %>% nrow,
             " SNPs making anti-parallel clines located in the inversions"))
