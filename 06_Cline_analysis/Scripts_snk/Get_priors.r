# Import libraries
libraries <- c("tidyverse", "adegenet", "vcfR", "readxl", "argparse")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)


################## Import arguments  ##################
parser <- ArgumentParser(description = "This program calculates the priors for the cline analysis")

# Add the arguments that are used
parser$add_argument("--metadata", "-m", help = "The metadata path")
parser$add_argument("--output", "-o", help = "The output path")
parser$add_argument("--genetic_data", "-g", help = "The path to the vcf file")

xargs <- parser$parse_args()

metadata_path <- xargs$metadata
output_file <- xargs$output
vcf_file <- xargs$genetic_data
################## Useful functions  ##################
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

################## Import the vcf file  ##################
data <- read.vcfR(vcf_file) %>% 
  vcfR2genind()

# We add the poputations to the vcf object
data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor
data@other$exposition <- (data@tab %>% rownames %>% str_split_fixed(., "_", 6))[, 5] %>% as.factor
data@other$exposition[which(data@other$exposition == "TRANSI")] <- "TRANS" %>% as.factor
data@other$exposition <- data@other$exposition %>% droplevels()
print("Got data")

# Run the PCA
pca <- scaleGen(data, NA.method="mean", scale=FALSE, center=TRUE) %>%
    dudi.pca(scale=TRUE, nf=5, scannf=FALSE)

################## Import the metadata  ##################
metadata <- read_excel(path = metadata_path,
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
                                                  Extreme_values = pca$li,
                                                  var = "Axis2",
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

High_delta_freqs %>%
 write.table(paste0(output_file, "High_delta_freqs.tsv"),
             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

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

Priors_high_delta_freqs %>%
 write.table(paste0(output_file, "Priors_high_delta_freqs.tsv"),
             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
