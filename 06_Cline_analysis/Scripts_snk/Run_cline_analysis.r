# Import libraries
libraries <- c("tidyverse", "adegenet", "vcfR", "readxl", "argparse")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)


################## Import arguments  ##################
parser <- ArgumentParser(description = "This program calculates the priors for the cline analysis")

# Add the arguments that are used
parser$add_argument("--priors", "-p", help = "Priors")
parser$add_argument("--vcf_file", "-v", help = "Vcf file")
parser$add_argument("--metadata", "-m", help = "Metadata")
parser$add_argument("--output", "-o", help = "The output path")
parser$add_argument("--high_delta_freqs", "-d", help = "The high delta freq values")

xargs <- parser$parse_args()

metadata_path <- xargs$metadata
vcf_file <- xargs$vcf_file
prior_path <- xargs$priors
output_file <- xargs$output
high_freq_path <- xargs$high_delta_freqs
print("Imported arguments")
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
print("Got metadata")
################## Import the priors  ##################
Priors_high_delta_freqs <- read.table(prior_path, sep = "\t", header = TRUE)

################## Filter SNPs for which delta freqs are bigger than 30%  ##################
High_delta_freqs <- read.table(high_freq_path, sep = "\t", header = TRUE) %>%
    filter(Position %in% Priors_high_delta_freqs$Position)

################## Fit the cline models  ##################
print("Starting optimisation")
optimise_clines(Priors = Priors_high_delta_freqs,
                logarithm = TRUE,
                batch_size = 100,
                genetic_data = data,
                SNP_subset = High_delta_freqs,
                meta_data = metadata,
                write_output = output_file)
