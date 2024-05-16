# Import libraries
require("anyLib")
#devtools::install_github("thomasp85/patchwork")
anyLib(c("tidyverse", "adegenet", "vcfR", "readxl", "statgenGWAS", "ggforce", "ggh4x", "patchwork"))


################## Useful functions  ##################
"%!in%" <- function(x, y){!(x %in% y)}

source("/shared/projects/pacobar/finalresult/bpajot/scripts/Phenotypic_analysis/Cline_functions.R")
source("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/Functions_optimise_plot_clines.r")

################################ Useful variables ################################
# Color palette to be reused everywhere with the shell size
size_palette = c("#4e79a7", "grey75", "#f28e2b")
# Basic theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

################## Import the vcf file  ##################
data <- read.vcfR("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/filtering_vcf_files/Final_outputs/Fully_filtered_thinned_Hobs.vcf.gz") %>% 
  vcfR2genind()

# We use the summary function to calculate the Hobs of the data we imported
#data_info <- summary(data)
#Hobs <- data_info$Hobs
#hist(Hobs)

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
metadata <- read_excel(path = "/shared/projects/pacobar/finalresult/bpajot/Data/data_Fabalis_resequencing_Basile.xlsx",
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

################## Making the delta frequencies correlation plot  ##################
Delta_freqs_whole_genome <- get_delta_freqs_and_F4(
  genetic_data = data,
  Extreme_values = pca$li,
  var = "Axis2",
  metadata = metadata
) %>% 
  left_join(pca$co %>% 
              rownames_to_column("Position"),
            by = "Position")

# Delta freqs colored by PC2
delta_PC2 <- Delta_freqs_whole_genome %>%
  ggplot() +
  geom_point(aes(Delta_freq_France, Delta_freq_Sweden, color = Comp2 %>% abs), alpha = 0.1, size = 3) +
  scale_color_gradientn(name = "Contribution PC2\n(valeur absolue)",
                        colors = c("olivedrab4", "#f28e2b", "brown4"), 
                        limits = c(min(pca$co$Comp2 %>% abs, na.rm=TRUE) - 0.01, 
                                 max(pca$co$Comp2 %>% abs, na.rm=TRUE) + 0.01),
                        breaks = c(0.025, 0.5, 0.75)) +
  labs(x = expression(Delta * "freq_France"),
       y = expression(Delta * "freq_Sweden"),
       tag = "(A)") +
  my_theme

# Delta freqs colored by PC3
delta_PC3 <- Delta_freqs_whole_genome %>%
  ggplot() +
  geom_point(aes(Delta_freq_France, Delta_freq_Sweden, color = Comp3 %>% abs), alpha = 0.1, size = 3) +
  scale_color_gradientn(name = "Contribution PC3\n(valeur absolue)",
                        colors = c("deepskyblue", "deeppink", "darkred"), 
                        limits = c(min(pca$co$Comp3 %>% abs, na.rm=TRUE) - 0.01, 
                                   max(pca$co$Comp3 %>% abs, na.rm=TRUE) + 0.01),
                        breaks = c(0.05, 0.5, 1)) +
  labs(x = expression(Delta * "freq_France"),
       y = expression(Delta * "freq_Sweden"),
       tag = "(B)") +
  my_theme


delta_PC2 + delta_PC3 +
  plot_layout(guides = "collect")


