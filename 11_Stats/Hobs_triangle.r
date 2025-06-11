# Import libraries
libraries <- c("tidyverse", "adegenet", "vcfR", "readxl", "ggforce", "ggh4x")
if (!require("pacman")) install.packages("pacman")
for (lib in libraries){
  pacman::p_load(lib, character.only = TRUE)
}
rm(lib)
if (!require("patchwork")) devtools::install_github("thomasp85/pacthwork")


################## Useful functions  ##################
source("../General_scripts/Functions_optimise_plot_clines.r")

################################ Useful variables ################################
# Color palette to be reused everywhere with the shell size
size_palette <- c("#4e79a7", "grey75", "#f28e2b")
# Basic theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

################## Import the vcf file  ##################
data <- read.vcfR("../../Output/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz") %>% 
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
# Extract the percentage of explained variance of interesting axis
summary(pca)
var_ax1 <- ((pca$eig[1] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax2 <- ((pca$eig[2] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax3 <- ((pca$eig[3] / sum(pca$eig)) * 100) %>% round(digits = 2)

################## Import the metadata  ##################
metadata <- read_excel(path = "../../Input_Data/Data/data_Fabalis_resequencing_Basile.xlsx",
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


################## Calculate the heterozygosity of each individual  ##################
# Calculate the differences in allelic frequencies
Delta_freqs_whole_genome <- get_delta_freqs_and_F4(
  genetic_data = data,
  Extreme_values = pca$li,
  var = "Axis2",
  meta_data = metadata
) %>% 
  left_join(pca$co %>% 
              rownames_to_column("Position"),
            by = "Position")

# Get the names of the distinctive SNPs
Distinctive_SNPs <- Delta_freqs_whole_genome %>% 
  filter((Delta_freq_Sweden > 0.8) | (Delta_freq_France > 0.8)) %>% 
  select(Position) %>% 
  mutate(pos = Position) %>% 
  separate(pos, c("SNP_name", "allele"), "\\.") %>% 
  select(-allele)
# See where the individuals are hetarozygote
Hobs <- rep(0, nrow(data@tab))
for (i in 1:nrow(data@tab)){
  Hobs[i] <- data[loc = Distinctive_SNPs$SNP_name]@tab[i, ] %>% 
    as.data.frame() %>% 
    rename(Indiv = ".") %>% 
    filter(Indiv == 1) %>% 
    drop_na %>% 
    colSums(na.rm = TRUE) / ((data[loc = Distinctive_SNPs$SNP_name]@tab[i, ] %>% length) - (data[loc = Distinctive_SNPs$SNP_name]@tab[i, ] %>%  
                                                                                              as.data.frame() %>% 
                                                                                              is.na() %>% colSums(na.rm = TRUE)))
}
names(Hobs) <- rownames(data@tab)
#Hobs <- Hobs/2
polygons <- data.frame(
  "x" = c(Hobs %>% 
          as.data.frame() %>% 
          rename(Hobs = ".") %>% 
          rownames_to_column("Sample_Name") %>% 
          left_join(metadata, by = "Sample_Name") %>% 
            filter(Population == "France") %>% 
          left_join(pca$li %>% 
                      rownames_to_column("Sample_Name"),
                    by = "Sample_Name") %>% select(Axis2) %>% min,
        Hobs %>% 
          as.data.frame() %>% 
          rename(Hobs = ".") %>% 
          rownames_to_column("Sample_Name") %>% 
          left_join(metadata, by = "Sample_Name") %>% 
          filter(Population == "France") %>% 
          left_join(pca$li %>% 
                      rownames_to_column("Sample_Name"),
                    by = "Sample_Name") %>% select(Axis2) %>% max,
        Hobs %>% 
          as.data.frame() %>% 
          rename(Hobs = ".") %>% 
          rownames_to_column("Sample_Name") %>% 
          left_join(metadata, by = "Sample_Name") %>%
          filter(Population == "France") %>% 
          left_join(pca$li %>% 
                      rownames_to_column("Sample_Name"),
                    by = "Sample_Name") %>% select(Axis2) %>% colMeans,
        Hobs %>% 
          as.data.frame() %>% 
          rename(Hobs = ".") %>% 
          rownames_to_column("Sample_Name") %>% 
          left_join(metadata, by = "Sample_Name") %>% 
          filter(Population == "Sweden") %>% 
          left_join(pca$li %>% 
                      rownames_to_column("Sample_Name"),
                    by = "Sample_Name") %>% select(Axis2) %>% min,
        Hobs %>% 
          as.data.frame() %>% 
          rename(Hobs = ".") %>% 
          rownames_to_column("Sample_Name") %>% 
          left_join(metadata, by = "Sample_Name") %>% 
          filter(Population == "Sweden") %>% 
          left_join(pca$li %>% 
                      rownames_to_column("Sample_Name"),
                    by = "Sample_Name") %>% select(Axis2) %>% max,
        Hobs %>% 
          as.data.frame() %>% 
          rename(Hobs = ".") %>% 
          rownames_to_column("Sample_Name") %>% 
          left_join(metadata, by = "Sample_Name") %>%
          filter(Population == "Sweden") %>% 
          left_join(pca$li %>% 
                      rownames_to_column("Sample_Name"),
                    by = "Sample_Name") %>% select(Axis2) %>% colMeans),
  "y" = c(0, 0, 1),
  "Group" = c(1, 1, 1),
  Population = c(rep("France", 3), rep("Sweden", 3))
) %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France")))

Hobs %>% 
  as.data.frame() %>% 
  rename(Hobs = ".") %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  left_join(pca$li %>% 
              rownames_to_column("Sample_Name"),
            by = "Sample_Name") %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France")),
         Habitat = Habitat %>% 
           factor(levels = c("Exposed", "Transition", "Sheltered"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis2, y = Hobs, color = Habitat), size = 3, alpha = 0.7) +
  scale_color_manual(name = "Habitat",
                     values = c("Exposed" = "orange2", "Transition" = "deeppink", "Sheltered" = "dodgerblue3")) +
  geom_polygon(data = polygons, aes(x = x, y = y, group = Group), fill = NA, colour = "black", lwd = 1.2) +
  my_theme +
  facet_row(vars(Population)) +
  ylim(0, 1) +
  labs(x = paste0("Axis 2 (", var_ax2, " %)"))
