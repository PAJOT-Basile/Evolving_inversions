# Import libraries
libraries <- c("tidyverse", "adegenet", "vcfR", "readxl", "ggforce", "ggh4x", "ggnewscale")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(charaters = libraries, character.only = TRUE)
rm(libraries)
if (!require("patchwork")) devtools::install_github("thomasp85/patchwork")

################## Useful functions  ##################
source("../General_scripts/Functions_optimise_plot_clines.r")
################################ Useful variables ################################
# Color palette to be reused everywhere with the shell size
size_palette <-  c("#4e79a7", "grey75", "#f28e2b")
# Basic theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))
# Palette to use on the chromosomes
chromosome_palette <- c("grey71", "turquoise4")

# Position of the inversions
position_inversions <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
                                  sep = "\t", header = TRUE)

grouped_inversions <- position_inversions %>% 
  filter(Inversion_grouped %in% c("Inv_4.1", "Inv_14.1", "Inv_3.1")) %>% 
  mutate(Inversion = Inversion_grouped)
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
pca <- scaleGen(data, NA.method="mean", scale=FALSE, center=TRUE) %>% 
  dudi.pca(scale=TRUE, nf=5, scannf=FALSE)

# Extract the percentage of explained variance of interesting axis
var_ax1 <- ((pca$eig[1] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax2 <- ((pca$eig[2] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax3 <- ((pca$eig[3] / sum(pca$eig)) * 100) %>% round(digits = 2)
################## Import the metadata  ##################
metadata <- read_excel(path = "../Input_Data/Data/data_Fabalis_resequencing_Basile.xlsx",
                       sheet = 1,
                       col_names = TRUE,
                       trim_ws = TRUE) %>%
  
  # Convert to the correct formats
  mutate(Population = as.factor(Population),
         Shell_colour = factor(Shell.colour %>% str_to_title, levels = c("Black", "Black/Square", "Brown", "Brown/Square", "Dark", "Yellow", "Yellow/Brown", "Yellow/Square", "Grey", "White", "Banded", NA)),
         LCmeanDist = as.numeric(LCmeanDist),
         Length = as.numeric(length),
         Habitat = ifelse(Habitat %in% c("EXPOS"), "Exposed", Habitat),
         Habitat = ifelse(Habitat %in% c("HARB", "SHELT"), "Sheltered", Habitat),
         Habitat = ifelse(Habitat %in% c("TRANS", "TRANSI"), "Transition", Habitat),
         Habitat = as.factor(Habitat)
  ) %>%
  
  # Select only the necessary columns for the analysis
  select(-length) %>% 

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


################## Plot PCA1 vs PCA2  ##################
ax1_v_ax2 <- pca$li %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  mutate(Habitat = Habitat %>% 
           factor(levels = c("Exposed", "Transition", "Sheltered"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis2, y = Axis1, color = Habitat), alpha = 0.7, size = 4) +
  scale_color_manual(name = "Habitat",
                     values = c("Exposed" = "orange2", "Transition" = "deeppink", "Sheltered" = "dodgerblue3")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis2, y = Axis1, color = Population, label = Population),
                    lwd = 1.2, label.fontsize = 15) +
  scale_color_manual(values = c("Sweden" = "forestgreen", "France" = "darkorchid4"),
                     guide = "none") +
  labs(x = paste0("Axis 2 (", var_ax2, " %)"),
       y = paste0("Axis 1 (", var_ax1, " %)"),
       tag = "(A)") +
  my_theme +
  xlim(-550, 400) +
  ylim(-475, 400)

################## Plot PCA1 vs PCA3  ##################
ax1_v_ax3 <- pca$li %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  mutate(Habitat = Habitat %>% 
           factor(levels = c("Exposed", "Transition", "Sheltered"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis3, y = Axis1, color = Habitat), alpha = 0.7, size = 4) +
  scale_color_manual(name = "Habitat",
                     values = c("Exposed" = "orange2", "Transition" = "deeppink", "Sheltered" = "dodgerblue3")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis3, y = Axis1, color = Population, label = Population),
                    lwd = 1.2, label.fontsize = 15) +
  scale_color_manual(values = c("Sweden" = "forestgreen", "France" = "darkorchid4"),
                     guide = "none") +
  labs(x = paste0("Axis 3 (", var_ax3, " %)"),
       y = paste0("Axis 1 (", var_ax1, " %)"),
       tag = "(B)") +
  my_theme +
  xlim(-300, 200) +
  ylim(-450, 350)

################## Plot PC2 along the transect in Sweden and France  ##################
pc2_trans <- pca$li %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  ggplot() +
  geom_point(aes(x = LCmeanDist, y = Axis2, color = Length), size = 4) +
  scale_colour_gradientn(name = "Shell size\n(mm)",
                         colors=size_palette, 
                         limits=c(min(metadata$Length, na.rm=TRUE) - 1, 
                                  max(metadata$Length, na.rm=TRUE) + 1),
                         breaks=c(5.6, 10.5, 13.8)) +
  facet_col(facets = vars(Population), scales = "free") +
  labs(x = "Position along the transect (m)",
       y = "Individual's PC 2 score",
       tag = "(C)") +
  my_theme

################## Plot PC3 along the transect in Sweden and France  ##################
pc3_trans <- pca$li %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  ggplot() +
  geom_point(aes(x = LCmeanDist, y = - Axis3, color = Length), size = 4) +
  scale_colour_gradientn(name = "Shell size\n(mm)",
                         colors=size_palette, 
                         limits=c(min(metadata$Length, na.rm=TRUE) - 1, 
                                  max(metadata$Length, na.rm=TRUE) + 1),
                         breaks=c(5.6, 10.5, 13.8)) +
  facet_col(facets = vars(Population), scales = "free") +
  labs(x = "Position along the transect (m)",
       y = "Individual's PC 3 score",
       tag = "(D)") +
  my_theme

################## Manhattan plot of contributions to PC2 and PC3  ##################
# Get a table of the maximum number of positions along each chromosome
Cumulative_position_per_chromosome <- data@loc.fac %>% 
  as_tibble %>% 
  rename(Position = value) %>% 
  transform_position_ade2tidy() %>% 
  group_by(Chromosome) %>% 
  summarize(Pos_max_chromosome = max(Position)) %>% 
  arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
  mutate(Pos_max_chromosome_cumul = lag(cumsum(Pos_max_chromosome), default = 0))

# Plot the contribution of each SNP to the PCA
contribs <- (pca$co %>% 
               rownames_to_column("Position") %>% 
               mutate(Comp2 = abs(Comp2)) %>% 
               select(Position, Comp2) %>%
               select_good_SNPs(Comp2) %>% 
               rbind(pca$co %>% 
                       rownames_to_column("Position") %>%
                       select_good_SNPs(Comp3) %>% 
                       mutate(Comp2 = -abs(Comp3)) %>% 
                       select(Position, Comp2)) %>%
               geom_manhattan(aes(y = Comp2), absolute = FALSE, palette = chromosome_palette)) +
  geom_hline(yintercept = 0, color = "black", lwd = 0.9) +
  geom_polygon(data = position_inversions %>% 
                 filter(Population == "France") %>% 
                 positions_to_polygons(values = c(-0.05, 0)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion), fill = "black") +
  geom_polygon(data = position_inversions %>% 
                 filter(Population == "Sweden") %>% 
                 positions_to_polygons(values = c(0, 0.05)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion), fill = "black") +
  geom_polygon(data = grouped_inversions %>% 
                 filter(Population == "Sweden") %>% 
                 positions_to_polygons(values = c(0, 0.05)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion),
               color = "red", fill = NA) +
  geom_polygon(data = grouped_inversions %>% 
                 filter(Population == "France") %>% 
                 positions_to_polygons(values = c(-0.05, 0)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion),
               color = "red", fill = NA) +
  labs(y = "SNP's eigenvalue to\nAxis 3 | Axis 2",
       x = "Chromosomes",
       tag = "(E)") +
  my_theme


################## Merge graphs together  ##################
layout_pca <- "
AAAAAABBBBBB
AAAAAABBBBBB
"

pca_plots <- ax1_v_ax2 + ax1_v_ax3 +
  plot_layout(design = layout_pca,
              guides = "collect")
pca_transect_plots <- pc2_trans + pc3_trans +
  plot_layout(design = layout_pca,
              guides = "collect")


pca_plots / pca_transect_plots / contribs
