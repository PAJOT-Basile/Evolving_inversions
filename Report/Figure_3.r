# Import libraries
require("anyLib")
#devtools::install_github("thomasp85/patchwork")
anyLib(c("tidyverse", "adegenet", "vcfR", "readxl", "statgenGWAS", "ggforce", "ggh4x", "patchwork"))


################## Useful functions  ##################
"%!in%" <- function(x, y){!(x %in% y)}

source("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/Functions_optimise_plot_clines.r")

################################ Useful variables ################################
# Color palette to be reused everywhere with the shell size
size_palette <-  c("#4e79a7", "grey75", "#f28e2b")
# Basic theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 30))
# Palette to use on the chromosomes
chromosome_palette <- c("grey71", "turquoise4")

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
# Extract the percentage of explained variance of interesting axis
summary(pca)
var_ax1 <- 11.09
var_ax2 <- 5.19
var_ax3 <- 0.97
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


################## Plot PCA1 vs PCA2  ##################
ax1_v_ax2 <- pca$li %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  mutate(Habitat = ifelse(Habitat == "Exposed", "Exposé", ifelse(Habitat == "Transition", "Transition", "Abrité")) %>% 
           factor(levels = c("Exposé", "Transition", "Abrité")),
         Population = ifelse(Population == "Sweden", "Suède", "France")) %>% 
  ggplot() +
  geom_point(aes(x = Axis2, y = Axis1, color = Habitat), alpha = 0.7, size = 4) +
  scale_color_manual(name = "Habitat",
                     values = c("Exposé" = "orange2", "Transition" = "deeppink", "Abrité" = "dodgerblue3")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis2, y = Axis1, color = Population, label = Population),
                    lwd = 1.2, label.fontsize = 20) +
  scale_color_manual(values = c("Suède" = "forestgreen", "France" = "darkorchid4"),
                     guide = "none") +
  labs(x = paste0("Axe 2 (", var_ax2, " %)"),
       y = paste0("Axe 1 (", var_ax1, " %)"),
       tag = "(A)") +
  xlim(-425, 310) +
  ylim(-410, 310) +
  my_theme

################## Plot PCA1 vs PCA3  ##################
ax1_v_ax3 <- pca$li %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  mutate(Habitat = ifelse(Habitat == "Exposed", "Exposé", ifelse(Habitat == "Transition", "Transition", "Abrité")) %>% 
           factor(levels = c("Exposé", "Transition", "Abrité")),
         Population = ifelse(Population == "Sweden", "Suède", "France")) %>% 
  ggplot() +
  geom_point(aes(x = Axis3, y = Axis1, color = Habitat), alpha = 0.7, size = 4) +
  scale_color_manual(name = "Habitat",
                     values = c("Exposé" = "orange2", "Transition" = "deeppink", "Abrité" = "dodgerblue3")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis3, y = Axis1, color = Population, label = Population),
                    lwd = 1.2, label.fontsize = 20) +
  scale_color_manual(values = c("Suède" = "forestgreen", "France" = "darkorchid4"),
                     guide = "none") +
  labs(x = paste0("Axe 3 (", var_ax3, " %)"),
       y = paste0("Axe 1 (", var_ax1, " %)"),
       tag = "(B)") +
  xlim(-250, 130) +
  ylim(-410, 310) +
  my_theme

################## Plot PC2 along the transect in Sweden and France  ##################
pc2_trans <- pca$li %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France"))) %>% 
  ggplot() +
  geom_point(aes(x = LCmeanDist, y = Axis2, color = Length), size = 4) +
  scale_colour_gradientn(name = "Taille\n(mm)",
                         colors=size_palette, 
                         limits=c(min(metadata$Length, na.rm=TRUE) - 0.01, 
                                  max(metadata$Length, na.rm=TRUE) + 0.01),
                         breaks=c(8, 13, 18)) +
  facet_col(facets = vars(Population), scales = "free") +
  labs(x = "Position le long du transect (m)",
       y = paste0("Axe 2 (", var_ax2, " %)"),
       tag = "(C)") +
  my_theme

################## Plot PC3 along the transect in Sweden and France  ##################
pc3_trans <- pca$li %>% 
  rownames_to_column("Sample_Name") %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France"))) %>% 
  ggplot() +
  geom_point(aes(x = LCmeanDist, y = - Axis3, color = Length), size = 4) +
  scale_colour_gradientn(name = "Taille\n(mm)",
                         colors=size_palette, 
                         limits=c(min(metadata$Length, na.rm=TRUE) - 0.01, 
                                  max(metadata$Length, na.rm=TRUE) + 0.01),
                         breaks=c(8, 13, 18)) +
  facet_col(facets = vars(Population), scales = "free") +
  labs(x = "Position le long du transect (m)",
       y = paste0("Axe 3 (", var_ax3, " %)"),
       tag = "(D)") +
  my_theme

################## Manhattan plot of contributions to PC2 and PC3  ##################
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
  plot_manhattan(Comp2, absolute = FALSE, palette = chromosome_palette)) +
  labs(y = "Contributions des SNPs à\nAxe 3 | Axe 2 ",
       x = "Groupes de liaison",
       tag = "(E)") +
  geom_hline(yintercept = 0, color = "black", lwd = 0.9) +
  theme(text = element_text(size = 30))


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
