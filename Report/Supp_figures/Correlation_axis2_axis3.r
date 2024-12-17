# Import libraries
libraries <- c("tidyverse", "adegenet", "vcfR", "readxl", "ggforce", "ggh4x", "ggnewscale")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)
rm(lib)
################## Useful functions  ##################
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

# Basic theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 30))

################## Import the vcf file  ##################
data <- read.vcfR("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/outputs/02_Filter_VCF_Files/11_LAM_and_LOK/VCF_File_Hobs.vcf.gz") %>% 
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
var_ax1 <- ((pca$eig[1] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax2 <- ((pca$eig[2] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax3 <- ((pca$eig[3] / sum(pca$eig)) * 100) %>% round(digits = 2)


############################ Correlation between contributions ##########################
mean_contributions <- pca$co %>%
  rownames_to_column("Position") %>% 
  select_good_SNPs(Comp1) %>% 
  transform_position_ade2tidy() %>% 
  mutate(Position = Position / 1e5) %>% 
  group_by(Position) %>% 
  summarize(Comp1 = mean(Comp1, na.rm = TRUE),
            Comp2 = mean(Comp2, na.rm = TRUE),
            Comp3 = mean(Comp3, na.rm = TRUE),
            Comp4 = mean(Comp4, na.rm = TRUE),
            Comp5 = mean(Comp5, na.rm = TRUE))

R2 <- (lm(Comp2 ~ Comp3, data = mean_contributions) %>% 
         summary)$r.squared

mean_contributions %>% 
  ggplot() +
  geom_point(aes(x = abs(Comp2), y = abs(Comp3))) +
  geom_line(data = data.frame(
    xpos = seq(0, 1, length.out = nrow(mean_contributions)),
    ypos = seq(0, 1, length.out = nrow(mean_contributions))
  ), aes(x = xpos, y = ypos), lwd = 1.2, color = "blue") +
  my_theme +
  labs(x = "Mean contribution of SNPs to Axis 2",
       y = "Mean contribution of SNPs to Axis 3") +
  xlim(0, 1) +
  ylim(0, 1)

lm(Comp2 ~ Comp3, data = mean_contributions) %>% 
  summary
