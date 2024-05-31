# Libraries
require("anyLib")
anyLib(c("tidyverse", "adegenet", "vcfR", "readxl", "ggpubr"))


################## Useful functions  ##################
source("../General_scripts/Functions_optimise_plot_clines.r")

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


################## Import GWAS results  ##################
## France
GWAS_France <- read.csv("../Data/output_france.csv", header = TRUE, sep = "\t")

## Sweden
GWAS_sweden <- read.csv("../Data/output_sweden.csv", header = TRUE, sep = "\t")

################## Import inversion pca results  ##################
clust_groups <- read.csv("../Data/clust_groups.csv", header = TRUE, sep = "\t")
local_inv_indivs <- read.csv("../Data/local_inversions_indivs.csv", header = TRUE, sep = "\t")

################## Identify region for size QTL  ##################
GWAS_sweden %>% 
  filter(trait == "Length",
         grepl("LG14_", snp)) %>% 
  mutate(log_pval = -log10(pValue)) %>% 
  rename(Position = snp) %>% 
  geom_manhattan(aes(y = log_pval), palette = c("grey71", "turquoise4")) +
# We add the delimitatins of inversions in chromosome 14
  geom_vline(xintercept = 200000, color = "black") +
  geom_vline(xintercept = 23600000, color = "black") +
  geom_vline(xintercept = 35900000, color = "red") +
  geom_vline(xintercept = 45000000, color = "red")

print("We are looking at inversion 14.1")  

################## Identify region for color QTL  ##################
GWAS_France %>% 
  filter(trait == "Shell_color_naive",
         grepl("LG6_", snp)) %>% 
  mutate(log_pval = -log10(pValue)) %>% 
  rename(Position = snp) %>% 
  geom_manhattan(aes(y = log_pval), palette = c("grey71", "turquoise4")) +
# We add the delimitatins of inversions in chromosome 6
  geom_vline(xintercept = 200000, color = "black") +
  geom_vline(xintercept = 24800000, color = "black") +
  geom_vline(xintercept = 42700000, color = "red") +
  geom_vline(xintercept = 55200000, color = "red")

print("We are looking at inversion 6.1")  


################## Plot the local inversion PCA for 6.1 and 14.1  ##################
D <- local_inv_indivs %>% 
# We add some columns of the metadata to have the color of the individuals
  left_join(metadata, by = "Sample_Name", relationship = "many-to-one") %>% 
# We keep just the inversion we identified earlier as the one in which is located the color QTL
  filter(Inversion == "Inv_6.1") %>% 
# We reorganise the levels of the shell color and the population to have them in the legend
  mutate(Shell_color_naive = ifelse(Shell_color_naive == "Brown", "Marron", "Jaune") %>% 
           factor(levels = c("Marron", "Jaune")),
         Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis2, color = Shell_color_naive), size = 3, alpha = 0.7) +
  scale_color_manual(name = "Couleur\ncoquille",
                     values = c("Marron" = "brown", "Jaune" = "orange")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis1, y = Axis2, color = Population), lwd = 1.2) +
  scale_color_manual(name = "Population",
                     values = c("Suède" = "forestgreen", "France" = "darkorchid4")) +
  labs(x = "Axe 1",
       y = "Axe 2",
       tag = "(D)") +
  theme_bw() + 
  theme(text = element_text(size = 30))

A <- local_inv_indivs %>% 
# In the same way, we keep the inversion we identified earlier that has the size QTL and we 
# reorganise the levels for the legend
  left_join(metadata, by = "Sample_Name", relationship = "many-to-one") %>% 
  filter(Inversion == "Inv_14.1") %>% 
  mutate(Shell_color_naive = ifelse(Shell_color_naive == "Brown", "Marron", "Jaune") %>% 
           factor(levels = c("Marron", "Jaune")),
         Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis2, color = Length), size = 3, alpha = 0.7) +
  scale_colour_gradientn(name = "Taille\n(mm)",
                         colors=c("#4e79a7", "grey75", "#f28e2b"), 
                         limits=c(min(metadata$Length, na.rm=TRUE) - 0.01, 
                                  max(metadata$Length, na.rm=TRUE) + 0.01),
                         breaks=c(8, 13, 18)) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis1, y = Axis2, color = Population), lwd = 1.2) +
  scale_color_manual(name = "Population",
                     values = c("Suède" = "forestgreen", "France" = "darkorchid4")) +
  labs(x = "Axe 1",
       y = "Axe 2",
       tag = "(A)") +
  theme_bw() + 
  theme(text = element_text(size = 30))

################## Make cline of QTLs for size ##################
# In the size QTL, as we do not have many individuals and that size is a very polygenic
# trait, we choose to lower the threshold value usually used to identify QTLs (usually 5)
# to 4.
size_snps <- GWAS_sweden %>% 
  filter(trait == "Length",
         grepl("LG14_", snp)) %>% 
  mutate(log_pval = -log10(pValue)) %>% 
  rename(Position = snp) %>% 
  filter(log_pval >= 4) %>% 
  mutate(pos = Position) %>% 
  separate(pos, c("SNP_name", "allele"), "\\.") %>% 
  select(Position, SNP_name)

# Prepare priors to fit allelic frequency variation models
Priors_clines <- get_allelic_frequencies(genetic_data = data,
                                             SNP_subset = size_snps,
                                             Extreme_values = pca$li,
                                             var = "Axis2",
                                             meta_data = metadata) %>% 
  mutate(Centre_prior = ifelse(Population == "France", 200, 70),
         Width_prior = ifelse(Population == "France", 50, 40),
         Centre_max = ifelse(Population == "France", 300, 150),
         Centre_min = ifelse(Population == "France", 100, 0),
         Width_max = ifelse(Population == "France", 700, 360),
         Width_min = 2)
# Optimise allelic frequency variation models
cline_fitting_size <- optimise_clines(Priors = Priors_clines,
                                      logarithm = FALSE,
                                      batch_size = 100,
                                      genetic_data = data,
                                      SNP_subset = size_snps,
                                      meta_data = metadata) %>% 
  select_clinal_SNPs()

# Get the plotting values of these models along the transect
plotting_clines_size <- plot_clines(cline_fitting_size,
                                    genetic_data = data,
                                    SNP_subset = size_snps,
                                    meta_data = metadata)

# Plotting them
plotting_clines_size %>% 
  left_join(cline_fitting_size, by = c("Position", "Population")) %>% 
  filter(Best_model == "Clinal") %>% 
  ggplot() +
  geom_line(aes(x = LCmeanDist, y = Frequency, group = Position), lty = "dashed") +
  facet_row(vars(Population), scales = "free")

################## Make cline of QTLs for color ##################
# Candidate SNPs
color_snps <- GWAS_France %>% 
  filter(trait == "Shell_color_naive",
         grepl("LG6_", snp)) %>% 
  mutate(log_pval = -log10(pValue)) %>% 
  rename(Position = snp) %>% 
  filter(log_pval >= 5) %>% 
  mutate(pos = Position) %>% 
  separate(pos, c("SNP_name", "allele"), "\\.") %>% 
  select(Position, SNP_name)

# Prepare priors to fit allelic frequency variation models
Priors_clines_color <- get_allelic_frequencies(genetic_data = data,
                                             SNP_subset = color_snps,
                                             Extreme_values = pca$li,
                                             var = "Axis2",
                                             meta_data = metadata) %>% 
  mutate(Centre_prior = ifelse(Population == "France", 200, 70),
         Width_prior = ifelse(Population == "France", 50, 40),
         Centre_max = ifelse(Population == "France", 300, 150),
         Centre_min = ifelse(Population == "France", 100, 0),
         Width_max = ifelse(Population == "France", 700, 360),
         Width_min = 2)
# Optimise allelic frequency variation models
cline_fitting_color <- optimise_clines(Priors = Priors_clines_color,
                                      logarithm = FALSE,
                                      batch_size = 1000,
                                      genetic_data = data,
                                      SNP_subset = color_snps,
                                      meta_data = metadata) %>% 
  select_clinal_SNPs()
# Get the plotting values of these models along the transect
plotting_clines_color <- plot_clines(cline_fitting_color,
                                    genetic_data = data,
                                    SNP_subset = color_snps,
                                    meta_data = metadata)

# Plotting them
plotting_clines_color %>% 
  left_join(cline_fitting_color, by = c("Position", "Population")) %>% 
  filter(Best_model == "Clinal") %>% 
  ggplot() +
  geom_line(aes(x = LCmeanDist, y = Frequency, group = Position), lty = "dashed", alpha = 0.7) +
  facet_row(vars(Population), scales = "free")


################## Get the genotype of inversion for size  ##################
# Find the frequencies of the inversion in the sheltered part in sweden.
# To do so, we use the frequency of individuals presenting one version of the rearrangement in the transect ends.
# Therefore, we first must take the individuals that are on the left or on the right side of the transect.
# This is what is done here. For each of the 4 repetitions, we subset a group of individuals in the left 
# (= sheltered) and the right (= exposed) part of the transect
indivs_shelt_freq_sweden_14.1 <- clust_groups %>% 
  filter(Inversion == "Inv_14.1") %>% 
  mutate(modif_gr = ifelse(Group == 2, 3, ifelse(Group == 3, 2, 1)),
         Group = modif_gr) %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  filter(Population == "Sweden") %>% 
  left_join(pca$li %>% 
              rownames_to_column("Sample_Name"), by = "Sample_Name") %>% 
  arrange(LCmeanDist) %>% 
  dplyr::slice(1:30) %>% 
  arrange(Axis2) %>% 
  dplyr::slice(1:20) %>% 
  select(Sample_Name, Group)

# Once we did it for one end of one transect, we repeat it on the other ends of the two transects
indivs_exp_freq_sweden_14.1 <- clust_groups %>% 
  filter(Inversion == "Inv_14.1") %>% 
  mutate(modif_gr = ifelse(Group == 2, 3, ifelse(Group == 3, 2, 1)),
         Group = modif_gr) %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  filter(Population == "Sweden") %>% 
  left_join(pca$li %>% 
              rownames_to_column("Sample_Name"), by = "Sample_Name") %>% 
  arrange(LCmeanDist %>% desc) %>% 
  dplyr::slice(1:30) %>% 
  arrange(Axis2 %>% desc) %>% 
  dplyr::slice(1:20) %>% 
  select(Sample_Name, Group)
# Find the frequencies of the inversion in the sheltered part in france
indivs_shelt_freq_france_14.1 <- clust_groups %>% 
  filter(Inversion == "Inv_14.1") %>% 
  mutate(modif_gr = ifelse(Group == 2, 3, ifelse(Group == 3, 2, 1)),
         Group = modif_gr) %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  filter(Population == "France") %>% 
  left_join(pca$li %>% 
              rownames_to_column("Sample_Name"), by = "Sample_Name") %>% 
  arrange(LCmeanDist) %>% 
  dplyr::slice(1:30) %>% 
  arrange(Axis2) %>% 
  dplyr::slice(1:20) %>% 
  select(Sample_Name, Group)

indivs_exp_freq_france_14.1 <- clust_groups %>% 
  filter(Inversion == "Inv_14.1") %>% 
  mutate(modif_gr = ifelse(Group == 2, 3, ifelse(Group == 3, 2, 1)),
         Group = modif_gr) %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  filter(Population == "France") %>% 
  left_join(pca$li %>% 
              rownames_to_column("Sample_Name"), by = "Sample_Name") %>% 
  arrange(LCmeanDist %>% desc) %>% 
  dplyr::slice(1:30) %>% 
  arrange(Axis2 %>% desc) %>% 
  dplyr::slice(1:20) %>% 
  select(Sample_Name, Group)

# Here, we make a function to calculate the allelic frequencies in one end of the transect
calculate_freq_df <- function(df){
  sum_inv <- 0
  for (i in 1:nrow(df)){
    if (df$Group[i] == 3){
      sum_inv <- sum_inv + 2
    }else if (df$Group[i] == 2){
      sum_inv <- sum_inv + 1
    }
  }
  return(sum_inv/(2*nrow(df)))
}

# Then, we use the created function to calculate the allelic frequencies on the transect
shelt_freq_sweden_14.1 <- indivs_shelt_freq_sweden_14.1 %>% 
  calculate_freq_df()

exp_freq_sweden_14.1 <- indivs_exp_freq_sweden_14.1 %>% 
  calculate_freq_df()

shelt_freq_france_14.1 <- indivs_shelt_freq_france_14.1 %>% 
  calculate_freq_df()

exp_freq_france_14.1 <- indivs_exp_freq_france_14.1 %>% 
  calculate_freq_df()


# We use these frequencies as priors to fit the allelic frequency variations modles
priors_clines_inv_14.1 <- data.frame(
  "p_left_shelt" = c(shelt_freq_sweden_14.1, shelt_freq_france_14.1),
  "p_right_expos" = c(exp_freq_sweden_14.1, exp_freq_france_14.1),
  "Population" = c("Sweden", "France")
) %>% 
  mutate(Centre_prior = ifelse(Population == "France", 200, 70),
         Width_prior = ifelse(Population == "France", 50, 40) %>% log,
         Centre_max = ifelse(Population == "France", 300, 150),
         Centre_min = ifelse(Population == "France", 100, 0),
         Width_max = ifelse(Population == "France", 700, 360) %>% log,
         Width_min = 2 %>% log,
         p_left_max = p_left_shelt + 0.1 %>% logit,
         p_left_min = p_left_shelt - 0.1 %>% logit,
         p_right_max = p_right_expos + 0.1 %>% logit,
         p_right_min = p_right_expos - 0.1 %>% logit,
         p_left_shelt = max(p_left_shelt, 0.000001) %>% logit,
         p_right_expos = p_right_expos %>% logit
  )
# Fit the models
cline_inv_14.1 <- priors_clines_inv_14.1 %>% 
  rename(Pop = Population) %>% 
  group_by(Pop) %>% 
  bow(tie(Centre, Width, Left, Right) := mle2(clineflog,
                                              list(centre = Centre_prior,
                                                   width = Width_prior ,
                                                   left = p_left_shelt,
                                                   right = p_right_expos),
                                              data = list(x = clust_groups %>% 
                                                            filter(Inversion == "Inv_14.1") %>% 
                                                            left_join(metadata, by = "Sample_Name") %>% 
                                                            filter(Population == Pop) %>% 
                                                            select(LCmeanDist) %>% 
                                                            drop_na %>% 
                                                            as.vector %>%
                                                            unlist %>%
                                                            unname,
                                                          g = clust_groups %>% 
                                                            filter(Inversion == "Inv_14.1") %>% 
                                                            left_join(metadata, by = "Sample_Name") %>% 
                                                            filter(Population == Pop) %>%
                                                            select(Group) %>% 
                                                            mutate(Group = Group - 1) %>%
                                                            as.vector %>%
                                                            unlist %>%
                                                            unname,
                                                          n = 2),
                                              method = "L-BFGS-B",
                                              upper = list(centre = Centre_max,
                                                           width = Width_max,
                                                           left = p_left_max,
                                                           right = p_right_max),
                                              lower = list(centre = Centre_min,
                                                           width = Width_min,
                                                           left = p_left_min,
                                                           right = p_right_min)) %>%
        coef() %>% 
        round(digits = 3)) %>% 
  mutate(Width = exp(Width),
         Left = invlogit(Left),
         Right = invlogit(Right))


# And then, we get the plotting values of the optimised model
plot_inv_14.1 <- clinef(optimisation = FALSE,
                        x = seq((metadata %>% filter(Population == "France"))$LCmeanDist %>% min, (metadata %>% filter(Population == "France"))$LCmeanDist %>% max, 1),
                        g = (clust_groups %>% filter(Inversion == "Inv_14.1") %>% left_join(metadata, by = "Sample_Name") %>% mutate(Group = Group - 1) %>% 
                               filter(Population == "France"))$Group,
                        centre = (cline_inv_14.1 %>% filter(Pop == "France"))$Centre,
                        width = (cline_inv_14.1 %>% filter(Pop == "France"))$Width,
                        left = (cline_inv_14.1 %>% filter(Pop == "France"))$Left,
                        right = (cline_inv_14.1 %>% filter(Pop == "France"))$Right) %>% 
  mutate(Population = "France") %>% 
  rbind(clinef(optimisation = FALSE,
               x = seq((metadata %>% filter(Population == "Sweden"))$LCmeanDist %>% min, (metadata %>% filter(Population == "Sweden"))$LCmeanDist %>% max, 1),
               g = (clust_groups %>% filter(Inversion == "Inv_14.1") %>% left_join(metadata, by = "Sample_Name") %>% mutate(Group = Group - 1) %>% 
                      filter(Population == "Sweden"))$Group,
               centre = (cline_inv_14.1 %>% filter(Pop == "Sweden"))$Centre,
               width = (cline_inv_14.1 %>% filter(Pop == "Sweden"))$Width,
               left = (cline_inv_14.1 %>% filter(Pop == "Sweden"))$Left,
               right = (cline_inv_14.1 %>% filter(Pop == "Sweden"))$Right) %>% 
          mutate(Population = "Sweden"))


# That we plot with the previous candidate SNPs for the QTL in the identified inversion for size
B <- plotting_clines_size %>% 
  left_join(cline_fitting_size, by = c("Position", "Population")) %>% 
  filter(Best_model == "Clinal") %>% 
  mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France"))) %>% 
  ggplot() +
  geom_line(aes(x = LCmeanDist, y = Frequency, group = Position), lty = "dashed", lwd = 1.2) +
  geom_line(data = plot_inv_14.1 %>% 
              mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
                       factor(levels = c("Suède", "France"))),
            aes(x = position, y = phen_cline), lwd = 1.8, color = "aquamarine3") +
  facet_row(vars(Population), scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 30)) +
  labs(x = "Position le long du transect (m)",
       y = "Fréquence",
       tag = "(B)")


################## Get the genotype of inversion for color  ##################
# Then, we repeat this whole process for the inversion identified as having a color QTL
# Find the frequencies of the inversion in the sheltered part in sweden
indivs_shelt_freq_sweden_6.1 <- clust_groups %>% 
  filter(Inversion == "Inv_6.1") %>% 
  mutate(modif_gr = ifelse(Group == 2, 3, ifelse(Group == 3, 2, 1)),
         Group = modif_gr) %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  filter(Population == "Sweden") %>% 
  left_join(pca$li %>% 
              rownames_to_column("Sample_Name"), by = "Sample_Name") %>% 
  arrange(LCmeanDist) %>% 
  dplyr::slice(1:30) %>% 
  arrange(Axis2) %>% 
  dplyr::slice(1:20) %>% 
  select(Sample_Name, Group)
indivs_exp_freq_sweden_6.1 <- clust_groups %>% 
  filter(Inversion == "Inv_6.1") %>% 
  mutate(modif_gr = ifelse(Group == 2, 3, ifelse(Group == 3, 2, 1)),
         Group = modif_gr) %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  filter(Population == "Sweden") %>% 
  left_join(pca$li %>% 
              rownames_to_column("Sample_Name"), by = "Sample_Name") %>% 
  arrange(LCmeanDist %>% desc) %>% 
  dplyr::slice(1:30) %>% 
  arrange(Axis2 %>% desc) %>% 
  dplyr::slice(1:20) %>% 
  select(Sample_Name, Group)
# Find the frequencies of the inversion in the sheltered part in france
indivs_shelt_freq_france_6.1 <- clust_groups %>% 
  filter(Inversion == "Inv_6.1") %>% 
  mutate(modif_gr = ifelse(Group == 2, 3, ifelse(Group == 3, 2, 1)),
         Group = modif_gr) %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  filter(Population == "France") %>% 
  left_join(pca$li %>% 
              rownames_to_column("Sample_Name"), by = "Sample_Name") %>% 
  arrange(LCmeanDist) %>% 
  dplyr::slice(1:30) %>% 
  arrange(Axis2) %>% 
  dplyr::slice(1:20) %>% 
  select(Sample_Name, Group)
indivs_exp_freq_france_6.1 <- clust_groups %>% 
  filter(Inversion == "Inv_6.1") %>% 
  mutate(modif_gr = ifelse(Group == 2, 3, ifelse(Group == 3, 2, 1)),
         Group = modif_gr) %>% 
  left_join(metadata, by = "Sample_Name") %>% 
  filter(Population == "France") %>% 
  left_join(pca$li %>% 
              rownames_to_column("Sample_Name"), by = "Sample_Name") %>% 
  arrange(LCmeanDist %>% desc) %>% 
  dplyr::slice(1:30) %>% 
  arrange(Axis2 %>% desc) %>% 
  dplyr::slice(1:20) %>% 
  select(Sample_Name, Group)


shelt_freq_sweden_6.1 <- indivs_shelt_freq_sweden_6.1 %>% 
  calculate_freq_df()

exp_freq_sweden_6.1 <- indivs_exp_freq_sweden_6.1 %>% 
  calculate_freq_df()

shelt_freq_france_6.1 <- indivs_shelt_freq_france_6.1 %>% 
  calculate_freq_df()

exp_freq_france_6.1 <- indivs_exp_freq_france_6.1 %>% 
  calculate_freq_df()


# Prepare cline fitting
priors_clines_inv_6.1 <- data.frame(
  "p_left_shelt" = c(shelt_freq_sweden_6.1, shelt_freq_france_6.1),
  "p_right_expos" = c(exp_freq_sweden_6.1, exp_freq_france_6.1),
  "Population" = c("Sweden", "France")
) %>% 
  mutate(Centre_prior = ifelse(Population == "France", 200, 70),
         Width_prior = ifelse(Population == "France", 50, 40) %>% log,
         Centre_max = ifelse(Population == "France", 300, 150),
         Centre_min = ifelse(Population == "France", 100, 0),
         Width_max = ifelse(Population == "France", 700, 360) %>% log,
         Width_min = 2 %>% log,
         p_left_max = p_left_shelt + 0.1 %>% logit,
         p_left_min = p_left_shelt - 0.1 %>% logit,
         p_right_max = p_right_expos + 0.1 %>% logit,
         p_right_min = p_right_expos - 0.1 %>% logit,
         p_left_shelt = max(p_left_shelt, 0.000001) %>% logit,
         p_right_expos = min(p_right_expos, 0.999999) %>% logit
)
# Fit the clines
cline_inv_6.1 <- priors_clines_inv_6.1 %>% 
  rename(Pop = Population) %>% 
  group_by(Pop) %>% 
  bow(tie(Centre, Width, Left, Right) := mle2(clineflog,
                                              list(centre = Centre_prior,
                                                   width = Width_prior ,
                                                   left = p_left_shelt,
                                                   right = p_right_expos),
                                              data = list(x = clust_groups %>% 
                                                            filter(Inversion == "Inv_6.1") %>% 
                                                            left_join(metadata, by = "Sample_Name") %>% 
                                                            filter(Population == Pop) %>% 
                                                            select(LCmeanDist) %>% 
                                                            drop_na %>% 
                                                            as.vector %>%
                                                            unlist %>%
                                                            unname,
                                                          g = clust_groups %>% 
                                                            filter(Inversion == "Inv_6.1") %>% 
                                                            left_join(metadata, by = "Sample_Name") %>% 
                                                            filter(Population == Pop) %>%
                                                            select(Group) %>% 
                                                            mutate(Group = Group - 1) %>%
                                                            as.vector %>%
                                                            unlist %>%
                                                            unname,
                                                          n = 2),
                                              method = "L-BFGS-B",
                                              upper = list(centre = Centre_max,
                                                           width = Width_max,
                                                           left = p_left_max,
                                                           right = p_right_max),
                                              lower = list(centre = Centre_min,
                                                           width = Width_min,
                                                           left = p_left_min,
                                                           right = p_right_min)) %>%
        coef() %>% 
        round(digits = 3)) %>% 
  mutate(Width = exp(Width),
         Left = invlogit(Left),
         Right = invlogit(Right))



plot_inv_6.1 <- clinef(optimisation = FALSE,
                        x = seq((metadata %>% filter(Population == "France"))$LCmeanDist %>% min, (metadata %>% filter(Population == "France"))$LCmeanDist %>% max, 1),
                        g = (clust_groups %>% filter(Inversion == "Inv_6.1") %>% left_join(metadata, by = "Sample_Name") %>% mutate(Group = Group - 1) %>% 
                               filter(Population == "France"))$Group,
                        centre = (cline_inv_6.1 %>% filter(Pop == "France"))$Centre,
                        width = (cline_inv_6.1 %>% filter(Pop == "France"))$Width,
                        left = (cline_inv_6.1 %>% filter(Pop == "France"))$Left,
                        right = (cline_inv_6.1 %>% filter(Pop == "France"))$Right) %>% 
  mutate(Population = "France") %>% 
  rbind(clinef(optimisation = FALSE,
               x = seq((metadata %>% filter(Population == "Sweden"))$LCmeanDist %>% min, (metadata %>% filter(Population == "Sweden"))$LCmeanDist %>% max, 1),
               g = (clust_groups %>% filter(Inversion == "Inv_6.1") %>% left_join(metadata, by = "Sample_Name") %>% mutate(Group = Group - 1) %>% 
                      filter(Population == "Sweden"))$Group,
               centre = (cline_inv_6.1 %>% filter(Pop == "Sweden"))$Centre,
               width = (cline_inv_6.1 %>% filter(Pop == "Sweden"))$Width,
               left = (cline_inv_6.1 %>% filter(Pop == "Sweden"))$Left,
               right = (cline_inv_6.1 %>% filter(Pop == "Sweden"))$Right) %>% 
          mutate(Population = "Sweden"))


# Put the two together
E <- plotting_clines_color %>% 
  left_join(cline_fitting_color, by = c("Position", "Population")) %>% 
  filter(Best_model == "Clinal") %>% 
  mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France"))) %>% 
  ggplot() +
  geom_line(aes(x = LCmeanDist, y = Frequency, group = Position), lty = "dashed", lwd = 1.2, alpha = 0.5) +
  geom_line(data = plot_inv_6.1 %>% 
              mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
                       factor(levels = c("Suède", "France"))),
            aes(x = position, y = phen_cline), lwd = 1.8, color = "orange") +
  facet_row(vars(Population), scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 30)) +
  labs(x = "Position le long du transect (m)",
       y = "Fréquence",
       tag = "(E)")


################## Arrange the plots together  ##################
ggarrange(A, B, D, E, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right", widths = c(1, 2))

