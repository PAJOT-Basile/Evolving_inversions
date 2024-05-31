# Libraries
require("anyLib")
anyLib(c("tidyverse", "vcfR", "adegenet", "readxl"))

############################ Import vcf file ##########################
# Import vcf file
data <- read.vcfR("C:/Documents/Stage_M2_2024/Project/Genetic data/Fully_filtered_thinned_Hobs.vcf.gz") %>% 
  vcfR2genind()

# We add the poputations to the vcf object
data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor
data@other$exposition <- (data@tab %>% rownames %>% str_split_fixed(., "_", 6))[, 5] %>% as.factor
data@other$exposition[which(data@other$exposition == "TRANSI")] <- "TRANS" %>% as.factor
data@other$exposition <- data@other$exposition %>% droplevels()

############################ Import metadata ##########################
metadata <- read_excel(path = "../../Data/data_Fabalis_resequencing_Basile.xlsx",
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

########################### Import delimitations of the inversion ##########################
Delim_inversions <- read.csv("../../Data/Delimitation_inversions.csv", header = TRUE, sep = "\t")
############################ Make a local pca on the inversions ########################## Separate the inversions from the rest of the genome
local_pca_inversions_indivs <- data.frame()
local_pca_inversions_contribs <- data.frame()
clust_groups <- data.frame()
# Here, we iterate over each inversion and run a local pca on it
for (inversion in Delim_inversions$Inversion %>% as.character){
  print(inversion)
  min_pos_inv <- Delim_inversions %>% 
    ungroup %>% 
    filter(Inversion == inversion) %>% 
    select(Pos_min_inv) %>% 
    as.vector %>% unname %>% unlist
  max_pos_inv <- Delim_inversions %>% 
    ungroup %>% 
    filter(Inversion == inversion) %>% 
    select(Pos_max_inv) %>% 
    as.vector %>% unname %>% unlist
  chrom <- Delim_inversions %>% 
    ungroup %>% 
    filter(Inversion == inversion) %>% 
    select(Chromosome) %>% 
    as.vector %>% unname %>% unlist
  
  snp_names <- data@loc.fac %>% 
    as_tibble %>% 
    filter(grepl(paste0("^", chrom), value)) %>% 
    mutate(value = value %>% as.character) %>% 
    mutate(Position = str_split_fixed(value, "_", 4)[, 4] %>% as.numeric) %>% 
    filter((Position > min_pos_inv) & (Position < max_pos_inv)) %>% 
    as.vector %>% unlist %>% unname
  print("    Made list of SNP names")
  # Select the SNPs that are located in the inversions
  genetic_data <- data[loc = snp_names] %>% suppressWarnings
  print("    Sampled data")
  
  X <- scaleGen(genetic_data, NA.method="mean", scale=FALSE, center=TRUE)
  print("    Scaled")
  
  # Run the PCA
  pca <- dudi.pca(X, scale=TRUE, nf=5, scannf=FALSE)
  print("    PCA done")
  
  clust_fr <- pca$li %>% 
    rownames_to_column("Sample_Name") %>% 
    left_join(metadata, by = "Sample_Name") %>% 
    filter(Population == "France") %>% 
    select(Sample_Name, contains("Axis")) %>% 
    column_to_rownames("Sample_Name") %>% 
    find.clusters(x = .,
                  n.pca = 5,
                  n.clust = 3)
  
  clust_sw <- pca$li %>% 
    rownames_to_column("Sample_Name") %>% 
    left_join(metadata, by = "Sample_Name") %>% 
    filter(Population == "Sweden") %>% 
    select(Sample_Name, contains("Axis")) %>% 
    column_to_rownames("Sample_Name") %>% 
    find.clusters(x = .,
                  n.pca = 5,
                  n.clust = 3)
  print("    Found clusters")
  
  local_pca_inversions_indivs <- add_table_to_df_in_iteration(local_pca_inversions_indivs, pca$li %>% 
                                 mutate(Inversion = inversion,
                                        Chromosome = chrom) %>% 
                                   rownames_to_column("Sample_Name"))
  local_pca_inversions_contribs <- add_table_to_df_in_iteration(local_pca_inversions_contribs, pca$co %>% 
                                                                  mutate(Inversion = inversion,
                                                                         Chromosome = chrom) %>% 
                                                                  rownames_to_column("Position"))
  clust_groups <- clust_fr$grp %>% 
    as.data.frame() %>% 
    rename(Group = ".") %>% 
    rownames_to_column("Sample_Name") %>% 
    mutate(Inversion = inversion) %>% 
    rbind(clust_sw$grp %>% 
            as.data.frame %>% 
            rename(Group = ".") %>% 
            rownames_to_column("Sample_Name") %>% 
            mutate(Inversion = inversion)) %>% 
    add_table_to_df_in_iteration(clust_groups, .)
  print(nrow(clust_groups))
  
}

# Clean environment
rm(chrom, max_pos_inv, min_pos_inv, inversion, snp_names, pca, X, clust_fr, clust_sw)

# Plot the results
fig_sup_3 <- local_pca_inversions_indivs %>%
  mutate(Inversion = Inversion %>% 
           factor(levels = Delim_inversions %>%
                    ungroup %>% 
                    select(Inversion) %>% 
                    unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist)) %>% 
  left_join(metadata, by = "Sample_Name", relationship = "many-to-one") %>% 
  mutate(Habitat = ifelse(Habitat == "Exposed", "Exposé", ifelse(Habitat == "Sheltered", "Abrité", "Transition")) %>% 
           factor(levels = c("Exposé", "Transition", "Abrité")),
         Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France")),
         Shell_color_naive = ifelse(Shell_color_naive == "Brown", "Marron", 'Jaune') %>% 
           factor(levels = c("Marron", "Jaune"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis2, color = Habitat), size = 3, alpha = 0.5) +
  scale_color_manual(name = "Habitat",
                     values = c("Exposé" = "orange2", "Transition" = "deeppink", "Abrité" = "dodgerblue3")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis1, y = Axis2, color = Population),
                    lwd = 1.2, label.fontsize = 20) +
  scale_color_manual(values = c("Suède" = "forestgreen", "France" = "darkorchid4")) +
  facet_wrap(vars(Inversion)) +
  theme_bw() +
  theme(text = element_text(size = 30)) +
  labs(x = "Axe 1",
       y = "Axe 2")

fig_sup_4 <- local_pca_inversions_indivs %>%
  mutate(Inversion = Inversion %>% 
           factor(levels = Delim_inversions %>%
                    ungroup %>% 
                    select(Inversion) %>% 
                    unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist)) %>% 
  left_join(metadata, by = "Sample_Name", relationship = "many-to-one") %>% 
  mutate(Habitat = ifelse(Habitat == "Exposed", "Exposé", ifelse(Habitat == "Sheltered", "Abrité", "Transition")) %>% 
           factor(levels = c("Exposé", "Transition", "Abrité")),
         Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France")),
         Shell_color_naive = ifelse(Shell_color_naive == "Brown", "Marron", 'Jaune') %>% 
           factor(levels = c("Marron", "Jaune"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis3, color = Habitat), size = 3, alpha = 0.5) +
  scale_color_manual(name = "Habitat",
                     values = c("Exposé" = "orange2", "Transition" = "deeppink", "Abrité" = "dodgerblue3")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis1, y = Axis3, color = Population),
                    lwd = 1.2, label.fontsize = 20) +
  scale_color_manual(values = c("Suède" = "forestgreen", "France" = "darkorchid4")) +
  facet_wrap(vars(Inversion)) +
  theme_bw() +
  theme(text = element_text(size = 30)) +
  labs(x = "Axe 1",
       y = "Axe 3")
