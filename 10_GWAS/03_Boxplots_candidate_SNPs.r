########## Libraries ##########
libraries <- c("tidyverse", "readxl", "vcfR", "adegenet", "ggh4x")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(libraries, character.only = TRUE)
rm(libraries)


########## Useful functions ##########
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")
########## Useful variables ##########
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

########## Import data ##########
GWAS <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/GWAS/GWAS_output.tsv",
                   sep = "\t", header = TRUE) %>% 
  rename(log_pval = logp_val)

# Import metadata
metadata <- read_excel(path = "/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Phenotypic/data_Fabalis_resequencing_Basile.xlsx",
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


# Import the genetic data
data <- read.vcfR("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz") %>% 
  vcfR2genind()

# We add the poputations and expositions to the vcf object
data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor
data@other$exposition <- (data@tab %>% rownames %>% str_split_fixed(., "_", 6))[, 5] %>% as.factor
data@other$exposition[which(data@other$exposition == "TRANSI")] <- "TRANS" %>% as.factor
data@other$exposition <- data@other$exposition %>% droplevels()
########## Boxplots of sizes ##########
# Candidate SNPs for size determinism in Sweden
candidate_SNPs_Sweden <- GWAS %>%
    filter(Trait == "Length",
         Population == "Sweden",
         log_pval >= 5) %>% 
  select(Position) %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  unite(SNP_Name, Chromosome, Position, sep = "_") %>% 
  rename(Position = pos)

# Candidate SNPs for size determinism in France
candidate_SNPs_France <- GWAS %>% 
  filter(Trait == "Length",
         Population == "France",
         log_pval >= 5) %>% 
  select(Position) %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  unite(SNP_Name, Chromosome, Position, sep = "_") %>% 
  rename(Position = pos)


# Get the genotypes of the SNPs both in France and Sweden
genotypes_candidate_SNPs <- get_genotype_transect(genetic_data = data,
                                                  SNP_subset = candidate_SNPs_Sweden %>% 
                                                    rbind(candidate_SNPs_France),
                                                  meta_data = metadata) %>% 
  rename_all(list(~ str_replace_all(., ".0$|.1$", "")))

# Prepare the order of the positions to be plotted
order_snps <- c(
  candidate_SNPs_Sweden %>% 
    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", candidate_SNPs_Sweden$SNP_Name))) %>% 
    pull(SNP_Name),
  candidate_SNPs_France %>% 
    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", candidate_SNPs_France$SNP_Name))) %>% 
    pull(SNP_Name)
)

# Plot the results
(genotypes_candidate_SNPs %>% 
    pivot_longer(starts_with("LG"), values_to = "Genotype", names_to = "Position") %>% 
    mutate(Genotype = (Genotype / 2) %>% 
             as.character(),
           Position = Position %>% 
             factor(levels = order_snps)) %>% 
    drop_na %>% 
    left_join(metadata,
              by = c("Sample_Name", "Population"),
              relationship = "many-to-one") %>% 
    mutate(Coloration = case_when(
      Population == "Sweden" & Position %in% candidate_SNPs_Sweden$SNP_Name ~ "Dark",
      Population == "France" & Position %in% candidate_SNPs_France$SNP_Name ~ "Dark",
      TRUE ~ "Light"
    )) %>% 
    ggplot() +
    geom_boxplot(aes(x = Genotype, y = Length, fill = Coloration)) +
    scale_fill_manual(values = c("Dark" = "grey30", "Light" = "white"),
                      guide = "none") +
    my_theme +
    labs(y = "Shell size (mm)") +
    facet_grid2(Population ~ Position)) %>% 
  ggsave(filename = "/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/GWAS/Size_QTL_boxplots.png",
         height = 800, width = 1800, units = "px", scale = 5)

########## Boxplots of colors ##########
# Candidate SNPs for color determinism in France
candidate_SNPs_France <- GWAS %>% 
  filter(Trait == "Shell_color_naive",
         Population == "France",
         log_pval >= 18) %>% 
  select(Position) %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  unite(SNP_Name, Chromosome, Position, sep = "_") %>% 
  rename(Position = pos)

# Candidate SNPs for color determinism in Sweden
candidate_SNPs_Sweden <- GWAS %>% 
  filter(Trait == "Shell_color_naive",
         Population == "Sweden",
         (log_pval >= 10) | (log_pval > 6 & grepl("LG6_", Position))) %>% 
  select(Position) %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  unite(SNP_Name, Chromosome, Position, sep = "_") %>% 
  rename(Position = pos)

# Get the genotypes of the SNPs both in France and Sweden
genotypes_candidate_SNPs <- get_genotype_transect(genetic_data = data,
                                                  SNP_subset = candidate_SNPs_Sweden %>% 
                                                    rbind(candidate_SNPs_France),
                                                  meta_data = metadata) %>% 
  rename_all(list(~ str_replace_all(., ".0$|.1$", "")))

# Prepare the order of the positions to be plotted
order_snps <- c(
  candidate_SNPs_Sweden %>% 
    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", candidate_SNPs_Sweden$SNP_Name))) %>% 
    pull(SNP_Name),
  candidate_SNPs_France %>% 
    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", candidate_SNPs_France$SNP_Name))) %>% 
    pull(SNP_Name)
)

# Plot the results
# As there are a lot of SNPs, we separate the plotting into two graphs
## First Sweden
(genotypes_candidate_SNPs %>% 
    pivot_longer(starts_with("LG"), values_to = "Genotype", names_to = "Position") %>% 
    mutate(Genotype = (Genotype / 2),
           Position = Position %>% 
             factor(levels = order_snps)) %>% 
    drop_na %>% 
    inner_join(metadata %>% 
                 filter(Population == "Sweden"),
               by = c("Sample_Name", "Population"),
               relationship = "many-to-one") %>% 
    mutate(Coloration = case_when(
      Population == "Sweden" & Position %in% candidate_SNPs_Sweden$SNP_Name ~ "Dark",
      TRUE ~ "Light"
    )) %>% 
    ggplot() +
    geom_boxplot(aes(x = Shell_color_naive, y = Genotype, fill = Coloration, color = Shell_color_naive),
                 lwd = 1.1) +
    scale_fill_manual(values = c("Dark" = "grey30", "Light" = "white"),
                      guide = "none") +
    scale_color_manual(values = c("Brown" = "brown", "Yellow" = "orange"),
                       guide = "none") +
    my_theme +
    labs(x = "Shell color") +
    facet_wrap(vars(Position))) %>% 
  ggsave(filename = "/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/GWAS/Color_QTL_boxplots_Sweden.png",
         height = 1200, width = 1800, units = "px", scale = 5)

## Then France
(genotypes_candidate_SNPs %>% 
    pivot_longer(starts_with("LG"), values_to = "Genotype", names_to = "Position") %>% 
    mutate(Genotype = (Genotype / 2),
           Position = Position %>% 
             factor(levels = order_snps)) %>% 
    drop_na %>% 
    inner_join(metadata %>% 
                 filter(Population == "France"),
               by = c("Sample_Name", "Population"),
               relationship = "many-to-one") %>% 
    mutate(Coloration = case_when(
      Population == "France" & Position %in% candidate_SNPs_France$SNP_Name ~ "Dark",
      TRUE ~ "Light"
    )) %>% 
    ggplot() +
    geom_boxplot(aes(x = Shell_color_naive, y = Genotype, fill = Coloration, color = Shell_color_naive),
                 lwd = 1.1) +
    scale_fill_manual(values = c("Dark" = "grey30", "Light" = "white"),
                      guide = "none") +
    scale_color_manual(values = c("Brown" = "brown", "Yellow" = "orange"),
                       guide = "none") +
    my_theme +
    labs(x = "Shell color") +
    facet_wrap(vars(Position))) %>% 
  ggsave(filename = "/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/GWAS/Color_QTL_boxplots_France.png",
         height = 1200, width = 1800, units = "px", scale = 5)
