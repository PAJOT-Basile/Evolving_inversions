# Import libraries
libraries <- c("tidyverse", "adegenet", "vcfR", "readxl", "statgenGWAS", "ggforce", "ggh4x", "ggpubr")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(characters = libraries, character.only = TRUE)
rm(libraries)


################## Useful variables  ##################
# Palette to use for the chromosome colors
chromosome_palette <- c("grey71", "turquoise4")

# Theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

# Position of the inversions
position_inversions <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
                                  sep = "\t", header = TRUE)

grouped_inversions <- position_inversions %>% 
  filter(Inversion_grouped %in% c("Inv_4.1", "Inv_14.1", "Inv_3.1")) %>% 
  mutate(Inversion = Inversion_grouped)

################## Useful functions  ##################
source("../General_scripts/Functions_optimise_plot_clines.r")

################## Import the vcf file  ##################
data <- read.vcfR("../../Output/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz") %>% 
  vcfR2genind()

# We add the poputations to the vcf object
data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor
data@other$exposition <- (data@tab %>% rownames %>% str_split_fixed(., "_", 6))[, 5] %>% as.factor
data@other$exposition[which(data@other$exposition == "TRANSI")] <- "TRANS" %>% as.factor
data@other$exposition <- data@other$exposition %>% droplevels()


############################ Scale the genome ##########################

# Scale the genome to get rid of missing data
X <- scaleGen(data, NA.method="mean", scale=FALSE, center=TRUE)

pca <- dudi.pca(X, scale=TRUE, nf=5, scannf=FALSE)

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

################## Make a map of the position of each snp on the chromosomes  ##################
Delta_freqs_whole_genome <- get_delta_freqs_and_F4(genetic_data = data,
                                                   SNP_subset = NULL,
                                                   nb_extreme_indivs = 30,
                                                   nb_indivs_to_keep = 20,
                                                   meta_data = metadata) %>% 
  select(-c(starts_with("p_"))) %>% 
  mutate(across(contains("Delta_"), .fns = abs, .names = "{.col}_abs")) %>% 
  # Calculate the mean delta freq
  mutate(Mean_delta_freq = rowMeans(select(., ends_with("_abs")), na.rm = TRUE)) %>% 
  select_good_SNPs(Delta_freq_Sweden)


map_chromosomes <- Delta_freqs_whole_genome %>% 
  select(Position) %>% 
  mutate(pos = Position) %>% 
  separate(pos, c("LG", "SUPER", "SUPER_frag", "pos"), "_") %>% 
  separate(pos, c("pos", "allele"), "\\.") %>% 
  select(-allele) %>% 
  unite(chr, LG, SUPER, SUPER_frag) %>% 
  mutate(pos = pos %>% as.numeric) %>% 
  column_to_rownames("Position")

################## Matrix of genotypes for the french population  ##################
# We use the scaled genotype to get rid of the missing values in the genotype matrix
## France
geno_france <- X %>% 
  as.data.frame %>% 
  rownames_to_column("Sample_Name") %>% 
  inner_join(metadata %>% 
               filter(Population == "France"), by = "Sample_Name") %>%
  select(Sample_Name, rownames(map_chromosomes)) %>% 
  column_to_rownames("Sample_Name")

## Sweden
geno_sweden <- X %>% 
  as.data.frame %>% 
  rownames_to_column("Sample_Name") %>% 
  inner_join(metadata %>% 
               filter(Population == "Sweden"), by = "Sample_Name") %>%
  select(Sample_Name, rownames(map_chromosomes)) %>% 
  column_to_rownames("Sample_Name")

################## Phenotypes  ##################
# The gData object needs a dataframe containing the trait(s) to analyse and the name
# of the individuals in the first column. This column must be called "genotype".
## France
Pheno_france <- metadata %>% 
  filter(Population == "France") %>% 
  select(Sample_Name, Shell_color_naive, Length) %>% 
  mutate(Shell_color_naive = ifelse(Shell_color_naive == "Yellow", 0, 1),
         indivs = Sample_Name) %>% 
  rename(genotype = indivs) %>% 
  relocate(genotype) %>% 
  column_to_rownames("Sample_Name")


## Sweden
Pheno_sweden <- metadata %>% 
  filter(Population == "Sweden") %>% 
  select(Sample_Name, Shell_color_naive, Length) %>% 
  mutate(Shell_color_naive = ifelse(Shell_color_naive == "Yellow", 0, 1),
         indivs = Sample_Name) %>% 
  rename(genotype = indivs) %>% 
  relocate(genotype) %>% 
  column_to_rownames("Sample_Name")

################## Scaled relatedness matrix (calculated from vcftools)  ##################
# First, we need to select the reference individuals needed to calculate the relatedness and
# save them in text files
metadata %>%
 filter(Population == "Sweden") %>%
 select(Sample_Name) %>%
 write.table("../../Output/Data/Reference_indivs_expos/Swedish_pop.txt",
             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

metadata %>%
 filter(Population == "France") %>%
 select(Sample_Name) %>%
 write.table("../../Output/Data/Reference_indivs_expos/French_pop.txt",
             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
## First, we import the relatedness calculated with vcftools with this command lines
## France
kin_france <- read.table("../../Output/Sweden_France_parallelism/10_GWAS/Relatedness/French_scaled_relatedness.tsv", header = FALSE)

## Sweden
kin_sweden <- read.table("../../Output/Sweden_France_parallelism/10_GWAS/Relatedness/Swedish_scaled_relatedness.tsv", header = FALSE)

################## Make the gData objects  ##################
## France
gData_france <- createGData(geno = geno_france, map = map_chromosomes, pheno = Pheno_france, kin = kin_france)

## Sweden
gData_sweden <- createGData(geno = geno_sweden, map = map_chromosomes, pheno = Pheno_sweden, kin = kin_sweden)

################## Run the GWAS  ##################
## France
GWAS_france <- runSingleTraitGwas(gData_france)

## Sweden
GWAS_sweden <- runSingleTraitGwas(gData_sweden)

# Save the GWAS output
GWAS_france$GWAResult$Pheno_france %>% 
  select(trait, snp, chr, pos, pValue) %>% 
  rename(Trait = trait,
         Position = snp,
         Chromosome = chr) %>% 
  mutate(pValue = -log10(pValue)) %>% 
  rename(logp_val = pValue) %>% 
  mutate(Population = "France") %>% 
  rbind(GWAS_sweden$GWAResult$Pheno_sweden %>% 
          select(trait, snp, chr, pos, pValue) %>% 
          rename(Trait = trait,
                 Position = snp,
                 Chromosome = chr) %>% 
          mutate(pValue = -log10(pValue)) %>% 
          rename(logp_val = pValue) %>% 
          mutate(Population = "Sweden")) %>% 
  write.table("../../Output/Sweden_France_parallelism/10_GWAS/GWAS_output.tsv",
              col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
