# Import libraries
require("anyLib")
anyLib(c("tidyverse", "adegenet", "vcfR", "readxl", "statgenGWAS", "ggforce", "ggh4x", "ggpubr"))


################## Useful variables  ##################
# Palette to use for the chromosome colors
chromosome_palette <- c("grey71", "turquoise4")
################## Useful functions  ##################
"%!in%" <- function(x, y){!(x %in% y)}

source("../General_scripts/Functions_optimise_plot_clines.r")

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


############################ Scale the genome ##########################

# Scale the genome to get rid of missing data
X <- scaleGen(data, NA.method="mean", scale=FALSE, center=TRUE)

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

################## Import the inversion delimitations  ##################
# First, we read the table for the inversion delimitations
Delim_inversions <- read.csv("../Data/Delim_inversions.csv", header = TRUE, sep = "\t")
# Then, we make a table containing the cumulative position of each chromosome
x <- Delta_freqs_whole_genome %>%
  select(Position, F4_stat) %>%
  transform_position_ade2tidy() %>% 
  separate(Chromosome, c("Chromosome", "SUPER", "SUPER_frag"), "_") %>% 
  mutate(Chromosome = Chromosome %>% 
           factor(levels = Delta_freqs_whole_genome %>%
                    select(Position, F4_stat) %>%
                    transform_position_ade2tidy() %>% 
                    separate(Chromosome, c("Chromosome", "SUPER", "SUPER_frag"), "_") %>%
                    select(Chromosome) %>% 
                    unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                    as.vector %>% unname %>% unlist)) %>% 
  group_by(Chromosome) %>% 
  summarise(bc_cum = max(Position)) %>% 
  mutate(bp_add = lag(cumsum(bc_cum), default = 0))

# We add this to the inversion delimitations to know where to place the inversions in the manhattan plot and prepare 
# the table to be used with the geom_polygon function, requiring that the coordinates of the vertexes be specified.
x1 <- Delim_inversions %>%
  left_join(x, relationship = "many-to-one") %>%
  mutate(bp_Pos_min_inv = Pos_min_inv + bp_add,
         bp_Pos_max_inv = Pos_max_inv + bp_add) %>% 
  select(Chromosome, bp_Pos_min_inv, bp_Pos_max_inv, Inversion) %>% 
  pivot_longer(starts_with("bp_"), names_to = "Position", values_to = "x_value") %>% 
  group_by(Chromosome, Inversion) %>% 
  rbind(., .) %>% 
  arrange(Chromosome, x_value) %>% 
  cbind("y_value" = rep(c(0, 0.5, 0.5, 0), nrow(.)/4))
################## Make a map of the position of each snp on the chromosomes  ##################
Delta_freqs_whole_genome <- get_delta_freqs_and_F4(genetic_data = data,
                                                   SNP_subset = NULL,
                                                   nb_extreme_indivs = 30,
                                                   nb_indivs_to_keep = 20,
                                                   metadata = metadata) %>% 
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
# First, we import the relatedness calculated with vcftools with this command lines
## France
#vcf_file <- "/shared/home/bpajot/littorina/finalresult/bpajot/genomic_analysis/filtering_vcf_files/Final_outputs/Fully_filtered_thinned_Hobs.vcf.gz"
#vcftools <- "/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools"
#system2(vcftools,args =c(paste0("--gzvcf ",
#                                vcf_file,
#                                " --relatedness",
#                                " --keep /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/FST/French_pop.txt",
#                                " --out /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/GWAS/French_relatedness")))
relatedness_france <- read.table("../Data/French_relatedness.relatedness", header = TRUE) %>% 
  rename(Relatedness = RELATEDNESS_AJK) %>% 
  pivot_wider(names_from = INDV2, values_from = Relatedness) %>% 
  column_to_rownames("INDV1")

## Sweden
#system2(vcftools,args =c(paste0("--gzvcf ",
#                                vcf_file,
#                                " --relatedness",
#                                " --keep /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/FST/Swedish_pop.txt",
#                                " --out /shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/GWAS/Swedish_relatedness")))
relatedness_sweden <- read.table("../Data/Swedish_relatedness.relatedness", header = TRUE) %>% 
  rename(Relatedness = RELATEDNESS_AJK) %>% 
  pivot_wider(names_from = INDV2, values_from = Relatedness) %>% 
  column_to_rownames("INDV1")

# The matrix is not symmetrical but triangular, so, we make it symmetrical by hand
## France
for (i in 1:nrow(relatedness_france)){
  for (j in 1:ncol(relatedness_france)){
    if (i == j){
      relatedness_france[i, j] <- 1
    }else if (relatedness_france[i, j] %>% is.na){
      relatedness_france[i, j] <- relatedness_france[j, i]
    }
  }
}

## Sweden
for (i in 1:nrow(relatedness_sweden)){
  for (j in 1:ncol(relatedness_sweden)){
    if (i == j){
      relatedness_sweden[i, j] <- 1
    }else if (relatedness_sweden[i, j] %>% is.na){
      relatedness_sweden[i, j] <- relatedness_sweden[j, i]
    }
  }
}

# Now, we scale the matrix to only have a positive relatedness between 0 and 1
## France
min_rel_france <- relatedness_france %>% min
max_rel_france <- relatedness_france %>% max
kin_france <- relatedness_france %>% 
  mutate(across(everything(), .fns = function(x) (x + abs(min_rel_france)) / (max_rel_france - min_rel_france))) %>% 
  as.matrix

## Sweden
min_rel_sweden <- relatedness_sweden %>% min
max_rel_sweden <- relatedness_sweden %>% max
kin_sweden <- relatedness_sweden %>% 
  mutate(across(everything(), .fns = function(x) (x + abs(min_rel_sweden)) / (max_rel_sweden - min_rel_sweden))) %>% 
  as.matrix

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

################## Plot the results  ##################
### Size
Size_plot <- GWAS_france %>% 
                filter(trait  == "Length") %>% 
                select(snp, pValue) %>% 
                rename(Position = snp) %>% 
                mutate(logp_val = -log10(pValue),
                       Population = "France") %>% 
                rbind(GWAS_sweden %>%
                        filter(trait == "Length") %>% 
                        select(snp, pValue) %>% 
                        rename(Position = snp) %>% 
                        mutate(logp_val = -log10(pValue),
                               Population = "Sweden")) %>% 
  mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France"))) %>% 
  geom_manhattan(aes(y = logp_val, supp = Population), palette = chromosome_palette) +
  geom_hline(yintercept = 5, color = "black", lwd = 0.9, lty = "dashed") +
  geom_polygon(data = x1, aes(x=x_value, y=y_value, group = Inversion), color = NA, fill = "grey10") +
  facet_col(vars(Population)) +
  labs(x = "Chromosome",
       y = expression("log"[10]* "(p-value)"),
       tag = "(A)") +
  ylim(0, 7) +
  theme(text = element_text(size = 20))

### Color
Color_plot <- GWAS_france %>% 
                 filter(trait  == "Shell_color_naive") %>% 
                 select(snp, pValue) %>% 
                 rename(Position = snp) %>% 
                 mutate(logp_val = -log10(pValue),
                        Population = "France") %>% 
                 rbind(GWAS_sweden %>%
                         filter(trait == "Shell_color_naive") %>% 
                         select(snp, pValue) %>% 
                         rename(Position = snp) %>% 
                         mutate(logp_val = -log10(pValue),
                                Population = "Sweden")) %>% 
                 mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
                          factor(levels = c("Suède", "France"))) %>% 
                 geom_manhattan(aes(y = logp_val, supp = Population), palette = chromosome_palette) +
  geom_hline(yintercept = 5, color = "black", lwd = 0.9, lty = "dashed") +
  geom_polygon(data = x1, aes(x=x_value, y=y_value, group = Inversion), color = NA, fill = "grey10") +
  facet_col(vars(Population)) +
  labs(x = "Chromosome",
       y = expression("log"[10]* "(p-value)"),
       tag = "(B)") +
  ylim(0, 25) +
  theme(text = element_text(size = 20))


# Merge them together
ggarrange(Size_plot,
          Color_plot,
          ncol = 1)

