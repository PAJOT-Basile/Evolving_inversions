################## Libraries  ##################
libraries <- c("tidyverse", "statgenGWAS", "vcfR", "adegenet", "readxl")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)
if (!require("patchwork")) devtools::install_github("thomasp85/patchwork")


################## Useful variables  ##################
# Palette to use for the chromosome colors
chromosome_palette <- c("grey71", "turquoise4")

# Palette to use to color the points depending on the size of the individuals
color_palette <-  c("Brown" = "brown", "Yellow" = "orange")

# Theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

################## Useful functions  ##################
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

################## Import data  ##################
# Positions of the inversions
position_inversions <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
                                  sep = "\t", header = TRUE) %>% 
  # Prepare the variables to trace the trees
  mutate(Is_inversion = TRUE,
         Is_split = case_when(
           Inversion_grouped %in% c("Inv_3.1", "Inv_4.1", "Inv_14.1") ~ TRUE,
           TRUE ~ FALSE
         ))

# Isolate the names of the grouped inversions
grouped_inversions <- position_inversions %>% 
  filter(Inversion_grouped %in% c("Inv_4.1", "Inv_14.1", "Inv_3.1")) %>% 
  mutate(Inversion = Inversion_grouped)

# Import the vcf file
vcf_file <- "/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz"
data <- read.vcfR(vcf_file) %>% 
  vcfR2genind()

# We add the poputations to the vcf object
data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor
# And the exposition of the individuals
data@other$exposition <- (data@tab %>% rownames %>% str_split_fixed(., "_", 6))[, 5] %>% as.factor
data@other$exposition[which(data@other$exposition == "TRANSI")] <- "TRANS" %>% as.factor
data@other$exposition <- data@other$exposition %>% droplevels()

# Import metadata on the individuals
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

# Import the delta freqs for the whole genome
Delta_freqs_whole_genome <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/Delta_freqs_whole_genome.tsv",
                                       sep = "\t", header = TRUE) %>% 
  select_good_SNPs(Delta_freq_Sweden)

# Import the local PCA per inversion
groups_pca <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
                         sep = "\t", header = TRUE) %>% 
  # Filter to get only the inversions that are kept after the 4 step filtering of inversions
  inner_join(position_inversions,
             by = c("Chromosome", "Inversion", "Population"),
             relationship = "many-to-many")

# Import the reference individuals
ref_individuals <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/Reference_indivs/France_exposed.txt",
                              sep = "\t", header = FALSE) %>%
  mutate(Population = "France",
         Exposition = "Exposed") %>%
  rbind(read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/Reference_indivs/France_sheltered.txt",
                   sep = "\t", header = FALSE) %>%
          mutate(Population = "France",
                 Exposition = "Sheltered")) %>%
  rbind(read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/Reference_indivs/Sweden_exposed.txt",
                   sep = "\t", header = FALSE) %>%
          mutate(Population = "Sweden",
                 Exposition = "Exposed")) %>%
  rbind(read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/Reference_indivs/Sweden_sheltered.txt",
                   sep = "\t", header = FALSE) %>%
          mutate(Population = "Sweden",
                 Exposition = "Sheltered")) %>%
  rename(Sample_Name = V1)
############################ Scale the genome and run PCA ##########################
# Scale the genome to get rid of missing data
X <- scaleGen(data, NA.method="mean", scale=FALSE, center=TRUE)

# Run the pca
pca <- dudi.pca(X, scale=TRUE, nf=5, scannf=FALSE)

################## Make a map of the position of each snp on the chromosomes  ##################
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
  select(Sample_Name, Shell_color_naive) %>% 
  mutate(Shell_color_naive = ifelse(Shell_color_naive == "Yellow", 0, 1),
         indivs = Sample_Name) %>% 
  rename(genotype = indivs) %>% 
  relocate(genotype) %>% 
  column_to_rownames("Sample_Name")


## Sweden
Pheno_sweden <- metadata %>% 
  filter(Population == "Sweden") %>% 
  select(Sample_Name, Shell_color_naive) %>% 
  mutate(Shell_color_naive = ifelse(Shell_color_naive == "Yellow", 0, 1),
         indivs = Sample_Name) %>% 
  rename(genotype = indivs) %>% 
  relocate(genotype) %>% 
  column_to_rownames("Sample_Name")

################## Scaled relatedness matrix (calculated from vcftools)  ##################
# First, we need to select the reference individuals needed to calculate the relatedness and
# save them in text files
# metadata %>%
#  filter(Population == "Sweden") %>%
#  select(Sample_Name) %>%
#  write.table("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/Swedish_pop.txt",
#              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# 
# metadata %>%
#  filter(Population == "France") %>%
#  select(Sample_Name) %>%
#  write.table("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/French_pop.txt",
#              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
## First, we import the relatedness calculated with vcftools with this command lines
### France
# vcf_file <- "/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz"
# vcftools <- "/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools"
# system2(vcftools,args =c(paste0("--gzvcf ",
#                                vcf_file,
#                                " --relatedness",
#                                " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/French_pop.txt",
#                                " --out /shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/Relatedness/French_relatedness")))
relatedness_france <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/Relatedness/French_relatedness.relatedness", header = TRUE) %>%
  rename(Relatedness = RELATEDNESS_AJK) %>%
  pivot_wider(names_from = INDV2, values_from = Relatedness) %>%
  column_to_rownames("INDV1")

## Sweden
# system2(vcftools,args =c(paste0("--gzvcf ",
#                                vcf_file,
#                                " --relatedness",
#                                " --keep /shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/Data/Reference_indivs_expos/Swedish_pop.txt",
#                                " --out /shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/Relatedness/Swedish_relatedness")))
relatedness_sweden <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/Relatedness/Swedish_relatedness.relatedness", header = TRUE) %>%
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

# Save the GWAS output
# GWAS_france$GWAResult$Pheno_france %>% 
#   select(trait, snp, chr, pos, pValue) %>% 
#   rename(Trait = trait,
#          Position = snp,
#          Chromosome = chr) %>% 
#   mutate(pValue = -log10(pValue)) %>% 
#   rename(logp_val = pValue) %>% 
#   mutate(Population = "France") %>% 
#   rbind(GWAS_sweden$GWAResult$Pheno_sweden %>% 
#           select(trait, snp, chr, pos, pValue) %>% 
#           rename(Trait = trait,
#                  Position = snp,
#                  Chromosome = chr) %>% 
#           mutate(pValue = -log10(pValue)) %>% 
#           rename(logp_val = pValue) %>% 
#           mutate(Population = "Sweden")) %>% 
#   write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/GWAS/GWAS_output.tsv",
#               col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

################## Plot the results  ##################
# Make a summarising table of the cumulative positions of the chromosomes along the genome
cumulative_positions <- data@loc.fac %>% 
  as_tibble %>% 
  rename(Position = value) %>% 
  transform_position_ade2tidy() %>% 
  group_by(Chromosome) %>% 
  summarize(Min_pos = min(Position, na.rm = TRUE),
            Max_pos = max(Position, na.rm = TRUE)) %>% 
  ungroup %>% 
  arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
  mutate(Position_start = cumsum(lag(Max_pos, default = 0))) %>% 
  select(Chromosome, Position_start)

# Know where the inversions are in this cumulative position
inversion_positions_color <- position_inversions %>% 
  group_by(Population) %>% 
  positions_to_polygons(values = c(0, 2)) %>% 
  left_join(cumulative_positions, by = "Chromosome") %>% 
  mutate(Position = Position + Position_start)


# Plot the results
### Size
GWAS_Color <- (GWAS_france$GWAResult$Pheno_france %>%
                filter(trait  == "Shell_color_naive") %>% 
                select(snp, pValue) %>% 
                rename(Position = snp) %>% 
                mutate(logp_val = -log10(pValue),
                       Population = "France") %>% 
                rbind(GWAS_sweden$GWAResult$Pheno_sweden %>%
                        filter(trait == "Shell_color_naive") %>% 
                        select(snp, pValue) %>% 
                        rename(Position = snp) %>% 
                        mutate(logp_val = -log10(pValue),
                               Population = "Sweden")) %>% 
                mutate(Population = Population %>% factor(levels = c("Sweden", "France"))) %>% 
                sample_n(250000) %>% 
                thresholds_manhattan(aes(y = logp_val, facets = Population), values = 5, palette = chromosome_palette)) +
  geom_polygon(data = inversion_positions_color %>% 
                 mutate(Population = Population %>% 
                          factor(levels = c("Sweden", "France"))),
               aes(x = Position, y = Height_polygon, group = Inversion), fill = "black") +
  geom_polygon(data = grouped_inversions %>% 
                 group_by(Population) %>% 
                 positions_to_polygons(values = c(0, 2)) %>% 
                 left_join(cumulative_positions, by = "Chromosome") %>% 
                 mutate(Position = Position + Position_start,
                        Population = Population %>% 
                          factor(levels = c("Sweden", "France"))),
               aes(x = Position, y = Height_polygon, group = Inversion),
               color = "red", fill = NA) +
  facet_col(vars(Population), scale = "free_x") +
  labs(x = "Chromosomes",
       y = expression("-log"[10]* "(p-value)"),
       tag = "(A)") +
  ylim(0, 25) +
  my_theme
################## Clean environment  ##################
rm(chromosome_palette, X, map_chromosomes, geno_france, geno_sweden, Pheno_france, Pheno_sweden, relatedness_france, relatedness_sweden,
   min_rel_france, min_rel_sweden, max_rel_france, max_rel_sweden, gData_france, gData_sweden, cumulative_positions,
   inversion_positions_color, j, kin_france, kin_sweden, grouped_inversions, Delta_freqs_whole_genome, i)
################## Run the local PCA  ##################
clust_groups_together <- (position_inversions %>%
                            filter(Inversion_grouped == "Inv_6.1") %>%
                            mutate(Inversion = Inversion_grouped) %>%
                            run_loca_pca_inversions())$PC_scores

# Plot the results
Local_PCA <- clust_groups_together %>% 
  mutate(Population = ifelse(grepl("LOK", Sample_Name), "Sweden", "France") %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  # In the same way, we keep the inversion we identified earlier that has the size QTL and we 
  # reorganise the levels for the legend
  left_join(metadata, by = c("Sample_Name", "Population"), relationship = "many-to-one") %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France")),
         Shell_color_naive = Shell_color_naive %>% 
           factor(levels = c("Brown", "Yellow"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis2, color = Shell_color_naive), size = 3, alpha = 0.8) +
  scale_color_manual(name = "Color",
                     values = color_palette) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis1, y = Axis2, color = Population), lwd = 1.2) +
  scale_color_manual(name = "Population",
                     values = c("Sweden" = "forestgreen", "France" = "darkorchid4")) +
  labs(x = "Axis 1",
       y = "Axis 2",
       tag = "(C)") +
  my_theme +
  xlim(-150, 160) +
  ylim(-70, 70)

# Clean the environment
rm(clust_groups_together, color_palette)
################## Make a tree of the inversion content  ##################
# Path to write output
output_path <- "/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/08_Trees/"

# Change the table of the delimitation of inversions to be compatible with the functions
delim_invs <- position_inversions %>% 
  mutate(Inversion = Inversion_grouped) %>% 
  select(-Inversion_grouped)

# Variable to be used in the functions to discriminate an inversion from colinear genome
is_inversion <- delim_invs %>%
  filter(Inversion == "Inv_6.1") %>%
  pull(Is_inversion) %>%
  unique

# Subset the vcf file to keep only homokaryotypes for the inversion and the positions
# inside the inversions
list_expos_indivs <- Subset_genetic_data("Inv_6.1", delim_invs, groups_pca, vcf_file, output_path)

# Import the vcf file of the inversion
data_inv <- read.vcfR(paste0(output_path, "VCF_inversion/Inv_6.1.vcf.gz")) %>% 
  vcfR2genind()

# Prepare a haploid genome by randomly sampling one genotype at each locus
# to be able to draw the tree
Prepare_haploid_genome("Inv_6.1", output_path, genetic_data = data_inv, force = TRUE)

# Run Phylogeny
Phylogeny <- Run_and_trace_phylogeny("Inv_6.1", output_path, list_expos_indivs)

# Add a tag to the plot
Tree <- Phylogeny$Tree  +
  annotate("text", label = "Sweden EE", x = -0.15, y = 0.17, size = 4) +
  annotate("text", label = "Sweden SS", x = -0.26, y = 0.015, size = 4) +
  annotate("text", label = "France SS", x = -0.17, y = -0.31, size = 4) +
  annotate("text", label = "France EE", x = 0.04, y = -0.135, size = 4) +
  labs(tag = "(D)") +
  theme(legend.position = "none")

################## Clean environment  ##################
rm(vcf_file, output_path, delim_invs, is_inversion, list_expos_indivs,
   data_inv, Phylogeny)
################## Cline of inversion along transect  ##################
# Candidate SNPs
color_snps <- GWAS_france$GWAResult$Pheno_france %>% 
  rename(Trait = trait,
         Position = snp,
         Chromosome = chr) %>% 
  mutate(log_pval = - log10(pValue),
         Population = "France") %>% 
  filter(grepl("LG6_", Position)) %>% 
  filter(log_pval >= 5) %>% 
  inner_join(position_inversions,
             by = c("Chromosome", "Population"), relationship = "many-to-many") %>% 
  filter(pos >= Start & pos <= End) %>% 
  select(Position) %>% 
  mutate(pos = Position) %>% 
  separate(pos, c("SNP_name", "allele"), "\\.") %>% 
  select(Position, SNP_name)

color_snps_sw <- GWAS_sweden$GWAResult$Pheno_sweden %>% 
  rename(Trait = trait,
         Position = snp,
         Chromosome = chr) %>% 
  mutate(log_pval = - log10(pValue),
         Population = "Sweden") %>% 
  filter(grepl("LG6_", Position)) %>% 
  filter(log_pval >= 5) %>% 
  inner_join(position_inversions,
             by = c("Chromosome", "Population"), relationship = "many-to-many") %>% 
  filter(pos >= Start & pos <= End) %>% 
  select(Position) %>% 
  mutate(pos = Position) %>% 
  separate(pos, c("SNP_name", "allele"), "\\.") %>% 
  select(Position, SNP_name)


# Prepare priors to fit allelic frequency variation models
Priors_clines <- get_allelic_frequencies(genetic_data = data,
                                         SNP_subset = color_snps %>% 
                                           rbind(color_snps_sw),
                                         Extreme_values = pca$li,
                                         var = "Axis2",
                                         meta_data = metadata) %>% 
  mutate(Centre_prior = ifelse(Population == "France", 200, 70),
         Width_prior = ifelse(Population == "France", 50, 40),
         Centre_max = ifelse(Population == "France", 300, 150),
         Centre_min = ifelse(Population == "France", 100, 0),
         Width_max = ifelse(Population == "France", 300, 150),
         Width_min = 2)
# Optimise allelic frequency variation models
cline_fitting_color <- optimise_clines(Priors = Priors_clines %>% 
                                         filter(Position %!in% color_snps_sw$Position),
                                      logarithm = FALSE,
                                      batch_size = 100,
                                      genetic_data = data,
                                      SNP_subset = color_snps,
                                      meta_data = metadata) %>% 
  select_clinal_SNPs() %>% 
  rbind(optimise_clines(
    Priors = Priors_clines %>% 
      filter(Position %in% color_snps_sw$Position),
    logarithm = FALSE,
    genetic_data = data,
    SNP_subset = color_snps_sw,
    meta_data = metadata
  ))

# Select only the SNPs that are clinal in France, but not in Sweden
SNPs_cline_fr_other_sw <- cline_fitting_color %>%
  select(Population, Position, Best_model, Significant) %>%
  pivot_wider(names_from = Population,
              values_from = Significant) %>%
  filter(France == "Clinal" & Sweden != "Clinal") %>%
  pull(Position) %>% 
  c(color_snps_sw$Position)

# Get the plotting values of these models along the transect
plotting_clines_color <- plot_clines(cline_fitting_color %>% 
                                       mutate(Best_model = Significant),
                                    genetic_data = data,
                                    real_distance = TRUE,
                                    SNP_subset = color_snps %>% 
                                      rbind(color_snps_sw),
                                    meta_data = metadata)

# Find the frequencies of the inversion in the sheltered part in sweden.
# To do so, we use the frequency of individuals presenting one version of the rearrangement in the transect ends.
# Therefore, we first must take the individuals that are on the left or on the right side of the transect.
# This is what is done here. For each of the 4 repetitions, we subset a group of individuals in the left 
# (= sheltered) and the right (= exposed) part of the transect
indivs_shelt_freq_sweden_6.1 <- groups_pca %>% 
  mutate(Group = Group - 1) %>% 
  select(-Inversion) %>% 
  inner_join(position_inversions, by = c("Chromosome", "Population"),
             relationship = "many-to-many") %>% 
  filter(Inversion == "Inv_6.1") %>% 
  left_join(metadata, by = c("Sample_Name", "Population")) %>% 
  filter(Population == "Sweden") %>% 
  inner_join(ref_individuals, by = c("Sample_Name", "Population")) %>% 
  filter(Exposition == "Sheltered") %>% 
  select(Sample_Name, Group)

# Once we did it for one end of one transect, we repeat it on the other ends of the two transects
indivs_exp_freq_sweden_6.1 <- groups_pca %>% 
  mutate(Group = Group - 1) %>% 
  select(-Inversion) %>% 
  inner_join(position_inversions, by = c("Chromosome", "Population"),
             relationship = "many-to-many") %>% 
  filter(Inversion == "Inv_6.1") %>% 
  left_join(metadata, by = c("Sample_Name", "Population")) %>% 
  filter(Population == "Sweden") %>% 
  inner_join(ref_individuals, by = c("Sample_Name", "Population")) %>% 
  filter(Exposition == "Exposed") %>% 
  select(Sample_Name, Group)
# Find the frequencies of the inversion in the sheltered part in france
indivs_shelt_freq_france_6.1 <- groups_pca %>% 
  mutate(Group = Group - 1) %>% 
  select(-Inversion) %>% 
  inner_join(position_inversions, by = c("Chromosome", "Population"),
             relationship = "many-to-many") %>% 
  filter(Inversion == "Inv_6.1") %>% 
  left_join(metadata, by = c("Sample_Name", "Population")) %>% 
  filter(Population == "France") %>% 
  inner_join(ref_individuals, by = c("Sample_Name", "Population")) %>% 
  filter(Exposition == "Sheltered") %>% 
  select(Sample_Name, Group)

indivs_exp_freq_france_6.1 <- groups_pca %>% 
  mutate(Group = Group - 1) %>% 
  select(-Inversion) %>% 
  inner_join(position_inversions, by = c("Chromosome", "Population"),
             relationship = "many-to-many") %>% 
  filter(Inversion == "Inv_6.1") %>% 
  left_join(metadata, by = c("Sample_Name", "Population")) %>% 
  filter(Population == "France") %>% 
  inner_join(ref_individuals, by = c("Sample_Name", "Population")) %>% 
  filter(Exposition == "Exposed") %>% 
  select(Sample_Name, Group)

# Here, we make a function to calculate the allelic frequencies in one end of the transect
calculate_freq_df <- function(df){
  sum_inv <- df$Group %>% sum
  return(sum_inv/(2*nrow(df)))
}

# Then, we use the created function to calculate the allelic frequencies on the transect
shelt_freq_sweden_6.1 <- indivs_shelt_freq_sweden_6.1 %>% 
  calculate_freq_df()

exp_freq_sweden_6.1 <- indivs_exp_freq_sweden_6.1 %>% 
  calculate_freq_df()

shelt_freq_france_6.1 <- indivs_shelt_freq_france_6.1 %>% 
  calculate_freq_df()

exp_freq_france_6.1 <- indivs_exp_freq_france_6.1 %>% 
  calculate_freq_df()


# We use these frequencies as priors to fit the allelic frequency variations modles
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
         p_right_expos = max(p_right_expos, 0.000001) %>% logit
  )

# Fit the models
cline_inv_6.1 <- priors_clines_inv_6.1 %>% 
  rename(Pop = Population) %>% 
  group_by(Pop) %>% 
  bow(tie(Centre, Width, Left, Right) := mle2(clineflog,
                                              list(centre = Centre_prior,
                                                   width = Width_prior ,
                                                   left = p_left_shelt,
                                                   right = p_right_expos),
                                              data = list(x = groups_pca %>% 
                                                            filter(Inversion == "Inv_14.1") %>% 
                                                            left_join(metadata, by = c("Sample_Name", "Population")) %>% 
                                                            filter(Population == Pop) %>% 
                                                            select(LCmeanDist) %>% 
                                                            drop_na %>% 
                                                            as.vector %>%
                                                            unlist %>%
                                                            unname,
                                                          g = groups_pca %>% 
                                                            filter(Inversion == "Inv_14.1") %>% 
                                                            left_join(metadata, by = c("Sample_Name", "Population")) %>% 
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


# Plot the results
Clines_along_transect <- plotting_clines_color %>%
  filter(Position %in% SNPs_cline_fr_other_sw) %>%
  left_join(cline_fitting_color, by = c("Position", "Population")) %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France")),
         Kind_SNP = case_when(
           Position %in% color_snps$Position ~ "QTL SNPs in France",
           TRUE ~ "QTL SNP in Sweden"
         ) %>% 
           factor(levels = c("QTL SNPs in France", "QTL SNP in Sweden", "Inv_6.1"))) %>% 
  select(Population, Position, LCmeanDist, Frequency, Kind_SNP) %>% 
  rbind(plot_inv_6.1 %>% 
          mutate(Population = Population %>% 
                   factor(levels = c("Sweden", "France")),
                 Position = "Inv_6.1",
                 Kind_SNP = "Inv_6.1" %>% 
                   factor(levels = c("QTL SNPs in France", "QTL SNP in Sweden", "Inv_6.1")),
                 Frequency = 1 - phen_cline) %>% 
          rename(LCmeanDist = position) %>% 
          select(-phen_cline)) %>% 
  ggplot() +
  geom_line(aes(x = LCmeanDist, y = Frequency, group = Position, color = Kind_SNP, lty = Kind_SNP, lwd = Kind_SNP), alpha = 0.5) +
  scale_color_manual(name = "Curve",
                     values = c("QTL SNPs in France" = "black", "QTL SNP in Sweden" = "orange", "Inv_6.1" = "orange")) +
  scale_linetype_manual(name = "Curve",
                        values = c("QTL SNPs in France" = 2, "QTL SNP in Sweden" = 2, "Inv_6.1" = 1)) +
  scale_linewidth_manual(name = "Curve",
                         values = c("QTL SNPs in France" = 1.2, "QTL SNP in Sweden" = 1.2, "Inv_6.1" = 1.8)) +
  facet_row(vars(Population), scales = "free") +
  my_theme +
  theme(legend.position = "bottom") +
  labs(x = "Position along the transect (m)",
       y = "Frequency",
       tag = "(E)") +
  ylim(0, 1)

################## Remove everything that is not needed ##################
rm(position_inversions, data, metadata, groups_pca, pca, color_snps, Priors_clines,
   cline_fitting_color, plotting_clines_color, indivs_shelt_freq_france_6.1,
   indivs_shelt_freq_sweden_6.1, indivs_exp_freq_france_6.1,
   indivs_exp_freq_sweden_6.1, calculate_freq_df, shelt_freq_sweden_6.1,
   shelt_freq_france_6.1, exp_freq_france_6.1, exp_freq_sweden_6.1,
   priors_clines_inv_6.1, cline_inv_6.1, plot_inv_6.1, ref_individuals,
   my_theme, GWAS_france)
################## Group everything together ##################
layout_pca_tree <- "
AAAAAABBBBBB
AAAAAABBBBBB
"

Pca_tree <- Local_PCA + Tree +
  plot_layout(design = layout_pca_tree)

(GWAS_Color / Pca_tree / Clines_along_transect) %>% 
  ggsave(plot = ., filename = "/shared/home/bpajot/Report_presentations/Figures/Anglais/Figure_6.svg", device = "svg", units = "px",
         height = 1800, width = 1500, scale = 3.5, dpi = 450)
