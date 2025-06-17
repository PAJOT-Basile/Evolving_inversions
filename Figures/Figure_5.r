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
size_palette <-  c("#4e79a7", "grey75", "#f28e2b")

# Theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

################## Useful functions  ##################
source("../General_scripts/Functions_optimise_plot_clines.r")

################## Import data  ##################
# Positions of the inversions
position_inversions <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
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
vcf_file <- "../../Output/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz"
data <- read.vcfR(vcf_file) %>% 
  vcfR2genind()

# We add the poputations to the vcf object
data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor
# And the exposition of the individuals
data@other$exposition <- (data@tab %>% rownames %>% str_split_fixed(., "_", 6))[, 5] %>% as.factor
data@other$exposition[which(data@other$exposition == "TRANSI")] <- "TRANS" %>% as.factor
data@other$exposition <- data@other$exposition %>% droplevels()

# Import the metadata
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

# Import the delta freqs for the whole genome
Delta_freqs_whole_genome <- read.table("../../Output/Sweden_France_parallelism/Delta_freqs_whole_genome.tsv",
                                       sep = "\t", header = TRUE) %>% 
  select_good_SNPs(Delta_freq_Sweden)

# Import the local PCA per inversion
groups_pca <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
                         sep = "\t", header = TRUE) %>% 
  # Filter to get only the inversions that are kept after the 4 step filtering of inversions
  inner_join(position_inversions,
             by = c("Chromosome", "Inversion", "Population"),
             relationship = "many-to-many")

# Import the reference individuals
ref_individuals <- read.table("../../Output/Data/Reference_indivs/France_exposed.txt",
                              sep = "\t", header = FALSE) %>%
  mutate(Population = "France",
         Exposition = "Exposed") %>%
  rbind(read.table("../../Output/Data/Reference_indivs/France_sheltered.txt",
                   sep = "\t", header = FALSE) %>%
          mutate(Population = "France",
                 Exposition = "Sheltered")) %>%
  rbind(read.table("../../Output/Data/Reference_indivs/Sweden_exposed.txt",
                   sep = "\t", header = FALSE) %>%
          mutate(Population = "Sweden",
                 Exposition = "Exposed")) %>%
  rbind(read.table("../../Output/Data/Reference_indivs/Sweden_sheltered.txt",
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
  select(Sample_Name, Length) %>% 
  mutate(indivs = Sample_Name) %>% 
  rename(genotype = indivs) %>% 
  relocate(genotype) %>% 
  column_to_rownames("Sample_Name")


## Sweden
Pheno_sweden <- metadata %>% 
  filter(Population == "Sweden") %>% 
  select(Sample_Name, Length) %>% 
  mutate(indivs = Sample_Name) %>% 
  rename(genotype = indivs) %>% 
  relocate(genotype) %>% 
  column_to_rownames("Sample_Name")

################## Scaled relatedness matrix (calculated from vcftools)  ##################
relatedness_france <- read.table("../../Output/Sweden_France_parallelism/10_GWAS/Relatedness/French_scaled_relatedness.tsv", header = TRUE)
relatedness_sweden <- read.table("../../Output/Sweden_France_parallelism/10_GWAS/Relatedness/Swedish_scaled_relatedness.tsv", header = TRUE)
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
inversion_positions_size <- position_inversions %>% 
  group_by(Population) %>% 
  positions_to_polygons() %>% 
  left_join(cumulative_positions, by = "Chromosome") %>% 
  mutate(Position = Position + Position_start)


# Plot the results
### Size
GWAS_Size <- (GWAS_france$GWAResult$Pheno_france %>% 
   filter(trait  == "Length") %>% 
   select(snp, pValue) %>% 
   rename(Position = snp) %>% 
   mutate(logp_val = -log10(pValue),
          Population = "France") %>% 
   rbind(GWAS_sweden$GWAResult$Pheno_sweden %>%
           filter(trait == "Length") %>% 
           select(snp, pValue) %>% 
           rename(Position = snp) %>% 
           mutate(logp_val = -log10(pValue),
                  Population = "Sweden")) %>% 
   mutate(Population = Population %>% factor(levels = c("Sweden", "France"))) %>% 
   thresholds_manhattan(aes(y = logp_val, facets = Population), values = 5, palette = chromosome_palette)) +
  geom_polygon(data = inversion_positions_size %>% 
                 mutate(Population = Population %>% 
                          factor(levels = c("Sweden", "France"))),
               aes(x = Position, y = Height_polygon, group = Inversion), fill = "black") +
  geom_polygon(data = grouped_inversions %>% 
                 group_by(Population) %>% 
                 positions_to_polygons() %>% 
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
  ylim(0, 8) +
  my_theme

################## Clean environment  ##################
rm(chromosome_palette, X, map_chromosomes, geno_france, geno_sweden, Pheno_france, Pheno_sweden, relatedness_france, relatedness_sweden,
   min_rel_france, min_rel_sweden, max_rel_france, max_rel_sweden, gData_france, gData_sweden, GWAS_france, cumulative_positions,
   inversion_positions_size, j, kin_france, kin_sweden, i)
################## Run the local PCA  ##################
clust_groups_together <- (position_inversions %>%
  filter(Inversion_grouped == "Inv_14.1") %>%
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
           factor(levels = c("Sweden", "France"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis2, color = Length), size = 3, alpha = 0.8) +
  scale_colour_gradientn(name = "Shell size\n(mm)",
                         colors=size_palette, 
                         limits=c(min(metadata$Length, na.rm=TRUE) - 1, 
                                  max(metadata$Length, na.rm=TRUE) + 1),
                         breaks=c(5, 10, 14)) +
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
rm(clust_groups_together)
################## Make a tree of the inversion content  ##################
# Path to write output
output_path <- "../../Output/Sweden_France_parallelism/08_Trees/"

# Change the table of the delimitation of inversions to be compatible with the functions
delim_invs <- position_inversions %>% 
  mutate(Inversion = Inversion_grouped) %>% 
  select(-Inversion_grouped)

# Variable to be used in the functions to discriminate an inversion from colinear genome
is_inversion <- delim_invs %>%
  filter(Inversion == "Inv_14.1") %>%
  pull(Is_inversion) %>%
  unique

# Subset the vcf file to keep only homokaryotypes for the inversion and the positions
# inside the inversions
list_expos_indivs <- Subset_genetic_data("Inv_14.1", delim_invs, groups_pca, vcf_file, output_path)

# Import the vcf file of the inversion
data_inv <- read.vcfR(paste0(output_path, "VCF_inversion/Inv_14.1.vcf.gz")) %>% 
  vcfR2genind()

# Prepare a haploid genome by randomly sampling one genotype at each locus
# to be able to draw the tree
Prepare_haploid_genome("Inv_14.1", output_path, genetic_data = data_inv)

# Run Phylogeny
Phylogeny <- Run_and_trace_phylogeny("Inv_14.1", output_path, list_expos_indivs)

# Add a tag to the plot
Tree <- Phylogeny$Tree +
  labs(tag = "(D)") +
  theme(legend.position = "none")

################## Clean environment  ##################
rm(size_palette, vcf_file, output_path, delim_invs, is_inversion, list_expos_indivs,
   data_inv, Phylogeny)
################## Cline of inversion along transect  ##################
# In the size QTL, as we do not have many individuals and that size is a very polygenic
# trait, we choose to lower the threshold value usually used to identify QTLs (usually 5)
# to 4.
size_snps <- GWAS_sweden$GWAResult$Pheno_sweden %>% 
  rename(Trait = trait,
         Position = snp,
         Chromosome = chr) %>% 
  mutate(log_pval = - log10(pValue),
         Population = "Sweden") %>% 
  filter(Trait == "Length",
         grepl("LG14_", Position)) %>% 
  filter(log_pval >= 4) %>% 
  inner_join(grouped_inversions,
             by = c("Chromosome", "Population"), relationship = "many-to-many") %>% 
  filter(pos >= Start & pos <= End) %>% 
  select(Position) %>% 
  rbind(Delta_freqs_whole_genome %>% 
          select(Position, F4_stat) %>% 
          arrange(F4_stat) %>% 
          dplyr::slice(1:6) %>% 
          select(Position)) %>% 
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
                                      logarithm = TRUE,
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

# Find the frequencies of the inversion in the sheltered part in sweden.
# To do so, we use the frequency of individuals presenting one version of the rearrangement in the transect ends.
# Therefore, we first must take the individuals that are on the left or on the right side of the transect.
# This is what is done here. For each of the 4 repetitions, we subset a group of individuals in the left 
# (= sheltered) and the right (= exposed) part of the transect
indivs_shelt_freq_sweden_14.1 <- groups_pca %>% 
  mutate(Group = Group - 1) %>% 
  select(-Inversion) %>% 
  inner_join(grouped_inversions, by = c("Chromosome", "Population"),
             relationship = "many-to-many") %>% 
  filter(Inversion == "Inv_14.1") %>% 
  left_join(metadata, by = c("Sample_Name", "Population")) %>% 
  filter(Population == "Sweden") %>% 
  inner_join(ref_individuals, by = c("Sample_Name", "Population")) %>% 
  filter(Exposition == "Sheltered") %>% 
  select(Sample_Name, Group)

# Once we did it for one end of one transect, we repeat it on the other ends of the two transects
indivs_exp_freq_sweden_14.1 <- groups_pca %>% 
  mutate(Group = Group - 1) %>% 
  select(-Inversion) %>% 
  inner_join(grouped_inversions, by = c("Chromosome", "Population"),
             relationship = "many-to-many") %>% 
  filter(Inversion == "Inv_14.1") %>% 
  left_join(metadata, by = c("Sample_Name", "Population")) %>% 
  filter(Population == "Sweden") %>% 
  inner_join(ref_individuals, by = c("Sample_Name", "Population")) %>% 
  filter(Exposition == "Exposed") %>% 
  select(Sample_Name, Group)
# Find the frequencies of the inversion in the sheltered part in france
indivs_shelt_freq_france_14.1 <- groups_pca %>% 
  mutate(Group = Group - 1) %>% 
  select(-Inversion) %>% 
  inner_join(grouped_inversions, by = c("Chromosome", "Population"),
             relationship = "many-to-many") %>% 
  filter(Inversion == "Inv_14.1") %>% 
  left_join(metadata, by = c("Sample_Name", "Population")) %>% 
  filter(Population == "France") %>% 
  inner_join(ref_individuals, by = c("Sample_Name", "Population")) %>% 
  filter(Exposition == "Sheltered") %>% 
  select(Sample_Name, Group)

indivs_exp_freq_france_14.1 <- groups_pca %>% 
  mutate(Group = Group - 1) %>% 
  select(-Inversion) %>% 
  inner_join(grouped_inversions, by = c("Chromosome", "Population"),
             relationship = "many-to-many") %>% 
  filter(Inversion == "Inv_14.1") %>% 
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
         p_right_expos = max(p_right_expos, 0.000001) %>% logit
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
Snps_neg_f4 <- Delta_freqs_whole_genome %>% 
  select(Position, F4_stat) %>% 
  arrange(F4_stat) %>% 
  mutate(pos = Position) %>% 
  filter(grepl("LG14_", Position)) %>% 
  transform_position_ade2tidy() %>% 
  mutate(pos_lag = lag(Position, default = 0)) %>% 
  filter(Position != pos_lag) %>% 
  select(-pos_lag) %>% 
  inner_join(position_inversions %>%
               filter(Inversion_grouped == "Inv_14.1") %>% 
               select(-c(Inversion, Length)) %>% 
               mutate(toto = rep(c(1, 2), 2)) %>% 
               pivot_wider(names_from = "toto",
                           values_from = c("Start", "End")) %>% 
               rename(Inversion = Inversion_grouped) %>% 
               group_by(Chromosome, Inversion) %>%
               summarize(Start_1 = min(Start_1),
                         Start_2 = min(Start_2),
                         End_1 = max(End_1),
                         End_2 = max(End_2)),
             by = "Chromosome",
             relationship = "many-to-many") %>% 
  filter((Position <= End_1 & Position >= Start_1) | (Position <= End_2 & Position >= Start_2)) %>% 
  dplyr::slice(1:6) %>% 
  pull(pos)

# Filter to get the names of the positions that have a negative f4
Snps_to_keep_size <- plotting_clines_size %>% 
  left_join(cline_fitting_size, by = c("Position", "Population")) %>% 
  filter(Best_model == "Clinal",
         Position %!in% Snps_neg_f4) %>% 
  ungroup %>% 
  pull(Position) %>% 
  unique


# Plot the results
Clines_along_transect <- plotting_clines_size %>% 
  left_join(cline_fitting_size, by = c("Position", "Population")) %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France")),
         F4_sign = case_when(
           Position %in% Snps_neg_f4 ~ "Negative",
           Position %in% Snps_to_keep_size ~ "Positive",
           TRUE ~ NA
         )) %>% 
  ggplot() +
  geom_line(aes(x = LCmeanDist, y = Frequency, group = Position, colour = F4_sign), lty = "dashed", lwd = 1.2) +
  scale_colour_manual(name = "F4 sign",
                      values = c("Positive" = "black", "Negative" = "darkred"),
                      guide = "none") +
  geom_line(data = plot_inv_14.1 %>% 
              mutate(Population = Population %>% 
                       factor(levels = c("Sweden", "France"))),
            aes(x = position, y = 1 - phen_cline), lwd = 1.8, color = "aquamarine3") +
  facet_row(vars(Population), scales = "free") +
  my_theme +
  labs(x = "Position along the transect (m)",
       y = "Frequency",
       tag = "(E)") +
  ylim(0, 1)

################## Remove everything that is not needed ##################
rm(position_inversions, grouped_inversions, data, metadata, Delta_freqs_whole_genome,
   groups_pca, pca, GWAS_sweden, size_snps, Priors_clines, cline_fitting_size,
   plotting_clines_size, indivs_shelt_freq_france_14.1, indivs_shelt_freq_sweden_14.1,
   indivs_exp_freq_france_14.1, indivs_exp_freq_sweden_14.1, calculate_freq_df,
   shelt_freq_sweden_14.1, shelt_freq_france_14.1, exp_freq_france_14.1, exp_freq_sweden_14.1,
   priors_clines_inv_14.1, cline_inv_14.1, plot_inv_14.1, Snps_neg_f4, Snps_to_keep_size,
   ref_individuals)
################## Group everything together ##################
layout_pca_tree <- "
AAAAAABBBBBB
AAAAAABBBBBB
"

Pca_tree <- Local_PCA + Tree +
  plot_layout(design = layout_pca_tree)

(GWAS_Size / Pca_tree / Clines_along_transect) %>% 
  ggsave(plot = ., filename = "../../Output/Figures/Figure_5.png", device = "png", units = "px",
        height = 1800, width = 1500, scale = 3.5)
