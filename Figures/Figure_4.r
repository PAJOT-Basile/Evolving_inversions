# Import libraries
libraries <- c("tidyverse", "adegenet", "vcfR", "readxl", "ggforce", "ggh4x", "ggpubr", "ggtext")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(characters = libraries, character.only = TRUE)
rm(libraries)
if (!require("patchwork")) devtools::install_github("thomasp85/patchwork")


################## Useful functions  ##################
source("../General_scripts/Functions_optimise_plot_clines.r")

################################ Useful variables ################################
# Basic theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 30))

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
# Run the PCA
pca <- scaleGen(data, NA.method="mean", scale=FALSE, center=TRUE) %>% 
  dudi.pca(scale=TRUE, nf=5, scannf=FALSE)

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

################## Making the delta frequencies correlation plot  ##################
# Calculate the differences in allelic frequencies
Delta_freqs_whole_genome <- get_delta_freqs_and_F4(
  genetic_data = data,
  Extreme_values = pca$li,
  var = "Axis2",
  meta_data = metadata
) %>% 
  select(-c(starts_with("p_"))) %>% 
  mutate(across(contains("Delta_"), .fns = abs, .names = "{.col}_abs")) %>% 
  # Calculate the mean delta freq
  mutate(Mean_delta_freq = rowMeans(select(., ends_with("_abs")), na.rm = TRUE))

# Save the output
Delta_freqs_whole_genome %>% 
  write.table("../../Output/Sweden_France_parallelism/Delta_freqs_whole_genome.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

Delta_freqs_whole_genome <- Delta_freqs_whole_genome %>% 
  left_join(pca$co %>% 
              rownames_to_column("Position"),
            by = "Position")

################## Plot the differences in allelic frequencies along genome  ##################
# Get a table of the maximum number of positions along each chromosome
Cumulative_position_per_chromosome <- data@loc.fac %>% 
  as_tibble %>% 
  rename(Position = value) %>% 
  transform_position_ade2tidy() %>% 
  group_by(Chromosome) %>% 
  summarize(Pos_max_chromosome = max(Position)) %>% 
  arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
  mutate(Pos_max_chromosome_cumul = lag(cumsum(Pos_max_chromosome), default = 0))

delta_freqs <- Delta_freqs_whole_genome %>% 
  select(Position, contains("Delta")) %>% 
  pivot_longer(cols = contains("Delta"), names_to = "Population", values_to = "Delta_freqs") %>% 
  mutate(Population = ifelse(Population == "Delta_freq_France", "France", "Sweden") %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  thresholds_manhattan(aes(y = Delta_freqs, facets = Population), values = 0.3, size = 3, palette = c("grey71", "turquoise4")) +
  geom_polygon(data = position_inversions %>% 
                 group_by(Population) %>% 
                 positions_to_polygons(values = c(0, 0.05)) %>% 
                 ungroup %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul,
                        Population = Population %>% 
                          factor(levels = c("Sweden", "France"))),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion), fill = "black") +
  geom_polygon(data = grouped_inversions %>% 
                 group_by(Population) %>% 
                 positions_to_polygons(values = c(0, 0.05)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul) %>%
                 arrange(Population, Inversion, Cumulative_position),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion),
               color = "red", fill = NA) +
  my_theme +
  facet_col(vars(Population), scales = "free") +
  labs(x = "Chromosomes",
       y = expression(Delta * "freq"),
       tag = "(B)")

################## Plot the F4  ##################
F4_plot <- Delta_freqs_whole_genome %>% 
  select(Position, F4_stat) %>% 
  thresholds_manhattan(aes(y = F4_stat), values = c(-0.09, 0.09),
                       size = 3, absolute = FALSE,
                       palette = c("grey71", "turquoise4")) +
  geom_polygon(data = position_inversions %>% 
                 filter(Population == "France") %>% 
                 positions_to_polygons(values = c(-0.02, 0)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion),
               fill = "black") +
  geom_polygon(data = position_inversions %>% 
                 filter(Population == "Sweden") %>% 
                 positions_to_polygons(values = c(0, 0.02)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion),
               fill = "black") +
  geom_polygon(data = grouped_inversions %>% 
                 filter(Population == "Sweden") %>%
                 positions_to_polygons(values = c(0, 0.02)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion),
               color = "red", fill = NA) +
  geom_polygon(data = grouped_inversions %>% 
                 filter(Population == "France") %>%
                 positions_to_polygons(values = c(-0.02, 0)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion),
               color = "red", fill = NA) +
  labs(x = "Chromosomes",
       y = expression(~italic('f')~'4'),
       tag = "(C)") +
  my_theme +
  geom_hline(yintercept = 0, color = "black", lwd = 1.2)

################## Run the local PCA  ##################
bin_size <- 500000
SNP_positions <- data@loc.fac %>% 
  as_tibble %>% 
  rename(pos=value) %>% 
  separate(pos, c("LG", "SUPER", "SUPER_frag", "Position"), "_") %>% 
  unite("Chromosome", LG, SUPER, SUPER_frag, remove=TRUE) %>% 
  mutate(bin = (as.numeric(Position)/bin_size %>% trunc) %>% floor,
         chrom = Chromosome) %>%
  unite("Position", chrom, Position) %>% 
  mutate(bin = as.factor(bin))


pca_result_per_bin <- data.frame(
  Chromosome = c(),
  Bin = c(),
  Population = c(),
  Axis1 = c(),
  Axis2 = c(),
  Axis3 = c(),
  Axis4 = c(),
  Axis5 = c(),
  Sample_Name = c()
)
for (chromosome in (SNP_positions$Chromosome %>% unique)){
  print(chromosome)
  # First, we get all the bins per chromosome and iterate over them
  bins <- SNP_positions$bin[SNP_positions$Chromosome == chromosome]
  for (bin in (bins %>% unique)){
    # For each bin, we make a list of all the positions that are in this bin
    list_of_positions_to_sample_per_bin <- SNP_positions$Position[which((SNP_positions$Chromosome == chromosome) & (SNP_positions$bin == bin))]
    if (list_of_positions_to_sample_per_bin %>% length > 10){
      # Then, we extract the genind object of all SNPs contained in said bin
      data_per_bin <- data[loc=list_of_positions_to_sample_per_bin]
      # We separate the populations
      for (population in data@pop %>% levels){
        # First, we get all the individuals from each population
        data_per_bin_per_pop <- data_per_bin[which(data@pop == population)]
        
        # We scale it to get rid of missing data
        scaled_data_bin_per_pop <- scaleGen(data_per_bin_per_pop, NA.method="mean", scale=FALSE, center=TRUE)
        
        # Run the PCA
        pca_per_bin_per_pop <- dudi.pca(scaled_data_bin_per_pop, scale=TRUE, nf=5, scannf=FALSE)
        
        # Changing the sign of the calculations depending on the population and the transect
        pca_points <- pca_per_bin_per_pop$li %>% 
          rownames_to_column("indiv") %>% 
          mutate(Sample_Name = indiv) %>%
          separate(indiv, c("Species", "Cov", "Population", "Sex", "Exposition", "Size", "Sample_ID"), sep="_") %>% 
          mutate(Exposition = ifelse(Exposition == "TRANSI", "TRANS", Exposition))
        # Axis 1
        if ((pca_points$Axis1[which(pca_points$Exposition == "SHELT")] %>% mean) < 0){
          pca_points <- pca_points %>% 
            mutate(Axis1 = - Axis1)
        }
        # Axis 2
        if ((pca_points$Axis2[which(pca_points$Exposition == "SHELT")] %>% mean) < 0){
          pca_points <- pca_points %>% 
            mutate(Axis2 = - Axis2)
        }
        # Axis 3
        if ((pca_points$Axis3[which(pca_points$Exposition == "SHELT")] %>% mean) < 0){
          pca_points <- pca_points %>% 
            mutate(Axis3 = - Axis3)
        }
        # Axis 4
        if ((pca_points$Axis4[which(pca_points$Exposition == "SHELT")] %>% mean) < 0){
          pca_points <- pca_points %>% 
            mutate(Axis4 = - Axis4)
        }
        # Axis 5
        if ((pca_points$Axis5[which(pca_points$Exposition == "SHELT")] %>% mean) < 0){
          pca_points <- pca_points %>% 
            mutate(Axis5 = - Axis5)
        }
        
        # Add the pca data to the dataframe
        pca_result_per_bin <- rbind(pca_result_per_bin,
                                    cbind(
                                      "Chromosome" = rep(chromosome, pca_per_bin_per_pop$li %>% nrow),
                                      "Bin" = rep(bin, pca_per_bin_per_pop$li %>% nrow),
                                      pca_points %>% select(-c(Species, Cov, Sex, Exposition, Size, Sample_ID))
                                    ))
      }
    }
  }
}

pca_result_per_bin_modif <- pca_result_per_bin %>% 
  filter(abs(Axis1) <= 50) %>% 
  mutate(Chromosome = Chromosome %>% factor(levels = c("LG1_SUPER_1", "LG2_SUPER_2",
                                                       "LG3_SUPER_4", "LG4_SUPER_7",
                                                       "LG5_SUPER_9", "LG6_SUPER_11",
                                                       "LG7_SUPER_10", "LG8_SUPER_5",
                                                       "LG9_SUPER_8", "LG10_SUPER_14",
                                                       "LG11_SUPER_16", "LG12_SUPER_3",
                                                       "LG13_SUPER_17", "LG14_SUPER_15",
                                                       "LG15_SUPER_13", "LG16_SUPER_12",
                                                       "LG17_SUPER_6")),
         LG = ((Chromosome %>% str_split_fixed(., "_", 2))[, 1]) %>% 
           factor(levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7",
                             "LG8", "LG9", "LG10", "LG11", "LG12", "LG13",
                             "LG14", "LG15", "LG16", "LG17")),
         Population = ifelse(Population == "LOKn", "Sweden", "France") %>% 
           factor(levels = c("Sweden", "France")),
         Bin = (Bin %>% as.numeric) * bin_size) 
################## Plot results of local PCA  ##################
pca_data_cum <- pca_result_per_bin_modif %>% 
  group_by(Chromosome) %>% 
  summarise(bin_max = Bin %>% max) %>% 
  mutate(bin_add = lag(cumsum(bin_max), default=0)) %>% 
  select(Chromosome, bin_add)

pca_data_manhat <- pca_result_per_bin_modif %>% 
  inner_join(pca_data_cum, by="Chromosome") %>% 
  mutate(Bin = Bin %>% as.numeric,
         bin_cum = Bin + bin_add)


center_bin_chrom <- pca_data_manhat %>% 
  group_by(Chromosome) %>%
  summarise(center = bin_cum %>% mean) %>% 
  separate(Chromosome, c("Chromosome", "SUPER", "SUPER_frag"), "_")


color_chromosome_local_pca <- rep(c("grey71","turquoise4"), (pca_data_manhat$Chromosome %>% unique %>% length)/2) %>% 
  cbind(c("LG1_SUPER_1", "LG2_SUPER_2",
          "LG3_SUPER_4", "LG4_SUPER_7",
          "LG5_SUPER_9", "LG6_SUPER_11",
          "LG7_SUPER_10", "LG8_SUPER_5",
          "LG9_SUPER_8", "LG10_SUPER_14",
          "LG11_SUPER_16", "LG12_SUPER_3",
          "LG13_SUPER_17", "LG14_SUPER_15",
          "LG15_SUPER_13", "LG16_SUPER_12",
          "LG17_SUPER_6")) %>% 
  as_tibble %>% 
  rename(Color = ".",
         Chromosome = V2)

chromosome_delimitations_polygons <- pca_data_manhat %>% 
  group_by(Chromosome) %>% 
  mutate(bin_max = bin_cum %>% max) %>% 
  ungroup %>% 
  pivot_longer(cols = c(bin_add, bin_max),
               values_to = "bin_delim_x") %>% 
  select(-name) %>% 
  group_by(Population) %>% 
  mutate(max_bin_y = (Axis1 * 5) %>% max,
         min_bin_y = (Axis1 * 5) %>% min) %>% 
  pivot_longer(cols = c(max_bin_y, min_bin_y),
               values_to = "bin_delim_y") %>% 
  select(-name) %>% 
  ungroup() %>% 
  select(Chromosome, bin_delim_x, bin_delim_y, Population) %>% 
  unique() %>% 
  mutate(order_poly = rep(c(1, 2, 4, 3), length.out = nrow(.))) %>% 
  arrange(Chromosome, order_poly) %>%
  select(-order_poly) %>% 
  inner_join(color_chromosome_local_pca, by="Chromosome") %>% 
  mutate(Chromosome = Chromosome %>% factor(levels = c("LG1_SUPER_1", "LG2_SUPER_2",
                                                       "LG3_SUPER_4", "LG4_SUPER_7",
                                                       "LG5_SUPER_9", "LG6_SUPER_11",
                                                       "LG7_SUPER_10", "LG8_SUPER_5",
                                                       "LG9_SUPER_8", "LG10_SUPER_14",
                                                       "LG11_SUPER_16", "LG12_SUPER_3",
                                                       "LG13_SUPER_17", "LG14_SUPER_15",
                                                       "LG15_SUPER_13", "LG16_SUPER_12",
                                                       "LG17_SUPER_6")))

# Plotting the manhattan plot
local_pca <- pca_data_manhat %>% 
  select(-Chromosome) %>% 
  rename(Chromosome = LG) %>% 
  left_join(center_bin_chrom, by="Chromosome") %>%
  left_join(metadata, by=c("Sample_Name", "Population")) %>%
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France")),
         Habitat = Habitat %>% 
           factor(levels = c("Exposed", "Transition", "Sheltered"))) %>% 
  ggplot() +
  geom_polygon(data = chromosome_delimitations_polygons %>% 
                 mutate(Population = Population %>% 
                          factor(levels = c("Sweden", "France"))),
               aes(x=bin_delim_x, y=bin_delim_y, group=Chromosome, fill = Color), alpha=0.3) +
  scale_fill_manual(values = chromosome_delimitations_polygons$Color %>% unique,
                    labels = chromosome_delimitations_polygons$Color %>% unique) +
  new_scale_color() +
  geom_line(aes(x = bin_cum, y = Axis1 * 5, color = Habitat, group = Sample_Name), lwd = 0.1, alpha = 0.5) +
  scale_color_manual(name = "Habitat",
                     values = c("Exposed" = "orange2", "Transition" = "deeppink", "Sheltered" = "dodgerblue3")) +
  geom_polygon(data = position_inversions %>% 
                 group_by(Population) %>% 
                 positions_to_polygons(values = c(-200, -150)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion),
               fill = "black") +
  geom_polygon(data = grouped_inversions %>% 
                 group_by(Population) %>% 
                 positions_to_polygons(values = c(-200, -150)) %>% 
                 left_join(Cumulative_position_per_chromosome, by = "Chromosome") %>% 
                 mutate(Cumulative_position = Position + Pos_max_chromosome_cumul),
               aes(x = Cumulative_position, y = Height_polygon, group = Inversion),
               color = "red", fill = NA) +
  scale_x_continuous(label = str_split_fixed(center_bin_chrom$Chromosome, "G", 2)[, 2],
                     breaks = center_bin_chrom$center) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30)) +
  labs(x = "Chromosomes",
       y = "Individual PC1 score",
       tag = "(A)") + 
  facet_col(facets = vars(Population), scales="free_y", space="free") +
  guides(fill = "none",
         color = guide_legend(override.aes = list(lwd = 2.5, alpha = 1)))


################## Merge everything together  ##################

ggarrange(local_pca,
          delta_freqs,
          F4_plot,
          nrow=3,
          common.legend = TRUE,
          legend = "top") %>% 
  ggsave(plot = ., filename = "../../Output/Figures/Figure_4.png", device = "png", units = "px",
         height = 2100, width = 1800, scale = 3.5)
