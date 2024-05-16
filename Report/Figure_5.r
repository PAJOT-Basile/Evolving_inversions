# Import libraries
require("anyLib")
#devtools::install_github("thomasp85/patchwork")
anyLib(c("tidyverse", "adegenet", "vcfR", "readxl", "statgenGWAS", "ggforce", "ggh4x", "patchwork", "ggpubr"))


################## Useful functions  ##################
"%!in%" <- function(x, y){!(x %in% y)}

source("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/Functions_optimise_plot_clines.r")

################################ Useful variables ################################
# Basic theme to use in the graphs
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

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

################## Making the delta frequencies correlation plot  ##################
# Calculate the differences in allelic frequencies
Delta_freqs_whole_genome <- get_delta_freqs_and_F4(
  genetic_data = data,
  Extreme_values = pca$li,
  var = "Axis2",
  metadata = metadata
) %>% 
  left_join(pca$co %>% 
              rownames_to_column("Position"),
            by = "Position")

################## Plot the differences in allelic frequencies along genome  ##################
delta_freqs <- (Delta_freqs_whole_genome %>% 
                  select(Position, contains("Delta")) %>% 
                  pivot_longer(cols = contains("Delta"), names_to = "Population", values_to = "Delta_freqs") %>% 
                  mutate(Population = ifelse(Population == "Delta_freq_Sweden", "Suède", "France") %>% 
                           factor(levels = c("Suède", "France"))) %>% 
                  plot_manhattan(Delta_freqs, facets = Population, size = 3, palette = c("grey71", "turquoise4"))) +
  theme(text = element_text(size = 30)) +
  labs(x = "Groupes de liaison",
       y = expression(Delta * "freq"),
       tag = "(B)")

################## Plot the F4  ##################
F4_plot <- (Delta_freqs_whole_genome %>% 
              select(Position, F4_stat) %>% 
              plot_manhattan(F4_stat, size = 3, absolute = FALSE, palette = c("grey71", "turquoise4"))) +
  labs(x = "Groupes de liaison",
       y = "F4",
       tag = "(C)") +
  theme(text = element_text(size = 30)) +
  geom_hline(yintercept = 0, color = "black", lwd = 1.2)

################## Run the local PCA  ##################
SNP_positions <- data@loc.fac %>% 
  as_tibble %>% 
  rename(pos=value) %>% 
  separate(pos, c("LG", "SUPER", "SUPER_frag", "Position"), "_") %>% 
  unite("Chromosome", LG, SUPER, SUPER_frag, remove=TRUE) %>% 
  mutate(bin = (as.numeric(Position)/500000 %>% trunc) %>% floor,
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

pca_result_per_bin <- read.csv("./Local_PCA_500k_bins.csv", header = TRUE, sep = "\t")
pca_result_per_bin_modif <- pca_result_per_bin %>% 
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
           factor(levels = c("Sweden", "France"))) %>% 
  filter(abs(Axis1) < 50)
################## Plot results of local PCA  ##################
pca_data_cum <- pca_result_per_bin_modif %>% 
  group_by(Chromosome) %>% 
  mutate(Bin = Bin %>% as.numeric) %>% 
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
  mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France")),
         Habitat = ifelse(Habitat == "Exposed", "Exposé", ifelse(Habitat == "Sheltered", "Abrité", "Transition")) %>% 
           factor(levels = c("Exposé", "Transition", "Abrité"))) %>% 
  ggplot() +
  geom_polygon(data = chromosome_delimitations_polygons %>% 
                 mutate(Population = ifelse(Population == "Sweden", "Suède", "France") %>% 
                          factor(levels = c("Suède", "France"))),
               aes(x=bin_delim_x, y=bin_delim_y, group=Chromosome, fill = Color), alpha=0.3) +
  scale_fill_manual(values = chromosome_delimitations_polygons$Color %>% unique,
                    labels = chromosome_delimitations_polygons$Color %>% unique) +
  new_scale_color() +
  geom_line(aes(x = bin_cum, y = Axis1 * 5, color = Habitat, group = Sample_Name), lwd = 0.1, alpha = 0.5) +
  #geom_line(aes(x = bin_cum, y = Axis1 * 5, color = Length, group = Sample_Name), lwd = 0.2) +
  #scale_color_gradientn(colours = c("#4e79a7","grey75","#f28e2b"),
  #                      limits = c((metadata$Length %>% min)-0.01,
  #                                 (metadata$Length %>% max)+0.01),
  #                      breaks = c(8, 13, 18)) +
  scale_color_manual(name = "Habitat",
                     values = c("Exposé" = "orange2", "Transition" = "deeppink", "Abrité" = "dodgerblue3")) +
  scale_x_continuous(label = str_split_fixed(center_bin_chrom$Chromosome, "G", 2)[, 2],
                     breaks = center_bin_chrom$center) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30)) +
  labs(x = "Groupes de liaison",
       y = "Variations des scores\nindividuels sur le PC1",
       tag = "(A)") + 
  facet_col(facets = vars(Population), scales="free_y", space="free")+
  guides(fill = "none",
         color = guide_legend(override.aes = list(lwd = 2.5, alpha = 1)))


################## Merge everything together  ##################

ggarrange(local_pca,
          delta_freqs,
          F4_plot,
          nrow=3,
          common.legend = TRUE,
          legend = "top")
