################## Libraries  ##################
libraries <- c("vcfR", "adegenet", "tidyverse", "readxl")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)

colfunc <- colorRampPalette(colors = c("#813d32","#ef852f","#fdc441" ,"#e3cfb4"))
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))
################## Useful functions  ##################
source("../General_scripts/Functions_optimise_plot_clines.r")

# Function to get the exposition of a group of individuals using the DAPC cluster outputs
# This function requires a Group (DAPC clustering) and a Sample_Name (name of the samples) column to work
find_exposition_of_individuals <- function(df){
  Exposed <- (df %>% 
                filter(Group == 1))$Sample_Name %>% 
    unique
  trans_indivs <- (df %>% 
                     filter(Group == 2))$Sample_Name %>% 
    unique
  sheltered_indivs <- (df %>% 
                         filter(Group == 3))$Sample_Name %>% 
    unique
  
  return(list(
    "Exposed" = Exposed,
    "Transition" = trans_indivs,
    "Sheltered" = sheltered_indivs
  ))
}

# Function to compute the heterozygozity of a group of indiviudals
calculate_heterozygosity <- function(genetic_data, exposition){
  # Add the sample names to the genetic table
  df <- genetic_data@tab %>% 
    as_tibble %>% 
    mutate(Sample_Name = rownames(genetic_data@tab))
  
  # Initialise the loop to treat every level of the list separatedly
  Hobs <- vector(mode = "list", length = length(exposition))
  # For each karyotype, calculate the mean heterozygosity
  for (level in names(exposition)){
    # Get the position of the karyotype in the list
    nb_in_list <- which(names(exposition) == level)
    
    # Get only the individuals with one karyotype
    # Get the number of missing data
    Hobs[[nb_in_list]] <- df %>% 
      filter(Sample_Name %in% exposition[[level]]) %>% 
      Compute_Hobs()
  }
  # Rename the list to get the same names as in entry
  names(Hobs) <- names(exposition)
  
  return(Hobs)
}

################## Import metadata  ##################
metadata <- read_excel(path = "../../Input_Data/Data/Phenotypic/data_Fabalis_resequencing_Basile.xlsx",
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
  droplevels

################## Import the pca ouptut data  ##################
output_pca <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
                         sep = "\t", header = TRUE)

################## Import the inversion delimitations  ##################
delimitation_inversions <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Candidate_inversion_post_pca.tsv",
                                      sep = "\t", header = TRUE) %>% 
  inner_join(output_pca, by = c("Inversion", "Population", "Chromosome")) %>%
  pivot_wider(names_from = "Population", values_from = c("Start", "End", "Length")) %>% 
  relocate(ends_with("Sweden"), .after = Inversion) %>% 
  arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
  group_by(Inversion) %>% 
  mutate(Start_large_delimitations = min(Start_France, Start_Sweden, na.rm = TRUE),
         End_large_delimitations = max(End_France, End_Sweden, na.rm = TRUE))


################## Import the genetic data  ##################
data <- read.vcfR("../../Output/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz") %>% 
  vcfR2genind()


################## Loop over inversions  ##################
SV_stats <- data.frame()
for (inversion in delimitation_inversions$Inversion %>% unique){
  cat("\n", inversion)
  cat(inversion, "\t", progress_bar(0, 3))
  chromosome <- (delimitation_inversions %>% 
                   filter(Inversion == inversion))$Chromosome %>% 
    unique
  Start_inv <- (delimitation_inversions %>% 
              filter(Inversion == inversion))$Start_large_delimitations %>% 
    unique
  
  End_inv <- (delimitation_inversions %>% 
                filter(Inversion == inversion))$End_large_delimitations %>% 
    unique
  
  Length_inv <- (delimitation_inversions %>% 
                   filter(Inversion == inversion) %>% 
                   mutate(Length = End_large_delimitations - Start_large_delimitations))$Length %>% 
    unique
  # Take the large delimitations for the inversions
  ## Prepare the data
  ### Make a list of all the positions in the inversion
  list_positions_in_inversion <- seq(
    from = Start_inv,
    to = End_inv,
    by = 1
  )
  
  ### Isolate all the positions in the candidate inversion
  name_positions_in_inversion <- data@tab %>%
    colnames %>% 
    as_tibble %>% 
    rename(Position = value) %>% 
    filter(grepl(chromosome, Position)) %>% 
    transform_position_ade2tidy() %>% 
    filter(Position %in% list_positions_in_inversion) %>% 
    unite(Position, Chromosome, Position) %>% 
    unique %>% 
    as.vector %>% unname %>% unlist
  
  ### Subset the data for the positions in the table
  data_inv <- data[loc = name_positions_in_inversion]
  cat(inversion, "\t", progress_bar(1, 3))
  ### Make a list of the exposed individuals on this inversion
  summary_exposition_individuals <- delimitation_inversions %>% 
    filter(Inversion == inversion) %>% 
    find_exposition_of_individuals()
  
  ## Calculate the heterozygosity on each arrangement
  Hobs_distrib <- calculate_heterozygosity(data_inv, summary_exposition_individuals)
  Hobs <- Hobs_distrib %>% 
    lapply(mean)
  
  Hobs_France_distrib <- calculate_heterozygosity(
    data_inv[which(grepl("LAM", indNames(data_inv)))],
    lapply(summary_exposition_individuals, function(x){
      x[which(grepl("LAM", x))]
    })
  )
  
  Hobs_France <- Hobs_France_distrib %>% 
    lapply(mean)
  
  Hobs_Sweden_distrib <- calculate_heterozygosity(
    data_inv[which(grepl("LOK", indNames(data_inv)))],
    lapply(summary_exposition_individuals, function(x){
      x[which(grepl("LOK", x))]
    })
  )
  
  Hobs_Sweden <- Hobs_Sweden_distrib %>% 
    lapply(mean)
  cat(inversion, "\t", progress_bar(2, 3))

  ## Calculate the mean LD
  LD_table <- read.table(paste0("../../Output/Sweden_France_parallelism/05_LD_computations/France/LD_values/", chromosome, ".geno.ld"),
                         sep = "\t", header = TRUE) %>% 
    rename(LD = "R.2") %>% 
    mutate(Population = "France") %>% 
    rbind(read.table(paste0("../../Output/Sweden_France_parallelism/05_LD_computations/Sweden/LD_values/", chromosome, ".geno.ld"),
                     sep = "\t", header = TRUE) %>% 
            rename(LD = "R.2") %>% 
            mutate(Population = "Sweden")) %>% 
    group_by(Population) %>% 
    filter(POS1 %in% list_positions_in_inversion,
           POS2 %in% list_positions_in_inversion) %>% 
    summarize(Mean_LD = mean(LD, na.rm = TRUE))
  
  cat(inversion, "\t", progress_bar(3, 3))
  df_add <- data.frame(
    "Chromosome" = chromosome,
    "Inversion" = inversion,
    "Start" = Start_inv,
    "End" = End_inv,
    "Length" = Length_inv,
    "Nb_Snps" = (ncol(data_inv@tab) / 2),
    "Sweden_LD" = ((LD_table %>% filter(Population == "Sweden"))$Mean_LD) %>% round(digits = 4),
    "France_LD" = ((LD_table %>% filter(Population == "France"))$Mean_LD) %>% round(digits = 4),
    "Hobs_EE" = (Hobs$Exposed) %>% round(digits = 4),
    "Hobs_ES" = (Hobs$Transition) %>% round(digits = 4),
    "Hobs_SS" = (Hobs$Sheltered) %>% round(digits = 4),
    "Hobs_Sweden_EE" = (Hobs_Sweden$Exposed) %>% round(digits = 4),
    "Hobs_Sweden_ES" = (Hobs_Sweden$Transition) %>% round(digits = 4),
    "Hobs_Sweden_SS" = (Hobs_Sweden$Sheltered) %>% round(digits = 4),
    "Hobs_France_EE" = (Hobs_France$Exposed) %>% round(digits = 4),
    "Hobs_France_ES" = (Hobs_France$Transition) %>% round(digits = 4),
    "Hobs_France_SS" = (Hobs_France$Sheltered) %>% round(digits = 4))
    
  if (!file.exists("../../Output/Sweden_France_parallelism/04_Inversions/Stats_inversions.tsv")){
    df_add %>%
      write.table("../../Output/Sweden_France_parallelism/04_Inversions/Stats_inversions.tsv", append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }else{
    df_add %>%
      write.table("../../Output/Sweden_France_parallelism/04_Inversions/Stats_inversions.tsv", append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  SV_stats <- add_table_to_df_in_iteration(SV_stats,
                                           df_add)
}
rm(inversion, chromosome, list_positions_in_inversion, df_add, LD_table, Distinctive_SNPs_France, Distinctive_SNPs_Sweden,
   Poly_France, Poly_Sweden, Poly_two_pops, Hobs, summary_exposition_individuals, data_inv)

################## Select the inversions that fit the description  ##################
SV_stats <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Stats/Stats_inversions.tsv",
                       header = TRUE, sep = "\t") %>% 
  select(-c(Nb_Snps, Hobs_EE, Hobs_ES, Hobs_SS, contains("LD"), Start, End, Length)) %>% 
  pivot_longer(cols = contains("Hobs"), names_to = "Genotype", values_to = "Hobs") %>% 
  mutate(Population = str_split_fixed(Genotype, "_", 3)[, 2] %>% 
           factor(levels = c("Sweden", "France")),
         Genotype = str_split_fixed(Genotype, "_", 3)[, 3]) %>% 
  pivot_wider(names_from = Genotype, values_from = Hobs)

# Filter to keep only putative inversions with a higher heterozygosity in the heterokaryote genotype of the inversions
SV_stats %>% 
  drop_na %>% 
  filter(ES > EE | ES > SS) %>% 
  arrange(Population,
          as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
  select(-c(EE, ES, SS)) %>% 
  left_join(read.table("../../Output/Sweden_France_parallelism/04_Inversions/Candidate_inversion_post_pca.tsv",
                       sep = "\t", header = TRUE),
            by = c("Chromosome", "Inversion", "Population")) %>% 
  # And save the outuput
  write.table("../../Output/Sweden_France_parallelism/04_Inversions/Candidate_inversions_post_stats.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
