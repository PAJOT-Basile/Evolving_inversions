################## Libraries  ##################
libraries <- c("vcfR", "adegenet", "tidyverse", "readxl")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)

colfunc <- colorRampPalette(colors = c("#813d32","#ef852f","#fdc441" ,"#e3cfb4"))
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))
################## Useful functions  ##################
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

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

compute_polymorphism <- function(genetic_data, exposition){
  # Add the sample names to the genetic table and filter the individuals we want
  df <- genetic_data@tab %>% 
    as_tibble %>% 
    mutate(Sample_Name = rownames(genetic_data@tab))
  
  # Initialise the loop
  Nb_poly_markers <- vector(mode = "list", length = length(exposition))
  for (level in names(summary_exposition_individuals)){
    # Get the position of the karyotype in the list
    nb_in_list <- which(names(exposition) == level)
    
    # Get only the individuals with one karyotype
    df_level <- df %>% 
      filter(Sample_Name %in% exposition[[level]])
    
    # Get the number of polymporhic loci
    nb_poly_loci <- df_level %>% 
      select(ends_with("0")) %>% 
      pivot_longer(cols = everything(), names_to = "Locus", values_to = "Genotype") %>%
      group_by(Locus) %>%
      summarize(distinct = n_distinct(Genotype)) %>% 
      filter(distinct > 1) %>% 
      nrow
    
    # Get the proportion of polymorphic loci
    Nb_poly_markers[nb_in_list] <- nb_poly_loci
  }
  names(Nb_poly_markers) <- names(exposition)
  return(Nb_poly_markers)
}

calculate_freq_genotype <- function(df){
  # Missing data
  Na_count <- df %>% 
    is.na %>% 
    colSums(na.rm = TRUE)
  
  # Sum of the p allele
  Nb_p <- (2 * (df == 0) + (df == 1)) %>% 
    colSums(na.rm = TRUE)
  
  # Compute p allele frequency
  p_freq <- Nb_p / ((nrow(df) - Na_count) * 2)
  
  # Compute q freq
  q_freq <- 1 - p_freq
  
  # Return the frequencies
  return(list("p" = p_freq, "q" = q_freq))
  
}

reorder_freqs <- function(freq_exp, freq_shelt){
  df <- freq_exp$p %>% 
    cbind(freq_shelt$p) %>% 
    as.data.frame %>% 
    rename(Expos = ".",
           Shelt = V2) %>% 
    mutate(Expos_freq = ifelse(Expos > Shelt, 1 - Expos, Expos),
           Shelt_freq = ifelse(Expos > Shelt, 1 - Shelt, Shelt))
  
  
  p_expos <- df$Expos_freq
  names(p_expos) <- names(freq_exp$p)
  
  p_shelt <- df$Shelt_freq
  names(p_shelt) <- names(freq_exp$p)
  return(list(
    "Exposed" = list("p" = p_expos, "q" = 1 - p_expos),
    "Sheltered" = list("p" = p_shelt, "q" = 1 - p_shelt)
  ))
}

compute_discriminant_snps <- function(genetic_data, exposition, diff_value = 1){
  # Add the sample names to the genetic table and filter the individuals we want
  df <- genetic_data@tab %>%
    as_tibble %>%
    select(ends_with("0")) %>%
    mutate(Sample_Name = rownames(genetic_data@tab),
           Karyotype = ifelse(Sample_Name %in% exposition$Exposed, "EE",
                              ifelse(Sample_Name %in% exposition$Sheltered, "SS", "ES")))
  
  # Calculate allele frequencies in the exposed individuals
  exposed_freqs <- df %>% 
    filter(Karyotype == "EE") %>% 
    select(-c(Sample_Name, Karyotype)) %>%
    calculate_freq_genotype()
  
  # Compute allele frequencies in the sheltered individuals
  sheltered_freqs <- df %>% 
    filter(Karyotype == "SS") %>% 
    select(-c(Sample_Name, Karyotype)) %>% 
    calculate_freq_genotype()
  
  
  # Reorder the allele frequencies so they can be comparable
  corrected_freqs <- reorder_freqs(exposed_freqs, sheltered_freqs)
  
  # Get the number of distinctive SNPs
  nb_discriminant_snps <- corrected_freqs$Exposed$p %>% 
    cbind(corrected_freqs$Sheltered$p) %>% 
    as.data.frame %>% 
    rename(Expos = ".",
           Shelt = V2) %>% 
    rownames_to_column("Position") %>% 
    filter(Shelt >= diff_value & Expos <= 1 - diff_value) %>% 
    nrow
  
  return(nb_discriminant_snps)
}

compute_heterozygote_deficit <- function(df){
  
  # Compute observed heterozygosity
  Hobs <- Compute_Hobs(df)
  
  # Compute expected heterozygosity
  ## Allele frequencies
  freqs <- calculate_freq_genotype(df)
  
  ## Expected heterozygosity
  Hexp <- 2 * freqs$p * freqs$q
  
  # Compute deficit of heterozygotes
  FIS <- (Hexp - Hobs) / Hexp
  return(FIS)
}

################## Import metadata  ##################
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
  droplevels

################## Import the cline parameters  ##################
cline_params <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/06_Cline_analysis/Cline_parameters.tsv",
                           sep = "\t", header = TRUE) %>% 
  pivot_wider(names_from = Population, values_from = c(Centre, Width, Left, Right))
################## Import the pca ouptut data  ##################
output_pca <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
                         sep = "\t", header = TRUE)

################## Import the inversion delimitations  ##################
delimitation_inversions <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Candidate_inversion_post_pca.tsv",
                                      sep = "\t", header = TRUE) %>% 
  inner_join(output_pca, by = c("Inversion", "Population", "Chromosome")) %>%
  pivot_wider(names_from = "Population", values_from = c("Start", "End", "Length")) %>% 
  relocate(ends_with("Sweden"), .after = Inversion) %>% 
  arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
  group_by(Inversion) %>% 
  mutate(Start_large_delimitations = min(Start_France, Start_Sweden, na.rm = TRUE),
         End_large_delimitations = max(End_France, End_Sweden, na.rm = TRUE))


################## Import the genetic data  ##################
data <- read.vcfR("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz") %>% 
  vcfR2genind()


################## Loop over inversions  ##################
SV_stats <- data.frame()
for (inversion in delimitation_inversions$Inversion %>% unique){
  cat("\n", inversion)
  cat(inversion, "\t", progress_bar(0, 5))
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
  cat(inversion, "\t", progress_bar(1, 5))
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
  cat(inversion, "\t", progress_bar(2, 5))
  # Create dataframe of heterozygosity
  (Hobs_distrib %>% 
    Transform_Hobs2df() %>% 
    mutate(Exposition = Names,
           Population = "Both") %>% 
    rbind(Hobs_France_distrib %>% 
            Transform_Hobs2df() %>% 
            mutate(Exposition = Names,
                   Population = "France"))  %>% 
    rbind(Hobs_Sweden_distrib %>% 
            Transform_Hobs2df() %>% 
            mutate(Exposition = Names,
                   Population = "Sweden")) %>% 
    mutate(Population = Population %>% 
             factor(levels = c("Sweden", "Both", "France")),
           Genotype = case_when(
             Exposition == "Exposed" ~ "EE",
             Exposition == "Sheltered" ~ "SS",
             TRUE ~ "ES"
           ) %>% 
             factor(levels = c("SS", "ES", "EE"))) %>% 
    ggplot() +
    geom_violin(aes(x = Genotype, y = Hobs, fill = Genotype), alpha = 0.6) +
    scale_fill_manual(values = c("EE" = "dodgerblue", "ES" = "grey41", "SS" = "orange2")) +
    facet_wrap(vars(Population)) +
    my_theme) %>% 
    ggsave(plot = ., filename = paste0("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Stats/", inversion, ".png"),
           device = "png", units = "px", width = 1200, height = 800, scale = 3.5)
  
  ## Calculate proportion of polymorphic SNPs with two pops
  # Poly_two_pops <- compute_polymorphism(data_inv, summary_exposition_individuals)
  # Poly_Sweden <- compute_polymorphism(
  #   data_inv[which(grepl("LOK", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LOK", x))]
  #   })
  # )
  
  # Poly_France <- compute_polymorphism(
  #   data_inv[which(grepl("LAM", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LAM", x))]
  #   })
  # )
  cat(inversion, "\t", progress_bar(3, 5))
  
  ## Find proportion of distinctive SNPs between exposed and Sheltered
  # Distinctive_SNPs_France <- compute_discriminant_snps(
  #   data_inv[which(grepl("LAM", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LAM", x))]
  #   })
  # )
  
  # Distinctive_SNPs_Sweden <- compute_discriminant_snps(
  #   data_inv[which(grepl("LOK", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LOK", x))]
  #   })
  # )
  
  # Diagnostic_Snps_both <- compute_discriminant_snps(data_inv, summary_exposition_individuals)
  
  # Distinctive_SNPs_France095 <- compute_discriminant_snps(
  #   data_inv[which(grepl("LAM", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LAM", x))]
  #   }), diff_value = 0.95
  # )
  
  # Distinctive_SNPs_Sweden095 <- compute_discriminant_snps(
  #   data_inv[which(grepl("LOK", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LOK", x))]
  #   }), diff_value = 0.95
  # )
  
  # Diagnostic_Snps_both095 <- compute_discriminant_snps(data_inv, summary_exposition_individuals, diff_value = 0.95)
  
  # Distinctive_SNPs_France090 <- compute_discriminant_snps(
  #   data_inv[which(grepl("LAM", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LAM", x))]
  #   }), diff_value = 0.90
  # )
  
  # Distinctive_SNPs_Sweden090 <- compute_discriminant_snps(
  #   data_inv[which(grepl("LOK", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LOK", x))]
  #   }), diff_value = 0.90
  # )
  
  # Diagnostic_Snps_both090 <- compute_discriminant_snps(data_inv, summary_exposition_individuals, diff_value = 0.90)
  
  # Distinctive_SNPs_France080 <- compute_discriminant_snps(
  #   data_inv[which(grepl("LAM", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LAM", x))]
  #   }), diff_value = 0.80
  # )
  
  # Distinctive_SNPs_Sweden080 <- compute_discriminant_snps(
  #   data_inv[which(grepl("LOK", indNames(data_inv)))],
  #   lapply(summary_exposition_individuals, function(x){
  #     x[which(grepl("LOK", x))]
  #   }), diff_value = 0.80
  # )
  
  # Diagnostic_Snps_both080 <- compute_discriminant_snps(data_inv, summary_exposition_individuals, diff_value = 0.80)
  
  # Heterozygote deficit 
  ## France
  # FIS_France <- delimitation_inversions %>% 
  #   filter(Inversion == inversion) %>% 
  #   left_join(metadata, by = "Sample_Name") %>% 
  #   left_join(cline_params, by = "Inversion") %>% 
  #   select(-contains("Sweden")) %>% 
  #   filter(Population == "France",
  #          LCmeanDist <= Centre_France + Width_France/2,
  #          LCmeanDist >= Centre_France - Width_France/2) %>% 
  #   mutate(Group = Group - 1) %>% 
  #   ungroup %>% 
  #   select(Group) %>% 
  #   compute_heterozygote_deficit()
  
  # FIS_France_left <- delimitation_inversions %>% 
  #   filter(Inversion == inversion) %>% 
  #   left_join(metadata, by = "Sample_Name") %>% 
  #   left_join(cline_params, by = "Inversion") %>% 
  #   select(-contains("Sweden")) %>% 
  #   filter(Population == "France",
  #          LCmeanDist < Centre_France - Width_France/2) %>% 
  #   mutate(Group = Group - 1) %>% 
  #   ungroup %>% 
  #   select(Group) %>% 
  #   compute_heterozygote_deficit()
  
  # FIS_France_right <- delimitation_inversions %>% 
  #   filter(Inversion == inversion) %>% 
  #   left_join(metadata, by = "Sample_Name") %>% 
  #   left_join(cline_params, by = "Inversion") %>% 
  #   select(-contains("Sweden")) %>% 
  #   filter(Population == "France",
  #          LCmeanDist > Centre_France + Width_France/2) %>% 
  #   mutate(Group = Group - 1) %>% 
  #   ungroup %>% 
  #   select(Group) %>% 
  #   compute_heterozygote_deficit()

  ## Sweden
  # FIS_Sweden <- delimitation_inversions %>% 
  #   filter(Inversion == inversion) %>% 
  #   left_join(metadata, by = "Sample_Name") %>% 
  #   left_join(cline_params, by = "Inversion") %>% 
  #   select(-contains("France")) %>% 
  #   filter(Population == "Sweden",
  #          LCmeanDist <= Centre_Sweden + Width_Sweden/2,
  #          LCmeanDist >= Centre_Sweden - Width_Sweden/2) %>% 
  #   mutate(Group = Group - 1) %>% 
  #   ungroup %>% 
  #   select(Group) %>% 
  #   compute_heterozygote_deficit()
  
  # FIS_Sweden_left <- delimitation_inversions %>% 
  #   filter(Inversion == inversion) %>% 
  #   left_join(metadata, by = "Sample_Name") %>% 
  #   left_join(cline_params, by = "Inversion") %>% 
  #   select(-contains("France")) %>% 
  #   filter(Population == "Sweden",
  #          LCmeanDist < Centre_Sweden - Width_Sweden/2) %>% 
  #   mutate(Group = Group - 1) %>% 
  #   ungroup %>% 
  #   select(Group) %>% 
  #   compute_heterozygote_deficit()
  
  
  # FIS_Sweden_right <- delimitation_inversions %>% 
  #   filter(Inversion == inversion) %>% 
  #   left_join(metadata, by = "Sample_Name") %>% 
  #   left_join(cline_params, by = "Inversion") %>% 
  #   select(-contains("France")) %>% 
  #   filter(Population == "Sweden",
  #          LCmeanDist > Centre_Sweden + Width_Sweden/2) %>% 
  #   mutate(Group = Group - 1) %>% 
  #   ungroup %>% 
  #   select(Group) %>% 
  #   compute_heterozygote_deficit()
  
    cat(inversion, "\t", progress_bar(4, 5))
  ## Calculate the mean LD
  LD_table <- read.table(paste0("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/05_LD_computations/France/LD_values/", chromosome, ".geno.ld"),
                         sep = "\t", header = TRUE) %>% 
    rename(LD = "R.2") %>% 
    mutate(Population = "France") %>% 
    rbind(read.table(paste0("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/05_LD_computations/Sweden/LD_values/", chromosome, ".geno.ld"),
                     sep = "\t", header = TRUE) %>% 
            rename(LD = "R.2") %>% 
            mutate(Population = "Sweden")) %>% 
    group_by(Population) %>% 
    filter(POS1 %in% list_positions_in_inversion,
           POS2 %in% list_positions_in_inversion) %>% 
    summarize(Mean_LD = mean(LD, na.rm = TRUE))
  
  cat(inversion, "\t", progress_bar(5, 5))
  df_add <- data.frame(
    "Chromosome" = chromosome,
    "Inversion" = inversion,
    "Start" = Start_inv,
    "End" = End_inv,
    "Length" = Length_inv,
    "Nb_Snps" = (ncol(data_inv@tab) / 2),
    # "FIS_Sweden" = FIS_Sweden %>% round(digits = 4),
    # "FIS_Sweden_left" = FIS_Sweden_left %>% round(digits = 4),
    # "FIS_Sweden_right" = FIS_Sweden_right %>% round(digits = 4),
    # "FIS_France" = FIS_France %>% round(digits = 4),
    # "FIS_France_left" = FIS_France_left %>% round(digits = 4),
    # "FIS_France_right" = FIS_France_right %>% round(digits = 4),
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
    "Hobs_France_SS" = (Hobs_France$Sheltered) %>% round(digits = 4)#,
    # "N_Poly_both_EE" = (Poly_two_pops$Exposed) %>% round(digits = 4),
    # "N_Poly_both_ES" = (Poly_two_pops$Transition) %>% round(digits = 4),
    # "N_Poly_both_SS" = (Poly_two_pops$Sheltered) %>% round(digits = 4),
    # "N_Poly_France_EE" = (Poly_France$Exposed) %>% round(digits = 4),
    # "N_Poly_France_ES" = (Poly_France$Transition) %>% round(digits = 4),
    # "N_Poly_France_SS" = (Poly_France$Sheltered) %>% round(digits = 4),
    # "N_Poly_Sweden_EE" = (Poly_Sweden$Exposed) %>% round(digits = 4),
    # "N_Poly_Sweden_ES" = (Poly_Sweden$Transition) %>% round(digits = 4),
    # "N_Poly_Sweden_SS" = (Poly_Sweden$Sheltered) %>% round(digits = 4),
    # "N_Diagnostic_Snps_Fr_Sw" = (Diagnostic_Snps_both) %>% round(digits = 4),
    # "N_Discriminant_Prop_France" = (Distinctive_SNPs_France) %>% round(digits = 4),
    # "N_Discriminant_Prop_Sweden" = (Distinctive_SNPs_Sweden) %>% round(digits = 4),
    # "N_Diagnostic_Snps_Fr_Sw0.95" = (Diagnostic_Snps_both095) %>% round(digits = 4),
    # "N_Discriminant_Prop_France0.95" = (Distinctive_SNPs_France095) %>% round(digits = 4),
    # "N_Discriminant_Prop_Sweden0.95" = (Distinctive_SNPs_Sweden095) %>% round(digits = 4),
    # "N_Diagnostic_Snps_Fr_Sw0.90" = (Diagnostic_Snps_both090) %>% round(digits = 4),
    # "N_Discriminant_Prop_France0.90" = (Distinctive_SNPs_France090) %>% round(digits = 4),
    # "N_Discriminant_Prop_Sweden0.90" = (Distinctive_SNPs_Sweden090) %>% round(digits = 4),
    # "N_Diagnostic_Snps_Fr_Sw0.80" = (Diagnostic_Snps_both080) %>% round(digits = 4),
    # "N_Discriminant_Prop_France0.80" = (Distinctive_SNPs_France080) %>% round(digits = 4),
    # "N_Discriminant_Prop_Sweden0.80" = (Distinctive_SNPs_Sweden080) %>% round(digits = 4)
  ) #%>% 
    # mutate(across(starts_with("N_"), ~ round(.x / Nb_Snps, digits = 4), .names = "p_{.col}")) %>% 
    # rename_all(
    #   .funs = list(
    #     function(x){
    #       ifelse(grepl("p_N_", make.names(names(.))),
    #              gsub("p_N_", "p_", make.names(names(.))),
    #              make.names(names(.)))
    #     }
    #   )
    # )
  if (!file.exists("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Stats_inversions.tsv")){
    df_add %>%
      write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Stats_inversions.tsv", append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }else{
    df_add %>%
      write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Stats_inversions.tsv", append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  SV_stats <- add_table_to_df_in_iteration(SV_stats,
                                           df_add)
}
rm(inversion, chromosome, list_positions_in_inversion, df_add, LD_table, Distinctive_SNPs_France, Distinctive_SNPs_Sweden,
   Poly_France, Poly_Sweden, Poly_two_pops, Hobs, summary_exposition_individuals, data_inv)

################## Select the inversions that fit the description  ##################
SV_stats <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Stats/Stats_inversions.tsv",
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
  left_join(read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Candidate_inversion_post_pca.tsv",
                       sep = "\t", header = TRUE),
            by = c("Chromosome", "Inversion", "Population")) %>% 
  # And save the outuput
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Candidate_inversions_post_stats.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
