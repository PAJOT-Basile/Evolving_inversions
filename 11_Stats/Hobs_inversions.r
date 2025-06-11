########## Libraries ##########
libraries <- c("tidyverse", "vcfR", "adegenet", "readxl", "see", "ggh4x")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(libraries, character.only = TRUE)
rm(libraries)

########## Useful variables ##########
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))
########## Useful functions ##########
source("../General_scripts/Functions_optimise_plot_clines.r")

calculate_heterozygosity <- function(genetic_data, exposition, ...){
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
      Compute_Hobs(...)
  }
  # Rename the list to get the same names as in entry
  names(Hobs) <- names(exposition)
  
  return(Hobs)
}

Bind_rows <- function(xlist, colName = "Hobs", add_names = TRUE){
  
  colName <- as.name(substitute(colName))
  
  df2return <- data.frame()
  Exposition_indivs <- c()
  for (list_element in names(xlist)){
    if (dim(df2return)[1] == 0){
      df2return <- tibble(
        column_name = xlist[[list_element]]
      ) %>% 
        rename(!!quo_name(colName) := column_name)
    }else{
      df2return <- df2return %>% 
        rbind(xlist[[list_element]] %>% 
                as_tibble %>% 
                rename(!!quo_name(colName) := value))
    }
    
    if (add_names){
      Exposition_indivs <- c(
        Exposition_indivs,
        rep(list_element, length(xlist[[list_element]]))
      )
    }
  }
  if (add_names){
    df2return <- df2return %>% 
      mutate(Exposition = Exposition_indivs)
  }
  
  return(df2return)
  
}

Find_positions_in_split_invs <- function(delim_invs, genetic_data, .chromosome = chromosome){
  
  delims_inv_i <- delim_invs %>% 
    filter(Inversion == inversion) 
  
  Number_delims <- delims_inv_i %>% 
    group_by(Population) %>% 
    summarize(count = n()) %>% 
    pull(count) %>% 
    unique
  
  pos2keep <- data.frame()
  for (i in 1:Number_delims){
    # Get the start of the inversion
    Starts_split_inversion <- delims_inv_i %>% 
      arrange(Start) %>% 
      dplyr::slice(Number_delims * (i - 1) + 1) %>% 
      pull(Start)
    
    # Get the end of the inversion
    if (Number_delims == 1) i <- dim(delims_inv_i)[1]
    Ends_split_inversion <- delims_inv_i %>% 
      arrange(End) %>% 
      dplyr::slice(Number_delims * i) %>% 
      pull(End)
    
    # Subset the genetic data
    pos2keep <- pos2keep %>% 
      add_table_to_df_in_iteration(
        genetic_data@loc.fac[which(grepl(chromosome, genetic_data@loc.fac))] %>% 
          as_tibble() %>% 
          rename(Position = value) %>% 
          transform_position_ade2tidy() %>% 
          unique %>% 
          filter(Position <= Ends_split_inversion,
                 Position >= Starts_split_inversion) %>% 
          unite(SNP_Name, Chromosome, Position, sep = "_")
      )
  }
  return(pos2keep)
}



#
########## Import data ##########
# Import the genetic data
data <- read.vcfR("../../Output/Sweden_France_parallelism/02_Filter_VCF/09_Maf_thin/VCF_File.vcf.gz") %>% 
  vcfR2genind()

# Import the metadata of the individuals
metadata <- read_excel(path = "../../Input_Data/Data/data_Fabalis_resequencing_Basile.xlsx",
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

# Import the pca outputs
output_pca <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
                         sep = "\t", header = TRUE)

# Import the delimitation of the inversions
delim_invs <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
                        sep = "\t", header = TRUE) %>% 
  # Prepare the variables to trace the trees
  mutate(Is_split = case_when(
           Inversion_grouped %in% c("Inv_3.1", "Inv_4.1", "Inv_14.1") ~ TRUE,
           TRUE ~ FALSE
         ),
         Inversion = Inversion_grouped,
         Inversion = ifelse(Inversion == "Inv_4.3", "Inv_4.2", Inversion))

########## Compute the heterozygosity for one inversion ##########
summary_df <- data.frame()
for (inversion in delim_invs$Inversion %>% unique){
  print(inversion)
  
  delim_inv_i <- delim_invs %>% 
    filter(Inversion == inversion)
  
  chromosome <- delim_inv_i$Chromosome %>% 
    unique
  
  list_positions_in_inversion <- Find_positions_in_split_invs(delim_inv_i, data, chromosome) %>% 
    pull(SNP_Name)
  print("    Subsetted")
  
  # Subset the genetic data to get the positions in the inversion
  data_inv <- data[loc = list_positions_in_inversion]
  # Get the exposition of each individual so we can classify them
  summary_exposition_individuals <- delim_inv_i %>% 
    inner_join(output_pca %>% 
                 filter(Inversion != "Inv_4.2") %>% 
                 mutate(Inversion = ifelse(Inversion == "Inv_4.3", "Inv_4.2", Inversion)), by = c("Population", "Inversion", "Chromosome")) %>% 
    filter(Inversion == inversion) %>% 
    find_exposition_of_individuals()
  
  # COmpute the individual heterozygosity
  Heterozygosity <- calculate_heterozygosity(data_inv, summary_exposition_individuals, MARGIN = 2)
  print("    Heterozygotes")
  # Transform into a dataframe to be usable in the tidy format
  summary_df <- summary_df %>%
    add_table_to_df_in_iteration(
      Heterozygosity %>%
            Bind_rows() %>%
            cbind(summary_exposition_individuals %>%
                    Bind_rows(colName = "Sample_Name", add_names = FALSE)) %>%
            mutate(Genotype = case_when(
              Exposition == "Exposed_indivs" ~ "EE",
              Exposition == "Transition" ~ "ES",
              TRUE ~ "SS"
            ),
            Chromosome = chromosome,
            Inversion = inversion)
      )
  
}


# Plot the results
(summary_df %>% 
    left_join(metadata, by = c("Sample_Name")) %>% 
    mutate(Population = Population %>% 
             factor(levels = c("Sweden" ,"France")),
           Genotype = Genotype %>% 
             factor(levels = c("EE", "ES", "SS")),
           Inversion = Inversion %>% 
             factor(levels = summary_df %>% 
                      select(Inversion) %>% 
                      unique %>% 
                      arrange(as.numeric(gsub("\\D*(\\d+).*\\.", "\\1", Inversion))) %>% 
                      pull(Inversion))) %>% 
    filter(Inversion %in% c("Inv_1.1", "Inv_2.1", "Inv_3.1",
                            "Inv_4.1", "Inv_4.2", "Inv_6.1",
                            "Inv_6.2", "Inv_6.3")) %>% 
    ggplot(aes(x = Genotype, y = Hobs)) +
    geom_point(aes(color = Genotype), position = position_jitterdodge(), size = 2, alpha = 0.7) +
    geom_violinhalf(aes(fill = Genotype), position = position_nudge(x = 0.2), alpha = 0.5, trim = FALSE) +
    scale_color_manual(values = c("EE" = "#4e79a7", "ES" = "grey75", "SS" = "#f28e2b")) +
    scale_fill_manual(values = c("EE" = "#4e79a7", "ES" = "grey75", "SS" = "#f28e2b")) +
    facet_grid2(Population ~ Inversion) +
    labs(x = "PCA clusters",
         y = "Individual heterozygosity") +
    my_theme) %>% 
  ggsave(filename = "../../Output/Sweden_France_parallelism/11_Stats/Heterozygosity_per_inv_1.png",
         scale = 4, height = 500, width = 1800, units = "px")

(summary_df %>% 
    left_join(metadata, by = c("Sample_Name")) %>% 
    mutate(Population = Population %>% 
             factor(levels = c("Sweden" ,"France")),
           Genotype = Genotype %>% 
             factor(levels = c("EE", "ES", "SS")),
           Inversion = Inversion %>% 
             factor(levels = summary_df %>% 
                      select(Inversion) %>% 
                      unique %>% 
                      arrange(as.numeric(gsub("\\D*(\\d+).*\\.", "\\1", Inversion))) %>% 
                      pull(Inversion))) %>% 
    filter(Inversion %!in% c("Inv_1.1", "Inv_2.1", "Inv_3.1",
                            "Inv_4.1", "Inv_4.2", "Inv_6.1",
                            "Inv_6.2", "Inv_6.3")) %>% 
    ggplot(aes(x = Genotype, y = Hobs)) +
    geom_point(aes(color = Genotype), position = position_jitterdodge(), size = 2, alpha = 0.7) +
    geom_violinhalf(aes(fill = Genotype), position = position_nudge(x = 0.2), alpha = 0.5, trim = FALSE) +
    scale_color_manual(values = c("EE" = "#4e79a7", "ES" = "grey75", "SS" = "#f28e2b")) +
    scale_fill_manual(values = c("EE" = "#4e79a7", "ES" = "grey75", "SS" = "#f28e2b")) +
    facet_grid2(Population ~ Inversion) +
    labs(x = "PCA clusters",
         y = "Individual heterozygosity") +
    my_theme) %>% 
  ggsave(filename = "../../Output/Sweden_France_parallelism/11_Stats/Heterozygosity_per_inv_2.png",
         scale = 4, height = 500, width = 1800, units = "px")
