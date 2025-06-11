#Script by Alan Le Moan : alan.le.moan@gmail.com
#Made for polar cod inversion structure and phylogeny on the 10/02/2023
# Modified by Basile Pajot for Littorina fabalis inversion structure and phylogeny on the 06/03/2025

################## Libraries  ##################
library(vcfR)
library(adegenet)
library(tidyverse)
library(ggpubr)
library(phangorn)
library(reshape2)
library(ggtree)
if (!require("MetBrewer")) install.packages("MetBrewer")
library(MetBrewer)

################## Import data  ##################
# Import the inversion delimitations
delim_invs <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
                         sep = "\t", header = TRUE) %>% 
  mutate(Is_inversion = TRUE,
         Inversion = Inversion_grouped) %>% 
  select(-Inversion_grouped) %>% 
  rbind(read.table("./Not_inversions.tsv",
                   sep = "\t", header = TRUE) %>%
          mutate(Is_inversion = FALSE,
                 Length = End - Start) %>%
          relocate(Length, .after = End)) %>% 
  mutate(Is_split = case_when(
    Inversion %in% c("Inv_3.1", "Inv_4.1", "Inv_14.1") ~ TRUE,
    TRUE ~ FALSE
  ))

# Import the output of the local PCA
groups_pca <- read.table("../../Output/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
                         sep = "\t", header = TRUE) %>% 
  # Filter to get only the inversions that are kept after the 4 step filtering of inversions
  inner_join(delim_invs,
             by = c("Chromosome", "Inversion"),
             relationship = "many-to-many")


# Import the reference individuals for the cases where we are working with colinear genome
ref_individuals <- read.table("../../Output/Sweden_France_parallelism/Reference_indivs/France_exposed.txt",
                sep = "\t", header = FALSE) %>%
  mutate(Population = "France",
         Exposition = "Exposed") %>%
  rbind(read.table("../../Output/Sweden_France_parallelism/Reference_indivs/France_sheltered.txt",
                  sep = "\t", header = FALSE) %>%
  mutate(Population = "France",
         Exposition = "Sheltered")) %>%
  rbind(read.table("../../Output/Sweden_France_parallelism/Reference_indivs/Sweden_exposed.txt",
                  sep = "\t", header = FALSE) %>%
  mutate(Population = "Sweden",
         Exposition = "Exposed")) %>%
    rbind(read.table("../../Output/Sweden_France_parallelism/Reference_indivs/Sweden_sheltered.txt",
                    sep = "\t", header = FALSE) %>%
  mutate(Population = "Sweden",
         Exposition = "Sheltered")) %>%
  rename(Sample_Name = V1)


# Define paths to use for softwares and where to find the vcf file
vcftools <- "/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools"
bcftools <- "/shared/software/miniconda/envs/bcftools-1.9/bin/bcftools"
vcf_file <- "../../Output/Sweden_France_parallelism/02_Filter_VCF/08_Hobs/VCF_File.vcf.gz"

# Path to write output
output_path <- "../../Output/Sweden_France_parallelism/08_Trees/"

print("Imported data")

################## Useful functions  ##################
# Function to get the exposition of the individuals depending on the Group they are given in the local PCA
find_exposition_of_individuals <- function(df){
  exposed_indivs <- (df %>% 
                       filter(Group == 1))$Sample_Name %>% 
    unique
  trans_indivs <- (df %>% 
                     filter(Group == 2))$Sample_Name %>% 
    unique
  sheltered_indivs <- (df %>% 
                         filter(Group == 3))$Sample_Name %>% 
    unique
  
  return(list(
    "Exposed_indivs" = exposed_indivs,
    "Transition" = trans_indivs,
    "Sheltered" = sheltered_indivs
  ))
}

# Function not in
'%!in%' <- function(x,y)!('%in%'(x,y))

# Create a directory if it does not exist
Create_dir_not_exist <- function(paths){
  for (path in paths){
    if (!dir.exists(path)){
      dir.create(path, recursive = TRUE)
    }
  }
}

Subset_genetic_data <- function(inversion, delim_invs, groups_pca, vcf_file, output_path, .is_inversion = is_inversion, .vcftools = vcftools, .bcftools = "/shared/software/miniconda/envs/bcftools-1.9/bin/bcftools",
                                force = FALSE){
  # First, create the directories to use in the function
  Create_dir_not_exist(c(paste0(output_path, "List_pos_per_inversion/"),
                         paste0(output_path, "List_indivs_inversion/"),
                         paste0(output_path, "VCF_inversion/")))
  
  
  if (force | !file.exists(paste0(output_path, "List_pos_per_inversion/", inversion, ".tsv"))){
    # Verify if we have a split inversion or not
    is_split <- delim_invs %>% 
      filter(Inversion == inversion) %>% 
      pull(Is_split) %>% 
      unique
    
    # Get the name of the chromosome we are working with
    chromosome <- (delim_invs %>% 
                     filter(Inversion == inversion))$Chromosome %>% 
      unique
    
    # Select the delimitations of the inversion to consider
    if (is_split){
      
      delims_inv_i <- delim_invs %>% 
        filter(Inversion == inversion) 
      
      Number_delims <- delims_inv_i %>% 
        group_by(Population) %>% 
        summarize(count = n()) %>% 
        pull(count) %>% 
        unique
      
      Starts_split_inversion <- c()
      Ends_split_inversion <- c()
      for (i in 1:Number_delims){
        Starts_split_inversion <- c(
          Starts_split_inversion,
          delims_inv_i %>% 
            arrange(Start) %>% 
            dplyr::slice(Number_delims * (i - 1) + 1) %>% 
            pull(Start)
        )
        Ends_split_inversion <- c(
          Ends_split_inversion,
          delims_inv_i %>% 
            arrange(End) %>% 
            dplyr::slice(Number_delims * i) %>% 
            pull(End)
        )
      }
      
      # Make a list of all the positions to keep
      pos2keep <- mapply(function(x, y) seq(x, y, 1), Starts_split_inversion, Ends_split_inversion) %>% 
        unlist
      
    }else{
      delims_inv_i <- delim_invs %>% 
        filter(Inversion == inversion) %>% 
        group_by(Inversion) %>% 
        summarize(Start = min(Start),
                  End = max(End)) %>% 
        ungroup
      
      # Make alist of all the positions to keep
      pos2keep <- seq(
        from = delims_inv_i$Start,
        to = delims_inv_i$End,
        by = 1
      )
    }
    
    pos2keep %>%
      as.data.frame %>% 
      rename(Pos = ".") %>% 
      mutate(Chr = chromosome) %>% 
      relocate(Pos, .after = Chr) %>% 
      write.table(paste0(output_path, "List_pos_per_inversion/", inversion, ".tsv"),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  # Make a list of the invdividuals to keep
  if (.is_inversion){
    list_expos_indivs <- groups_pca %>% 
      filter(Inversion == inversion) %>% 
      find_exposition_of_individuals()
  }else{
    list_expos_indivs <-  list(
      "Exposed_indivs" = (ref_individuals %>%
                            filter(Exposition == "Exposed"))$Sample_Name,
      "Sheltered" = (ref_individuals %>%
                       filter(Exposition == "Sheltered"))$Sample_Name
    )
  }
  
  if (force | !file.exists(paste0(output_path, "List_indivs_inversion/", inversion, ".txt"))){
    # Save the names of the individiuals that are homozygous for the two inversion genotypes
    list_expos_indivs$Exposed_indivs %>% 
      as.data.frame %>% 
      rbind(list_expos_indivs$Sheltered %>% 
              as.data.frame) %>% 
      write.table(paste0(output_path, "List_indivs_inversion/", inversion, ".txt"),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  if (force | !file.exists(paste0(output_path, "VCF_inversion/", inversion, ".vcf.gz"))){
    system2(.vcftools,
            args = paste0(" --gzvcf ", vcf_file,
                          " --positions ", output_path, "List_pos_per_inversion/", inversion, ".tsv",
                          " --keep ", output_path, "List_indivs_inversion/", inversion, ".txt",
                          " --recode --stdout | ", .bcftools,
                          " view -c 2 -Oz",
                          " > ", output_path, "VCF_inversion/", inversion, ".vcf.gz"))
  }
  return(list_expos_indivs)
}

# Sample genotypes effectively
Sample_genotypes <- function(df){
  
  # First, we need all the unique combinations of locus and sample name and their count in the dataframe
  unique_comb_count <- df %>% 
    group_by(locus, Sample_Name) %>% 
    summarize(count = n()) %>% 
    ungroup
  
  # Here we initialise the output data frame and iterate over the unique counts that we have
  Sampled_genotypes <- data.frame()
  for (unique_count in unique_comb_count$count %>% unique){
    
    # Filter the unique combinations that contain the right count
    unique_comb_i <- unique_comb_count %>% 
      filter(count == unique_count) %>% 
      arrange(locus, Sample_Name)
    
    # Select these unique combinations in the raw dataframe
    df_i <- df %>% 
      filter(locus %in% unique_comb_i$locus,
             Sample_Name %in% unique_comb_i$Sample_Name) %>% 
      arrange(locus, Sample_Name)
    
    # Sample a number between 1 and the number of time the combination is used in the original dataframe
    Random_positions_to_sample <- sample(1:unique_count, size = dim(unique_comb_i)[1], replace = TRUE)
    # Here we transform the random number into the actual position of the row in the dataframe
    Random_rows_to_sample_in_df <- Random_positions_to_sample + c(0, rep(unique_count, length(Random_positions_to_sample) - 1)) %>% cumsum
    # And sample these rows in the dataframe
    Positions_sampled_count <- df_i[Random_rows_to_sample_in_df, ]
    
    # Then, add these rows to the data frame to return
    if (dim(Sampled_genotypes)[1] == 0){
      Sampled_genotypes <- Positions_sampled_count
    }else{
      Sampled_genotypes <- Sampled_genotypes %>% 
        rbind(Positions_sampled_count)
    }
  }
  
  return(Sampled_genotypes)
}

# Function to prepare a haploid genome to trace the trees
Prepare_haploid_genome <- function(inversion, output_path, genetic_data = data, .vcftools = vcftools, force = FALSE){
  # First create the directories we need to use
  Create_dir_not_exist(c(paste0(output_path, "Fastas_inversions/"),
                         paste0(output_path, "Ped_inversions/")))
  
  if (force | !file.exists(paste0(output_path, "Ped_inversions/", inversion, ".ped"))){
    # Run this command to create a ped file (genotype of the individuals for each position)
    system2(.vcftools,args =c(paste0("--gzvcf ", output_path, "VCF_inversion/", inversion, ".vcf.gz",
                                     " --plink --out ", output_path, "Ped_inversions/", inversion)))
  }
  
  if (force | !file.exists(paste0(output_path, "Fastas_inversions/", inversion, "sequences_random.txt"))){
    # Import the newly created ped file
    pedfile <- read.table(paste0(output_path, "Ped_inversions/", inversion, ".ped"),
                          header = FALSE)
    
    # Transform the table to get a table with the genotype of each individual at each position
    data_fasta <- genetic_data@loc.fac %>% 
      rbind(pedfile[, 7:dim(pedfile)[2]]) %>% 
      t %>% 
      as.data.frame
    Sample_Name <- genetic_data@tab %>% rownames
    colnames(data_fasta) <- c("locus", Sample_Name)
    
    # Make a table with the genotype of every individual at every locus
    reshape_fasta <- data_fasta %>% 
      pivot_longer(!locus,
                   names_to = "Sample_Name",
                   values_to = "Genotype") %>% 
      mutate(Genotype = str_trim(Genotype))
    
    
    # select one random SNP per individuals and marker
    data_fasta_sub <- reshape_fasta %>%
      Sample_genotypes()
    
    # Replace missing data or weird data in the subset fasta
    data_fasta_sub$Genotype[which(data_fasta_sub$Genotype %!in% c("A","T","C","G","N"))] <- "N"
    
    # Transform 3 colomn dataframe into a classic haploid genotype matrix
    data_fasta_random <- data_fasta_sub %>% 
      pivot_wider(names_from = locus, values_from = Genotype)
    
    data_fasta_random %>% 
      select(-Sample_Name) %>% 
      write.table(paste0(output_path, "Fastas_inversions/", inversion, "sequences_random.txt"),
                sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  if (force | !file.exists(paste0(output_path, "Fastas_inversions/", inversion, "sequences_random.fasta"))){
    # Create fasta with the sampled random allele
    data_fasta_random %>% 
      select(Sample_Name) %>% 
      write.table(paste0(output_path, "Ped_inversions/", inversion, "Sample_Names.txt"), sep = "\t",
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    system(command = paste0("sed 's/$/_sequenceR/' ", output_path, "Ped_inversions/", inversion, "Sample_Names.txt",
                            " | sed 's/^/>/' > ", output_path, "Ped_inversions/", inversion, "Corrected.txt"))
    
    system(command = paste0("paste ", output_path, "Ped_inversions/", inversion, "Corrected.txt ", output_path, "Fastas_inversions/", inversion, "sequences_random.txt",
                            " | tr '\t' '\n' > ", output_path, "Fastas_inversions/", inversion, "sequences_random.fasta"))
    
    # Clean up the environment to only keep files that are useful for the next steps
    system(command = paste0("rm ",
                            output_path, "Ped_inversions/", inversion, "Sample_Names.txt ",
                            output_path, "Ped_inversions/", inversion, "Corrected.txt ",
                            output_path, "Fastas_inversions/", inversion, "sequences_random.txt ",
                            output_path, "VCF_inversion/", inversion, ".vcf.gz"))
  }
}

Get_Karyotype_name <- function(df, .list_expos_indivs = list_expos_indivs){
  df %>%
    mutate(Karyotype = ifelse((Sample_Name %in% .list_expos_indivs$Exposed) & (grepl("LOK", Sample_Name)), "SW_EE",
                              ifelse((Sample_Name %in% .list_expos_indivs$Sheltered) & (grepl("LOK", Sample_Name)), "SW_SS",
                                     ifelse((Sample_Name %in% .list_expos_indivs$Exposed) & (grepl("LAM", Sample_Name)), "FR_EE",
                                            ifelse((Sample_Name %in% .list_expos_indivs$Sheltered) & (grepl("LAM", Sample_Name)), "FR_SS", NA)))) %>%
             factor(levels = c("SW_EE", "SW_SS", "FR_EE", "FR_SS"))) %>%
    return()
}

# Function that loads a fasta file, computes the genetic distance between the sequences and draws a NJ tree
Run_and_trace_phylogeny <- function(inversion, .output_path = output_path, .list_expos_indivs = list_expos_indivs){
  
  # Load the fasta file
  data_phylo <- read.FASTA(paste0(.output_path, "Fastas_inversions/", inversion, "sequences_random.fasta"),
                           type = "DNA")
  
  # transform data
  phyDat <- phyDat(data_phylo, type = "DNA")
  fabalis <- as.phyDat(phyDat)
  
  # compute genetic distance based on polymorphic sites
  dna_dist <- dist.ml(fabalis, model = "F81")
  
  # draw NJ tree
  fabalis_NJ  <- NJ(dna_dist)
  fabalis_NJ$tip.label <- fabalis_NJ$tip.label %>% str_remove_all("_sequenceR")
  
  # Get the names of the individuals and classify them by karyotype
  Sample_Name <- fabalis_NJ$tip.label %>%
    as.data.frame %>% 
    rename(Sample_Name = ".") %>% 
    Get_Karyotype_name() %>%
    mutate(Population = ifelse(grepl("LOK", Sample_Name), "Sweden", "France"),
           Color_ID = case_when(
             Karyotype == "SW_EE" ~ met.brewer("Peru1", 7)[6],
             Karyotype == "SW_SS" ~ met.brewer("Peru1", 7)[5],
             Karyotype == "FR_EE" ~ met.brewer("Peru1", 7)[2],
             Karyotype == "FR_SS" ~ met.brewer("Peru1", 7)[3],
             TRUE ~ met.brewer("Peru1", 7)[4]
           ))
  
  # Regroup the individuals in the tree depending on their karyotype
  group_Info <- list(
    "SW_EE" = Sample_Name %>% filter(Karyotype == "SW_EE") %>% pull(Sample_Name),
    "SW_SS" = Sample_Name %>% filter(Karyotype == "SW_SS") %>% pull(Sample_Name),
    "FR_EE" = Sample_Name %>% filter(Karyotype == "FR_EE") %>% pull(Sample_Name),
    "FR_SS" = Sample_Name %>% filter(Karyotype == "FR_SS") %>% pull(Sample_Name)
  )
  
  fabalis_NJ_ <- groupOTU(fabalis_NJ, group_Info)
  
  # Plot the tree
  tree <- fabalis_NJ_ %>% 
    ggtree(aes(colour = group), layout = "unrooted", lwd = 1.05) +
    scale_color_manual(name = "Genotype",
                       values = c("FR_EE" = met.brewer("Peru1", 6)[1],
                                  "FR_SS" = met.brewer("Peru1", 6)[2],
                                  "SW_EE" = met.brewer("Peru1", 6)[6],
                                  "SW_SS" = met.brewer("Peru1", 6)[5])) +
    theme(text = element_text(size = 20))
  
  # Compute the distances between the ecotypes in the two countries
  #extract phylogenetic distance  
  phylo_dist <- cophenetic.phylo(fabalis_NJ)
  karyo_info_1 <- data.frame(cbind(fabalis_NJ$tip.label, Sample_Name$Karyotype %>% as.character))
  colnames(karyo_info_1) <- c("Var1", "Karyotype_ID1")
  karyo_info_2 <- data.frame(cbind(fabalis_NJ$tip.label, Sample_Name$Karyotype %>% as.character))
  colnames(karyo_info_2) <- c("Var2", "Karyotype_ID2")
  
  
  # #reshape phylogenetic distance matrix
  test <- as.matrix(phylo_dist)
  test[upper.tri(test, diag = TRUE)] <- NA
  test2 <- melt(test) %>% 
    filter(!is.na(value)) %>% 
    mutate(Population = ifelse((grepl("LOK", Var1) & grepl("LOK", Var2)), "Sweden", ifelse(
      (grepl("LAM", Var1) & grepl("LAM", Var2)), "France", NA)
    )) %>% 
    filter(!is.na(Population))
  test4 <- test2 %>% 
    left_join(karyo_info_2, by = "Var2") %>% 
    left_join(karyo_info_1, by = "Var1") %>% 
    rowwise() %>% 
    mutate(comp = paste0(sort(c(Karyotype_ID1, Karyotype_ID2)), collapse = "_"))
  
  #calculate mean net distance between group
  mean_divergence_karyotype <- aggregate(value ~ comp, data = test4, FUN = "mean") %>% 
    rename(Genotype = comp,
           Divergence = value) %>% 
    filter(Genotype %in% c("FR_EE_FR_SS", "SW_EE_SW_SS")) %>% 
    mutate(Divergence = Divergence %>% round(digits = 4))
  
  Divergence_karyotypes_countries <- data.frame(
    "Population" = c("Sweden", "France"),
    "Divergence" = c(
      (mean_divergence_karyotype %>% filter(grepl("SW", Genotype)))$Divergence,
      (mean_divergence_karyotype %>% filter(grepl("FR", Genotype)))$Divergence
    )
  )
  
  
  return(list(
    "Tree" = tree,
    "Divergence" = Divergence_karyotypes_countries
  ))
}

# Function that loads a genind object, scales the data, runs a pca and plots it
Run_and_trace_local_pca <- function(inversion, .output_path = output_path, .groups_pca = groups_pca, .ref_individuals = ref_individuals, .is_inversion = is_inversion){

  # Import the vcf file
  data <- read.vcfR(paste0(.output_path, "VCF_inversion/", inversion, ".vcf.gz")) %>% 
    vcfR2genind()
  
  # Scale the data to remove missing data and run a PCA on this filtered data
  pca <- scaleGen(data, NA.method="mean",scale=FALSE) %>% 
    dudi.pca(scale=FALSE, nf = 3,scannf = F)
  
  # Get the percentage of variance of axes
  pc1 <- round(pca$eig[1] / sum(pca$eig) * 100)
  pc2 <- round(pca$eig[2] / sum(pca$eig) * 100)

  if(.is_inversion){
    local_pca <- pca$li %>% 
      rownames_to_column("Sample_Name") %>%
      inner_join(.groups_pca %>% 
                  select(-starts_with("Axis")), by = "Sample_Name") %>% 
      Get_Karyotype_name() %>% 
      rename(Ecotype = Karyotype) %>%
      ggplot() +
      geom_point(aes(x = Axis1, y = Axis2, color = Ecotype), size = 3, alpha = 0.8) +
      scale_color_manual(name = "Genotype",
                        values = c("FR_EE" = met.brewer("Peru1", 6)[1],
                                    "FR_SS" = met.brewer("Peru1", 6)[2],
                                    "SW_EE" = met.brewer("Peru1", 6)[6],
                                    "SW_SS" = met.brewer("Peru1", 6)[5])) +
      labs(x = paste0("PC1 (",round(pca$eig[1]/sum(pca$eig)*100),"%)"),
          y = paste0("PC2 (",round(pca$eig[2]/sum(pca$eig)*100),"%)")) +
      theme_classic() +
      theme(text = element_text(size = 20))

  }else{
  local_pca <- pca$li %>% 
    rownames_to_column("Sample_Name") %>%
    inner_join(.ref_individuals, by = "Sample_Name") %>%
    Get_Karyotype_name () %>%
    rename(Ecotype = Karyotype) %>%
    ggplot() +
    geom_point(aes(x = Axis1, y = Axis2, color = Ecotype), size = 3, alpha = 0.8) +
    scale_color_manual(name = "Genotype",
                       values = c("FR_EE" = met.brewer("Peru1", 6)[1],
                                  "FR_SS" = met.brewer("Peru1", 6)[2],
                                  "SW_EE" = met.brewer("Peru1", 6)[6],
                                  "SW_SS" = met.brewer("Peru1", 6)[5])) +
    labs(x = paste0("PC1 (",round(pca$eig[1]/sum(pca$eig)*100),"%)"),
         y = paste0("PC2 (",round(pca$eig[2]/sum(pca$eig)*100),"%)")) +
    theme_classic() +
    theme(text = element_text(size = 20))
  }
return(list(
  "Plot" = local_pca,
  "Genetic_data" = data,
  "PCA" = pca
))
}

################## Run the functions for each inversion  ##################
# Subset the vcf file for the positions in the inversion
for (inversion in delim_invs$Inversion %>% unique){
  print(inversion)
  
  is_inversion <- delim_invs %>%
    filter(Inversion == inversion) %>%
    pull(Is_inversion) %>%
    unique
  
  chromosome <- delim_invs %>% 
    filter(Inversion == inversion) %>% 
    pull(Chromosome) %>% 
    unique
  
  # Subset the vcf file to keep only homokaryotypes for the inversion and the positions
  # inside the inversions
  list_expos_indivs <- Subset_genetic_data(inversion, delim_invs, groups_pca, vcf_file, output_path)
  print("    Subsetting done for the data")

  # Run the local PCA
  pca <- Run_and_trace_local_pca(inversion, output_path, groups_pca, ref_individuals)
  data <- pca$Genetic_data
  print("    Ran PCA")


  # Prepare a haploid genome by randomly sampling one genotype at each locus
  # to be able to draw the tree
  Prepare_haploid_genome(output_path, genetic_data = data)
  print("    Prepared haploid genome")
  

  # Run Phylogeny
  Phylogeny <- Run_and_trace_phylogeny(inversion, output_path, list_expos_indivs)
  print("    Ran Phylogeny")
  Divergence <- Phylogeny$Divergence
  
  # Create the directory if it does not exist
  Create_dir_not_exist(paste0(output_path, "Plots_inversions/"))

      ggarrange(pca$Plot,
            NA,
            Phylogeny$Tree,
            nrow=1,
            common.legend = TRUE,
            legend = "right",
            labels = inversion,
            hjust = -9.75,
            widths = c(2, 0.5, 2)) %>% 
    annotate_figure(top = text_grob(paste0("Inversion: ", inversion,
                        "\nDivergence Sweden (blue) = ", (Divergence %>% filter(Population == "Sweden"))$Divergence,
                        "\nDivergence France (red) = ", (Divergence %>% filter(Population == "France"))$Divergence))) %>%
    ggsave(plot = ., filename = paste0(output_path, "Plots_inversions/", inversion, ".png"), device = "png", units = "px",
           height = 850, width = 1500, scale = 3.5)

  
  
  # Add the divergence to the table
  if (!file.exists(paste0(output_path, "Divergence_karyotypes_countries.tsv"))){
    Divergence %>%
      mutate(Inversion = inversion,
             Chromosome = chromosome) %>% 
      relocate(Chromosome) %>% 
      relocate(Inversion, .after = Chromosome) %>% 
      write.table(paste0(output_path, "Divergence_karyotypes_countries.tsv"), append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }else{
    Divergence %>%
      mutate(Inversion = inversion,
             Chromosome = chromosome) %>% 
      relocate(Chromosome) %>% 
      relocate(Inversion, .after = Chromosome) %>% 
      write.table(paste0(output_path, "Divergence_karyotypes_countries.tsv"), append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
    print("    Done inversion")
  
}
