
############### Import libraries ###############
library(tidyverse)

############### Useful functions ###############
source("../General_scripts/Functions_optimise_plot_clines.r")

############### Useful variables ###############
# Make a function to sample colours in a user-defined scale
colfunc <- colorRampPalette(colors = c("#813d32","#ef852f","#fdc441" ,"#e3cfb4") %>% rev, bias = 2)

############### Plotting the LD per population ###############
# Loop over the populations
for (population in c("Sweden", "France")){
  print(population)

  # Get the name of the directory where the LD values are stored 
  Storage_directory <- paste0("../../Output/Sweden_France_parallelism/05_LD_computations/", population, "/LD_values/")

  # Initialise the output data frame containig the LD and iterate over the files in the storage directory to load them successively
  LD_df <- data.frame()
  for (file in list.files(Storage_directory)){
    print(paste0("    ", file))

    # Load the table of the considered chromosome
    temp_df <- read.table(paste0(Storage_directory, file), sep = "\t", header = TRUE) %>% 
      rename(Chromosome = CHR,
             LD = R.2) %>% 
      mutate(Position1 = (POS1 / 1e5) %>% round,
             Position2 = (POS2 / 1e5) %>% round) %>% 
    # Calculate the mean LD per bin (1e5 bp) to make the plotting easier
      group_by(Position1, Position2) %>% 
      summarize(mean_LD = mean(LD), .groups = "drop_last") %>% 
      ungroup

    # Get the name of the chromosome
    chromosome <- (file %>% 
                     str_split_fixed("\\.", 2))[, 1]

    # Add the loaded LD info for the chromosome to the complete dataframe
    LD_df <- add_table_to_df_in_iteration(LD_df, temp_df %>% 
                                                  mutate(Chromosome = chromosome))
  }
  
  # Clear the environment a bit
  rm(temp_df, chromosome, file)
  
  # Plot the output and save it
  print("    Plotting")
  (LD_df %>% 
   # Reorder the chromosome order for the facets to be in the correct order
      mutate(Chromosome = Chromosome %>% 
               factor(levels = LD_df %>% select(Chromosome) %>% 
                        unique %>%
                        arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                        as.vector %>% unname %>% unlist)) %>% 
   # Plot the output
      ggplot() +
      geom_tile(aes(x = Position1, y = Position2, fill = mean_LD)) +
      scale_fill_gradientn(name = "Mean LD\n(R²)",
                           colors = colfunc(20)) +
      facet_wrap(facets = vars(Chromosome), scales = "free") +
      theme_bw() +
      theme(text = element_text(size = 40),
            legend.key.size = unit(1.5, "cm")) +
      labs(x = "Position along the chromosome (*1e5)",
           y = "Position along the chromosome (*1e5)")) %>% 
  # Save the graph into a file
    ggsave(plot = ., filename = paste0("../../Output/Sweden_France_parallelism/05_LD_computations/", population, "/LD_heatmaps/LD_heatmaps.png"),
           device = "png", width = 2000, height = 1800, units = "px", scale = 5, limitsize = FALSE)
}
print("Done")
