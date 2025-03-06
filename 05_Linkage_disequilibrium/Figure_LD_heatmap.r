
#libraries <- c("tidyverse")
#pacman::p_load(char = libraries, character.only = TRUE)
#rm(libraries)

library(tidyverse)
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")
# Make a function to sample colours in a user-defined scale
colfunc <- colorRampPalette(colors = c("#813d32","#ef852f","#fdc441" ,"#e3cfb4") %>% rev, bias = 2)

# Import and bind the files that are in the same directory
for (population in c("Sweden", "France")){
  print(population)
  
  Storage_directory <- paste0("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/05_LD_computations/", population, "/LD_values/")
  
  complete_df <- data.frame()
  for (file in list.files(Storage_directory)){
    print(paste0("    ", file))
    temp_df <- read.table(paste0(Storage_directory, file), sep = "\t", header = TRUE) %>% 
      rename(Chromosome = CHR,
             LD = R.2) %>% 
      mutate(Position1 = (POS1 / 1e5) %>% round,
             Position2 = (POS2 / 1e5) %>% round) %>% 
      group_by(Position1, Position2) %>% 
      summarize(mean_LD = mean(LD), .groups = "drop_last") %>% 
      ungroup
    chromosome <- (file %>% 
                     str_split_fixed("\\.", 2))[, 1]
    
    complete_df <- add_table_to_df_in_iteration(complete_df, temp_df %>% 
                                                  mutate(Chromosome = chromosome))
  }
  rm(temp_df, chromosome, file)
  # Plot the output and save it
  print("    Plotting")
  (complete_df %>% 
      mutate(Chromosome = Chromosome %>% 
               factor(levels = complete_df %>% select(Chromosome) %>% 
                        unique %>%
                        arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                        as.vector %>% unname %>% unlist)) %>% 
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
    ggsave(plot = ., filename = paste0("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/05_LD_computations/", population, "/LD_heatmaps/LD_heatmaps.png"),
           device = "png", width = 2000, height = 1800, units = "px", scale = 5, limitsize = FALSE)
}
print("Done")
