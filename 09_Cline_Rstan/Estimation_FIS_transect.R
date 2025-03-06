### Script for Stan fits of allele frequency clines with Fis
### 
###
### Roger Butlin - March 2022
###

# recommended rstan start-up options
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')

# add loo and bayesplot
library("loo")
library(readxl)
library("bayesplot")
library(dplyr)
library(MetBrewer)
austria <- met.brewer("Austria", 7)
archambault <- met.brewer("Archambault", 7)
peru <- met.brewer("Peru2", 6)
colorpalette <- c("Inv_1.1" = austria[1],
                  "Inv_2.1" = austria[2],
                  "Inv_3.1" = austria[3],
                  "Inv_4.1" = austria[5],
                  "Inv_4.2" = austria[6],
                  "Inv_6.1" = archambault[2],
                  "Inv_6.2" = archambault[3],
                  "Inv_6.3" = archambault[4],
                  "Inv_7.1" = archambault[5],
                  "Inv_8.1" = archambault[6],
                  "Inv_11.1" = austria[4],
                  "Inv_13.1" = archambault[1],
                  "Inv_14.1" = peru[1],
                  "Inv_16.1" = peru[5],
                  "Inv_16.2" = peru[6])

overall_data_allLG<-NULL
data_fited_allLG<-NULL
All_model_AF_allLG<-NULL
model_comp_allLG<-NULL

######################################################################################################################
## input data
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

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
# Import the delimiation of the inversions
delim_inversions <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
                               sep = "\t", header = TRUE) %>% 
  filter(Inversion %!in% c("Inv_3.2", "Inv_4.2", "Inv_14.2"))
# Import output pca
output_pca <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
                         sep = "\t", header = TRUE) %>% 
  inner_join(delim_inversions %>%
               select(-c(Start, End, Length, Inversion)) %>% 
               rename(Inversion = Inversion_grouped) %>% 
               unique,
             by = c("Chromosome", "Population", "Inversion"),
             relationship = "many-to-many") %>% 
  mutate(Inversion = Inversion %>% 
           factor(levels = read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
                                      sep = "\t", header = TRUE) %>%
                    select(Inversion) %>% 
                    unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist))


fab <- metadata %>% 
  left_join(output_pca, by = c("Population", "Sample_Name")) %>% 
  select(Inversion, Sample_Name, Population, Group, LCmeanDist) %>% 
  unique %>% 
  mutate(Group = Group - 1) %>% 
  pivot_wider(names_from = "Inversion", values_from = "Group")



params_clines <- data.frame()
cline_plots <- data.frame()
for (population in fab$Population %>% unique){
  print(population)
  
  input <- fab %>% 
    filter(Population == population) %>% 
    arrange(LCmeanDist)
  
  if (population == "Sweden"){
    input <- input %>% 
      select(-c("Inv_2.1"))
  }
  
  Inversion <- input %>% 
    select(starts_with("Inv_")) %>% 
    colnames()
  
  
  for(i in 1:length(Inversion)){
    inv_name <- Inversion[i]
    print(paste0("    ", inv_name))
    inv <- input[[inv_name]]
    
    
    nchains = 1
    
    ## fit AF model
    
    cline_data <- list(N = length(input$LCmeanDist), d_x = input$LCmeanDist, gen = inv, 
                       upc=max(input$LCmeanDist),loc=min(input$LCmeanDist), 
                       uplw=log(max(input$LCmeanDist)/1.5))
    
    AF <- stan(
      file = "/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/09_Cline_Rstan/AFcline_Fis.stan",  # Stan program - simple cline
      data = cline_data,      # named list of data
      chains = nchains,             # number of Markov chains
      warmup = 2000,          # number of warmup iterations per chain
      iter = 3000,            # total number of iterations per chain
      cores = 4,              # number of cores (could use one per chain)
      refresh = 500,          # progress shown
      control = list(adapt_delta = 0.95)
    )
    
    AF_summary <- summary(AF,probs=c(0.025,0.1,0.25,0.5,0.75,0.9,0.975)) 
    params_clines_i <- AF_summary$summary[c("centre", "w", "lp", "rp","f"),] %>% 
      as.data.frame %>% 
      rownames_to_column("parameters") %>% 
      mutate(Inversion = inv_name,
             Population = population) %>% 
      select(Population, Inversion, parameters, mean, "2.5%", "97.5%") %>% 
      pivot_longer(cols = c(mean, "2.5%", "97.5%"), names_to = "Params", values_to = "values") %>% 
      pivot_wider(names_from = "parameters", values_from = "values") %>% 
      rename(Centre = centre,
             Width = w,
             Left = lp,
             Right = rp,
             Fis = f)
    
    params_clines <- params_clines %>% 
      add_table_to_df_in_iteration(params_clines_i)
    
    params_clines_i <- params_clines_i %>% 
      filter(Params == "mean")
    
    # Plot the clines
    plot_cline_i <- clinef(
      x = input$LCmeanDist,
      g = inv,
      n = 2,
      centre = params_clines_i$Centre,
      width = params_clines_i$Width,
      left = params_clines_i$Left,
      right = params_clines_i$Right,
      optimisation=FALSE
    ) %>% 
      rename(LCmeanDist = position,
             Freq = phen_cline) %>% 
      mutate(Population = population,
             Inversion = inv_name) %>% 
      cbind(AF_summary$summary[grepl("f_", rownames(AF_summary$summary)),] %>%
              as.data.frame %>%
              select(mean)) %>% 
      rename(Fis = mean)
    rownames(plot_cline_i) <- NULL
    
    
    cline_plots <- add_table_to_df_in_iteration(cline_plots,
                                                plot_cline_i)
  }
}

params_clines %>% 
  select(Centre, Width, Params, Population) %>% 
  filter(Params == "mean") %>% 
  pivot_longer(cols = c("Centre", "Width"),
               names_to = "Parameters",
               values_to = "values") %>% 
  ggplot(aes(x = Parameters, y = values)) +
  geom_violin() +
  facet_col(vars(Population), scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 20))

params_clines %>% 
  select(Left, Right, Fis, Params, Population) %>% 
  filter(Params == "mean") %>% 
  pivot_longer(cols = c("Left", "Right", "Fis"),
               names_to = "Parameters",
               values_to = "values") %>% 
  ggplot(aes(x = Parameters, y = values)) +
  geom_violin() +
  facet_col(vars(Population), scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 20))

params_clines %>% 
  write.table("/shared/home/bpajot/fabalis/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/09_Cline_Rstan/Params_clines_inversions.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

cline_plots %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France")),
         Inversion = Inversion %>% 
           factor(levels = cline_plots %>% 
                    select(Inversion) %>% 
                    unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist)) %>% 
  ggplot() +
  geom_line(aes(x = LCmeanDist, y = 1 - Freq, colour = Inversion), lwd = 1.2) +
  scale_color_manual(values = colorpalette) +
  labs(x = "Position along the transect (m)",
       y = "Inversion frequency") + 
  facet_col(vars(Population), scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 20))

cline_plots %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France")),
         Inversion = Inversion %>% 
           factor(levels = cline_plots %>% 
                    select(Inversion) %>% 
                    unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist)) %>% 
  ggplot() +
  geom_line(aes(x = LCmeanDist, y = Fis, colour = Inversion), lwd = 1.2) +
  scale_color_manual(values = colorpalette) +
  labs(x = "Position along the transect (m)",
       y = "Inversion frequency") + 
  facet_col(vars(Population), scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 20))
