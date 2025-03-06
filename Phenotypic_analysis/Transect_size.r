# Libraries
if (!require("anyLib")) install.packages("anyLib")
anyLib::anyLib(c("tidyverse", "readxl", "ggh4x", "bbmle", "tie", "ggforce", "ggpubr", "ggnewscale", "rmarkdown"))

################################ Useful functions ################################
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

zip <- function(...){
  mapply(list, ..., SIMPLIFY=FALSE)
}

enumerate <- function(...){
  zip(ix=seq_along(..1), ...)
}
################################ Useful variables ################################
# Color palette to be reused everywhere with the shell size
size_palette = c("#4e79a7", "grey75", "#f28e2b")
################################ Import metadata ################################
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
         Shell.colour = factor(Shell.colour %>% str_to_title, levels = c("Black", "Black/Square", "Brown", "Brown/Square", "Dark", "Yellow", "Yellow/Brown", "Yellow/Square", "Grey", "White", "Banded", NA)),
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
  
  # Rename the column to get rid of the space (transformed into a dot in the importation)
  rename(Shell_color = Shell.colour) %>% 
  
  # Select only the necessary columns for the analysis
  select(-c(Penial.glands, NGI_ID, target, TRANSECT, ID, length, biSIZE, Dist, Angle)) %>% 
  
  # select only the metadata we need (the one on fabalis)
  filter(Species == "FAB",
         Population != "BRE") %>% 
  
  # Add a new single location column
  mutate(Single_location = paste(Country, Population, Transect, sep="_")) %>% 
  filter(Single_location %in% c("FRA_LAM_n", "SWE_LOK_n")) %>% 
  mutate(Single_location = ifelse(Single_location == "FRA_LAM_n", "France", "Sweden") %>% factor(levels = c("Sweden", "France")))


# Deal with missing data (drop lines that have less than the calculated threshold of 5% of missing data)
#threshold <- ((metadata %>% nrow) * .05) %>% floor
#metadata <- metadata %>% drop_na(((metadata %>% is.na %>% colSums) < threshold) %>% names)

## Visualising variation in shell size depending on the position along the transect
metadata %>% 
  ggplot(aes(LCmeanDist, Length, color=Length)) +
  geom_point(size=3, alpha=0.7) +
  scale_colour_gradientn(name = "Shell size\n(mm)",
                         colors=size_palette, 
                         limits=c(min(metadata$Length, na.rm=TRUE) - 0.01, 
                                  max(metadata$Length, na.rm=TRUE) + 0.01),
                         breaks=c(8, 13, 18)) +
  labs(x="Position along the transect (m)",
       y="Shell size (mm)",
       title = "Distribution of shell size along the transect") +
  theme_bw() +
  facet_col(facets = vars(Single_location), scales="free", space="free") +
  theme(text = element_text(size = 20))

################################ Fit the cline model to the data ################################
# As the clines are not all the same, we use some  different initial values for each one in the cline model
Priors_values_cline_pheno <- data.frame(Centre_size = c(200, 100),
                                        Width_size = c(50, 40),
                                        Left_size = c(11, 7.5),
                                        Right_size = c(8, 12.5),
                                        Sd_left_size = c(1, 4),
                                        Sd_right_size = c(2, 4),
                                        Sd_centre_size = c(2, 8),
                                        Single_location = c("France", "Sweden"))

# We fit the cline model to the data
Clinal_model_size <- metadata %>% 
  merge(Priors_values_cline_pheno, by="Single_location") %>% 
  group_by(Single_location) %>% 
  bow(tie(Centre_size, Width_size, Left_size, Right_size,
          Std_left, Std_centre, Std_right) :=
        mle2(cline_phen, list(centre = Centre_size %>% unique(),
                              width = Width_size %>% unique(),
                              left = Left_size %>% unique(),
                              right = Right_size %>% unique(),
                              sl = Sd_left_size %>% unique(),
                              sr = Sd_right_size %>% unique(),
                              sc = Sd_centre_size %>% unique()),
             data = list(phen = Length,
                       x = LCmeanDist,
                       n = 1),
             method = "L-BFGS-B",
             upper = list(centre = max(LCmeanDist),
                          width = max(LCmeanDist)*2,
                          left = 15,
                          right = 15,
                          sl = 10,
                          sc = 20,
                          sr = 10),
             lower = list(centre = 1,
                          width = 1,
                          left = 1,
                          right = 1,
                          sl = 0.1,
                          sc = 0.1,
                          sr = 0.1)) %>%
        
        coef() %>% 
        round(digits = 3))

################################ Confidence intervals ################################
# We calculate the confidence intervals of the cline parameters
Confidence_interval_cline_size <- metadata %>% 
  left_join(Clinal_model_size, by="Single_location") %>% 
  group_by(Single_location) %>% 
  bow(tie(Centre_2.5, Width_2.5, Left_2.5, Right_2.5, Sd_left_2.5, Sd_centre_2.5, Sd_right_2.5, Centre_97.5, Width_97.5, Left_97.5, Right_97.5, Sd_left_97.5, Sd_centre_97.5, Sd_right_97.5) :=
        mle2(cline_phen, list(centre = Centre_size %>% unique(),
                              width = Width_size %>% unique(),
                              left = Left_size %>% unique(),
                              right = Right_size %>% unique(),
                              sl = Std_left %>% unique(),
                              sr = Std_right %>% unique(),
                              sc = Std_centre %>% unique()),
             data = list(phen = Length,
                       x = LCmeanDist,
                       n = 1),
             method = "L-BFGS-B",
             upper = list(centre = max(LCmeanDist),
                          width = max(LCmeanDist)*2,
                          left = 15,
                          right = 15,
                          sl = 10,
                          sc = 20,
                          sr = 10),
             lower = list(centre = 1,
                          width = 1,
                          left = 1,
                          right = 1,
                          sl = 0.1,
                          sc = 0.1,
                          sr = 0.1)) %>%  
        confint()) %>% 
  group_by(Single_location) %>% 
  # Here, we reorganise the calculated confidence interval values to get the 
  # smallest value in the 2.5th percentile and the 97.5th percentile as the largest
  # number
  mutate(Width_ci = max(Width_2.5, Width_97.5),
         Left_2.5i = min(Left_2.5, Left_97.5),
         Right_2.5i = min(Right_2.5, Right_97.5),
         Centre_2.5i = ifelse(Left_2.5i > Right_2.5i, min(Centre_2.5, Centre_97.5), max(Centre_2.5, Centre_97.5)),
         Left_97.5 = max(Left_2.5, Left_97.5),
         Right_97.5 = max(Right_2.5, Right_97.5),
         Centre_97.5 = ifelse(Left_2.5i > Right_2.5i, max(Centre_2.5, Centre_97.5), min(Centre_2.5, Centre_97.5)),
         Sd_left_ci = max(Sd_left_2.5, Sd_left_97.5, na.rm=TRUE),
         Sd_right_ci = max(Sd_right_2.5, Sd_right_97.5, na.rm=TRUE),
         Sd_centre_ci = max(Sd_centre_2.5, Sd_centre_97.5, na.rm=TRUE)) %>% 
  select(-c(Left_2.5, Right_2.5, Centre_2.5, Sd_left_2.5, Sd_left_97.5, Sd_right_2.5, Sd_right_97.5, Sd_centre_2.5, Sd_centre_97.5)) %>% 
  rename(Left_2.5 = Left_2.5i,
         Right_2.5 = Right_2.5i,
         Centre_2.5 = Centre_2.5i)

# We merge the data with the cline model parameters estimated above
data1 <- metadata %>% 
  left_join(Clinal_model_size, by = "Single_location") %>% 
  left_join(Confidence_interval_cline_size, by = "Single_location")

################################ Plot the estimated models ################################
# Prepare the plotting of results for the clinal model
LAMn <- data1 %>% filter(Single_location == "France")
LOKn <- data1 %>% filter(Single_location == "Sweden")
dfs <- list(LAMn, LOKn)
cline_fit_plots_size <- list()
for (i in enumerate(dfs)){
  df <- i[2] %>% as.data.frame()
  # We return the plot for the model
  cline_fit_plots_size <- append(cline_fit_plots_size, list(cline_phen(x = seq(from = min(df$LCmeanDist),
                                                                               to = max(df$LCmeanDist),
                                                                               by = 1),
                                                                       centre = df$Centre_size %>% unique(),
                                                                       w = df$Width_size %>% unique(),
                                                                       left = df$Left_size %>% unique(),
                                                                       right = df$Right_size %>% unique(),
                                                                       sl = df$Std_left %>% unique(),
                                                                       sc = df$Std_centre %>% unique(),
                                                                       sr = df$Std_right %>% unique(),
                                                                       optimisation=FALSE)))
}
# We recombine the data together for the two locations
LAMn_plot_pheno <- cline_fit_plots_size[1] %>% 
  as.data.frame() %>% 
  mutate(Single_location = "France")
LOKn_plot_pheno <- cline_fit_plots_size[2] %>% 
  as.data.frame() %>% 
  mutate(Single_location = "Sweden")

# We bind the three locations together and add this to the whole dataset
Plotting_data_pheno <- rbind(LAMn_plot_pheno, LOKn_plot_pheno) %>% 
  mutate(Single_location = Single_location %>% 
           factor(levels = c("Sweden", "France")))
# We plot the results
data1 %>% 
  ggplot(aes(LCmeanDist, Length)) +
  geom_point(aes(color = Length), size=2) +
  scale_colour_gradientn(name = "Shell size\n(mm)",
                         colors = size_palette, 
                         limits = c(min(metadata$Length, na.rm=TRUE) - 0.01, 
                                  max(metadata$Length, na.rm=TRUE) + 0.01),
                         breaks = c(8, 13, 18)) +
  geom_line(data = Plotting_data_pheno, aes(x = position, y = phen_cline),
            linewidth = 2, color = "orange2") +
  facet_col(facets = vars(Single_location),
            scales = "free_x",
            space = "free") +
  labs(x = "Position along the transect (m)",
       y = "Shell size (mm)",
       title = "Distribution of shell size along the transect") +
  theme(text = element_text(size = 20))


# Recapitulation table
size_table <- data1 %>% 
  group_by(Single_location) %>% 
  # Here, we modify the table so it can be understood easily. We merge some 
  # columns together
  summarise(Centre_size = Centre_size %>% unique %>% round(digits=1),
            Centre_size_ci = paste0("[", min(Centre_2.5, Centre_97.5) %>% unique %>% round(digits=1), ", ", max(Centre_97.5, Centre_2.5) %>% unique %>% round(digits=1), "]"),
            Left_size = Left_size %>% unique %>% round(digits=1),
            Left_size_ci = paste0("[", Left_2.5 %>% unique %>% round(digits=1), ", ", Left_97.5 %>% unique %>% round(digits=1), "]"),
            Right_size = Right_size %>% unique %>% round(digits=1),
            Right_size_ci = paste0("[", Right_2.5 %>% unique %>% round(digits=1), ", ", Right_97.5 %>% unique %>% round(digits=1), "]"),
            Width = Width_size %>% unique %>% round(digits=1),
            Width_ci = paste0("[", Width_2.5 %>% unique %>% round(digits=1), ", ", Width_97.5 %>% unique %>% round(digits=1), "]")) %>%
  mutate(Pays = Single_location %>% 
           factor(levels = c("Sweden", "France")),
         Centre_size = paste(Centre_size, Centre_size_ci),
         Left_size = paste(Left_size, Left_size_ci),
         Right_size = paste(Right_size, Right_size_ci),
         Width_size = paste(Width, Width_ci)) %>% 
  rename("Centre du cline" = Centre_size,
         "Zone abritée" = Left_size,
         "Zone exposée" = Right_size,
         "Largeur du cline" = Width_size) %>% 
  select(Pays, "Zone abritée", "Zone exposée", "Centre du cline", "Largeur du cline") %>% 
  paged_table()
size_table  
