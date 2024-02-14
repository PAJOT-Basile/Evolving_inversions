# Libraries
require("anyLib")
anyLib(c("tidyverse", "readxl", "forcats", "ggh4x", "bbmle", "tie", "ggforce", "ggpubr", "ggnewscale", "rmarkdown"))

# Create useful functions
"%!in%" <- function(x, y)!("%in%"(x, y))
zip <- function(...){
  mapply(list, ..., SIMPLIFY=FALSE)
}

enumerate <- function(...){
  zip(ix=seq_along(..1), ...)
}

# Import the data
data <- read_excel(path = "../../Data/data_Fabalis_resequencing_Basile.xlsx",
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
  
  # select only the data we need (the one on fabalis)
  filter(Species == "FAB",
         Population != "BRE") %>% 
  
  # Add a new single location column
  mutate(Single_location = paste(Country, Population, Transect, sep="_"))

# Color palette to be reused everywhere with the shell size
size_palette = c("#4e79a7", "grey75", "#f28e2b")


# Deal with missing data (drop lines that have less than the calculated threshold of 5% of missing data)
threshold <- ((data %>% nrow) * .05) %>% floor
data <- data %>% drop_na(((data %>% is.na %>% colSums) < threshold) %>% names)

## Visualising variation in shell size depending on the position along the transect
data %>% 
  ggplot(aes(LCmeanDist, Length, color=Length)) +
  geom_point(size=3, alpha=0.7) +
  scale_colour_gradientn(name = "Shell size\n(mm)",
                         colors=size_palette, 
                         limits=c(min(data$Length, na.rm=TRUE) - 0.01, 
                                  max(data$Length, na.rm=TRUE) + 0.01),
                         breaks=c(8, 13, 18)) +
  labs(x="Position along the transect (m)",
       y="Shell size (mm)",
       title = "Distribution of shell size along the transect") +
  theme_bw() +
  facet_col(facets = vars(Single_location), scales="free", space="free") +
  theme(text = element_text(size = 20))

# Apply the three models to the colors (stable model, linear model and clinal model)
source("Cline_functions.R")

# As the clines are not all the same, we use some  different initial values for each one in the cline model
Priors_values_cline_pheno <- data.frame(Centre_size = c(200, 100, 100),
                                        Width_size = c(50, 40, 40),
                                        Left_size = c(11, 7.5, 8.8),
                                        Right_size = c(8, 12.5, 11),
                                        Sd_left_size = c(1, 4, 2),
                                        Sd_right_size = c(2, 4, 7),
                                        Sd_centre_size = c(2, 8, 8),
                                        Single_location = c("FRA_LAM_n", "SWE_LOK_n", "SWE_LOK_s"))

# We fit the cline model to the data
Clinal_model_size <- data %>% 
  merge(Priors_values_cline_pheno, by="Single_location") %>% 
  group_by(Single_location) %>% 
  bow(tie(Centre_size, Width_size, Left_size, Right_size, Std_left, Std_centre, Std_right) := mle2(cline_phen, list(centre=Centre_size %>% unique(),
                                                                                                                    w=Width_size %>% unique(),
                                                                                                                    left=Left_size %>% unique(),
                                                                                                                    right=Right_size %>% unique(),
                                                                                                                    sl=Sd_left_size %>% unique(),
                                                                                                                    sr=Sd_right_size %>% unique(),
                                                                                                                    sc=Sd_centre_size %>% unique()),
                                                                                                   data=list(phen=Length,
                                                                                                             position=LCmeanDist,
                                                                                                             n=1),
                                                                                                   method="L-BFGS-B",
                                                                                                   upper=list(centre=max(LCmeanDist),w=max(LCmeanDist)*2,
                                                                                                              left=15,right=15,sl=10,sc=20,sr=10),
                                                                                                   lower=list(centre=1,width=1,
                                                                                                              left=1,right=1,sl=0.1,sc=0.1,sr=0.1)) %>%  
        coef() %>% 
        round(digits = 3))

# We calculate the confidence intervals of the cline parameters
Confidence_interval_cline_size <- data %>% 
  merge(Priors_values_cline_pheno, by="Single_location") %>% 
  group_by(Single_location) %>% 
  bow(tie(Centre_2.5, Width_2.5, Left_2.5, Right_2.5, Sd_left_2.5, Sd_centre_2.5, Sd_right_2.5, Centre_97.5, Width_97.5, Left_97.5, Right_97.5, Sd_left_97.5, Sd_centre_97.5, Sd_right_97.5) :=
        mle2(cline_phen, list(centre=Centre_size %>% unique(),
                              w=Width_size %>% unique(),
                              left=Left_size %>% unique(),
                              right=Right_size %>% unique(),
                              sl=Sd_left_size %>% unique(),
                              sr=Sd_right_size %>% unique(),
                              sc=Sd_centre_size %>% unique()),
             data=list(phen=Length,
                       position=LCmeanDist,
                       n=1),
             method="L-BFGS-B",
             upper=list(centre=max(LCmeanDist),w=max(LCmeanDist)*2,
                        left=15,right=15,sl=10,sc=20,sr=10),
             lower=list(centre=1,width=1,
                        left=1,right=1,sl=0.1,sc=0.1,sr=0.1)) %>%  
        confint()) %>% 
  group_by(Single_location) %>% 
  # We use the calculated values of the confidence interval to visualize them on the plot
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
data1 <- data %>% 
  merge(Clinal_model_size, by="Single_location") %>% 
  merge(Confidence_interval_cline_size, by="Single_location")

# Prepare the plotting of results for the clinal model
LAMn <- data1 %>% filter(Single_location == "FRA_LAM_n")
LOKn <- data1 %>% filter(Single_location == "SWE_LOK_n")
LOKs <- data1 %>% filter(Single_location == "SWE_LOK_s")
dfs <- list(LAMn, LOKn, LOKs)
cline_fit_plots_size <- list()
cline_fit_plots_size_2.5 <- list()
cline_fit_plots_size_97.5 <- list()
for (i in enumerate(dfs)){
  df <- i[2] %>% as.data.frame()
  # We return the plot for the model
  cline_fit_plots_size <- append(cline_fit_plots_size, list(cline_phen(position = seq(from=min(df$LCmeanDist), to=max(df$LCmeanDist), by=1),
                                                                       centre = df$Centre_size %>% unique(),
                                                                       w = df$Width_size %>% unique(),
                                                                       left = df$Left_size %>% unique(),
                                                                       right = df$Right_size %>% unique(),
                                                                       sl = df$Std_left %>% unique(),
                                                                       sc = df$Std_centre %>% unique(),
                                                                       sr = df$Std_right %>% unique(),
                                                                       optimisation=FALSE,
                                                                       plotting=TRUE)))
  # Here, we prepare the plotting for the smallest estimate of the cline model
  cline_fit_plots_size_2.5 <- append(cline_fit_plots_size_2.5, list(cline_phen(position = seq(from=min(df$LCmeanDist), to=max(df$LCmeanDist), by=1),
                                                                               centre = df$Centre_2.5 %>% unique(),
                                                                               w = df$Width_ci %>% unique(),
                                                                               left = df$Left_2.5 %>% unique(),
                                                                               right = df$Right_2.5 %>% unique(),
                                                                               sl = df$Sd_left_ci %>% unique(),
                                                                               sc = df$Sd_centre_ci %>% unique(),
                                                                               sr = df$Sd_right_ci %>% unique(),
                                                                               optimisation=FALSE,
                                                                               plotting=TRUE)))
  # Here, we prepare the plotting for the biggest estimate of the cline model
  cline_fit_plots_size_97.5 <- append(cline_fit_plots_size_97.5, list(cline_phen(position = seq(from=min(df$LCmeanDist), to=max(df$LCmeanDist), by=1),
                                                                                 centre = df$Centre_97.5 %>% unique(),
                                                                                 w = df$Width_ci %>% unique(),
                                                                                 left = (df$Left_97.5 %>% unique()),
                                                                                 right = (df$Right_97.5 %>% unique()),
                                                                                 sl = df$Sd_left_ci %>% unique(),
                                                                                 sc = df$Sd_centre_ci %>% unique(),
                                                                                 sr = df$Sd_right_ci %>% unique(),
                                                                                 optimisation=FALSE,
                                                                                 plotting=TRUE)))
}
# We recombine the data together for the three locations
LAMn_plot_pheno <- cline_fit_plots_size_2.5[1] %>%
  as.data.frame() %>% 
  rename(phen_cline_2.5 = phen_cline,
         sd_cline_2.5 = sd_cline) %>% 
  merge(cline_fit_plots_size_97.5[1], by="position") %>% 
  rename(phen_cline_97.5 = phen_cline,
         sd_cline_97.5 = sd_cline) %>% 
  merge(cline_fit_plots_size[1], by="position") %>% 
  as.data.frame() %>% 
  mutate(Single_location = "FRA_LAM_n")
LOKn_plot_pheno <- cline_fit_plots_size_2.5[2] %>%
  as.data.frame() %>% 
  rename(phen_cline_2.5 = phen_cline,
         sd_cline_2.5 = sd_cline) %>% 
  merge(cline_fit_plots_size_97.5[2], by="position") %>% 
  rename(phen_cline_97.5 = phen_cline,
         sd_cline_97.5 = sd_cline) %>% 
  merge(cline_fit_plots_size[2], by="position") %>% 
  as.data.frame() %>% 
  mutate(Single_location = "SWE_LOK_n")
LOKs_plot_pheno <- cline_fit_plots_size_2.5[3] %>%
  as.data.frame() %>% 
  rename(phen_cline_2.5 = phen_cline,
         sd_cline_2.5 = sd_cline) %>% 
  merge(cline_fit_plots_size_97.5[3], by="position") %>% 
  rename(phen_cline_97.5 = phen_cline,
         sd_cline_97.5 = sd_cline) %>% 
  merge(cline_fit_plots_size[3], by="position") %>% 
  as.data.frame() %>% 
  mutate(Single_location = "SWE_LOK_s")

# We bind the three locations together and add this to the whole dataset
Plotting_data_pheno <- rbind(LAMn_plot_pheno, LOKn_plot_pheno, LOKs_plot_pheno)

# We plot the results
data1 %>% 
  ggplot(aes(LCmeanDist, Length)) +
  geom_point(aes(color=Length), size=2) +
  facet_col(facets = vars(Single_location),
            scales="free_x",
            space="free") +
  scale_colour_gradientn(name = "Shell size\n(mm)",
                         colors=size_palette, 
                         limits=c(min(data$Length, na.rm=TRUE) - 0.01, 
                                  max(data$Length, na.rm=TRUE) + 0.01),
                         breaks=c(8, 13, 18)) +
  new_scale_color() +
  geom_ribbon(data = subset(Plotting_data_pheno, phen_cline < phen_cline_97.5), 
              aes(x=position, y=phen_cline, ymin=phen_cline, ymax=phen_cline_97.5), 
              fill = "blue", 
              alpha=0.2) +
  geom_ribbon(data = subset(Plotting_data_pheno, phen_cline > phen_cline_2.5), 
              aes(x=position, y=phen_cline, ymin=phen_cline_2.5, ymax=phen_cline), 
              fill = "blue", 
              alpha=0.2) +
  geom_line(data=Plotting_data_pheno, aes(x=position, y=phen_cline, color="Regression", linetype="Regression"), linewidth=2) +
  geom_line(data=Plotting_data_pheno, aes(x=position, y=phen_cline_2.5, color="Confidence interval", linetype = "Confidence interval")) +
  geom_line(data=Plotting_data_pheno, aes(x=position, y=phen_cline_97.5, color="Confidence interval", linetype = "Confidence interval")) +
  scale_linetype_manual(name = "Curve",
                        guide = "legend",
                        values = c("Regression" = 1, "Confidence interval" = 2)) +
  scale_color_manual(name = "Curve",
                     guide = "legend",
                     values = c("Regression"="orange", "Confidence interval"="lightblue4"),
                     ) +
  guides(colour = guide_legend(keywidth = 3, keyheight = 1), linetype = guide_legend(keywidth = 3, keyheight = 1)) +
  labs(x="Position along the transect (m)",
       y="Shell size (mm)",
       title = "Distribution of shell size along the transect",
       linetype = "Courbe",
       color = "Courbe") +
  theme(text = element_text(size = 20))

# TODO: check for model suitability (not necessary ?)

# Recapitulation table
size_table <- data1 %>% 
  group_by(Single_location) %>% 
  summarise(Centre_size = Centre_size %>% unique %>% round(digits=1),
            Centre_size_ci = paste0("[", Centre_2.5 %>% unique %>% round(digits=1), ", ", Centre_97.5 %>% unique %>% round(digits=1), "]"),
            Left_size = Left_size %>% unique %>% round(digits=1),
            Left_size_ci = paste0("[", Left_2.5 %>% unique %>% round(digits=1), ", ", Left_97.5 %>% unique %>% round(digits=1), "]"),
            Right_size = Right_size %>% unique %>% round(digits=1),
            Right_size_ci = paste0("[", Right_2.5 %>% unique %>% round(digits=1), ", ", Right_97.5 %>% unique %>% round(digits=1), "]"),
            Width = Width_size %>% unique %>% round(digits=1),
            Width_ci = paste0("[", Width_2.5 %>% unique %>% round(digits=1), ", ", Width_97.5 %>% unique %>% round(digits=1), "]")) %>% paged_table()
size_table  
