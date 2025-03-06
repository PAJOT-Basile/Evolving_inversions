# Libraries
libraries <- c("tidyverse", "readxl", "ggh4x", "bbmle", "tie", "ggforce", "ggpubr", "ggnewscale", "LaplacesDemon", "rmarkdown", "ggnewscale")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(characters = libraries, character.only = TRUE)
rm(libraries)
################################ Useful functions ################################
zip <- function(...){
  mapply(list, ..., SIMPLIFY=FALSE)
}

enumerate <- function(...){
  zip(ix=seq_along(..1), ...)
}

# Apply the three models to the colors (stable model, linear model and clinal model)
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

################################ Useful variables ################################
# Color palette to be reused everywhere with the shell colors
shell_palette = c("black", "grey28", "brown", "saddlebrown", "wheat4", "gold", "goldenrod", "khaki2", "grey", "blue", "rosybrown4")

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



# Deal with missing metadata (drop lines that have less than the calculated threshold of 5% of missing metadata)
#threshold <- ((metadata %>% nrow) * .05) %>% floor
#metadata <- metadata %>% drop_na(((metadata %>% is.na %>% colSums) < threshold) %>% names)

################################ Color part of the importation ################################

# Modify the colors to a bi-allelic gene
Shell_color <- metadata$Shell_color


metadata$Shell_color_naive_color <- str_split_fixed(string=Shell_color, pattern = "/", n=2)[, 1]
metadata$Shell_color_morphology <- str_split_fixed(string=Shell_color, pattern = "/", n=2)[, 2]

metadata <- metadata %>% 
  mutate(Shell_color_naive_color = ifelse(Shell_color_naive_color %in% c("Yellow", "White", "Grey"), "Yellow", Shell_color_naive_color),
         Shell_color_morphology = ifelse(Shell_color_naive_color == "Banded", "Banded", Shell_color_morphology),
         Shell_color_naive_color = ifelse(Shell_color_naive_color %in% c("Black", "Brown", "Dark"), "Brown", Shell_color_naive_color),
         Shell_color_naive_color = ifelse(Shell_color_naive_color == "Banded", "Brown", Shell_color_naive_color),
         Shell_color_naive_color = Shell_color_naive_color %>% factor(levels = c("Yellow", "Brown")),
         Shell_color_morphology = ifelse(! Shell_color_morphology %in% c("Banded", "Square"), "Uniform", "Banded")) %>% 
  mutate(Genotype_shell_color_naive = ifelse(Shell_color_naive_color == "Yellow", 0, 1)) #%>% 
#mutate(Single_location = ifelse(Single_location == "France", "France", "Suède") %>% factor(levels = c("Suède", "France")),
#Shell_color_naive_color = ifelse(Shell_color_naive_color == "Yellow", "Jaune", "Marron"))

# First, we need the priors. As the values are not the same for each location, we separate them
Priors_values_cline_color <- data.frame(Centre_color_Prior = c(200, 70),
                                        Width_color_Prior = c(20, 8e-6),
                                        Left_color_Prior = c(0.25, 0.1),
                                        Right_color_Prior = c(0.8, 0.05),
                                        Min_centre_color_Prior = c(200, 50),
                                        Min_width_color_Prior = c(0, 0),
                                        Single_location = c("France", "Sweden"))

# Then, we calculate the initial color frequency to use in the stable model. We also merge with the chosen priors
metadata <- metadata %>% 
  group_by(Single_location) %>% 
  summarise(p_yellow_location = mean(Shell_color_naive_color == "Yellow"),
            p_brown_location = 1 - p_yellow_location) %>% 
  left_join(metadata, by="Single_location") %>% 
  left_join(Priors_values_cline_color, by="Single_location")

# We run the cline model 
Clinal_model_color_france <- metadata %>%
  filter(Single_location == "France") %>% 
  group_by(Single_location) %>%
  bow(tie(Center, Width, Left, Right) := mle2(clinef,list(centre=Centre_color_Prior %>% unique,
                                                          width=Width_color_Prior %>% unique,
                                                          left=Left_color_Prior %>% unique,
                                                          right=Right_color_Prior %>% unique),
                                              data=list(x=LCmeanDist,
                                                        g=Genotype_shell_color_naive,
                                                        n=1),
                                              method="L-BFGS-B",
                                              upper=list(centre=max(LCmeanDist),
                                                         width=max(LCmeanDist)*2,
                                                         left=0.999,
                                                         right=0.999),
                                              lower=list(centre=Min_centre_color_Prior %>% unique,
                                                         width=Min_width_color_Prior %>% unique,
                                                         left=0.001, right=0.001)) %>%
        coef() %>%
        round(digits=3))

Stable_model_color_sweden <- metadata %>% 
  filter(Single_location == "Sweden") %>% 
  group_by(Single_location) %>% 
  summarise(p_brown_all = mle2(stable, list(p_all = p_brown_location %>% unique),
                               data=list(x=LCmeanDist, g=Genotype_shell_color_naive, n=1),
                               method="L-BFGS-B",
                               upper=list(p_all=0.999),
                               lower=list(p_all=0.001)) %>% 
              coef() %>% round(digits = 3))

# Merge the three dataframes and calculate the AIC for each model
data1_color <- metadata %>% 
  left_join(Clinal_model_color_france, by="Single_location") %>% 
  left_join(Stable_model_color_sweden, by = "Single_location")

# Prepare the plotting of results for the clinal model
LAMn_color <- data1_color %>% filter(Single_location == "France")
LOKn_color <- data1_color %>% filter(Single_location == "Sweden")

LAM_color_curve <- clinef(x = seq(min(LAMn_color$LCmeanDist),
                                  max(LAMn_color$LCmeanDist),
                                  1),
                          centre = LAMn_color$Center %>% unique,
                          width = LAMn_color$Width %>% unique,
                          left = LAMn_color$Left %>% unique,
                          right = LAMn_color$Right %>% unique,
                          optimisation = FALSE) %>% 
  mutate(Single_location = "France" %>% factor(levels = c("Sweden", "France")))

LOKn_color_curve <- data.frame(
  phen_cline = LOKn_color$p_brown_all %>% unique,
  position = seq(min(LOKn_color$LCmeanDist),
                 max(LOKn_color$LCmeanDist),
                 1)
) %>% 
  mutate(Single_location = "Sweden" %>% factor(levels = c("Sweden", "France")))


Plotting_data_color <- rbind(LAM_color_curve, LOKn_color_curve) %>%
  # We add a correction for a dominance effect on the brown allele
  mutate(pheno_cline_p = sqrt(1 - phen_cline),
         pheno_cline_q = 1 - pheno_cline_p,
         pheno_cline_homo_brown = pheno_cline_q **2,
         pheno_cline_homo_yellow = pheno_cline_p **2,
         pheno_cline_hetero = 2 * pheno_cline_p * pheno_cline_q)


################################ Size part of the importation ################################
# As the clines are not all the same, we use some  different initial values for each one in the cline model
Priors_values_cline_pheno <- data.frame(Centre_size_Priors = c(200, 100),
                                        Width_size_Priors = c(50, 40),
                                        Left_size_Priors = c(11, 7.5),
                                        Right_size_Priors = c(8, 12.5),
                                        Sd_left_size_Priors = c(1, 4),
                                        Sd_right_size_Priors = c(2, 4),
                                        Sd_centre_size_Priors = c(2, 8),
                                        Single_location = c("France", "Sweden"))

# We fit the cline model to the metadata
Clinal_model_size <- metadata %>% 
  merge(Priors_values_cline_pheno, by="Single_location") %>% 
  group_by(Single_location) %>% 
  bow(tie(Centre_size, Width_size, Left_size, Right_size, Std_left, Std_centre, Std_right) :=
        mle2(cline_phen, list(centre = Centre_size_Priors %>% unique(),
                              width = Width_size_Priors %>% unique(),
                              left = Left_size_Priors %>% unique(),
                              right = Right_size_Priors %>% unique(),
                              sl = Sd_left_size_Priors %>% unique(),
                              sr = Sd_right_size_Priors %>% unique(),
                              sc = Sd_centre_size_Priors %>% unique()),
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

# We merge the metadata with the cline model parameters estimated above
data1_size <- metadata %>% 
  merge(Clinal_model_size, by="Single_location")


# Prepare the plotting of results for the clinal model
LAMn_size <- data1_size %>% filter(Single_location == "France")
LOKn_size <- data1_size %>% filter(Single_location == "Sweden")
dfs <- list(LAMn_size, LOKn_size)
cline_fit_plots_size <- list()
for (i in enumerate(dfs)){
  df <- i[2] %>% as.data.frame()
  # We return the plot for the model
  cline_fit_plots_size <- append(cline_fit_plots_size, list(cline_phen(x = seq(from=min(df$LCmeanDist), to=max(df$LCmeanDist), by=1),
                                                                       centre = df$Centre_size %>% unique(),
                                                                       width = df$Width_size %>% unique(),
                                                                       left = df$Left_size %>% unique(),
                                                                       right = df$Right_size %>% unique(),
                                                                       sl = df$Std_left %>% unique(),
                                                                       sc = df$Std_centre %>% unique(),
                                                                       sr = df$Std_right %>% unique(),
                                                                       optimisation=FALSE)))
}
# We recombine the metadata together for the three locations
LAMn_plot_pheno <- cline_fit_plots_size[1] %>% 
  as.data.frame() %>% 
  mutate(Single_location = "France")
LOKn_plot_pheno <- cline_fit_plots_size[2] %>% 
  as.data.frame() %>% 
  mutate(Single_location = "Sweden")

# We bind the three locations together and add this to the whole dataset
Plotting_data_size <- rbind(LAMn_plot_pheno, LOKn_plot_pheno) %>% 
  mutate(Cline = "Size")

################################ Merge the two analysis ################################
data1_size <- data1_size %>%
  select(-contains("Prior"))

data1_color <- data1_color %>%
  select(-contains("Prior"))


data2 <- data1_size %>%
  full_join(data1_color, by = (data1_size %>% names)[which((data1_size %>% names) %in% (data1_color %>% names))]) %>% 
  select(-c(Sample_Name, Species, ID_number, Country, Population, Sex, Mreads, Gbp, Q30, x, y, Transect, Id, Bi_size)) %>%
  select(-contains(c("p_"))) %>%
  rename(Centre_color = Center,
         Width_color = Width,
         Left_color = Left,
         Right_color = Right)


################################ Confidence intervals ################################
# We calculate the confidence intervals of the cline parameters
Confidence_interval_cline_size <- data2 %>% 
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

Conf_interval_cline_color <- data2 %>% 
  filter(Single_location == "France") %>% 
  group_by(Single_location) %>% 
  tie::bow(tie(Centre_2.5, Width_2.5, Left_2.5, Right_2.5, Centre_97.5, Width_97.5, Left_97.5, Right_97.5) :=
             mle2(clineflog, list(centre = Centre_color %>% unique,
                                  width = Width_color %>% unique %>% log,
                                  left = Left_color %>% unique %>% logit,
                                  right = Right_color %>% unique %>% logit),
                  data = list(x = LCmeanDist,
                              g = Genotype_shell_color_naive,
                              n = 1),
                  method ="L-BFGS-B",
                  upper = list(centre = max(LCmeanDist),
                               width = (max(LCmeanDist) * 2) %>% log,
                               left = 0.999 %>% logit,
                               right = 0.999 %>% logit),
                  lower = list(centre = 1,
                               width = 1 %>% log,
                               left = 0.001 %>% logit,
                               right = 0.001 %>% logit)) %>%
             confint()) %>%
  as.data.frame() %>% 
  mutate(Centre_2.5i = ifelse(min(Centre_2.5, Centre_97.5, na.rm = TRUE) != Inf, min(Centre_2.5, Centre_97.5, na.rm = TRUE), NA),
         Centre_97.5 = ifelse(min(Centre_2.5, Centre_97.5, na.rm = TRUE) != Inf, max(Centre_2.5, Centre_97.5, na.rm = TRUE), NA),
         Width_2.5i = ifelse(min(Width_2.5, Width_97.5, na.rm = TRUE) != Inf, min(Width_2.5, Width_97.5, na.rm = TRUE) %>% exp, NA),
         Width_97.5 = ifelse(min(Width_2.5, Width_97.5, na.rm = TRUE) != Inf, max(Width_2.5, Width_97.5, na.rm = TRUE) %>% exp, NA),
         Left_2.5i = ifelse(min(Left_2.5, Left_97.5, na.rm = TRUE) != Inf, min(Left_2.5, Left_97.5, na.rm = TRUE) %>% invlogit, NA),
         Left_97.5 = ifelse(min(Left_2.5, Left_97.5, na.rm = TRUE) != Inf, max(Left_2.5, Left_97.5, na.rm = TRUE) %>% invlogit, NA),
         Right_2.5i = ifelse(min(Right_2.5, Right_97.5, na.rm = TRUE) != Inf, min(Right_2.5, Right_97.5, na.rm = TRUE) %>% invlogit, NA),
         Right_97.5 = ifelse(min(Right_2.5, Right_97.5, na.rm = TRUE) != Inf, max(Right_2.5, Right_97.5, na.rm = TRUE) %>% invlogit, NA)) %>%
  select(Single_location, contains(".5i"), contains("97.")) %>% 
  rename(Left_2.5 = Left_2.5i,
         Right_2.5 = Right_2.5i,
         Centre_2.5 = Centre_2.5i,
         Width_2.5 = Width_2.5i)

################################ Make the figures ################################
# In english
data2 %>%
  mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
           factor(levels = c("Sweden", "France"))) %>% 
  ggplot(aes(x = LCmeanDist, y = Length)) +
  # The size cline
  geom_point(data = data2 %>%
               mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                        factor(levels = c("Sweden", "France")),
                      Cline = "Size" %>% 
                        factor(levels = c("Size", "Color", "Off")),
                      Shell_color_naive_color = Shell_color_naive_color %>% 
                        factor(levels = c("Brown", "Yellow"))),
             aes(colour=Shell_color_naive_color), size=2, alpha = 0.7) +
  geom_line(data=Plotting_data_size %>% 
              mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                       factor(levels = c("Sweden", "France")),
                     Cline = "Size" %>% 
                       factor(levels = c("Size", "Color", "Off"))),
            aes(x=position, y=phen_cline), color="black", linewidth=2) +
  geom_vline(data = data2 %>%
               mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                        factor(levels = c("Sweden", "France")),
                      Cline = "Size" %>% 
                        factor(levels = c("Size", "Color", "Off"))) %>% 
               filter(!is.na(Cline)) %>% 
               select(Centre_size, Cline, Single_location) %>% unique,
             aes(xintercept = Centre_size),
             linetype = "dashed", linewidth = 1.2) +
  geom_vline(data = Confidence_interval_cline_size %>% 
               mutate(Single_location = Single_location %>% 
                        factor(levels = c("Sweden", "France")),
                      Cline = "Size" %>% 
                        factor(levels = c("Size", "Color", "Off"))) %>% 
               select(Cline, Centre_97.5, Centre_2.5) %>% 
               pivot_longer(contains("Centre"), names_to = "Centile", values_to = "Centre"),
             aes(xintercept = Centre),
             linetype = "dashed", linewidth = 1, color = "grey40") +
  # Color scale for the shell size
  geom_point(data = data2 %>%
               mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                        factor(levels = c("Sweden", "France")),
                      Cline = "Color" %>% 
                        factor(levels = c("Size", "Color", "Off")),
                      Shell_color_naive_color = ifelse(Shell_color_naive_color == "Yellow", "Yellow", "Brown") %>% 
                        factor(levels = c("Brown", "Yellow"))),
             aes(y=Genotype_shell_color_naive, color = Shell_color_naive_color), size=2, alpha=0.7) +
  geom_line(data= Plotting_data_color %>% 
              mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>% 
                       factor(levels = c("Sweden", "France")),
                     Cline = "Color" %>% 
                       factor(levels = c("Size", "Color", "Off"))),
            aes(x=position, y=phen_cline), color = "black", linewidth=2) +
  geom_vline(data = data2 %>%
               mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                        factor(levels = c("Sweden", "France")),
                      Cline = "Color" %>% 
                        factor(levels = c("Size", "Color", "Off"))) %>% 
               filter(!is.na(Cline)) %>% 
               select(Centre_color, Cline, Single_location) %>% unique,
             aes(xintercept = Centre_color),
             linetype = "dashed", linewidth = 1.25) +
  geom_vline(data = Conf_interval_cline_color %>% 
               mutate(Single_location = Single_location %>% 
                        factor(levels = c("Sweden", "France")),
                      Cline = "Color" %>% 
                        factor(levels = c("Size", "Color", "Off"))) %>% 
               select(Single_location, Cline, Centre_2.5, Centre_97.5) %>% 
               pivot_longer(contains("Centre"), values_to = "Centre", names_to = "Centile"),
             aes(xintercept = Centre),
             linetype = "dashed", linewidth = 1, color = "grey40") +
  scale_color_manual(name = "Shell color",
                     values=c("brown", "orange"),
                     guide = "none") +
  # new_scale_color() +
  # geom_point(data = data2 %>% 
  #              mutate(Single_location = Single_location %>% 
  #                       factor(levels = c("Sweden", "France")),
  #                     Cline = "Off" %>% 
  #                       factor(levels = c("Size", "Color", "Off")),
  #                     Habitat = Habitat %>% 
  #                       factor(levels = c("Sheltered", "Transition", "Exposed"))),
  #            aes(x = LCmeanDist, y = Length, color = Habitat)) +
  # scale_color_manual(name = "Habitat",
  #                 values = c("Sheltered" = "dodgerblue2", "Transition" = "deeppink", "Exposed" = "orange2")) +
  # Parameters of the graph
  facet_grid2(vars(Cline), vars(Single_location),
              scales="free",
              switch = "y",
              labeller = as_labeller(c("Size" = "Shell size (mm)", "Color" = "Frequency of color brown",
                                       "Sweden" = "Sweden", "France" = "France"))) +
  labs(x = "Position along the transect (m)",
       y = "") +
  theme_bw() +
  theme(strip.placement = "outside")


# In french
data2 %>%
  mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
           factor(levels = c("Sweden", "France"))) %>% 
  ggplot(aes(x = LCmeanDist, y = Length)) +
  # The size cline
  geom_point(data = data2 %>%
               mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                        factor(levels = c("Sweden", "France")),
                      Cline = "Size" %>% 
                        factor(levels = c("Size", "Color")),
                      Shell_color_naive_color = ifelse(Shell_color_naive_color == "Yellow", "Jaune", "Marron") %>% 
                        factor(levels = c("Marron", "Jaune"))),
             aes(colour=Shell_color_naive_color), size=2, alpha = 0.7) +
  geom_line(data=Plotting_data_size %>% 
              mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                       factor(levels = c("Sweden", "France")),
                     Cline = "Size" %>% 
                       factor(levels = c("Size", "Color"))),
            aes(x=position, y=phen_cline), color="black", linewidth=2) +
  geom_vline(data = data2 %>%
               mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                        factor(levels = c("Sweden", "France")),
                      Cline = "Size" %>% 
                        factor(levels = c("Size", "Color"))) %>% 
               filter(!is.na(Cline)) %>% 
               select(Centre_size, Cline, Single_location) %>% unique,
             aes(xintercept = Centre_size),
             linetype = "dashed", linewidth = 1) +
  # Color scale for the shell size
  geom_point(data = data2 %>%
               mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                        factor(levels = c("Sweden", "France")),
                      Cline = "Color" %>% 
                        factor(levels = c("Size", "Color")),
                      Shell_color_naive_color = ifelse(Shell_color_naive_color == "Yellow", "Jaune", "Marron") %>% 
                        factor(levels = c("Marron", "Jaune"))),
             aes(y=Genotype_shell_color_naive, color = Shell_color_naive_color), size=2, alpha=0.7) +
  geom_line(data= Plotting_data_color %>% 
              mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>% 
                       factor(levels = c("Sweden", "France")),
                     Cline = "Color" %>% 
                       factor(levels = c("Size", "Color"))),
            aes(x=position, y=phen_cline), color = "black", linewidth=2) +
  geom_vline(data = data2 %>%
               mutate(Single_location = ifelse(Single_location == "Sweden", "Sweden", "France") %>%
                        factor(levels = c("Sweden", "France")),
                      Cline = "Color" %>% 
                        factor(levels = c("Size", "Color"))) %>% 
               filter(!is.na(Cline)) %>% 
               select(Centre_color, Cline, Single_location) %>% unique,
             aes(xintercept = Centre_color),
             linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Couleur\nCoquille",
                     values=c("brown", "orange"),
                     guide = guide_legend(order = 2)) +
  # Parameters of the graph
  facet_grid2(vars(Cline), vars(Single_location),
              scales="free",
              switch = "y",
              labeller = as_labeller(c("Size" = "Taille \n(mm)", "Color" = "Fréquence de couleur\nmarron",
                                       "Sweden" = "Suède", "France" = "France"))) +
  labs(x = "Position le long du transect (m)",
       y = "") +
  theme_bw() +
  theme(text = element_text(size=20),
        strip.placement = "outside")


################################ Confidence intervals for color ################################

Conf_interval_cline_color_fr <- data1_color %>% 
  left_join(Priors_values_cline_color, by = "Single_location") %>% 
  filter(Single_location == "France") %>% 
  group_by(Single_location) %>% 
  bow(tie(Centre_2.5, Width_2.5, Left_2.5, Right_2.5, Centre_97.5, Width_97.5, Left_97.5, Right_97.5) :=
        mle2(clineflog, list(centre = Center %>% unique,
                          width = Width %>% unique %>% log,
                          left = Left %>% unique %>% logit,
                          right = Right %>% unique %>% logit),
             data = list(x = LCmeanDist,
                         g = Genotype_shell_color_naive,
                         n = 1),
             method ="L-BFGS-B",
             upper = list(centre = max(LCmeanDist),
                          width = (max(LCmeanDist) * 2) %>% log,
                          left = 0.999 %>% logit,
                          right = 0.999 %>% logit),
             lower = list(centre = Min_centre_color_Prior %>%
                            unique,
                          width = Min_width_color_Prior %>%
                            unique %>% log,
                          left = 0.001 %>% logit,
                          right = 0.001 %>% logit)) %>%
        confint()) %>%
  as.data.frame() %>% 
  mutate(Centre_2.5i = ifelse(min(Centre_2.5, Centre_97.5, na.rm = TRUE) != Inf, min(Centre_2.5, Centre_97.5, na.rm = TRUE), NA),
         Centre_97.5 = ifelse(min(Centre_2.5, Centre_97.5, na.rm = TRUE) != Inf, max(Centre_2.5, Centre_97.5, na.rm = TRUE), NA),
         Width_2.5i = ifelse(min(Width_2.5, Width_97.5, na.rm = TRUE) != Inf, min(Width_2.5, Width_97.5, na.rm = TRUE) %>% exp, NA),
         Width_97.5 = ifelse(min(Width_2.5, Width_97.5, na.rm = TRUE) != Inf, max(Width_2.5, Width_97.5, na.rm = TRUE) %>% exp, NA),
         Left_2.5i = ifelse(min(Left_2.5, Left_97.5, na.rm = TRUE) != Inf, min(Left_2.5, Left_97.5, na.rm = TRUE) %>% invlogit, NA),
         Left_97.5 = ifelse(min(Left_2.5, Left_97.5, na.rm = TRUE) != Inf, max(Left_2.5, Left_97.5, na.rm = TRUE) %>% invlogit, NA),
         Right_2.5i = ifelse(min(Right_2.5, Right_97.5, na.rm = TRUE) != Inf, min(Right_2.5, Right_97.5, na.rm = TRUE) %>% invlogit, NA),
         Right_97.5 = ifelse(min(Right_2.5, Right_97.5, na.rm = TRUE) != Inf, max(Right_2.5, Right_97.5, na.rm = TRUE) %>% invlogit, NA)) %>%
  select(Single_location, contains(".5i"), contains("97.")) %>% 
  rename(Left_2.5 = Left_2.5i,
         Right_2.5 = Right_2.5i,
         Centre_2.5 = Centre_2.5i,
         Width_2.5 = Width_2.5i)

Conf_interval_cline_color_sw <- data1_color %>% 
  left_join(Priors_values_cline_color, by = "Single_location") %>% 
  filter(Single_location == "Sweden") %>% 
  group_by(Single_location) %>% 
  bow(tie(Right_2.5, Right_97.5) :=
        mle2(stable, list(p_all = p_brown_all %>%
                            unique),
             data = list(x=LCmeanDist,
                         g = Genotype_shell_color_naive,
                         n = 1),
             method = "L-BFGS-B",
             upper = list(p_all = 0.999),
             lower = list(p_all = 0.001)) %>% 
        confint()) %>% 
  as.data.frame %>% 
  mutate(Width_2.5 = NA, Centre_97.5 = NA, Width_97.5 = NA, Left_97.5 = Right_97.5, Width_ci = NA, Left_2.5 = Right_2.5, Centre_2.5 = NA)


color_table <- data1_color %>%
  group_by(Single_location) %>%
  full_join(Conf_interval_cline_color_fr, by = "Single_location") %>%
  full_join(Conf_interval_cline_color_sw, by = c(names(Conf_interval_cline_color_fr))) %>% 
  select(Single_location, contains(".5"), Center, Left, Right, Width, p_brown_all) %>% 
  unique() %>% 
  mutate(Right = ifelse(Single_location == "Sweden", data1_color %>%
                          filter(Single_location == "Sweden") %>%
                          select(p_brown_all) %>% 
                          unique %>% as.vector %>% 
                          unname %>% unlist, Right),
         Left = ifelse(Single_location == "Sweden", Right, Left),
         Right_2.5 = ifelse(Right_2.5 == Right_97.5 & Right_2.5 > Right, NA, Right_2.5),
         Right_97.5 = ifelse(Right_2.5 == Right_97.5 & Right_97.5 < Right, NA, Right_97.5),
         Width_2.5 = ifelse(Width_2.5 == Width_97.5 & Width_2.5 > Width, NA, Width_2.5),
         Width_97.5 = ifelse(Width_2.5 == Width_97.5 & Width_97.5 < Width, NA, Width_97.5)) %>% 
  filter(!is.na(Left_97.5)) %>% 
  mutate(Centre_ci = ifelse(Single_location == "Sweden", "/", paste0("[", Centre_2.5 %>% round(digits=1), ", ", Centre_97.5 %>% round(digits=1), "]")),
         Left_ci = paste0("[", Left_2.5 %>% round(digits = 1), ", ", Left_97.5 %>% round(digits = 1), "]"),
         Right_ci = paste0("[", ifelse(is.na(Right_2.5), Right_97.5, Right_2.5) %>% round(digits = 1), ", ", ifelse(is.na(Right_2.5), Right_2.5, Right_97.5) %>% round(digits = 1), "]"),
         Width_ci = ifelse(Single_location == "Sweden", "/", paste0("[", Width_2.5 %>% round(digits = 1), ", ", Width_97.5 %>% round(digits = 1), "]")),
         Left = ifelse(is.na(Left), NA, Left %>% round(digits = 1)),
         Right = ifelse(is.na(Right), NA, Right %>% round(digits = 1)),
         Centre = ifelse(Single_location == "Sweden", "/", Center %>% round(digits = 1) %>% as.character),
         Width = ifelse(Single_location == "Sweden", "/", Width %>% round(digits = 1) %>% as.character),
         "Sheltered Zone" = paste(Left, Left_ci),
         "Exposed Zone" = paste(Right, Right_ci),
         "Cline Centre (m)" = ifelse(Single_location == "Sweden", "/", paste(Centre, Centre_ci)),
         "Cline Width (m)" = ifelse(Single_location == "Sweden", "/", paste(Width, Width_ci)),
         Country = Single_location %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  ungroup %>% 
  select(Country, "Sheltered Zone", "Exposed Zone", "Cline Centre (m)", "Cline Width (m)") %>% 
  paged_table()

color_table

################################ Confidence intervals for size ################################
Confidence_interval_cline_size <- data1_size %>% 
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
                       n=1),
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


size_table <- data1_size %>%
  group_by(Single_location) %>%
  full_join(Confidence_interval_cline_size, by = "Single_location") %>%
  select(Single_location, contains(".5"), Centre_size, Left_size, Right_size, Width_size) %>% 
  unique() %>% 
  mutate(Centre_ci = paste0("[", min(Centre_2.5, Centre_97.5) %>% round(digits=1), ", ", max(Centre_97.5, Centre_2.5) %>% round(digits=1), "]"),
         Left_ci = paste0("[", min(Left_2.5, Left_97.5) %>% round(digits = 1), ", ", max(Left_97.5, Left_2.5) %>% round(digits = 1), "]"),
         Right_ci = paste0("[", min(Right_2.5, Right_97.5) %>% round(digits = 1), ", ", max(Right_97.5, Right_2.5) %>% round(digits = 1), "]"),
         Width_ci = paste0("[", min(Width_2.5, Width_97.5) %>% round(digits = 1), ", ", max(Width_97.5, Width_2.5) %>% round(digits = 1), "]"),
         Left = Left_size %>% round(digits = 1),
         Right = Right_size %>% round(digits = 1),
         Centre = Centre_size %>% round(digits = 1),
         Width = Width_size %>% round(digits = 1),
         Country = Single_location %>% 
           factor(levels = c("Sweden", "France")),
         "Sheltered Zone" = paste(Left, Left_ci),
         "Exposed Zone" = paste(Right, Right_ci),
         "Cline Centre (m)" = paste(Centre, Centre_ci),
         "Cline Width (m)" = paste(Width, Width_ci)) %>% 
  ungroup %>% 
  select(Country, "Sheltered Zone", "Exposed Zone", "Cline Centre (m)", "Cline Width (m)") %>% 
  paged_table()

size_table


