# Libraries
libraries <- c("tidyverse", "readxl", "rmarkdown")
if (!require("pacman")) install.packages("pacman")
for (lib in libraries){
  pacman::p_load(lib, character.only = TRUE)
}
################################ Useful functions ################################
source("../General_scripts/Functions_optimise_plot_clines.r")

################## Import the metadata  ##################
metadata <- read_excel(path = "../Input_Data/Data/data_Fabalis_resequencing_Basile.xlsx",
                       sheet = 1,
                       col_names = TRUE,
                       trim_ws = TRUE) %>%
  
  # Convert to the correct formats
  mutate(Population = as.factor(Population),
         Shell_colour = factor(Shell.colour %>% str_to_title, levels = c("Black", "Black/Square", "Brown", "Brown/Square", "Dark", "Yellow", "Yellow/Brown", "Yellow/Square", "Grey", "White", "Banded", NA)),
         LCmeanDist = as.numeric(LCmeanDist),
         Length = as.numeric(length),
         Habitat = ifelse(Habitat %in% c("EXPOS"), "Exposed", Habitat),
         Habitat = ifelse(Habitat %in% c("HARB", "SHELT"), "Sheltered", Habitat),
         Habitat = ifelse(Habitat %in% c("TRANS", "TRANSI"), "Transition", Habitat),
         Habitat = as.factor(Habitat)
  ) %>%
  
  # Select only the necessary columns for the analysis
  select(-length) %>% 

  # Modify the population column to get only the name of the country
  mutate(Single_location = ifelse(Population == "LOK", "Sweden", "France") %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  # Drop unused levels
  droplevels

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
  mutate(Genotype_shell_color_naive = ifelse(Shell_color_naive_color == "Yellow", 0, 1))

################################ Fit the three models to the data ################################
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
  merge(metadata, by="Single_location") %>% 
  merge(Priors_values_cline_color, by="Single_location")

# We run the cline model 
Clinal_model_color <- metadata %>%
  group_by(Single_location) %>%
  tie::bow(tie(Center, Width, Left, Right) := mle2(clinef,list(centre = Centre_color_Prior %>% unique,
                                                          width = Width_color_Prior %>% unique,
                                                          left = Left_color_Prior %>% unique,
                                                          right = Right_color_Prior %>% unique),
                                              data = list(x = LCmeanDist,
                                                        g = Genotype_shell_color_naive,
                                                        n = 1),
                                              method = "L-BFGS-B",
                                              upper = list(centre = max(LCmeanDist),
                                                         width=max(LCmeanDist)*2,
                                                         left =  0.999,
                                                         right = 0.999),
                                              lower = list(centre = Min_centre_color_Prior %>% unique,
                                                         width = Min_width_color_Prior %>% unique,
                                                         left = 0.001, right = 0.001)) %>%
        coef() %>%
        round(digits = 3))

# We run the stable model
Stable_model_color <- metadata %>% 
  group_by(Single_location) %>% 
  summarise(p_brown_all = mle2(stable, list(p_all = p_brown_location %>% unique),
                               data = list(x = LCmeanDist,
                                           g = Genotype_shell_color_naive,
                                           n = 1),
                               method = "L-BFGS-B",
                               upper = list(p_all = 0.999),
                               lower = list(p_all = 0.001)) %>% 
              coef() %>% round(digits = 3))

# We run the linear model
Linear_model_color <- metadata %>%
  group_by(Single_location) %>%
  tie::bow(tie(p_left, p_right) := mle2(linear, list(p_left = Left_color_Prior %>% unique,
                                               p_right = Right_color_Prior %>% unique),
                                              data=list(x = LCmeanDist,
                                                        g = Genotype_shell_color_naive,
                                                        n = 1),
                                              method = "L-BFGS-B",
                                              upper = list(p_left = 0.999,
                                                         p_right = 0.999),
                                              lower = list(p_left = 0.001,
                                                         p_right = 0.001)) %>%
        coef() %>%
        round(digits=3))

# Merge the three dataframes and calculate the AIC for each model
data_AIC <- metadata %>% 
  left_join(Linear_model_color, by = "Single_location") %>% 
  left_join(Clinal_model_color, by = "Single_location") %>% 
  left_join(Stable_model_color, by = "Single_location") %>% 
  group_by(Single_location) %>% 
  mutate(AIC_stable = mle2(stable, list(p_all = p_brown_all %>% unique),
                           data = list(x = LCmeanDist,
                                       g = Genotype_shell_color_naive,
                                       n = 1),
                           method = "L-BFGS-B",
                           upper = list(p_all = 0.999999),
                           lower = list(p_all = 0.000001)) %>% 
           AIC,
         AIC_linear = mle2(linear, list(p_left = p_left %>% unique,
                                        p_right = p_right %>% unique),
                           data = list(x = LCmeanDist,
                                       g = Genotype_shell_color_naive,
                                       n = 1),
                           method = "L-BFGS-B",
                           upper = list(p_left = 0.999999,
                                        p_right = 0.999999),
                           lower = list(p_left = 0.000001,
                                        p_right = 0.00001)) %>% 
           AIC,
         AIC_clinal = mle2(clinef, list(centre = Center %>% unique,
                                        width = Width %>% unique,
                                        left = Left %>% unique,
                                        right = Right %>% unique),
                           data = list(x = LCmeanDist,
                                       g = Genotype_shell_color_naive,
                                       n = 1),
                           method = "L-BFGS-B",
                           upper = list(centre = max(LCmeanDist),
                                        width=max(LCmeanDist)*2,
                                        left =  0.999,
                                        right = 0.999),
                           lower = list(centre = Min_centre_color_Prior %>% unique,
                                        width = Min_width_color_Prior %>% unique,
                                        left = 0.001, right = 0.001)) %>% 
           AIC) %>% 
  select(Single_location, starts_with("AIC")) %>% 
  unique %>%
  # Here, we estimate which model is the best one using the differences with the
  # clinal model AIC.
  mutate(Delta_AIC_stable = AIC_clinal - AIC_stable, Delta_AIC_linear = AIC_clinal - AIC_linear,
                    Best_model = ifelse((Delta_AIC_stable > Delta_AIC_linear) & (Delta_AIC_stable > 0), "Stable",
                                        ifelse((Delta_AIC_stable > Delta_AIC_linear) & (Delta_AIC_stable < 0), "Clinal",
                                               ifelse((Delta_AIC_stable < Delta_AIC_linear) & (Delta_AIC_linear > 0), "Linear",
                                                      ifelse((Delta_AIC_stable < Delta_AIC_linear) & (Delta_AIC_linear < 0), "Clinal", NA)))),
         # Once the best model is selected, we select the according AIC
         Best_model_AIC = ifelse(Best_model == "Stable", AIC_stable, ifelse(Best_model == "Linear", AIC_linear, AIC_clinal)),
         # Then, we estimate which is the second best model using the AICs
         Second_best_model = ifelse(Best_model == "Stable" & Delta_AIC_linear > 0, "Linear",
                                    ifelse(Best_model == "Stable" & Delta_AIC_linear < 0, "Clinal",
                                           ifelse(Best_model == "Linear" & Delta_AIC_stable < 0, "Clinal",
                                                  ifelse(Best_model == "Linear" & Delta_AIC_stable > 0, "Stable",
                                                         ifelse(Best_model == "Clinal" & AIC_linear > AIC_stable, "Stable", "Linear"))))),
         # And we calculate the delta AIC between the two best models
         Delta_AIC_second_best_model = ifelse(Best_model == "Clinal" & Second_best_model == "Linear", AIC_clinal - AIC_linear,
                                              ifelse(Best_model == "Clinal" & Second_best_model == "Stable", AIC_clinal - AIC_stable,
                                                     ifelse(Best_model == "Linear" & Second_best_model == "Clinal", AIC_linear - AIC_clinal,
                                                            ifelse(Best_model == "Linear" & Second_best_model == "Stable", AIC_linear - AIC_stable,
                                                                   ifelse(Best_model == "Stable" & Second_best_model == "Clinal", AIC_stable - AIC_clinal, AIC_stable - AIC_linear)))))) %>% 
  select(-c(starts_with("AIC"), Delta_AIC_stable, Delta_AIC_linear))

# Now that we chose the best model for the two locations, we can use these to do the analysis
data1_color <- metadata %>% 
  left_join(Clinal_model_color %>% 
              filter(Single_location == "France"), by="Single_location") %>% 
  left_join(Stable_model_color %>% 
              filter(Single_location == "Sweden"), by = "Single_location")

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


################################ Confidence intervals for color ################################

Conf_interval_cline_color_fr <- data1_color %>% 
  filter(Single_location == "France") %>% 
  group_by(Single_location) %>% 
  tie::bow(tie(Centre_2.5, Width_2.5, Left_2.5, Right_2.5, Centre_97.5, Width_97.5, Left_97.5, Right_97.5) :=
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
  tie::bow(tie(Right_2.5, Right_97.5) :=
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
  mutate(Centre_ci = ifelse(Single_location == "Sweden", "/", paste0("[", Centre_2.5 %>% round(digits=3), ", ", Centre_97.5 %>% round(digits=3), "]")),
         Left_ci = paste0("[", Left_2.5 %>% round(digits = 3), ", ", Left_97.5 %>% round(digits = 3), "]"),
         Right_ci = paste0("[", ifelse(is.na(Right_2.5), Right_97.5, Right_2.5) %>% round(digits = 3), ", ", ifelse(is.na(Right_2.5), Right_2.5, Right_97.5) %>% round(digits = 3), "]"),
         Width_ci = ifelse(Single_location == "Sweden", "/", paste0("[", Width_2.5 %>% round(digits = 3), ", ", Width_97.5 %>% round(digits = 3), "]")),
         Left = ifelse(is.na(Left), NA, Left %>% round(digits = 3)),
         Right = ifelse(is.na(Right), NA, Right %>% round(digits = 3)),
         Centre = ifelse(Single_location == "Sweden", "/", Center %>% round(digits = 3) %>% as.character),
         Width = ifelse(Single_location == "Sweden", "/", Width %>% round(digits = 3) %>% as.character),
         "Zone abritée" = paste(Left, Left_ci),
         "Zone exposée" = paste(Right, Right_ci),
         "Centre du cline" = ifelse(Single_location == "Sweden", "/", paste(Centre, Centre_ci)),
         "Largeur du cline" = ifelse(Single_location == "Sweden", "/", paste(Width, Width_ci)),
         Pays = ifelse(Single_location == "Sweden", "Suède", "France") %>% 
           factor(levels = c("Suède", "France"))) %>% 
  ungroup %>% 
  select(Pays, "Zone abritée", "Zone exposée", "Centre du cline", "Largeur du cline") %>% 
  paged_table()

color_table

# Here, we also add the AIC table
data_AIC %>% paged_table()
print(data_AIC %>% as.data.frame)
