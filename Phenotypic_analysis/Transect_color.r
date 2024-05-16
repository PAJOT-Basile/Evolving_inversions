# Libraries
require("anyLib")
anyLib(c("tidyverse", "readxl", "ggh4x", "bbmle", "tie", "ggforce", "ggpubr", "ggnewscale", "LaplacesDemon"))

################################ Useful functions ################################
"%!in%" <- function(x, y){!("%in%"(x, y))}

zip <- function(...){
  mapply(list, ..., SIMPLIFY=FALSE)
}

enumerate <- function(...){
  zip(ix=seq_along(..1), ...)
}

# Apply the three models to the colors (stable model, linear model and clinal model)
source("/shared/projects/pacobar/finalresult/bpajot/genomic_analysis/scripts/01_Filtering_stats_vcf/Functions_optimise_plot_clines.r")

################################ Useful variables ################################
# Color palette to be reused everywhere with the shell colors
shell_palette = c("black", "grey28", "brown", "saddlebrown", "wheat4", "gold", "goldenrod", "khaki2", "grey", "blue", "rosybrown4")

################################ Import metadata ################################
metadata <- read_excel(path = "/shared/projects/pacobar/finalresult/bpajot/Data/data_Fabalis_resequencing_Basile.xlsx",
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
threshold <- ((metadata %>% nrow) * .05) %>% floor
metadata <- metadata %>% drop_na(((metadata %>% is.na %>% colSums) < threshold) %>% names)

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
  mutate(Genotype_shell_color_naive = ifelse(Shell_color_naive_color == "Yellow", 0, 1))

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
                          plotting = TRUE,
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


