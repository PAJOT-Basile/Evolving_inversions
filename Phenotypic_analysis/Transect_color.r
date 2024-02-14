# Libraries
require("anyLib")
anyLib(c("tidyverse", "readxl", "forcats", "ggh4x", "bbmle", "tie", "ggforce", "rmarkdown"))

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

# Color palette to be reused everywhere with the shell colors
shell_palette = c("black", "grey28", "brown", "saddlebrown", "wheat4", "gold", "goldenrod", "khaki2", "grey", "blue", "rosybrown4")

# Deal with missing data (drop lines that have less than the calculated threshold of 5% of missing data)
threshold <- ((data %>% nrow) * .05) %>% floor
data <- data %>% drop_na(((data %>% is.na %>% colSums) < threshold) %>% names)

# Shell color and pattern analysis
data %>% 
  ggplot(aes(x=LCmeanDist, y=Length)) +
  geom_point(aes(color=Shell_color)) +
  facet_col(facets = vars(Single_location), scales="free", space="free") +
  scale_color_manual(values=shell_palette) +
  labs(title = "Shell color distribution along the transects",
       x = "Distance along the transect (m)",
       y = "Shell size (mm)") +
  theme_bw() +
  theme(text = element_text(size = 20))

# Modify the colors to a bi-allelic gene
Shell_color <- data$Shell_color


data$Shell_color_naive_color <- str_split_fixed(string=Shell_color, pattern = "/", n=2)[, 1]
data$Shell_color_morphology <- str_split_fixed(string=Shell_color, pattern = "/", n=2)[, 2]

data <- data %>% 
  mutate(Shell_color_naive_color = ifelse(Shell_color_naive_color %in% c("Yellow", "White", "Grey"), "Yellow", Shell_color_naive_color),
         Shell_color_morphology = ifelse(Shell_color_naive_color == "Banded", "Banded", Shell_color_morphology),
         Shell_color_naive_color = ifelse(Shell_color_naive_color %in% c("Black", "Brown", "Dark"), "Brown", Shell_color_naive_color),
         Shell_color_naive_color = ifelse(Shell_color_naive_color == "Banded", "Brown", Shell_color_naive_color),
         Shell_color_naive_color = Shell_color_naive_color %>% factor(levels = c("Yellow", "Brown")),
         Shell_color_morphology = ifelse(! Shell_color_morphology %in% c("Banded", "Square"), "Uniform", "Banded")) %>% 
  mutate(Genotype_shell_color_naive = ifelse(Shell_color_naive_color == "Yellow", 0, 1))

data %>% 
  ggplot(aes(x=LCmeanDist, y=Length)) +
  geom_point(aes(color=Shell_color_naive_color), size=2) +
  facet_col(facets = vars(Single_location), scales="free", space="free") +
  scale_color_manual(name = "Shell color",
                     values=c("orange", "brown")) +
  labs(title = "Distribution of shell color in bi-color coding",
       x="Distance along the transect (m)",
       y = "Shell size (mm)") +
  theme_bw() +
  theme(text = element_text(size = 20))

# Apply the three models to the colors (stable model, linear model and clinal model)
source("Cline_functions.R")

# First, we need the priors. As the values are not the same for each location, we separate them
Priors_values_cline_color <- data.frame(Centre_color = c(200, 70, 70),
                                        Width_color = c(20, 8e-6, 10),
                                        Left_color = c(0.25, 0.1, 0.1),
                                        Right_color = c(0.8, 0.05, 0.05),
                                        Min_centre_color = c(200, 50, 50),
                                        Min_width_color = c(0, 0, 1e-7),
                                        Single_location = c("FRA_LAM_n", "SWE_LOK_n", "SWE_LOK_s"))

# Then, we calculate the initial color frequency to use in the stable model. We also merge with the chosen priors
data <- data %>% 
  group_by(Single_location) %>% 
  summarise(p_yellow_location = mean(Shell_color_naive_color == "Yellow"),
            p_brown_location = 1 - p_yellow_location) %>% 
  merge(data, by="Single_location") %>% 
  merge(Priors_values_cline_color, by="Single_location")

# We run the stable model
Stable_model_coef <- data %>%
  group_by(Single_location) %>%
  summarise(Stable_fit_coef = mle2(stable, list(p_all = p_brown_location %>% unique),
                                   data=list(x=LCmeanDist, g=Genotype_shell_color_naive, n=1),
                                   method="L-BFGS-B",
                                   upper=list(p_all=0.999),
                                   lower=list(p_all=0.001)) %>%
              coef() %>%
              round(digits=3))

# We run the linear model
Linear_model_coefs <- data %>%
  group_by(Single_location) %>%
  bow(tie(Linear_fit_coef_brown, Linear_fit_coef_yellow) := mle2(linear, list(p_left = Left_color %>% unique,
                                                                              p_right = Right_color %>% unique),
                                                                 data=list(x=LCmeanDist, g=Genotype_shell_color_naive, n=1),
                                                                 method="L-BFGS-B",
                                                                 upper=list(p_left=0.999, p_right=0.999),
                                                                 lower=list(p_left=0.001, p_right=0.001)) %>%
        coef() %>%
        round(digits=3))

# We run the cline model 
Clinal_model_color <- data %>%
  group_by(Single_location) %>%
  bow(tie(Center, Width, Left, Right) := mle2(clinef,list(centre=Centre_color %>% unique,
                                                          width=Width_color %>% unique,
                                                          left=Left_color %>% unique,
                                                          right=Right_color %>% unique),
                                              data=list(x=LCmeanDist,
                                                        g=Genotype_shell_color_naive,
                                                        n=1),
                                              method="L-BFGS-B",
                                              upper=list(centre=max(LCmeanDist),
                                                         width=max(LCmeanDist)*2,
                                                         left=0.999,
                                                         right=0.999),
                                              lower=list(centre=Min_centre_color %>% unique,
                                                         width=Min_width_color %>% unique,
                                                         left=0.001, right=0.001)) %>%
        coef() %>%
        round(digits=3))
  
# Merge the three dataframes and calculate the AIC for each model
data1 <- data %>% 
  merge(Stable_model_coef, by="Single_location") %>% 
  merge(Linear_model_coefs, by="Single_location") %>% 
  merge(Clinal_model_color, by="Single_location") %>% 
  group_by(Single_location) %>% 
  mutate(AIC_stable_color = mle2(stable, list(p_all = p_brown_location %>% unique),
                           data=list(x=LCmeanDist, g=Genotype_shell_color_naive, n=1),
                           method="L-BFGS-B",
                           upper=list(p_all=0.999),
                           lower=list(p_all=0.001)) %>% AIC,
         AIC_linear_color = mle2(linear, list(p_left = Left_color %>% unique,
                                           p_right = Right_color %>% unique),
                              data=list(x=LCmeanDist, g=Genotype_shell_color_naive, n=1),
                              method="L-BFGS-B",
                              upper=list(p_left=0.999, p_right=0.999),
                              lower=list(p_left=0.001, p_right=0.001)) %>% AIC,
         AIC_cline_color = mle2(clinef,list(centre=Centre_color %>% unique,
                                            width=Width_color %>% unique,
                                            left=Left_color %>% unique,
                                            right=Right_color %>% unique),
                                data=list(x=LCmeanDist,
                                          g=Genotype_shell_color_naive,
                                          n=1),
                                method="L-BFGS-B",
                                upper=list(centre=max(LCmeanDist),width=max(LCmeanDist)*2,
                                           left=0.999,right=0.999),
                                lower=list(centre=Min_centre_color %>% unique,width=Min_width_color %>% unique,left=0.001,right=0.001)) %>% AIC) %>% 
  ungroup()

# Getting confidence intervals for the colors
Confidence_interval_cline_size <- data1 %>% 
  merge(Priors_values_cline_pheno, by="Single_location") %>% 
  group_by(Single_location) %>% 
  bow(tie(Centre_2.5, Width_2.5, Left_2.5, Right_2.5, Centre_97.5, Width_97.5, Left_97.5, Right_97.5) :=
        mle2(clinef, list(centre=Centre_color %>% unique,
                          width =Width_color %>% unique,
                          left=Left_color %>% unique,
                          right=Right_color %>% unique),
             data=list(x=LCmeanDist,
                       g=Genotype_shell_color_naive,
                       n=1),
             method="L-BFGS-B",
             upper=list(centre=max(LCmeanDist), width=max(LCmeanDist)*2, left=0.999, right=0.999),
             lower=list(centre=Min_centre_color %>% unique, width=Min_width_color %>% unique, left=0.001, right=0.001)) %>%
        confint()) %>% as.data.frame() %>% 
  group_by(Single_location) %>% 
  mutate(Width_ci = ifelse(max(Width_2.5, Width_97.5, na.rm=TRUE) == -Inf, NA, max(Width_2.5, Width_97.5, na.rm=TRUE)),
         Left_2.5i = min(Left_2.5, Left_97.5, na.rm=TRUE),
         Right_2.5i = min(Right_2.5, Right_97.5, na.rm=TRUE),
         Centre_2.5i = ifelse(Left_2.5i > Right_2.5i, ifelse(min(Centre_2.5, Centre_97.5, na.rm=TRUE) == Inf, NA, min(Centre_2.5, Centre_97.5, na.rm=TRUE)), max(Centre_2.5, Centre_97.5, na.rm=TRUE)),
         Left_97.5 = max(Left_2.5, Left_97.5, na.rm=TRUE),
         Right_97.5 = max(Right_2.5, Right_97.5, na.rm=TRUE),
         Centre_97.5 = ifelse(Left_2.5i > Right_2.5i, ifelse(max(Centre_2.5, Centre_97.5, na.rm=TRUE) == -Inf, NA, max(Centre_2.5, Centre_97.5, na.rm=TRUE)), min(Centre_2.5, Centre_97.5, na.rm=TRUE))) %>% 
  select(-c(Left_2.5, Right_2.5, Centre_2.5)) %>% 
  rename(Left_2.5 = Left_2.5i,
         Right_2.5 = Right_2.5i,
         Centre_2.5 = Centre_2.5i)

# Prepare the plotting of results for the clinal model
LAMn <- data1 %>% filter(Single_location == "FRA_LAM_n")
LOKn <- data1 %>% filter(Single_location == "SWE_LOK_n")
LOKs <- data1 %>% filter(Single_location == "SWE_LOK_s")
dfs <- list(LAMn, LOKn, LOKs)
cline_fit_plots <- list()
for (i in enumerate(dfs)){
  df <- i[2] %>% as.data.frame()
  cline_fit_plots <- append(
    cline_fit_plots, list(
      clinef(x=seq(min(df$LCmeanDist), max(df$LCmeanDist), by=1),
             centre=df$Center %>% unique(),
             width=df$Width %>% unique(),
             left=df$Left %>% unique(),
             right=df$Right %>% unique(),
             plotting=TRUE,
             optimisation=FALSE)
    )
  )
}

# Get the maximum lengths of each group 
Max_lengths_per_location <- data1 %>% 
  group_by(Single_location) %>% 
  summarise(Max_Length = max(Length))
LAMn_plot <- cline_fit_plots[1] %>% 
  as.data.frame() %>% 
  mutate(Single_location = "FRA_LAM_n",
         Max_Length = Max_lengths_per_location$Max_Length[Max_lengths_per_location$Single_location == "FRA_LAM_n"])
LOKn_plot <- cline_fit_plots[2] %>% 
  as.data.frame() %>% 
  mutate(Single_location = "SWE_LOK_n",
         Max_Length = Max_lengths_per_location$Max_Length[Max_lengths_per_location$Single_location == "SWE_LOK_n"])
LOKs_plot <- cline_fit_plots[3] %>% 
  as.data.frame() %>% 
  mutate(Single_location = "SWE_LOK_s",
         Max_Length = Max_lengths_per_location$Max_Length[Max_lengths_per_location$Single_location == "SWE_LOK_s"])

Plotting_data <- rbind(LAMn_plot, LOKn_plot, LOKs_plot) %>%
  # We add a correction for a dominance effect on the brown allele
  mutate(pheno_cline_p = sqrt(1 - phen_cline),
    pheno_cline_q = 1 - pheno_cline_p,
    pheno_cline_homo_brown = pheno_cline_q **2,
    pheno_cline_homo_yellow = pheno_cline_p **2,
    pheno_cline_hetero = 2 * pheno_cline_p * pheno_cline_q)



data1 %>% 
  ggplot(aes(x=LCmenaDist)) +
  geom_point(data = Plotting_data, aes(x=position, y=sqrt(phen_cline)), color="purple") +
  geom_point(data=Plotting_data, aes(x=position, y=phen_cline), color="brown") +
  facet_col(facets = vars(Single_location), scales="free_x", space="free")
  

# Plot the results
data1 %>% 
  ggplot(aes(x=LCmeanDist)) +
  geom_point(aes(y=Genotype_shell_color_naive)) +
  geom_point(data= Plotting_data, aes(x=position, y=phen_cline, color = "Observed brown frequency")) +
  geom_point(data=Plotting_data, aes(x=position, y=(phen_cline + pheno_cline_hetero/2)/max(phen_cline + pheno_cline_hetero/2), color = "Corrected brown allele \nfrequency for a yellow \ndominant allele"))+
  facet_col(facets = vars(Single_location), scales="free", space="free") +
  labs(title = "Variation in color frequency along the transect",
       x = "Position along the transect (m)",
       y = "Frequency of color brown") +
  theme_bw() +
  scale_color_manual(name = "Curve",
                     breaks = c("Observed brown frequency", "Corrected brown allele \nfrequency for a yellow \ndominant allele"),
                     values = c("Observed brown frequency" = "brown", "Corrected brown allele \nfrequency for a yellow \ndominant allele" = "navyblue")) +
  theme(text = element_text(size=20))

#TODO: Find the confidence intervals of the clines

# Recapitulation tables
AIC_table_color <- data1 %>% 
  group_by(Single_location) %>% 
  summarise(AIC_stable = AIC_stable_color %>% unique %>% round(digits=1),
            AIC_linear = AIC_linear_color %>% unique %>% round(digits=1),
            AIC_cline = AIC_cline_color %>% unique %>% round(digits=1),
            Delta_AIC_Stable = (AIC_stable_color - AIC_cline_color) %>% unique %>% round(digits=1),
            Delta_AIC_linear = (AIC_linear_color - AIC_cline_color) %>% unique %>% round(digits=1)) %>% 
  paged_table()

AIC_table_color


color_table <- data1 %>%
  group_by(Single_location) %>%
  mutate(Delta_AIC_stable = AIC_stable_color - AIC_cline_color,
         Delta_AIC_linear = AIC_linear_color - AIC_cline_color) %>%
  summarise(Centre_color = Center %>% unique %>% round(digits=1),
            Left_color = Left %>% unique %>% round(digits=1),
            Right_color = Right %>% unique %>% round(digits=1),
            Width_color = Width %>% unique %>% round(digits=1),
            Delta_AIC_stable = Delta_AIC_stable %>% unique %>% round(digits=1),
            Delta_AIC_linear = Delta_AIC_linear %>% unique %>% round(digits=1))

color_table
