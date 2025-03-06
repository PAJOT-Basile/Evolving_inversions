############################ Libraries ##########################
# Import libraries
libraries <- c("tidyverse", "adegenet", "vcfR")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)

############################ Useful functions ##########################
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

# Import the regions where there are inversions
Confirmed_inversion <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
                                  sep = "\t", header = TRUE) %>% 
  group_by(Chromosome, Inversion) %>% 
  summarize(Start = min(Start, na.rm = TRUE),
            End = max(End, na.rm = TRUE),
            Length = End - Start) %>% 
  ungroup %>% 
  inner_join(read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
                        sep = "\t", header = TRUE) %>% 
               filter(Population == "France") %>% 
               select(Chromosome, Inversion, Inversion_grouped),
             by = c("Chromosome", "Inversion"))

# Import the vcf file
data <- read.vcfR("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/outputs/02_Filter_VCF_Files/11_LAM_and_LOK/VCF_File_Hobs.vcf.gz") %>% 
  vcfR2genind()

data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor


############### Run local PCA and clusterisation on candidate inversion localisations ###############
local_pca_inversions <- (Confirmed_inversion %>% 
                           filter(Inversion_grouped %!in% c("Inv_3.1", "Inv_4.1", "Inv_14.1")) %>% 
  run_loca_pca_inversions())$PC_scores %>% 
  rbind((Confirmed_inversion %>% 
           filter(Inversion_grouped %in% c("Inv_3.1", "Inv_4.1", "Inv_14.1")) %>% 
           mutate(Inversion = Inversion_grouped) %>% 
           run_loca_pca_inversions())$PC_scores)
  
# Represent all the potential inversions for PCA
local_pca_inversions %>% 
  mutate(Habitat = str_split_fixed(Sample_Name, "_", 6)[, 5],
         Habitat = ifelse(Habitat == "EXPOS", "Exposed", ifelse(Habitat == "SHELT", "Sheltered", "Transition")) %>% 
           factor(levels = c("Exposed", "Transition", "Sheltered")),
         Inversion = Inversion %>% 
           factor(levels = local_pca_inversions %>% 
                    select(Inversion) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist),
         Group = Group %>% factor(levels = c("1", "2", "3")),
         Population = ifelse(grepl("LOK", Sample_Name), "Sweden", "France") %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis2, colour = Habitat), size = 3, alpha = 0.7) +
  facet_wrap(vars(Inversion), scales = "free", ncol = 4) +
  scale_color_manual(values = c("Exposed" = "orange2", "Transition" = "deeppink", "Sheltered" = "dodgerblue2")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis1, y = Axis2, colour = Population), lwd = 1.2) +
  scale_color_manual(values = c("Sweden" = "darkorchid4", "France" = "forestgreen")) +
  my_theme

local_pca_inversions %>% 
  mutate(Habitat = str_split_fixed(Sample_Name, "_", 6)[, 5],
         Habitat = ifelse(Habitat == "EXPOS", "Exposed", ifelse(Habitat == "SHELT", "Sheltered", "Transition")) %>% 
           factor(levels = c("Exposed", "Transition", "Sheltered")),
         Inversion = Inversion %>% 
           factor(levels = local_pca_inversions %>% 
                    select(Inversion) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist),
         Group = Group %>% factor(levels = c("1", "2", "3")),
         Population = ifelse(grepl("LOK", Sample_Name), "Sweden", "France") %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis3, colour = Habitat), size = 3, alpha = 0.7) +
  facet_wrap(vars(Inversion), scales = "free", ncol = 4) +
  scale_color_manual(values = c("Exposed" = "orange2", "Transition" = "deeppink", "Sheltered" = "dodgerblue2")) +
  new_scale_color() +
  geom_mark_ellipse(aes(x = Axis1, y = Axis3, colour = Population), lwd = 1.2) +
  scale_color_manual(values = c("Sweden" = "darkorchid4", "France" = "forestgreen")) +
  my_theme


# Save the output
local_pca_inversions %>%
  mutate(Habitat = str_split_fixed(Sample_Name, "_", 6)[, 5],
         Habitat = ifelse(Habitat == "EXPOS", "Exposed", ifelse(Habitat == "SHELT", "Sheltered", "Transition")) %>% 
           factor(levels = c("Exposed", "Transition", "Sheltered")),
         Inversion = Inversion %>% 
           factor(levels = local_pca_inversions %>% 
                    select(Inversion) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    pull(Inversion)),
         Group = Group %>% factor(levels = c("1", "2", "3", "4", "5", "6")),
         Population = ifelse(grepl("LOK", Sample_Name), "Sweden", "France") %>% 
           factor(levels = c("Sweden", "France"))) %>%
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Pca_Together.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
