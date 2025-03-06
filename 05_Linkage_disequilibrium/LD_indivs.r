############### Libraries ###############
libraries <- c("tidyverse", "ggforce")
if(!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)

############### Useful variables ###############
colfunc <- colorRampPalette(colors = c("#813d32","#ef852f","#fdc441" ,"#e3cfb4") %>% rev)

############### Useful functions ###############
source("/shared/projects/pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Functions_optimise_plot_clines.r")

############### Import inversion delimitations ###############
delimitation_inversions <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Candidate_inversions_post_stats.tsv",
                                      header = TRUE)
############### Import local pca output ###############
pca_inversions <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Inversions_post_pca.tsv",
                             header = TRUE) 

############### Compute LD per population ###############
LD_between_inversions_Sweden <- pca_inversions %>% 
  filter(Population == "Sweden") %>% 
  inner_join(delimitation_inversions, by = c("Population", "Chromosome", "Inversion")) %>% 
  # And compute the inter-inversion LD
  Compute_LD_inversions(LD_version = "inter_chromosomic")

# Filter to only keep the inversions with the three genotypes
LD_between_inversions_France <- pca_inversions %>% 
  filter(Population == "France") %>% 
  inner_join(delimitation_inversions, by = c("Population", "Chromosome", "Inversion")) %>% 
  # And compute the inter-inversion LD
  Compute_LD_inversions(LD_version = "inter_chromosomic")

# Prepare output to save and to plot
LD_inversions_indivs <- LD_between_inversions_France %>% 
  mutate(Population = "France") %>% 
  rbind(LD_between_inversions_Sweden %>% 
          mutate(Population = "Sweden"))
LD_inversions_indivs <- LD_inversions_indivs %>% 
  mutate(Population = Population %>% 
           factor(levels = c("Sweden", "France")),
         Inv_1 = Inv_1 %>% 
           factor(levels = LD_inversions_indivs %>% 
                    select(Inv_1) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inv_1))) %>% 
                    as.vector %>% unname %>% unlist),
         Inv_2 = Inv_2 %>% 
           factor(levels = LD_inversions_indivs %>% 
                    select(Inv_2) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inv_2))) %>% 
                    as.vector %>% unname %>% unlist)) %>% 
  filter(paste0(Inv_1, Inv_2) != paste0(Inv_2, Inv_1))


LD_inversions_indivs %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/LD_inversions.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



LD_inversions_indivs %>% 
  ggplot() +
  geom_tile(aes(x = Inv_2, y = Inv_1, fill = Perc_indivs_same_cluster)) +
  scale_fill_gradientn(name = "% of individuals\nin common",
                       colors = colfunc(20)) +
  facet_row(facets = vars(Population), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(size = 20)) +
  labs(x = "Inversion",
       y = "Inversion")

############### Regroup inversions LD 1 ###############
delimitation_inversions %>% 
  mutate(Inversion_grouped = case_when(
    Inversion == "Inv_3.2" ~ "Inv_3.1",
    Inversion == "Inv_4.2" ~ "Inv_4.1",
    Inversion == "Inv_14.2" ~ "Inv_14.1",
    TRUE ~ Inversion
  )) %>% 
  # And save the output as the delimitation of the chromosomal inversions
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Delimitation_inversions.tsv",
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
