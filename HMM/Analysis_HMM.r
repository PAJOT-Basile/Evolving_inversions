############################ Libraries ##########################
require("tidyverse")

############################ Useful functions ##########################

source("../General_scripts/Functions_optimisation_visualisation.r")

############################ Load inputs to HMM ##########################
input_delta_freqs_france <- read.csv("./Whole_genome_france_freqs.csv",
                                     header = TRUE,
                                     sep = "\t") %>% 
  rownames_to_column("Position") %>% 
  rownames_to_column("correspondance")

input_delta_freqs_sweden <- read.csv("./Whole_genome_sweden_freqs.csv",
                                     header = TRUE,
                                     sep = "\t") %>% 
  rownames_to_column("Position") %>% 
  rownames_to_column("correspondance")
print("Finished importing input")
############################ Load HMM output ##########################
output_delta_freqs_france <- read.table("./Whole_genome_france_freqs_3state_HMMstates.txt",
                                        header = TRUE,
                                        sep = " ") %>% 
  rownames_to_column("correspondance") %>% 
  left_join(input_delta_freqs_france, by = "correspondance") %>% 
  select(-correspondance) %>% 
  rename(Group = x.x,
         Delta_freq = x.y) %>% 
  relocate(Position)

output_delta_freqs_sweden <- read.table("./Whole_genome_sweden_freqs_3state_HMMstates.txt",
                                        header = TRUE,
                                        sep = " ") %>% 
  rownames_to_column("correspondance") %>% 
  left_join(input_delta_freqs_sweden, by = "correspondance") %>% 
  select(-correspondance) %>% 
  rename(Group = x.x,
         Delta_freq = x.y) %>% 
  relocate(Position) %>% 
## WARNING ! The groups (1, 2, and 3) do not always correspond to the same values in the two populations 
# for example, 1 can be the group of the very differenciated SNPs in the swedish population, but the
# group of low differenciation in sweden. Therefore, here, we correct for this effect, but this line might need to be 
# changed according to the outputs that are given.
  mutate(Group = 4 - Group)

rm(list = apropos("input_delta_freqs"))
print("Finished importing output")

############################ Plot the outputs ##########################
output_delta_freqs_france %>% 
  mutate(Population = "France") %>% 
  rbind(output_delta_freqs_sweden %>% 
          mutate(Population = "Suède")) %>% 
  mutate(Population = Population %>% factor(levels = c("Suède", "France"))) %>% 
  filter(!grepl("unloc|unplaced", Position)) %>% 
  geom_manhattan(aes(y = Delta_freq, color = Group, facetting = Population), size = 2, alpha = 0.2, palette = c("grey71", "turquoise4")) +
  facet_col(vars(Population)) +
  theme(text = element_text(size = 30))


############################ Separate the inversions from the rest ##########################
# First, we transform the data to be able to use the positions in the genome as numeric values
output_delta_freqs_sweden_modif <- output_delta_freqs_sweden %>% 
  mutate(pos = Position) %>% 
  separate(pos, c("LG", "SUPER", "SUPER_frag", "pos"), "_") %>%
  separate(pos, c("pos", "allele"), "\\.") %>% 
  select(-allele) %>% 
  mutate(pos = pos %>% as.numeric) %>% 
  filter(!grepl("unloc|unplaced", Position))

# Then, we select the most common group (1, 2, or 3) on a fragment of the genome (the bin).
bin_size <- 1e5
# Then, we try indentifying the inversions. The numbers that have been chosen are taken by 
# examining the group values along each chromosome.
Position_inversions <- output_delta_freqs_sweden_modif %>% 
  mutate(bin = (pos / bin_size) %>% round) %>% 
  group_by(LG, bin, Group) %>% 
  summarize(counting = n()) %>% 
  filter(counting == max(counting)) %>% 
  select(-counting) %>% 
  ungroup %>% 
  rename(Chromosome = LG,
         Position = bin) %>% 
  filter(Group == 3) %>%
  group_by(Chromosome) %>% 
  mutate(pos_prev = lag(Position, default = 0),
         difference_pos_max = Position - pos_prev,
         Chromosome = Chromosome %>% factor(levels = paste0("LG", seq(1, 17, 1)))) %>% 
  ungroup %>% 
  mutate(Inversion = NA,
         Inversion = ifelse(Chromosome == "LG3" & Position != 76, "Inv_3.1", Inversion),
         Inversion = ifelse(Chromosome == "LG4" & Position %in% c(50:149), "Inv_4.1", Inversion),
         Inversion = ifelse(Chromosome == "LG4" & Position %in% c(341:467), "Inv_4.2", Inversion),
         Inversion = ifelse(Chromosome == "LG4" & Position %in% c(550:702), "Inv_4.3", Inversion),
         Inversion = ifelse(Chromosome == "LG6" & Position %in% c(0:248), "Inv_6.1", Inversion),
         Inversion = ifelse(Chromosome == "LG6" & Position %in% c(427:552), "Inv_6.2", Inversion),
         Inversion = ifelse(Chromosome == "LG7" & Position %in% c(25:98), "Inv_7.1", Inversion),
         Inversion = ifelse(Chromosome == "LG7" & Position %in% c(478:588), "Inv_7.2", Inversion),
         Inversion = ifelse(Chromosome == "LG8" & Position %in% c(484:717), "Inv_8.1", Inversion),
         Inversion = ifelse(Chromosome == "LG11" & Position %in% c(0:107), "Inv_11.1", Inversion),
         Inversion = ifelse(Chromosome == "LG13" & Position %in% c(249:449), "Inv_13.1", Inversion),
         Inversion = ifelse(Chromosome == "LG14" & Position %in% c(0:236), "Inv_14.1", Inversion),
         Inversion = ifelse(Chromosome == "LG14" & Position %in% c(359:450), "Inv_14.2", Inversion),
         Inversion = ifelse(Chromosome == "LG16" & Position %in% c(180:348), "Inv_16.1", Inversion)) %>% 
  filter(!is.na(Inversion)) %>% 
  select(Chromosome, Position, Inversion)

# We make a table that we then export
Position_inversions %>% 
  group_by(Chromosome, Inversion) %>% 
  mutate(Position = Position * bin_size,
         Inversion = Inversion %>% 
           factor(levels = Position_inversions %>%
                    select(Inversion) %>% 
                    unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Inversion))) %>% 
                    as.vector %>% unname %>% unlist)) %>% 
  summarise(Pos_min_inv = min(Position),
            Pos_max_inv = max(Position)) %>%
  write.table("C:/Documents/Stage_M2_2024/Project/Genetic data/Delimitation_inversions.csv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
