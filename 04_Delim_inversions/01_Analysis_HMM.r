############################ Libraries ##########################
# Import libraries
libraries <- c("tidyverse")

if (!require("pacman")) install.packages("pacman")
for (lib in libraries){
  pacman::p_load(lib, character.only = TRUE)
}
rm(libraries, lib)


############################ Useful functions ##########################
source("../General_scripts/Functions_optimise_plot_clines.r")

############################ Useful variables ##########################
bin_size <- 1e5


############################ Load inputs to HMM ##########################
input_delta_freqs_france <- read.csv("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/03_HMM/Data/France.weir.fst",
                                     header = TRUE,
                                     sep = "\t") %>% 
  rownames_to_column("correspondance") %>% 
  unite(Position, CHROM, POS)

input_delta_freqs_sweden <- read.csv("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/03_HMM/Data/Sweden.weir.fst",
                                     header = TRUE,
                                     sep = "\t") %>% 
  rownames_to_column("correspondance") %>% 
  unite(Position, CHROM, POS)
print("Finished importing input")
############################ Load HMM output ##########################
output_delta_freqs_france <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/03_HMM/Results_HMM/France.weir_2state_HMMstates.txt",
                                        header = TRUE,
                                        sep = " ") %>% 
  rownames_to_column("correspondance") %>% 
  # Merge with the input to get the position names
  left_join(input_delta_freqs_france, by = "correspondance") %>% 
  select(-correspondance) %>% 
  rename(Group = x,
         Fst = WEIR_AND_COCKERHAM_FST) %>% 
  relocate(Position) %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  mutate(Position = Position %>% as.numeric) %>% 
  filter(grepl("LG", Chromosome))

output_delta_freqs_sweden <- read.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/03_HMM/Results_HMM/Sweden.weir_2state_HMMstates.txt",
                                        header = TRUE,
                                        sep = " ") %>% 
  rownames_to_column("correspondance") %>% 
  # Merge with the input to get the position names
  left_join(input_delta_freqs_sweden, by = "correspondance") %>% 
  select(-correspondance) %>% 
  rename(Group = x,
         Fst = WEIR_AND_COCKERHAM_FST) %>% 
  relocate(Position) %>% 
  mutate(pos = Position) %>% 
  transform_position_ade2tidy() %>% 
  mutate(Position = Position %>% as.numeric) %>% 
  filter(grepl("LG", Chromosome))

rm(list = apropos("input_delta_freqs"))
print("Finished importing output")

############################ Plot the outputs ##########################
output_delta_freqs_sweden %>% 
  geom_manhattan(aes(y = Fst, color = Group), palette = c("grey70", "turquoise4"), filter_low_freqs = FALSE)

output_delta_freqs_france %>% 
  geom_manhattan(aes(y = Fst, color = Group), palette = c("grey70", "turquoise4"), filter_low_freqs = FALSE)

output_delta_freqs_france %>%
 mutate(Population = "France") %>%
 rbind(output_delta_freqs_sweden %>%
         mutate(Population = "Sweden")) %>%
 mutate(Population = Population %>% factor(levels = c("Sweden", "France"))) %>%
 geom_manhattan(aes(y = Fst, color = Group, facetting = Population), palette = c("grey70", "turquoise4"), filter_low_freqs = FALSE) +
  facet_col(vars(Population))


############################ Separate the inversions from the rest ##########################
# Sweden
## Plot the chromosomes
output_delta_freqs_sweden %>%
  mutate(bin = (Position /bin_size) %>% round) %>%
  group_by(Chromosome, bin, Group) %>%
  summarize(counting = n()) %>%
  filter(counting == max(counting)) %>%
  select(-counting) %>%
  ungroup %>%
  rename(Position = bin) %>%
  geom_manhattan(aes(y = Group, colour = Group, facets = Chromosome)) +
  facet_wrap(vars(Chromosome), scale = "free")


####################### Visualise and add the inversions to the list ##############
list_chroms <- output_delta_freqs_sweden %>% 
  select(Chromosome) %>% unique %>% as.vector %>% unlist %>% unname
output_delta_freqs_sweden %>% 
  filter(grepl("LG17_", Chromosome)) %>% 
  mutate(bin = (Position / bin_size) %>% round) %>% 
  group_by(Chromosome, bin, Group) %>% 
  summarize(counting = n()) %>% 
  filter(counting == max(counting)) %>% 
  select(-counting) %>% 
  ungroup %>% 
  rename(Position = bin) %>%
  # geom_manhattan(aes(y = Group, colour = Group))
  filter(Group == 2) %>%
  mutate(pos_prev = lag(Position, default = 0),
                                  difference_pos_max = Position - pos_prev) %>% 
  as.data.frame

Positions_candidate_inversions_sweden <- output_delta_freqs_sweden %>% 
  mutate(Position = (Position / bin_size) %>% round,
         Inversion = case_when(
           (grepl("LG1_", Chromosome) & Position %in% c(69:79)) ~ "Isl_diff_1.1",
           (grepl("LG1_", Chromosome) & Position %in% c(385:410)) ~ "Isl_diff_1.2",
           (grepl("LG1_", Chromosome) & Position %in% c(475:492)) ~ "Isl_diff_1.3",
           (grepl("LG1_", Chromosome) & Position %in% c(826:842)) ~ "Isl_diff_1.4",
           (grepl("LG2_", Chromosome) & Position %in% c(320:340)) ~ "Isl_diff_2.1",
           (grepl("LG2_", Chromosome) & Position %in% c(350:397)) ~ "Isl_diff_2.2",
           (grepl("LG2_", Chromosome) & Position %in% c(447:518)) ~ "Isl_diff_2.3",
           (grepl("LG2_", Chromosome) & Position %in% c(622:652)) ~ "Isl_diff_2.4",
           (grepl("LG2_", Chromosome) & Position %in% c(792:807)) ~ "Isl_diff_2.5",
           (grepl("LG3_", Chromosome) & Position %in% c(153:454)) ~ "Isl_diff_3.1",
           (grepl("LG3_", Chromosome) & Position %in% c(478:742)) ~ "Isl_diff_3.2",
           (grepl("LG4_", Chromosome) & Position %in% c(43:146)) ~ "Isl_diff_4.1",
           (grepl("LG4_", Chromosome) & Position %in% c(330:466)) ~ "Isl_diff_4.2",
           (grepl("LG4_", Chromosome) & Position %in% c(538:702)) ~ "Isl_diff_4.3",
           (grepl("LG6_", Chromosome) & Position %in% c(1:112)) ~ "Isl_diff_6.1",
           (grepl("LG6_", Chromosome) & Position %in% c(137:196)) ~ "Isl_diff_6.2",
           (grepl("LG6_", Chromosome) & Position %in% c(518:542)) ~ "Isl_diff_6.3",
           (grepl("LG7_", Chromosome) & Position %in% c(41:98)) ~ "Isl_diff_7.1",
           (grepl("LG7_", Chromosome) & Position %in% c(271:284)) ~ "Isl_diff_7.2",
           (grepl("LG7_", Chromosome) & Position %in% c(296:306)) ~ "Isl_diff_7.3",
           (grepl("LG7_", Chromosome) & Position %in% c(470:588)) ~ "Isl_diff_7.4",
           (grepl("LG8_", Chromosome) & Position %in% c(193:204)) ~ "Isl_diff_8.1",
           (grepl("LG8_", Chromosome) & Position %in% c(411:427)) ~ "Isl_diff_8.2",
           (grepl("LG8_", Chromosome) & Position %in% c(505:523)) ~ "Isl_diff_8.3",
           (grepl("LG8_", Chromosome) & Position %in% c(540:717)) ~ "Isl_diff_8.4",
           (grepl("LG10_", Chromosome) & Position %in% c(326:336)) ~ "Isl_diff_10.1",
           (grepl("LG10_", Chromosome) & Position %in% c(473:484)) ~ "Isl_diff_10.2",
           (grepl("LG11_", Chromosome) & Position %in% c(0:107)) ~ "Isl_diff_11.1",
           (grepl("LG11_", Chromosome) & Position %in% c(236:251)) ~ "Isl_diff_11.2",
           (grepl("LG12_", Chromosome) & Position %in% c(61:78)) ~ "Isl_diff_12.1",
           (grepl("LG12_", Chromosome) & Position %in% c(91:105)) ~ "Isl_diff_12.2",
           (grepl("LG12_", Chromosome) & Position %in% c(482:505)) ~ "Isl_diff_12.3",
           (grepl("LG13_", Chromosome) & Position %in% c(131:143)) ~ "Isl_diff_13.1",
           (grepl("LG13_", Chromosome) & Position %in% c(260:449)) ~ "Isl_diff_13.2",
           (grepl("LG14_", Chromosome) & Position %in% c(1:241)) ~ "Isl_diff_14.1",
           (grepl("LG14_", Chromosome) & Position %in% c(382:415)) ~ "Isl_diff_14.2",
           (grepl("LG14_", Chromosome) & Position %in% c(447:457)) ~ "Isl_diff_14.3",
           (grepl("LG15_", Chromosome) & Position %in% c(418:433)) ~ "Isl_diff_15.1",
           (grepl("LG16_", Chromosome) & Position %in% c(180:203)) ~ "Isl_diff_16.1",
           (grepl("LG16_", Chromosome) & Position %in% c(213:348)) ~ "Isl_diff_16.2",
           TRUE ~ NA
         )) %>%
  filter(!is.na(Inversion)) %>% 
  group_by(Chromosome, Inversion) %>% 
  summarize(Start = min(Position),
            End = max(Position),
            Length = End - Start) %>% 
  ungroup %>% 
  filter(Length >= 10) %>% 
  mutate(Length = Length * bin_size,
         End = End * bin_size,
         Start = Start * bin_size,
         Chromosome = Chromosome %>% 
           factor(levels = c(output_delta_freqs_france %>% 
                    select(Chromosome) %>% unique %>% 
                    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                    as.vector %>% unname %>% unlist))) %>% 
  arrange(Chromosome)


## France
list_chroms <- output_delta_freqs_france %>% 
  select(Chromosome) %>% unique %>% 
  arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
  as.vector %>% unname %>% unlist
output_delta_freqs_france %>% 
  filter(grepl("LG16_", Chromosome)) %>% 
  mutate(bin = (Position / bin_size) %>% round) %>% 
  group_by(Chromosome, bin, Group) %>% 
  summarize(counting = n()) %>% 
  filter(counting == max(counting)) %>% 
  select(-counting) %>% 
  ungroup %>% 
  rename(Position = bin) %>% 
  # geom_manhattan(aes(y = Group, colour = Group))
  filter(Group == 2) %>%
  mutate(pos_prev = lag(Position, default = 0),
         difference_pos_max = Position - pos_prev) %>% 
  as.data.frame

Positions_candidate_inversions_france <- output_delta_freqs_france %>% 
  mutate(Position = (Position / bin_size) %>% round,
         Inversion = case_when(
           (grepl("LG1_", Chromosome) & Position %in% c(388:398)) ~ "Isl_diff_1.1",
           (grepl("LG1_", Chromosome) & Position %in% c(475:507)) ~ "Isl_diff_1.2",
           (grepl("LG2_", Chromosome) & Position %in% c(327:339)) ~ "Isl_diff_2.1",
           (grepl("LG2_", Chromosome) & Position %in% c(372:417)) ~ "Isl_diff_2.2",
           (grepl("LG2_", Chromosome) & Position %in% c(452:510)) ~ "Isl_diff_2.3",
           (grepl("LG3_", Chromosome) & Position %in% c(155:453)) ~ "Isl_diff_3.1",
           (grepl("LG3_", Chromosome) & Position %in% c(484:742)) ~ "Isl_diff_3.2",
           (grepl("LG4_", Chromosome) & Position %in% c(43:109)) ~ "Isl_diff_4.1",
           (grepl("LG4_", Chromosome) & Position %in% c(330:466)) ~ "Isl_diff_4.2",
           (grepl("LG4_", Chromosome) & Position %in% c(549:702)) ~ "Isl_diff_4.3",
           (grepl("LG6_", Chromosome) & Position %in% c(1:104)) ~ "Isl_diff_6.1",
           (grepl("LG6_", Chromosome) & Position %in% c(155:197)) ~ "Isl_diff_6.2",
           (grepl("LG6_", Chromosome) & Position %in% c(518:536)) ~ "Isl_diff_6.3",
           (grepl("LG7_", Chromosome) & Position %in% c(47:98)) ~ "Isl_diff_7.1",
           (grepl("LG7_", Chromosome) & Position %in% c(272:300)) ~ "Isl_diff_7.2",
           (grepl("LG7_", Chromosome) & Position %in% c(478:577)) ~ "Isl_diff_7.3",
           (grepl("LG8_", Chromosome) & Position %in% c(518:530)) ~ "Isl_diff_8.1",
           (grepl("LG8_", Chromosome) & Position %in% c(540:589)) ~ "Isl_diff_8.2",
           (grepl("LG8_", Chromosome) & Position %in% c(600:684)) ~ "Isl_diff_8.3",
           (grepl("LG10_", Chromosome) & Position %in% c(311:325)) ~ "Isl_diff_10.1",
           (grepl("LG11_", Chromosome) & Position %in% c(0:106)) ~ "Isl_diff_11.1",
           (grepl("LG13_", Chromosome) & Position %in% c(265:445)) ~ "Isl_diff_13.1",
           (grepl("LG14_", Chromosome) & Position %in% c(1:241)) ~ "Isl_diff_14.1",
           (grepl("LG14_", Chromosome) & Position %in% c(380:415)) ~ "Isl_diff_14.2",
           (grepl("LG16_", Chromosome) & Position %in% c(14:25)) ~ "Isl_diff_16.1",
           (grepl("LG16_", Chromosome) & Position %in% c(180:196)) ~ "Isl_diff_16.2",
           (grepl("LG16_", Chromosome) & Position %in% c(216:348)) ~ "Isl_diff_16.3",
           TRUE ~ NA
         )) %>%
  filter(!is.na(Inversion)) %>% 
  group_by(Chromosome, Inversion) %>% 
  summarize(Start = min(Position),
            End = max(Position),
            Length = End - Start) %>% 
  ungroup %>% 
  filter(Length >= 10) %>% 
  mutate(Length = (End - Start) * bin_size,
         End = End * bin_size,
         Start = Start * bin_size,
         Chromosome = Chromosome %>% 
           factor(levels = c(output_delta_freqs_france %>% 
                               select(Chromosome) %>% unique %>% 
                               arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                               as.vector %>% unname %>% unlist))) %>% 
  arrange(Chromosome)


Positions_candidate_inversions_france %>% 
  mutate(Population = "France") %>% 
  rbind(Positions_candidate_inversions_sweden %>% 
          mutate(Population = "Sweden")) %>% 
  write.table("/shared/projects/pacobar/finalresult/Littorina_WGS_illumina/Sweden_France_parallelism/04_Inversions/Candidate_inversions.tsv", sep = "\t",
              col.names = TRUE, row.names = FALSE, quote = FALSE)
