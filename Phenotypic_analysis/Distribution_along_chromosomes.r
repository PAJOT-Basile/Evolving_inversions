# Libraries
require("anyLib")
library(anyLib)
anyLib(c("tidyverse", "ggh4x"))


# Import data
data <- read.csv("/shared/home/bpajot/littorina/finalresult/bpajot/outputs/Stats_VCF/Number_SNPs_per_region.csv", sep = ";", header=FALSE) %>% 
  rename(Chromosome = V1,
         Nb_SNPs = V2) %>% 
  separate(Chromosome, c("Chromosome", "Portion"), sep = ":") %>% 
  separate(Portion, c("Start", "End"), sep = "-") %>% 
  mutate(Chromosome = Chromosome %>% as.factor,
         Nb_SNPs = Nb_SNPs %>% as.numeric,
         Start = Start %>% as.numeric,
         End = End %>% as.numeric)


# Plot the distribution of the mapping regions
data %>% 
  ggplot(aes(x=Start, y=Nb_SNPs/(End-Start))) + 
  geom_point() +
  geom_line() +
  facet_wrap2(~Chromosome, scales="free_x") +
  labs(title = "Distribution of SNP when mapping L. fabalis on L. saxatilis reference genome")

data1 <- data
