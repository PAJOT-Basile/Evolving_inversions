# Import libraries
libraries <- c("argparse", "tidyverse")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

# Import the arguments
parser <- ArgumentParser(description = "This program is used to plot the LD heatmaps")

# Add the arguments that are used
parser$add_argument("--input", "-i", help = "The input file (LD output from vcftools)")
parser$add_argument("--output", "-o", help = "The output path and name.")
parser$add_argument("--bin_size", "-b", help = "The bin size to use")

xargs <- parser$parse_args()

input_path <- xargs$input
output_file <- xargs$output
bin_size <- xargs$bin_size %>% as.numeric


# Function to create a color palette
colfunc <- colorRampPalette(colors = c("#813d32","#ef852f","#fdc441" ,"#e3cfb4") %>% rev, bias = 2)

LG_heatmap <- read.table(input_path, header = TRUE, sep = "\t") %>% 
  rename(Chromosome = CHR,
         LD = R.2)
print("Imported")
# Compute a mean LD by the chosen bin
(LG_heatmap %>% 
  mutate(Position1 = (POS1 / bin_size) %>% round,
         Position2 = (POS2 / bin_size) %>% round) %>% 
  group_by(Position1, Position2) %>% 
  summarize(mean_LD = mean(LD)) %>% 
  ggplot() +
  geom_tile(aes(x = Position2, y = Position1, fill = mean_LD)) +
  scale_fill_gradientn(colors = colfunc(20)) +
  labs(x = paste0("Position along the chromosome (*", bin_size, ")"),
       y = paste0("Position along the chromosome (*", bin_size, ")")) +
  theme_bw() +
  theme(text = element_text(size = 20))) %>% 
  ggsave(plot = ., filename = output_file,
         device = "png", width = 800, height = 500, units = "px", scale = 5)
print("Done")
