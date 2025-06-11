
# Add a miror to download the libraries if needed
utils::setRepositories(ind = 0, addURLs = c(CRAN = "https://cran.irsn.fr/"))
# Install libraries if needed and load them
libraries <- c("tidyverse", "adegenet", "pegas", "vcfR", "readxl", "argparse")

if (!require("pacman")) install.packages("pacman")
for (lib in libraries){
  pacman::p_load(lib, character.only = TRUE)
}

############################ Parse and use arguments ##########################
# Take inputs from the snakemake program
parser <- ArgumentParser(description = "This program is used to make a list of individuals that are on the extreme parts of the transect")

# Add the arguments that are used
parser$add_argument("--vcf", "-v", help = "The path to the vcf file")
parser$add_argument("--metadata", "-m", help = "The path to the metadata file")
parser$add_argument("--outdir", "-o", help = "The path to the output")

xargs <- parser$parse_args()

vcf_file <- xargs$vcf
metadata_path <- xargs$metadata

outdir <- xargs$outdir
print("Imported arguments")

################## Get the argument given in input  ##################
source("../General_scripts/Functions_optimise_plot_clines.r")
################## Import the vcf file  ##################
data <- read.vcfR(vcf_file) %>% 
  vcfR2genind()

# We add the poputations to the vcf object
data@pop <- (data@tab %>% rownames %>% str_split_fixed(., "_", 4))[, 3] %>% as.factor
data@other$exposition <- (data@tab %>% rownames %>% str_split_fixed(., "_", 6))[, 5] %>% as.factor
data@other$exposition[which(data@other$exposition == "TRANSI")] <- "TRANS" %>% as.factor
data@other$exposition <- data@other$exposition %>% droplevels()
data@other$Sample_Name <- data@tab %>% rownames

print("Imported genetic data")

############################ PCA on whole genome ##########################

# Scale the genome to get rid of missing data
X <- scaleGen(data, NA.method="mean", scale=FALSE, center=TRUE)

# Run the PCA
pca <- dudi.pca(X, scale=TRUE, nf=5, scannf=FALSE)

rm(X)
print("Run PCA")

################## Import the metadata  ##################
metadata <- read_excel(path = metadata_path,
                       sheet = 1,
                       col_names = TRUE,
                       trim_ws = TRUE) %>%
  
  # Convert to the correct formats
  mutate(Species = as.factor(Species),
         ID_number = as.factor(ID_number),
         Population = as.factor(Population),
         Transect = as.factor(TRANSECT),
         Id = as.factor(ID),
         Shell_colour = factor(Shell.colour %>% str_to_title, levels = c("Black", "Black/Square", "Brown", "Brown/Square", "Dark", "Yellow", "Yellow/Brown", "Yellow/Square", "Grey", "White", "Banded", NA)),
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
  
  # Select only the necessary columns for the analysis
  select(-c(length, biSIZE, Shell.colour, ID, TRANSECT)) %>% 
  
  # select only the data we need (the one on fabalis just in the transects from LAM and LOKn)
  filter(Species == "FAB",
         Population != "BRE",
         Transect == "n") %>% 
  
  # Modify the population column to get only the name of the country
  mutate(Population = ifelse(Population == "LOK", "Sweden", "France") %>% 
           factor(levels = c("Sweden", "France"))) %>% 
  # Drop unused levels
  droplevels

print("Imported meta data")

################## Select the individuals that are on the extreme parts of the transect  ##################
# Get one random position to use as SNP subset
random_snp <- data@loc.fac[1] %>% 
       as.data.frame %>% 
       rename(SNP_name = ".") %>% 
       mutate(Position = paste0(SNP_name, ".1"))

# Get the extreme individuals on the ends of the transects
filtered_indivs_genotypes <- get_extreme_genotypes(
       genetic_data = data,
       SNP_subset = random_snp,
       Extreme_values = pca$li,
       var = "Axis2",
       meta_data = metadata
)

# Get the names for the individuals form the four populations
## We iterate over the four possible cases for the four populations that are considered here
for (pop in c("France", "Sweden")){
       for (expos in c("exposed", "sheltered")){
              ## For each population, we get the name of the table, 
              table_name <- paste0(outdir, pop, "_", expos, ".txt")
              ## We get the name of the population to use
              pop_name <- paste0(ifelse(pop == "France", "LAM", "LOK"), "_", ifelse(expos == "exposed", "EXPOS", "SHELT"))
              #And get the names of the samples to write in a table
              filtered_indivs_genotypes[pop_name] %>%
                     unname %>%
                     as.data.frame %>%
                     select(Sample_Name) %>%
                     write.table(table_name, col.names = FALSE, row.names = FALSE,
                                 sep = "\t", quote = FALSE)
       }
}

print("Done !")
