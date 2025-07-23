# Import libraries
libraries <- c("tidyverse", "ggforce", "viridis", "ggnewscale", "LaplacesDemon", "bbmle", "zeallot", "docstring", "phangorn", "MetBrewer", "ggtree", "reshape2")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(characters = libraries, character.only = TRUE)
if (!require("tie")) remotes::install_github("romainfrancois/tie")
rm(libraries)

#library(tidyverse); library(ggforce); library(viridis); library(ggnewscale); library(LaplacesDemon, lib.loc = "/shared/home/bpajot/R/x86_64-conda-linux-gnu-library/4.2/"); library(tie, lib.loc = "/shared/home/bpajot/R/x86_64-conda-linux-gnu-library/4.2/"); library(bbmle); library(zeallot); library(docstring, lib.loc = "/shared/home/bpajot/R/x86_64-conda-linux-gnu-library/4.2/")

# Import functions from the `Cline_functions.R` script
source("/shared/projects//pacobar/finalresult/bpajot/Stage_Roscoff/scripts/A_Genetic_analysis/General_scripts/Cline_functions.R")

################################################################################
######################### 1. Useful functions  #################################
################################################################################

is_numeric_in_character <- function(x){
  #' Checks for numerics
  #' 
  #' Checks if there is a numeric value contained in a character
  #' 
  #' @param x (string, factor, vector).
  #'      Object to check. It can be a character, a factor, a vector or anything
  #'      that can be converted to numeric.
  #' @returns (bool).
  #'      If there is a numeric value in the input, it returns TRUE, it returns
  #'      FALSE otherwise
  #' @export
  
  return(!is.na(x %>% as.numeric) %>% suppressWarnings)
}

is.continuous <- function(x){
  #' Checks for continuousity
  #' 
  #' This function checks if a variable is discrete or continuous.
  #' 
  #' @param x (numeric vector).
  #'      A vector containing values to see if it is continuous or discrete
  #' @returns (bool).
  #'      This function returns TRUE if the variable is continuous (more than 10
  #'      levels) or FALSE if the variable is discrete (less than 10 levels)
  #' @section Warning:
  #' Here a variable is considered discrete if it contains less than 10
  #' distinct levels. This approximation is done to simplify the distinction 
  #' between continuous or discrete variables, but it is not a real classification
  #' @export
  
  return(length(unique(x)) >= 10)
}

"%!in%" <- function(x, y){
  #' Not in
  #' 
  #' This function is the opposite of the "%in%" function. It checks if the 
  #' values contained in input x are NOT contained in y.
  #' 
  #'@param x (vector).
  #'    This input contains any type of object, as long as they are the same type
  #'    as the ones in the y vector
  #'@param y (vector).
  #'    This input contains any type of object, as long as they are the same type
  #'    as the ones in the x vector
  #'    
  #'@returns (bool)
  #'    This function returns TRUE for each value of x that is NOT in y and FALSE
  #'    for every value of x that is in y.
  #' @export
  
  return(!(x %in% y))
}

add_table_to_df_in_iteration <- function(df1, df2){
  #' Adds a data frame to another in an iteration
  #' 
  #' In an iteration, if we want to rbind a data frame to another, we have to 
  #' check if the data frame is empty and if it is not, we can rbind the two
  #' data frames together
  #' 
  #' @param df1 (data.frame).
  #'    This data frame is the name of the data frame the data will be added to
  #' @param df2 (data.frame).
  #'    This data frame contains the data to add to the df1.
  #'    
  #' @returns (data.frame).
  #'    This is the binding of the two data frames
  #'
  #' @export
  
  # First, we check if the first data frame exists. If if does not, we create the
  # variable
  #if (!exists(quote(df1), where = .GlobalEnv)){
  #  # If the data frame does not exist, we have to transform it to match what we
  #  # want
  #  df1_name <- as.name(substitute(df1)) %>% as.character
  #  # Then, we create an emtpty data fame
  #  x <- data.frame()
  #  assign(df1_name, x, pos = .GlobalEnv)
  #}
  
  # Then, we have to check if the first data frame is empty
  if (nrow(df1) == 0){
    # If the first data frame is empty, we simply attribute the value of df2 to df1
    df1 <- df2
  }else{
    # Otherwise we rbind the both of them
    df1 <- df1 %>% 
      rbind(df2)
  }
  return(df1)
}

progress_bar <- function(iteration_number, nb_iterations_to_do){
  #' Make a progress bar
  #' 
  #' This function allows to show a progress bar while a function with iterations
  #' is running. It is used in a for loop. 
  #' 
  #' @param iteration_number (numeric value).
  #'    This number indicates what iteration we are on
  #' @param nb_iteration_to_do (numeric value). 
  #'    This number indicates how many iterations are to be done.
  #' 
  #' @returns 
  #'  This function does not return anything. It simply prints a progress bar 
  #'  that is incremented for each iteration.
  #' @export
  
  # We transform the iteration into a percentage to see the advance
  percentage <- (iteration_number / nb_iterations_to_do * 100) %>% round(digits = 2)
  # We separate this into to bars. A filed bar and an empty bar to complete the
  # lines
  charged_bar <- ""
  empty_bar <- ""
  # We make the percentage smaller to represent it on the screen
  rep_perc <- (percentage/2) %>% round
  # We simply separate the case where the bars have to be empty
  if (rep_perc != 0){
    # We fill the bars by iterating over them
    for (j in 1:(rep_perc)){
      charged_bar <- paste0(charged_bar, "█")
    }
  }
  if (rep_perc != 50){
    for (j in 1:((50-rep_perc) %>% round)){
      empty_bar <- paste0(empty_bar, " ")
    }
  }
  # And we print the output.
  cat("\r", paste0(" | ", charged_bar, empty_bar, " | ", percentage, " % of positions                       "))
  flush.console()
}

Create_dir_not_exist <- function(paths){
  #' Create directory if it does not exist yet
  #' 
  #' @description
    #' This function checks if a directory exists. If it does, it continues,
    #' otherwise, it creates said directory
    #' 
    #' @param paths The path of the directory to check and create shall it not exist.
    #' It can also support a vector of several directories
    #' 
    #' @return There is no output of this function
    #' 
    #' @usage Create_dir_not_exist(path)
    #' 
    #' @examples
        #' # This checks if the directory "mytest_directory" exists and creates it if not
        #' Create_dir_not_exist("./mytest_directory")
  for (path in paths){
    if (!dir.exists(path)){
      dir.create(path, recursive = TRUE)
    }
  }}


################################################################################
############### 2. Functions for genomics with adegenet  #######################
################################################################################

transform_position_ade2tidy <- function(df){
  #' Transform adegenet format to tidy format
  #' 
  #' The position of each SNP in the adegenet format is as follows: 
  #' "CHROMOSOME_SNPposition.allele". 
  #' This function transforms the position column from the previous format into
  #' Two columns containing the chromosome in one column and the numeric 
  #' position in the other column.
  #' 
  #'@param df (data.frame).
  #'    The dataframe containing the position column to transform.
  #'  
  #'@returns (data.frame).
  #'    This function returns the input dataframe with two separate chromosome and
  #'    position columns.
  #' @export

  df %>%
    # First, we remove the allele at the end of the adegenet position
    separate(Position, c("Position", "Allele"), sep="\\.") %>% 
    select(-Allele) %>% 
    # Then, we separate the position from the name of the chromsome from the 
    # numeric position on said chromosome
    separate(Position, c("Chromosome", "Position"), "\\_(?!.*_)", extra = "merge") %>% 
    # Finaly, we transform the position to have a numeric value
    mutate(Position = Position %>% as.numeric) %>%
    return()
}

select_good_SNPs <- function(df, colName){
  #' Selects good SNPs
  #' 
  #' This function selects the allele version of SNPs that have a given value 
  #' greater or equal to 0 and if this condition is not sufficient to select one
  #' allele, it picks the first one in the table.
  #' 
  #'@param df (data.frame).
  #'    Dataframe containing a column called "Position" that can be duplicated 
  #'    and a column that is given to be filtered.
  #'@param colName (symbol or string).
  #'    It is the name of a column that is given to be filtered on. All the values
  #'    of this column that are greater or equal to 0 will be kept.
  #' 
  #'@returns (data.frame).
  #'    This dataframe contains only alleles of the first dataframe that have
  #'    values greater or equal to 0 in the selected column
  #' @export
  
  # First, we convert the name of the column to a name that we can use as column
  # name in the tidyverse
  colName <- as.name(substitute(colName))
  
  df %>% 
    # We arrange the data frame by increasing value of the Position column to be
    # sure that both alleles are in the same place
    arrange(Position) %>% 
    # We filter the alleles that have a value greater or equal to 0 for the
    # selected column
    filter(!!as.symbol(colName) >= 0) %>% 
    # We add a column that contains the position of the previous row
    mutate(pos_bis = lag(Position, default = "default.position")) %>% 
    # We compare this newly created column with the position column. If the two
    # physical positions are identical between following alleles, we remove the
    # second one
    filter(str_split_fixed(Position, "\\.", 2)[, 1] != str_split_fixed(pos_bis, "\\.", 2)[, 1]) %>%
    # Finally, we remove the duplicated position column created to filter
    select(-pos_bis) %>% 
    return()
}

select_clinal_SNPs <- function(df){
  #' Select clinal SNPs
  #' 
  #' This function selects only clinal SNPs in at least one location. It requires
  #' that the input to have a "Position" and a "Model_to_select" column
  #' 
  #'@param df (data.frame).
  #'    This dataframe must contain two columns at least with the "Position" and
  #'    the "Model_to_select" column. The "Position" column contains the position
  #'    of the SNP on the Genome. The "Model_to_select" column must contain the
  #'    name of the model to select: one of "Clinal", "Stable" or "Linear".
  #' 
  #'@returns (data.frame).
  #'    This dataframe contains all the SNPs that have at least one clinal model.
  #'    If the SNP is present several times in the dataframe, it will select all
  #'    the lines containing said SNP if at least one model is clinal.
  #' @export

  # We iterate over each unique SNP position in the dataframe
  for (position_i in df$Position %>% unique){
    # We count the number of rows for each SNP that use the "Clinal" model 
    is_clinal <- df %>% 
      filter(Position == position_i,
             Best_model == "Clinal" & abs(Delta_AIC_second_best_model) > 4) %>% 
      nrow
    
    # If the SNP does not fit the clinal model, no rows will be found, so we 
    # remove it from the dataframe
    if (is_clinal == 0){
      df <- df %>% 
        filter(Position != position_i)
    }
  }
  
  # Once we iterated over all the dataframe's SNPs, we return it
  return(df)
}

Compute_Hobs <- function(data, MARGIN = 1){
  #' Compute Hobs
  #' 
  #' @description
    #' Computes the observed heterozygosity of the genomic dataset
  #' @param data data.frame or matrix. It should contain some genetic data
  #'            with values of genotypes 0, 1, and 2 with the position in
  #'            columns and individuals in the rows.
  #' @param MARGIN (default = 1) Computes the Heterozygosity per row (1) or per column (2)
  #' @return Hobs: value of the observed heterozygosity for the given dataframe. 
  
  # Choose the function to use depending on the MARGIN value
  if (MARGIN == 1){
    func2use <- function(...){
      return(colSums(...))
    }
    
    size_data <- nrow(data)
    
  }else if (MARGIN == 2){
    func2use <- function(...){
      return(rowSums(...))
    }
    
    size_data <- ncol(data)
  }
  # Compute the missing data for each position
  missing_data <- data %>% 
    is.na %>% 
    func2use(na.rm = TRUE)
  
  # Compute the number of heterozygotes for each position
  Nb_hetero <- (data == 1) %>% 
    func2use(na.rm = TRUE)
  
  # Compute the heterozygosity
  Hobs <- Nb_hetero / (size_data - missing_data)
  
  # Return the computed Hobs
  return(Hobs)
}

Transform_Hobs2df <- function(x){
  x %>% 
    bind_rows() %>% 
    t %>% 
    as.data.frame() %>% 
    rownames_to_column("Position") %>% 
    filter(Position != "Sample_Name") %>% 
    rename(Exposed = V1,
           Transition = V2,
           Sheltered = V3) %>% 
    separate(Position, c("Position", "allele"), sep = "\\.") %>% 
    select(-allele) %>% 
    unique %>% 
    pivot_longer(!Position, names_to = "Names", values_to = "Hobs") %>% 
    return()
}


################################################################################
########## 3. Functions for allelic frequencies from adegent objects  ##########
################################################################################

#rm(genetic_data, SNP_subset)
get_genotype_transect <- function(genetic_data, SNP_subset = NULL, meta_data = metadata){
  #' Get genotypes along transect
  #' 
  #' This function uses an adegenet object and returns the genotype of all or
  #' part of the SNPs.
  #' 
  #'@param genetic_data (Genind object). 
  #'    It containing the SNPs from which to get the genotype
  #'@param SNP_subset (data.frame). (default = NULL)
  #'    This data frame contains columns "Position" and/or "SNP_name". If it is
  #'    given, it will only return the genotype from these SNPs. If this data
  #'    frame is used, it is required to have a column called "SNP_name".
  #'@param meta_data (data.frame).
  #'    This dataframe contains metadata on the individuals (for example the
  #'    position on the transect, the color, ...)
  #' 
  #'@returns (data.frame).
  #'    This data frame contains the genotypes for both alleles of the selected
  #'    SNPs as well as some metadata on the individuals (position along the 
  #'    transect)
  #' @export
  
  # We select the table of genotypes for all SNPs along the whole genome
  return_table <- genetic_data@tab %>%
    # We transform this into a data frame to keep the name of the columns
    as.data.frame %>% 
    # We add a column to the table to get the name of the individuals
    rownames_to_column("Sample_Name") %>% 
    # We merge this with the metadata using the name of the individuals
    inner_join(meta_data, by = "Sample_Name")
  
  # If the SNP subset is not selected, do this
  if ((SNP_subset %>% is.null) || ("Position" %!in% names(SNP_subset))){
    return_table %>% 
      # We select the columns we need to keep (position along the transect, 
      # name of the individuals, the population they are from and all the 
      # genotypes)
      select(Sample_Name, LCmeanDist, Population, starts_with("LG")) %>%
      return()
  }else{
    return_table %>% 
      # The only difference is in the following line, we select the columns
      # that contain the genotypes of only the alleles of the SNPs we want to
      # keep
      select(Sample_Name, LCmeanDist, Population, SNP_subset$Position) %>% 
      return()
    }
}

get_extreme_genotypes <- function(genetic_data, SNP_subset = NULL, Extreme_values = NULL, var = NULL, nb_extreme_indivs = 30, nb_indivs_to_keep = 20, meta_data = metadata){
  #' Get extreme genotypes
  #' 
  #' This function is used to get the genotypes of the individuals located on the
  #' extreme ends of the transect. If the selected individuals do not have the 
  #' right genotype, we can also use a supplementary dataset to re-filter the 
  #' individuals to select ones with more differences.
  #' 
  #'@param genetic_data (Genind object). 
  #'    It containing the SNPs from which to get the genotype
  #'@param SNP_subset (data.frame). (default = NULL)
  #'    This data frame contains columns "Position" and/or "SNP_name". If it is
  #'    given, it will only return the genotype from these SNPs. If this data
  #'    frame is used, it is required to have a column called "SNP_name".
  #'@param Extreme_values (data.frame). (default = NULL)
  #'    This data frame contains the values to use to re-filter the selected
  #'    individuals. This variable works in combination with the next one.
  #'@param var (string). (default = NULL)
  #'    This string is the name of the column from Extreme_values to use to 
  #'    refilter the individuals.
  #'@param nb_extreme_indivs (integer). (default = 30)
  #'    This number is the number of individuals to sample in the first samping
  #'    process (position on the transect).
  #'@param nb_indivs_to_keep (integer). (default = 20)
  #'    This number is the number of individuals to re-select in the
  #'    nb_extreme_indivs using the Extreme_values dataset. It must therefore be
  #'    smaller than the number indicated by nb_extreme_indivs
  #'@param meta_data (data.frame).
  #'    This dataframe contains metadata on the individuals (for example the
  #'    position on the transect, the color, ...)
  #'    
  #'@returns (list of 4 data.frames).
  #'    This list contains 4 dataframes corresponding to the 4 groups of individuals
  #'    (the two extremes of each transect). Each data frame contains the 
  #'    genotypes of the selected individuals.
  #' @export
  
  # First, we have to get the genotype of the considered SNPs on the whole transect
  genotype_transect <- get_genotype_transect(genetic_data = genetic_data,
                                             SNP_subset = SNP_subset,
                                             meta_data = meta_data)
  
  # If there is a second filtration step, the arguments "Extreme values" and
  # "var" will not be NULL. We separate this case from the rest. As the steps in
  # the four populations are the same, we will only comment the first one and then
  # call your attention to small modulations
  if (!is.null(Extreme_values) & !is.null(var)){
    # First, we take the example of the French population located in the exposed
    # part of the transect (maximal values of position)
    LAM_EXPOS <- genotype_transect %>% 
      # We select only the french individuals in the gentoypes dataframe
      filter(Population == "France") %>% 
      # We arrange the position along the transect in descending order to have
      # the most exposed individuals in the first rows 
      arrange(LCmeanDist %>% desc) %>% 
      # Then, we take only the first nb_extreme_indivs rows with the individuals
      # that are the most exposed
      dplyr::slice(1:nb_extreme_indivs) %>% 
      # Then, for the second filtration step, we have to merge the filtered data
      # frame to the "Extreme_values" data frame using the names of the individuals
      left_join(Extreme_values %>% 
                  select(var %>% matches) %>% 
                  rownames_to_column("Sample_Name"),
                by = "Sample_Name") %>% 
      # And again, we arrange the selected column by descending order to have 
      # the biggest values of said column in the first rows
      arrange(!!as.symbol(var) %>% desc) %>% 
      # And we keep only the first nb_indivs_to_keep rows
      dplyr::slice(1:nb_indivs_to_keep)
    
    # For the second end of the transect, the process is similar
    LAM_SHELT <- genotype_transect %>% 
      filter(Population == "France") %>% 
      # The only differences are that we arrange by increassing order the 
      # position along the transect to have the smallest values in the first rows
      arrange(LCmeanDist) %>% 
      dplyr::slice(1:nb_extreme_indivs) %>% 
      left_join(Extreme_values %>% 
                  select(var %>% matches) %>% 
                  rownames_to_column("Sample_Name"),
                by = "Sample_Name") %>% 
      # And in the same way, we order the selected column of the "Extreme_values"
      # data frame by increasing order.
      arrange(!!as.symbol(var)) %>% 
      dplyr::slice(1:nb_indivs_to_keep)
    
    # And then, we do the same thing for the swedish population by filtering only
    # the swedish individuals
    LOK_EXPOS <- genotype_transect %>% 
      filter(Population == "Sweden") %>% 
      arrange(LCmeanDist %>% desc) %>% 
      dplyr::slice(1:nb_extreme_indivs) %>% 
      left_join(Extreme_values %>% 
                  select(var %>% matches) %>% 
                  rownames_to_column("Sample_Name"),
                by = "Sample_Name") %>% 
      arrange(!!as.symbol(var) %>% desc) %>% 
      dplyr::slice(1:nb_indivs_to_keep)
    
    LOK_SHELT <- genotype_transect %>% 
      filter(Population == "Sweden") %>% 
      arrange(LCmeanDist) %>% 
      dplyr::slice(1:nb_extreme_indivs) %>% 
      left_join(Extreme_values %>% 
                  select(var %>% matches) %>% 
                  rownames_to_column("Sample_Name"),
                by = "Sample_Name") %>% 
      arrange(!!as.symbol(var)) %>% 
      dplyr::slice(1:nb_indivs_to_keep)
    
  }else{
    # If there is no second step of filtration, we simply do the beginnig.
    LAM_EXPOS <- genotype_transect %>% 
      # First, we filter by population
      filter(Population == "France") %>%
      # Then we arder the position along the transect by decresing order for the
      # exposed part and by increasing order for the sheltered part ...
      arrange(LCmeanDist %>% desc) %>% 
      # And keep only the first nb_indivs_to_keep rows to keep only the most 
      # exposed/sheltered individuals
      dplyr::slice(1:nb_indivs_to_keep)
    
    # Again, this is the same as before, the only difference is the direction by
    # which the position along the transect is ordered
    LAM_SHELT <- genotype_transect %>% 
      filter(Population == "France") %>% 
      arrange(LCmeanDist) %>% 
      dplyr::slice(1:nb_indivs_to_keep)
    
    # And we do the same thing for the swedish population
    LOK_EXPOS <- genotype_transect %>% 
      filter(Population == "Sweden") %>% 
      arrange(LCmeanDist %>% desc) %>% 
      dplyr::slice(1:nb_indivs_to_keep)
    
    LOK_SHELT <- genotype_transect %>% 
      filter(Population == "Sweden") %>% 
      arrange(LCmeanDist) %>% 
      dplyr::slice(1:nb_indivs_to_keep)
  }
  # Finally, we return a list containing all four data frames.
  return(list("LAM_EXPOS" = LAM_EXPOS, "LAM_SHELT" = LAM_SHELT, "LOK_EXPOS" = LOK_EXPOS, "LOK_SHELT" = LOK_SHELT))
}

get_allelic_frequencies <- function(genetic_data, SNP_subset = NULL, Extreme_values = NULL, var = NULL, nb_extreme_indivs = 30, nb_indivs_to_keep = 20, meta_data = metadata){
  #' Calculte allelic frequencies from adegenet objects
  #' 
  #' This function is used to calculate the allelic frequencies of the obtained
  #' genotypes for each extreme of each transect
  #' 
  #'@param genetic_data (Genind object). 
  #'    It containing the SNPs from which to get the genotype
  #'@param SNP_subset (data.frame). (default = NULL)
  #'    This data frame contains columns "Position" and/or "SNP_name". If it is
  #'    given, it will only return the genotype from these SNPs. If this data
  #'    frame is used, it is required to have a column called "SNP_name".
  #'@param Extreme_values (data.frame). (default = NULL)
  #'    This data frame contains the values to use to re-filter the selected
  #'    individuals. This variable works in combination with the next one.
  #'@param var (string). (default = NULL)
  #'    This string is the name of the column from Extreme_values to use to 
  #'    refilter the individuals.
  #'@param nb_extreme_indivs (integer). (default = 30)
  #'    This number is the number of individuals to sample in the first samping
  #'    process (position on the transect).
  #'@param nb_indivs_to_keep (integer). (default = 20)
  #'    This number is the number of individuals to re-select in the
  #'    nb_extreme_indivs using the Extreme_values dataset. It must therefore be
  #'    smaller than the number indicated by nb_extreme_indivs
  #'@param meta_data (data.frame).
  #'    This dataframe contains metadata on the individuals (for example the
  #'    position on the transect, the color, ...)
  #'    
  #'@returns (data.frame).
  #'      This data frame contains the allelic frequencies of the exposed and
  #'      sheltered parts of the transect for both populations. Each row is a SNP
  #'      and for each one, is given the population (here Sweden or France),
  #'      the allelic frequency in the sheltered and in the exposed part of the 
  #'      transect.
  #' @export
  
  # First, we get the genotypes of extreme individuals in the transect
  genotypes <- get_extreme_genotypes(genetic_data = genetic_data, 
                                     SNP_subset = SNP_subset, 
                                     Extreme_values = Extreme_values, 
                                     var = var, 
                                     nb_extreme_indivs = nb_extreme_indivs, 
                                     nb_indivs_to_keep = nb_indivs_to_keep,
                                     meta_data = meta_data)
  c(LAM_EXPOS, LAM_SHELT, LOK_EXPOS, LOK_SHELT) %<-% genotypes
  # Calculate the allelic frequencies on the selected individuals
  p_LAM_EXPOS <- (LAM_EXPOS %>% 
                    # We select all the genotype columns
                    select(starts_with("LG")) %>% 
                    # We sum the genotypes and divide them by the number of
                    # individuals, to which we remove the missing data if there
                    # is any
                    colSums(na.rm = TRUE)) / ((nrow(LAM_EXPOS) - (LAM_EXPOS %>% 
                                                                    select(starts_with("LG")) %>% 
                                                                    is.na %>% 
                                                                    colSums(na.rm = TRUE))) * 2)
  # The process is the same for the french sheltered part and for the following
  # swedish parts
  p_LAM_SHELT <- (LAM_SHELT %>% 
                    select(starts_with("LG")) %>% 
                    colSums(na.rm = TRUE)) / ((nrow(LAM_SHELT) - (LAM_SHELT %>% 
                                                                    select(starts_with("LG")) %>% 
                                                                    is.na %>% 
                                                                    colSums(na.rm = TRUE))) * 2)
  p_LOK_EXPOS <- (LOK_EXPOS %>% 
                    select(starts_with("LG")) %>% 
                    colSums(na.rm = TRUE)) / ((nrow(LOK_EXPOS) - (LOK_EXPOS %>% 
                                                                    select(starts_with("LG")) %>% 
                                                                    is.na %>% 
                                                                    colSums(na.rm = TRUE))) * 2)
  p_LOK_SHELT <- (LOK_SHELT %>% 
                    select(starts_with("LG")) %>% 
                    colSums(na.rm = TRUE)) / ((nrow(LOK_SHELT) - (LOK_SHELT %>% 
                                                                    select(starts_with("LG")) %>% 
                                                                    is.na %>% 
                                                                    colSums(na.rm = TRUE))) * 2)
  # Finally, we merge all these togeteher and return the genotypes table
  ## First, we bind the french extreme frequencies
  Allele_freqs <- rbind(p_LAM_SHELT, p_LAM_EXPOS) %>% 
    # We then reshape and transform this to a tibble
    matrix(nrow = p_LAM_SHELT %>% length, byrow=TRUE) %>% 
    as_tibble %>% 
    # We rename the columns to reecognise them
    rename(p_left_shelt = V1,
           p_right_expos = V2) %>% 
    # We add a population column and the position on the genome
    mutate(Population = "France",
           Position = p_LAM_SHELT %>% names) %>% 
    # And then, we repeat the process for the swedish population and bind the two
    # tibbles
    rbind(rbind(p_LOK_SHELT, p_LOK_EXPOS) %>% 
            matrix(nrow = p_LOK_SHELT %>% length, byrow=TRUE) %>% 
            as_tibble %>% 
            rename(p_left_shelt = V1,
                   p_right_expos = V2) %>% 
            mutate(Population = "Sweden",
                   Position = p_LOK_SHELT %>% names))
  
  return(Allele_freqs)
}

Calculate_allelic_frequency <- function(genotype){
  #' Calculate allelic frequency from a genotype
  #' 
  #' This function calculates the allelic frequency of one allele in a population
  #' using the genotype of this allele in the population.
  #' 
  #'@param genotype (vector of numbers).
  #'    This vector contains genotype values for all the individuals in the
  #'    population. It can contain missing values
  #'  
  #'@returns (numeric value)
  #'    This value is the allelic frequency of the considered allele in the 
  #'    population.
  #' @export
  
  (sum(genotype, na.rm = TRUE) / ((length(genotype) - (is.na(genotype) %>% sum(na.rm = TRUE))) * 2)) %>% 
    return()
}

get_delta_freqs_and_F4 <- function(genetic_data, SNP_subset = NULL, Extreme_values = NULL, var = NULL, nb_extreme_indivs = 30, nb_indivs_to_keep = 20, meta_data = metadata){
  #' Calculates Delat frequencies and F4
  #' 
  #' This function is used to calculate the differences in allelic frequencies
  #' between both ends of a transect. It also calculates the F4 statistic between
  #' two populations (the multiplication of differences in allelic frequencies).
  #' 
  #'@param genetic_data (Genind object). 
  #'    It containing the SNPs from which to get the genotype
  #'@param SNP_subset (data.frame). (default = NULL)
  #'    This data frame contains columns "Position" and/or "SNP_name". If it is
  #'    given, it will only return the genotype from these SNPs. If this data
  #'    frame is used, it is required to have a column called "SNP_name".
  #'@param Extreme_values (data.frame). (default = NULL)
  #'    This data frame contains the values to use to re-filter the selected
  #'    individuals. This variable works in combination with the next one.
  #'@param var (string). (default = NULL)
  #'    This string is the name of the column from Extreme_values to use to 
  #'    refilter the individuals.
  #'@param nb_extreme_indivs (integer). (default = 30)
  #'    This number is the number of individuals to sample in the first samping
  #'    process (position on the transect).
  #'@param nb_indivs_to_keep (integer). (default = 20)
  #'    This number is the number of individuals to re-select in the
  #'    nb_extreme_indivs using the Extreme_values dataset. It must therefore be
  #'    smaller than the number indicated by nb_extreme_indivs
  #'@param meta_data (data.frame).
  #'    This dataframe contains metadata on the individuals (for example the
  #'    position on the transect, the color, ...)
  #'    
  #'@returns (data.frame).
  #'      This data frame contains the allelic frequencies of the exposed and
  #'      sheltered parts of the transect for both populations. Each row is a SNP
  #'      and for each one is given the difference in allelic frequency between 
  #'      the exposed and the sheltered part of the transect in both populations,
  #'      and the F4 statistic.
  #' @export
  
  # First, we need the allelic frequencies of the SNPs on both ends of the transect
  get_allelic_frequencies(genetic_data = genetic_data,
                          SNP_subset = SNP_subset,
                          Extreme_values = Extreme_values,
                          var = var,
                          nb_extreme_indivs = nb_extreme_indivs,
                          nb_indivs_to_keep = nb_indivs_to_keep,
                          meta_data = meta_data) %>% 
    # Then, we calculate the differences in allelic frequencies
    mutate(Delta_freq = p_right_expos - p_left_shelt) %>% 
    # Then ,we transform the table to have the differences in allelic frequencies
    # in two columns (one per population)
    pivot_wider(names_from = Population, values_from = c(p_right_expos, p_left_shelt, Delta_freq)) %>% 
    # Finally, we calculate the F4 statistic
    mutate(F4_stat = Delta_freq_France * Delta_freq_Sweden) %>% 
    return()
}

get_extreme_values <- function(df, condition, slice = FALSE, threshold = NULL){
  #' Get extreme values from a dataframe
  #' 
  #' This function is used to get the extreme values of a data frame with a
  #' specific condition and return a table containing the "Position" and "SNP_name"
  #' columns that can be used in the previous allelic frequency functions.
  #' 
  #'@param df (data.frame).
  #'    This data frame is the data frame to filter, from which to extract
  #'    extreme values
  #'@param condition (string).
  #'    This string contains the condition to get the extreme values. For example
  #'    you could use '"Comp2 > 0"'. This would filter the values of column Comp2
  #'    to be greater than 0.
  #'@param slice (boolean). (default = FALSE)
  #'    This indicates if there is a need to get only the first n values of the
  #'    filtered data frame. If TRUE, the threshold argument is required and it
  #'    will keep only these values. If FALSE, it will simply filter.
  #'@param threshold (integer). (default = NULL)
  #'    This argument works in combination with the slice argument. It indicates
  #'    how many maximum values you need to keep in the final data frame.
  #' 
  #'@returns (data.frame).
  #'    This data frame contains the values that validate the given condition 
  #'    (and the most extreme values of the input data frame.)
  #'    
  #'@section Warning:
  #'  This function is not very replicable. It only works with some cases. It
  #'  could be perfected with more time and more rigor.
  #' @export
  
  # First, we consider wether the order in which the condition wants us to
  # filter the data
  var <- (condition %>% str_split_fixed(., " ", 3))[, 1]
  order_arranging <- grepl("> 0", condition)
  
  # Then, we use the condition to filter the data
  x <- df %>%  
    filter(!! rlang::parse_expr(condition))
  
  # If the name of the position is not already in the columns, it means that it
  # is in the rownames, so we create the "Position" column from the rownames
  if ("Position" %!in% (x %>% names)){
    x <- x %>% 
      rownames_to_column("Position")
  }
  
  # Then, we transform the table to get a column called "SNP_name" using the
  # adegenet naming format 
  x <- x %>% 
    mutate(SNP_name = str_split_fixed(Position, "\\.", 2)[, 1])
  
  # Finally, if there is a threshold and that we want to slice, ...
  if (!(threshold %>% is.null) & (slice)){
    # ... we first check if the specified condition is not too complicated
    if (grepl(" ", (condition %>% str_split_fixed(., " ", 3))[, 3])){
      # If it is too complicated, we present a warning
      warning("The condition may be too complex for the slice to be effective. Please use a simpler argument like 'Comp2 > 0'.")
    }
    # However, we select the column we are interested in and arrange it in 
    # increasing order if we want to get the smallest values or by decreasing
    # order if we want to get the highest values
    x %>%
      select(var %>% matches) %>%
      arrange(ifelse(rep(order_arranging, nrow(.)), desc(!!as.symbol(var)), !!as.symbol(var))) %>%
      # Finally, we keep only the selected value of rows
      dplyr::slice(1:(threshold * nrow(.))) %>% 
      return()
  }else{
    return(x)
  }
}

################################################################################
######################### 4. Visualise genome scans   ##########################
################################################################################

geom_box_background <- function(To_plot_manhat, colName, chromosome_centers){
  #' Make polygon delimitations for geom_polygon
  #' 
  #' This function creates a data frame that contains the polygon delimitations
  #' if it is required to separate chromosomes by using boxes in the visual
  #' representations
  #' 
  #'@param To_plot_manhat (data.frame).
  #'    This data frame contains the data on the plotting of the values in a 
  #'    manhattan plot. 
  #'@param colName (string or name).
  #'    This argument is the name of the column in "To_plot_manhat" that we want
  #'    to plot.
  #'@param chromosomes_centers (data.frame).
  #'    This data frame contains information on the chromosomes: their name and
  #'    the position of their center is the minimum required for this function to
  #'    work.
  #'    
  #'@returns (data.frame)
  #'    This data frame contains the information on delimitations of the polygons
  #'    (position of the angles on the x and y axis). 
  #' @export
  
  # First, we transform the colName to use it with tidyverse
  colName <- as.character(colName)
  colName <- as.name(substitute(colName))
  
  # Then, we make y delimitations (max and min) for the boxes. The height of the
  # boxes has to be bigger than the maximal value and minimal value of the
  # considered column so we add/remove 20% to its value
  max_col <- (To_plot_manhat %>% 
                select(as.symbol(colName)) %>% max(na.rm = TRUE)) * 1.2
  min_col <- To_plot_manhat %>% 
    select(as.symbol(colName)) %>% min(na.rm = TRUE)
  # Here, we just consider the case where the value of the minimum is positive
  # or negative
  min_col <- min_col * ifelse(min_col >= 0, 0.8, 1.2)
  
  
  # Then, we add this information as well as the x delimitations of the
  # chromosomes to the data frame to return 
  boxes <- To_plot_manhat  %>% 
    # First, we merge the data with the data frame containing the center of the
    # chromosomes to acess their center position
    right_join(chromosome_centers, by="Chromosome") %>% 
    # We group the data frame by the chromosome to treat each chromosome
    # separately
    group_by(Chromosome) %>% 
    # Then, we calculate the max and min positions on each chromosome
    summarise(Min = bp_cum %>% min,
              Max = bp_cum %>% max) %>% 
    # As the geom_polygon function requires to have all the x values in the same
    # column, we pivot the data frame to regroup the max and min value per 
    # chromosome
    pivot_longer(cols = c(Max, Min), values_to = "x_value", names_to = "operation") %>% 
    # The column "operation" contains the names max an min, so we get rid of it
    # because it does not bring any new information
    select(-operation) %>% 
    # As the geom_polygon function requires to have the positions of all the
    # polygon angles on both x and y axis, this means that we need to have the 
    # x/y positions of the first angle, followed by the one next to it, and so
    # on until you run out of angles. So, as we are building a rectangle, we need
    # 4 angles, so 4 positions along the x axis and 4 positions along the y axis.
    # Therefore, we duplicate the lines we already created for each chromosome
    rbind(., .) %>% 
    # We arrange the data frame by Chromosome and x value
    arrange(Chromosome, x_value) %>% 
    # And finally, we add the y values that need to have this order to be able
    # to build the polygons.
    cbind("y_value" = rep(c(min(c(0, min_col)),
                            max(c(1.02, max_col)),
                            max(c(1.02, max_col)),
                            min(c(0, min_col))),
                          nrow(.)/4))
  
  # Return the created variables
  return(boxes)
}

geom_manhattan <- function(df, mapping, thresholding = FALSE, absolute = TRUE, palette = c("grey71","orange2"), filter_low_freqs = TRUE, ...){
  #' Draw manhattan plot
  #' 
  #' This function traces a manhattan plot relying on the ggplot2 aesthetics.
  #' Some additional arguments have been added to be able to personalise the plots
  #' to the best.
  #' 
  #'@param df (data.frame).
  #'    This data frame contains the information to plot. It requires at least
  #'    two columns: the position in the genome and the column to plot along the
  #'    genome.
  #'@param mapping (mapping object).
  #'    This is the mapping of the graph just like in the ggplot2 library. The
  #'    mapping in this function does not require an x column because the position
  #'    along the genome is taken by default.
  #'@param thresholding (boolean). (default = FALSE)
  #'    This is an argument that changes the output of the function. It outputs
  #'    the graph and a maximum cumulative position along the whole genome. It is
  #'    to be used with the "thresholds_manhattan" function.
  #'@param absolute (boolean). (default = TRUE)
  #'    This arguments simply says if you want to plot the absolute values in the
  #'    selected column or if you plot the real value
  #'@param palette (vector). (default = c("grey71", "orange2"))
  #'    This argument is the color palette to use to distinguish between the 
  #'    successive chromosomes. 
  #'@param ...
  #'    In these arguments, you can add any arguments that you would give a 
  #'    ggplot2 graph outside of the aesthetics.
  #'    
  #'@returns (ggplot2 object)
  #'    This returns the manhattan plot of the required column values along the
  #'    genome
  #'@returns (numeric value)
  #'    If the "thresholding" argument is TRUE, then this function also returns
  #'    the maximum cumulative position in the whole genome.
  #'@section Warning:
  #' This function has been made to deal only with two-colored palettes.
  #' @export
  
  # First, we added some error checking to keep the function from running for
  # nothing
  if ("Position" %!in% names(df)){
    stop("This function needs a column called 'Position' that contains the position of each SNP on the genome")
  }
  if (is.numeric(df$Position[1])){
    if ("Chromosome" %!in% names(df)){
      stop("If the positions are already usable, please name a column 'Chromosome' to be used as a chromosome reference.")
    }
  }
  if ("y" %!in% names(mapping)){
    stop("This function requires a y aesthetic.")
  }
  
  # We store the original mapping because the mapping will be modified in this 
  # function, so we need to be able to compare the modified to the original one
  mapping_ori <- mapping
  
  # Then, we extract the name of the column that we have to represent along the
  # y axis.
  colName <- mapping$y[2] %>% as.character
  # And we transform it into a name so that we can use it with the tidyverse.
  colName <- as.name(substitute(colName))
  
  # Now that all the argument importation has been done, we are going to select
  # the columns in the input dataframe that we want to keep in the analysis. To 
  # do this, we make a list of all the parameters that were called to use in the 
  # "matches" function of select
  # We initialise the string "to_select" to an empty character that will be
  # complementeds
  to_select <- ""
  # We iterate over the names of the arguments passed to the function (except
  # for the first one which is the name of the datafame)
  for (i in 1:length(mapping)){
    # We differentiate the first one because there will be no "|" character
    # before it
    if (to_select == ""){
      to_select <- quo_name(mapping[[i]])
    }else{
      # Add the name of the parameters from the mapping to keep
      to_select <- paste0(to_select, "|", quo_name(mapping[[i]]))
    }
  }
  
  # Then, we separate cases where the "Position" column is already numeric
  # (position on one chromosome for example) from the case where the "Position"
  # column is under the adegenet format.
  if (is.numeric(df$Position[1])){
    # Now, we keep only the columns of interest in the dataframe
    To_plot <- df %>% 
      # We keep the columns of interest in the data frame
      select(Position, Chromosome, matches(to_select)) %>%
      # We re-order the levels of the "Chromosome" variable to have them in 
      # increasing order (from 1 to n)
      mutate(Chromosome = Chromosome %>% 
               factor(levels = df %>% select(Chromosome) %>% 
                        unique %>% 
                        arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                        as.vector %>% unname %>% unlist)) %>% 
      # We remove missing values if there are any
      drop_na()
  }else{
    # If the "Position" column is not numeric, it is likely they are at the 
    # adegent format
    To_plot <- df %>%
      # So, we use the created function to separate the position on the
      # chromosome and the chromosome name
      transform_position_ade2tidy() %>% 
      # Then, we re-order the levels of the "Chromosome" variable to have them in 
      # increasing order (from 1 to n)
      mutate(Chromosome = Chromosome %>% 
               factor(levels = df %>%
                        transform_position_ade2tidy() %>%
                        select(Chromosome) %>% 
                        unique %>% 
                        arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                        as.vector %>% unname %>% unlist)) %>% 
      # And we select the columns of interest
      select(Position, Chromosome, matches(to_select)) %>% 
      # Finally, we drop the missing values
      drop_na()
  }
  
  # Then, we make a cumulative data frame with the positions of each end and beginning of chromosome
  data_cum <- To_plot %>% 
    group_by(Chromosome) %>% 
    summarise(max_bp = Position %>% max) %>% 
    # Here, we make a new column that contains the maximum position of the
    # previous chromosome to add this value to the position of the positions of
    # said chromosome
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(Chromosome, bp_add)
  
  # We join it to the filtered dataframe
  To_plot_manhat <- To_plot %>% 
    inner_join(data_cum, by="Chromosome") %>%
    # And we add the cumulative position for each chromosome along the genome
    # to place each chromosome one after the other
    mutate(bp_cum = Position + bp_add)

  # Then, we find the center of the chromosomes
  chromosome_centers <- To_plot_manhat %>% 
    group_by(Chromosome) %>% 
    # The values of the centers of the chromosomes are approximated using the
    # mean function
    summarise(center = bp_cum %>% mean)
  
  # We select two colors that will be used to distinguish chromosomes
  color_chromosome <- rep(palette, To_plot$Chromosome %>% unique %>% length)
  
  # Once all this is done, we create a vector of boolean to use as indicator to
  # take or not the absolute value of the column to plot
  absolute_list <- rep(absolute, nrow(To_plot_manhat))
  
  # Finally, we mutate, if needed, the column to plot with the absolute value
  To_plot_manhat <- To_plot_manhat  %>% 
    right_join(chromosome_centers, by="Chromosome") %>% 
    mutate(!!as.symbol(colName) := ifelse(absolute_list == TRUE, !!as.symbol(colName) %>% abs, !!as.symbol(colName)))
  if (filter_low_freqs){
    # We also filter some values to lighten the plot a little bit
    To_plot_manhat <- To_plot_manhat %>% 
      filter(!!as.symbol(colName) %>% abs > 0.05)
  }

  
  # Once this is done, we create the architecture of the plot (x axis, name of
  # the axis and the theme to use)
  p <- ggplot(data = To_plot_manhat, aes(x = bp_cum)) +
    # The x scale we use is just to differenciate the chromosomes (indicated by
    # their number)
    scale_x_continuous(labels = ((chromosome_centers$Chromosome %>%
                                    str_split_fixed(., "_", 2))[, 1] %>%
                                   str_split_fixed(., "G", 2))[, 2],
                       breaks = chromosome_centers$center) +
    labs(x = "Chromosomes",
         y = colName) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          text = element_text(size = 20))
  
  # We are going to build a string containing the function to call in the ggplot.
  # To do this, we are going to make a list of the default values for aesthetic
  # parameters and we iterate over them.
  list_default_parameters <- list("y" = ".", "colour" = "Chromosome", "size" = 2)
  # We initialise the function to the geom_point function with the mapping we are
  # going to use
  function_to_call <- "geom_point(aes(mapping)"
  for (param in names(list_default_parameters)){
    default_param <- list_default_parameters[[param]]
    # For each parameter, if it is not in the mapping, we add it to the function
    # to call 
    if (param %!in% names(mapping)){
      if (param != "colour"){
        function_to_call <- paste0(function_to_call, ", ", param, " = ", default_param)
      }else{
        # The only exception to this is that if there is no colour mapping, we
        # add one to distinguish the chromosomes
        mapping <- c(mapping, aes(color = Chromosome %>% as.factor))
        class(mapping) <- "uneval"
      }
    }else{
      # If the parameters are in the mapping, we add their names to the labels
      # of the plot
      p$labels[[param]] <- mapping[[param]] %>% quo_name
    }
  }
  
  # If there are some supplementary ggplot arguments they are added to the plot 
  if (length(list(...))){
    function_to_call <- paste0(function_to_call, ", ...")
  }
  # The function to call variable is completed and closed here.
  function_to_call <- paste0(function_to_call, ", inherit.aes = TRUE)")
  if ("colour" %in% names(mapping_ori)){
    # If there is a colour mapping in the plot, we can not distinguish
    # chromosomes by coloring the points of the plot by chromosome, so we create
    # polygons in the background that will be coloured differently.
    # Fisrt, we isolate the name of the colour column
    color <- mapping$colour[2] %>% as.character
    # And we trasform it into a name to use it in the tidyverse
    color <- as.name(substitute(color))
    
    # We use the geom_box_background function to get the data frame of polygon
    # delimitations
    boxes <- geom_box_background(To_plot_manhat, colName, chromosome_centers)
    
    # We add the polygons to the graph
    p <- p +
      # First, the polygons so they are in the background
      geom_polygon(data = boxes, aes(x = x_value, y = y_value, fill = Chromosome), color = NA) +
      # Then, we chose how to color the chromosomes
      scale_fill_manual(values = color_chromosome, guide="none")
    
    # We transform the beginning of the funnction to call variable so that it does
    # not inherit the table boxes as argument
    function_to_call <- function_to_call %>% str_replace("geom_point\\(", "geom_point\\(data = To_plot_manhat, ")
    
    # Then, we use a continuous colour scale if we have a continuous variable and a 
    # discrete colour scale if we have discrete values.
    if (is.continuous(To_plot_manhat[[as.character(color)]])){
      p <- p +
        # We add the function to call to the graph
        eval(parse(text = function_to_call)) +
        # We use a continuous colour scale
        scale_color_gradientn(name = as.character(color),
                              colors = c("darkorchid4", "darkorchid", "mediumorchid1", "magenta"))
    }else{
      # In the case where the colour column is discrete, we have to separate the 
      # case where the colour column is the same as the one we want to plot the 
      # values of along the genome form the case where they are different.
      if (color == colName){
        # If the two columns are the same, we have to separate them to transform
        # the colour column into factors. So, we create a duplicate column
        To_plot_manhat <- To_plot_manhat %>% 
          mutate(color_column = !!as.symbol(color) %>% factor(levels = df %>% 
                                                                select(as.symbol(color)) %>% 
                                                                unique() %>% arrange(!!as.symbol(color)) %>% 
                                                                as.vector %>% unname %>% unlist))
        # In this case, we have to modify the mapping of the graph
        mapping$colour <- quo(color_column)
        # And select some colors to use
        color_polygon_chrom <- To_plot_manhat %>% 
          select(color_column) %>%
          n_distinct() %>% 
          plasma()
      }else{
        # If the columns are different, the process is the same, except we simply
        # use the colour column rather than creating a duplicate column.
        To_plot_manhat <- To_plot_manhat %>% 
          mutate(!!as.symbol(color) := !!as.symbol(color) %>% factor(levels = df %>% 
                                                                       select(as.symbol(color)) %>% 
                                                                       unique() %>% arrange(!!as.symbol(color)) %>% 
                                                                       as.vector %>% unname %>% unlist))
        
        color_polygon_chrom <- To_plot_manhat %>% 
          select(!!as.symbol(color)) %>%
          n_distinct() %>% 
          plasma()
      }
      # We add this to the plot and change the colour scale
      p <- p +
        eval(parse(text = function_to_call)) +
        scale_color_manual(name = as.character(color),
                           values = color_polygon_chrom,
                           drop = FALSE)
    }
    # Finally, we transform the layers of the plot so it can add the required 
    # mapping
    p$layers[[2]]$computed_mapping <- NULL
    p$layers[[2]]$mapping <- mapping
    
    
  }else{
    # If the colour argument is not in the mapping, we simply color using the 
    # chromosomes (one colour per chromosome)
    p <- p +
      # We add the points to the plot
      eval(parse(text = function_to_call)) +
      # We change the color scale
      scale_color_manual(values = color_chromosome, guide="none")
    
    # And in the same way, we modify the layers of the plot so as to use the 
    # required aesthetics
    p$layers[[1]]$computed_mapping <- NULL
    p$layers[[1]]$mapping <- mapping
    
  }
  # Finally, if there is a thresholding (i.e. if the function is embedded in the
  # thresholding function), then, we simply return a list of parameters
  if (thresholding) {
    return(list("plot" = p, "max_value" = To_plot_manhat %>% select(bp_cum) %>% max))
  }else{
    # Otherwise, we return the plot
    return(p)
  }
}

thresholds_manhattan <- function(df, mapping, percentages = NULL, values = NULL, add_text = TRUE, ...){
  #' Draw threshold(s) on a manhattan plot
  #' 
  #' This function re-uses the geom_manhattan function to get a manhattan plot
  #' and it adds threshold values to it (usually represented by horizontal lines)
  #' 
  #'@param df (data.frame).
  #'    This data frame contains the information to plot. It requires at least
  #'    two columns: the position in the genome and the column to plot along the
  #'    genome.
  #'@param mapping (mapping object).
  #'    This is the mapping of the graph just like in the ggplot2 library. The
  #'    mapping in this function does not require an x column because the position
  #'    along the genome is taken by default.
  #'@param percentages (numeric value or vector of numeric values). (default = NULL)
  #'    This contains values between 0 and 1. It will show the given percentages
  #'    of extreme values on the manhattan plot: if we want 20%, a line will 
  #'    appear on the plot showing where the cut of the 20% of highest values 
  #'    are. If this argument is used, do not use the values argument.
  #'@param values (numeric value or vector of numeric values). (default = NULL)
  #'    This contains numeric values in the range of values in the plot. It will
  #'    show the given values on the manhattan plot: if we want a cutoff value at
  #'    5, a vertical line will appear on the y axis at the value of 5.
  #'    If this argument is used, do not use the percentages argument.
  #'@param ...
  #'    supplementary arguments to pass to the geom_manhattan function
  #'    
  #'@returns (ggplot2 object).
  #'    It will return the same plot as the geom_manhattan function, except some
  #'    vertical lines will have appeared at the given values or percentages.
  #' @export
  
  # First, we make a security to be sure that some cutoff values are given.
  if (is_null(percentages) & is_null(values)){
    stop("At least one percentage or value threshold should be given")
  }
  
  # Then, depending on the used argument (values or percentages), we use define 
  # a variable containing these values
  if (is_null(percentages)){
    thresholds <- values
  }else{
    thresholds <- percentages
  }
  
  # We isolate and transform the column to represent on the y axis to be used in
  # the tidyverse
  colName <- mapping$y[2] %>% as.character
  colName <- as.name(substitute(colName))
  
  # We run the geom_manhattan function with the thresholding argument to TURE.
  # We also attribute the plot and the maximum cumulative position along the
  # genome to variables. 
  c(p, max_value) %<-% geom_manhattan(df, mapping, thresholding = TRUE, ...)
  
  # We get as many colours as there are threshold values
  colors_thresholds <- plasma(n = length(thresholds))
  
  # Once this is done, iterate over the threshold values to represent a
  # horizontal line for each one
  for (threshold in thresholds){
    # we get the index of the considered threshold value in this iteration
    i <- which(thresholds == threshold)
    # Select the y position of the line and the text on the graph
    if (is_null(percentages)){
      # If we are dealing with values that are not percentages, the y value of
      # the line is the given value
      y_value_line <- threshold
      # We also add a label to the line to differenciate them. It also has a y
      # position. We place it higher than the horizontal line.
      y_value_text <- threshold * 1.5
    }else{
      # If we are working with percentages, we have to get the top percentages
      # of the distribution of considered values
      y_value_line <- df %>% 
        # We select the column to work on in the input data frame
        select(colName) %>% 
        # We take the absolute value of the column
        mutate(Absolute = (!!as.symbol(colName)) %>% abs) %>% 
        # We arrange the absolute value column to have the biggest values up top
        arrange(desc(Absolute)) %>% 
        # We keep only the first percentage of the table
        dplyr::slice(1:(threshold * nrow(.))) %>% 
        # And we select the minimum value of the subseted ones
        select(Absolute) %>% min
      
      # We do the same thing for the text position, except we put it a little bit
      # higher.
      y_value_text <- df %>% 
        select(colName) %>% 
        mutate(Absolute = (!!as.symbol(colName)) %>% abs) %>% 
        arrange(desc(Absolute)) %>% 
        dplyr::slice(1:(threshold * nrow(.))) %>% 
        select(Absolute) %>% min * 1.5
    }
    # For each value, we append the plot with the new line and text
    p <- p +
      # We add the horizontal line
      geom_hline(yintercept = y_value_line, color = "black", lwd = 1.2, lty = "dashed")
    
    if (add_text){
      p <- p +
        # We add the text above
        annotate("text", x = max_value * 0.975, y = y_value_text,
                 label = ifelse(is_null(percentages), threshold %>% as.character,
                                paste0((threshold * 100) %>% as.character, "%")),
                 size = 20/.pt)
      
    }
  }
  # Once we iterated over all the possible cutoff values, we return the plot
  return(p)
}

positions_to_polygons <- function(df, columns = c("Start", "End"), values = c(0, 1)){
  #' Positions_to_polygons
  #' 
  #' This function transforms a dataframe with the positions of the inversions
  #' to be used to plot them as polygons on the manhattan plot
  #' 
  #'@param df (data.frame).
  #'    This dataframe must contain the columns Start and End indicating the
  #'    positions of the inversion breakpoints. It should also have a column "Inversion"
  #'    specifying the name of the inversion.
  #'    If the argument columns is specified, there is no need for the columns
  #'    Start and End.
  #'    
  #'@param columns (vector) (default = c("Start", "End")).
  #'    This vector is used to specify the names of the columns in which the
  #'    breakpoints of each inversion is specified.
  #'    
  #'@param values (vector) (default = c(0, 1)).
  #'    This vector is used to specify the height of the polygons. It can be adjusted
  #'    with one or two values. If only one value is specified, the other
  #'    default value is 0.
  #' 
  #'@returns (data.frame).
  #'    This dataframe contains all the SNPs that have at least one clinal model.
  #'    If the SNP is present several times in the dataframe, it will select all
  #'    the lines containing said SNP if at least one model is clinal.
  #' @export
  if (length(values) > 2){
    stop("The values argument has to be of length 1 or 2 to get ajustable delimitations for your manhattan plot")
  }
  # Check if the data frame is grouped or not
  is_grouped <- is_grouped_df(df)
  # Transform the dataframe to be able to use it
  df <- df %>% 
    # Summarise the positions for split inversions if needed
    group_by(Chromosome, Inversion) %>% 
    summarize(Start = min(Start, na.rm = TRUE),
              End = max(End, na.rm = TRUE)) %>% 
    ungroup %>% 
    inner_join(df %>% 
                 select(-c(Start, End, Length)),
               by = c("Chromosome", "Inversion")) %>% 
    arrange(Population) %>% 
    unique
  if (is_grouped){
    df <- df %>% 
      group_by(Population)
  }
  df <- df %>% 
    # Pivot longer to have the positions of the breakpoints in one column
    pivot_longer(cols = all_of(columns), names_to = "Breakpoints", values_to = "Position") %>% 
    # Double the dataframe to have twice the positions of the x value (necessary
    # to plot the polygons in ggplot2)
    rbind(., .) %>% 
    # Regroup every inversion and put the same positions consecutively
    arrange(Inversion, Position)
  # Add the y values that can be user-specified
  if (length(values) == 1){
    df <- df %>%
      mutate(Height_polygon = rep(c(0, values, values, 0), n()/4))
  }else{
    df <- df %>%
      mutate(Height_polygon = rep(c(values, rev(values)), n()/4))
  }
  # Return the dataframe
  return(df)
}

################################################################################
#################### 5. Cline fitting and visualisation   ######################
################################################################################

optimise_clines <- function(Priors, logarithm = FALSE, batch_size = 1000, write_output = NULL, ...){
  #' Optimise cline models
  #' 
  #' This function is made to optimise the best of three models (stable, linear
  #' and clinal) along a transect. It gives the best optimised parameters for
  #' all three models as well as the AIC of each model and which is better to use
  #' 
  #'@param Priors (data.frame).
  #'    This data frame contains the prior values of the parameters to optimise.
  #'    This data frame has to contain at least a "Population" and "Position" 
  #'    column as well as some columns with the allelic frequencies on both ends
  #'    of the transect, the center of a cline and the width of the cline.
  #'@param logarithm (boolean). (default = FALSE)
  #'    This arguemnt is used to know if the optimisation process should be done
  #'    using logarithm transformations for the cline model fitting (TRUE) or 
  #'    not (FALSE)
  #'@param batch_size (integer). (default = 1000)
  #'    How many rows to sub-sample from the Priors data frame to run in parallel.
  #'    If the batch is too big, it will be longer. If it is smaller, it will
  #'    do more iterations.
  #'@param write_output (string (defaul = NULL)
  #'    Where to write the output if it is a long process
  #'@param ...
  #'    Supplementary arguments to pass to the get_gentoype_transect function.
  #'    See this function to know which arguments to use.
  #' 
  #'@returns data.frame
  #'    This data frame contains the optimised values of the three models for 
  #'    each position in the input data frame. It also tests the models on their
  #'    AIC and selects the best model to use to represent it.
  #' @export
  
  # First, we need to get the genotypes of the individuals along the transect to
  # optimise and test the models. It is recommended to use a subset of SNPs as 
  # this function runs for a long time
  genotype_transect <- get_genotype_transect(...)
  
  # This determines which function to use for the optimisation of the cline model.
  # If we want to do a logarithm transformation, we use the clineflog function. 
  # Otherwise, we use the clinef function
  function_to_use <- ifelse(logarithm, clineflog, clinef)
  
  # Then, we slighly modify the priors to be used in the optimisation process
  Priors_func <- Priors %>% 
    # Then, we rename the population and position columns to not get confused
    rename(Pop = Population,
           Pos = Position) %>% 
    # We separate every SNP from every population to optimise it separately
    group_by(Pop, Pos) %>% 
    # Then, we modify parameters to be used, if necessary with the log transform
    # function
    mutate(
      # The three following lines do the same thing for different parameters: 
      # if we want to do a logarithm transformation, we transform the width to be
      # the log of the width. We do the same transformation for the max and min
      # priors
      Width_prior = ifelse(logarithm, Width_prior %>% log, Width_prior),
      Width_max = ifelse(logarithm, Width_max %>% log, Width_max),
      Width_min = ifelse(logarithm, Width_min %>% log, Width_min),
      # Then, as the cline optimisation functions work best with the allelic
      # frequency on the left being smaller than the one on the right side of
      # the transect, in the case where the allelic frequencies are reversed 
      # p_left > p_right), we want to take the opposite allele of frequency 1-p_left
      # and 1-p_right on both ends. So, to do this, first we check if the allelic
      # frequencies are reversed ...
      reversed = p_left_shelt > p_right_expos,
      # ... and if they are, we take the opposite allele's frequency
      p_left_shelt = ifelse(reversed, 1 - p_left_shelt, p_left_shelt),
      p_right_expos = ifelse(reversed, 1 - p_right_expos, p_right_expos),
      # Then, the logarithm transformation of the priors to do the log optimisation
      # requires that the allelic frequencies be transformed using the logit function.
      # So, in the same way as for the width, we take the logit of the allelic
      # frequencies if we want to do the logarithm transformation.
      # Wa also change fixed allele frequencies to be able to transform them with
      # the logit function
      p_left_shelt = ifelse(p_left_shelt == 1, 0.999999,
                            ifelse(p_left_shelt == 0, 0.000001, p_left_shelt)),
      p_left_max = ifelse(logarithm, min(p_left_shelt + 0.1, 0.999999) %>% logit,
                          min(p_left_shelt + 0.1, 0.999999)),
      p_left_min = ifelse(logarithm, max(p_left_shelt - 0.1, 0.000001) %>% logit,
                          max(p_left_shelt - 0.1, 0.000001)),
      p_left_shelt = ifelse(logarithm, p_left_shelt %>% logit, p_left_shelt),
      p_right_expos = ifelse(p_right_expos == 1, 0.999999,
                             ifelse(p_right_expos == 0, 0.000001, p_right_expos)),
      p_right_max = ifelse(logarithm, min(p_right_expos + 0.1, 0.999999) %>% logit,
                           min(p_right_expos + 0.1, 0.999999)),
      p_right_min = ifelse(logarithm, max(p_right_expos - 0.1, 0.000001) %>% logit,
                           max(p_right_expos - 0.1, 0.000001)),
      p_right_expos = ifelse(logarithm, p_right_expos %>% logit, p_right_expos)
    )
    
    # To win some time, we will iterate over all the lines in the table by 
  number_of_iterations <- ((nrow(Priors_func)/batch_size) %>% floor) + 1
  Comp_table <- data.frame()
  neutral_models <- c("Stable", "Linear")
  for (i in 1:number_of_iterations){
    # Print the progress bar to see where we are in the loop
    progress_bar(iteration_number = i, nb_iterations_to_do = number_of_iterations)
    
    # Then, we sample 10 values from the Priors_func for the clinal model
    subset_priors_func <- Priors_func %>% 
      head(n = batch_size * i) %>% 
      tail(n = batch_size)
    # And we do the same from the Priors for the sable and linear models
    subset_priors <- Priors %>% 
      head(n = batch_size * i) %>% 
      tail(n = batch_size)
    
    # Here, we make a vector of booleans to use in the tranformation of the priors
    # if we are using a logarithm transformation or not.
    list_size_fr <- subset_priors_func %>% 
      filter(Pop == "France") %>% 
      nrow()
    list_size_sw <- subset_priors_func %>% 
      filter(Pop == "Sweden") %>% 
      nrow()
    logarithm_list_fr <- rep(logarithm, list_size_fr)
    logarithm_list_sw <- rep(logarithm, list_size_sw)
    
      
    # Now that the priors are ready, we can start the optimisation process
    Clinal_model_i <- subset_priors_func %>% 
      # We separate every SNP from every population
        group_by(Pop, Pos) %>% 
        # And summarise the optimised parameters we are interested in
        bow(tie(Centre, Width, Left, Right) := mle2(function_to_use,
                                                    # Here, we give the priors to use
                                                    # for every parameter we want
                                                    list(centre = Centre_prior,
                                                         width = Width_prior ,
                                                         left = p_left_shelt,
                                                         right = p_right_expos),
                                                    # In the data part, we give the 
                                                    # position of each individual along
                                                    # the transect and the genotype
                                                    # of each individual for the chosen
                                                    # SNP in the given population
                                                    data = list(x = genotype_transect %>% 
                                                                  filter(Population == Pop) %>% 
                                                                  select(LCmeanDist) %>% 
                                                                  drop_na %>% 
                                                                  as.vector %>%
                                                                  unlist %>%
                                                                  unname,
                                                                g = genotype_transect %>% 
                                                                  filter(Population == Pop) %>%
                                                                  select(Pos) %>% 
                                                                  rename(Geno = starts_with("LG")) %>% 
                                                                  drop_na %>% 
                                                                  mutate(Geno = ifelse(rep(reversed, nrow(.)), 2 - Geno, Geno)) %>% 
                                                                  select(Geno) %>% 
                                                                  as.vector %>%
                                                                  unlist %>%
                                                                  unname,
                                                                n = 2),
                                                    # Method to use
                                                    method = "L-BFGS-B",
                                                    # Upper limits of the optimisation
                                                    # for the parameters we are trying
                                                    # to optimise
                                                    upper = list(centre = Centre_max,
                                                                 width = Width_max,
                                                                 left = p_left_max,
                                                                 right = p_right_max),
                                                    # lower limits of the optimisation
                                                    # for the parameters we are trying
                                                    # to optimise
                                                    lower = list(centre = Centre_min,
                                                                 width = Width_min,
                                                                 left = p_left_min,
                                                                 right = p_right_min)) %>%
              # We extract the coefficients 
              coef() %>% 
              # and round them to three digits
              round(digits = 3)) %>% 
        # Then, we merge this with the Priors table to know if we used a log 
        # transformation and if the allele frequencies are reversed or not
        left_join(subset_priors_func, by = c("Pop", "Pos")) %>% 
        # We select the columns we want to keep
        select(Pop, Pos, reversed, Centre, Width, Left, Right) %>% 
        # Backtransform the optimised parameters
        mutate(Width = ifelse(ifelse(Pop == "France", logarithm_list_fr, logarithm_list_sw), Width %>% exp, Width),
               Left = ifelse(ifelse(Pop == "France", logarithm_list_fr, logarithm_list_sw), Left %>% invlogit, Left),
               Right = ifelse(ifelse(Pop == "France", logarithm_list_fr, logarithm_list_sw), Right %>% invlogit, Right),
               # Backtransform with the reversed parameter
               Left = ifelse(reversed, 1 - Left, Left),
               Right = ifelse(reversed, 1 - Right, Right)) %>% 
        rename(Position = Pos,
               Population = Pop)
    
    # Once we estimated the clinal parameters, we are now going to estimate the 
    # parameter for the stable model
    Stable_model_i <- subset_priors %>% 
      # As the stable model describes a continuous allele frequency along the 
      # transect, we only need one parameter to optimise and it is the mean allelic
      # frequency in the transect
      mutate(Mean_freq = (p_left_shelt + p_right_expos)/2) %>% 
      # In the same way, we isolate each SNP from each population
      rename(Pop = Population,
             Pos = Position) %>% 
      group_by(Pop, Pos) %>% 
      # And summarise the optimised parameter for each SNP in each population
      summarise(Stable_model_fit = mle2(stable,
                                        # List of parameters to optimise
                                        list(p_all = Mean_freq),
                                        # In the data part, we give the 
                                        # position of each individual along
                                        # the transect and the genotype
                                        # of each individual for the chosen
                                        # SNP in the given population
                                        data = list(x = genotype_transect %>% 
                                                      filter(Population == Pop) %>% 
                                                      select(LCmeanDist) %>% 
                                                      drop_na %>% 
                                                      as.vector %>% unlist %>% unname,
                                                    g = genotype_transect %>% 
                                                      filter(Population == Pop) %>%
                                                      select(Pos) %>% 
                                                      drop_na %>% 
                                                      as.vector %>% unlist %>% unname,
                                                    n = 2),
                                        # Method to use
                                        method = "L-BFGS-B",
                                        # Upper bound of the estimated parameter
                                        upper = list(p_all = 0.999999),
                                        # lower bound of the estimated parameter
                                        lower = list(p_all = 0.000001)) %>% 
                  # We extract the coefficients and round them to three digits
                  coef() %>% round(digits = 3),
                .groups = "drop_last") %>% 
      rename(Position = Pos,
             Population = Pop)
    
    # Once we estimated the stable parameters, we are now going to estimate the 
    # parameters for the linear model
    Linear_model_i <- subset_priors %>% 
      # As the linear model fits best for p_left < p_right, we change the allelic
      # frequencies if they are reversed (p_left > p_right).
      mutate(reversed = p_left_shelt > p_right_expos,
             p_left_shelt = ifelse(reversed, 1 - p_left_shelt, p_left_shelt),
             p_right_expos = ifelse(reversed, 1 - p_right_expos, p_right_expos)) %>% 
      # Again, we separate each SNP from each Population
      rename(Pos = Position,
             Pop = Population) %>% 
      group_by(Pop, Pos) %>% 
      # And we summarise the parameters we want for the linear model
      bow(tie(p_left, p_right) := mle2(linear,
                                       # List of parameters to optimise
                                       list(p_left = p_left_shelt,
                                            p_right = p_right_expos),
                                       # In the data part, we give the 
                                       # position of each individual along
                                       # the transect and the genotype
                                       # of each individual for the chosen
                                       # SNP in the given population
                                       data = list(x = genotype_transect %>% 
                                                     filter(Population == Pop) %>% 
                                                     select(LCmeanDist) %>% 
                                                     drop_na %>% 
                                                     as.vector %>% unlist %>% unname,
                                                   g = genotype_transect %>% 
                                                     filter(Population == Pop) %>%
                                                     select(Pos) %>% 
                                                     drop_na %>% 
                                                     as.vector %>% unlist %>% unname,
                                                   n = 2),
                                       # Method to use
                                       method = "L-BFGS-B",
                                       # Upper bounds of the parameters to optimise
                                       upper = list(p_left = 0.999999,
                                                    p_right = 0.999999),
                                       # Lower bounds of the parameters to estimate
                                       lower = list(p_left = 0.000001,
                                                    p_right = 0.000001)) %>%
            # We get the coefficient of the parameters we want and round them to 
            # three digits
            coef() %>% round(digits = 3)) %>% 
      rename(Population = Pop,
             Position = Pos) %>% 
      # and we backtransform the reversion
      left_join(subset_priors %>% 
                  mutate(reversed = p_left_shelt > p_right_expos) %>% 
                  select(Population, Position, reversed),
                by = c("Population", "Position")) %>% 
      mutate(p_left = ifelse(reversed, 1 - p_left, p_left),
             p_right = ifelse(reversed, 1 - p_right, p_right))
    
    # The last step is to merger all the estimated parameters together and calculate
    # the AICs of the models using the estimated parameters.
    Comp_table_i <- Clinal_model_i %>% 
      # First, we merge the table of the estimated parameters together
      left_join(Stable_model_i, by = c("Population", "Position")) %>% 
      left_join(Linear_model_i, by = c("Population", "Position", "reversed")) %>% 
      left_join(subset_priors_func %>% 
                  rename(Population = Pop,
                         Position = Pos), by = c("Population", "Position", "reversed")) %>% 
      # Again, we separate each SNP of each population
      rename(Pop = Population,
             Pos = Position) %>% 
      group_by(Pop, Pos) %>% 
      # And we calculate the AIC of each model using the same method as in the 
      # optimisation process. For the clinal model, this suggests we need to do the
      # transformation to logarithm values beforehand if necessary
      mutate(
        # Width
        Width = ifelse(logarithm, Width %>% log, Width),
        Width_max = ifelse(logarithm, Width_max %>% log, Width_max),
        Width_min = ifelse(logarithm, Width_min %>% log, Width_min),
        # check if the allelic frequency are reversed
        Left = ifelse(reversed, 1 - Left, Left),
        Right = ifelse(reversed, 1 - Right, Right),
        p_left = ifelse(reversed, 1 - p_left, p_left),
        p_right = ifelse(reversed, 1 - p_right, p_right),
        # Then, the logarithm transformation of the priors to do the log optimisation
        # requires that the allelic frequencies be transformed using the logit function.
        # So, in the same way as for the width, we take the logit of the allelic
        # frequencies if we want to do the logarithm transformation.
        # Wa also change fixed allele frequencies to be able to transform them with
        # the logit function
        Left = ifelse(Left == 1, 0.999999,
                      ifelse(Left == 0, 0.000001, Left)),
        p_left_max = ifelse(logarithm, min(Left + 0.1, 0.999999) %>% logit,
                            min(Left + 0.1, 0.999999)),
        p_left_min = ifelse(logarithm, max(Left - 0.1, 0.000001) %>% logit,
                            max(Left - 0.1, 0.000001)),
        Left = ifelse(logarithm, Left %>% logit, Left),
        Right = ifelse(Right == 1, 0.999999,
                       ifelse(Right == 0, 0.000001, Right)),
        p_right_max = ifelse(logarithm, min(Right + 0.1, 0.999999) %>% logit,
                             min(Right + 0.1, 0.999999)),
        p_right_min = ifelse(logarithm, max(Right - 0.1, 0.000001) %>% logit,
                             max(Right - 0.1, 0.000001)),
        Right = ifelse(logarithm, Right %>% logit, Right),
      ) %>% 
      # Then, we can calculate the AICs
      ## First for the stable model
      mutate(AIC_stable = mle2(stable, list(p_all = Stable_model_fit),
                               data = list(x = genotype_transect %>% 
                                             filter(Population == Pop) %>% 
                                             select(LCmeanDist) %>% 
                                             drop_na %>% 
                                             as.vector %>% unlist %>% unname,
                                           g = genotype_transect %>% 
                                             filter(Population == Pop) %>%
                                             select(Pos) %>% 
                                             drop_na %>% 
                                             as.vector %>% unlist %>% unname,
                                           n = 2),
                               method = "L-BFGS-B",
                               upper = list(p_all = 0.999999),
                               lower = list(p_all = 0.000001)) %>% 
               AIC,
             ## Then for the linear model
             AIC_linear = mle2(linear, list(p_left = p_left,
                                            p_right = p_right),
                               data = list(x = genotype_transect %>% 
                                             filter(Population == Pop) %>% 
                                             select(LCmeanDist) %>% 
                                             drop_na %>% 
                                             as.vector %>% unlist %>% unname,
                                           g = genotype_transect %>% 
                                             filter(Population == Pop) %>%
                                             select(Pos) %>% 
                                             drop_na %>% 
                                             as.vector %>% unlist %>% unname,
                                           n = 2),
                               method = "L-BFGS-B",
                               upper = list(p_left = 0.999999,
                                            p_right = 0.999999),
                               lower = list(p_left = 0.000001,
                                            p_right = 0.000001)) %>% 
               AIC,
             ## And finally for the clinal model
             AIC_clinal = mle2(function_to_use,
                               list(centre = Centre,
                                    width = Width ,
                                    left = Left,
                                    right = Right),
                               data = list(x = genotype_transect %>% 
                                             filter(Population == Pop) %>% 
                                             select(LCmeanDist) %>% 
                                             drop_na %>% 
                                             as.vector %>% unlist %>% unname,
                                           g = genotype_transect %>% 
                                             filter(Population == Pop) %>%
                                             select(Pos) %>% 
                                             rename(Geno = starts_with("LG")) %>% 
                                             drop_na %>% 
                                             mutate(Geno = ifelse(rep(reversed, nrow(.)), 2 - Geno, Geno)) %>% 
                                             select(Geno) %>% 
                                             as.vector %>% unlist %>% unname,
                                           n = 2),
                               "L-BFGS-B",
                               upper = list(centre = Centre_max,
                                            width = Width_max,
                                            left = p_left_max,
                                            right = p_right_max),
                               lower = list(centre = Centre_min,
                                            width = Width_min,
                                            left = p_left_min,
                                            right = p_right_min)) %>%
               AIC) %>% 
      # Then, we need to backtransform the width and allelic frequencies
      mutate(Width = ifelse(logarithm, Width %>% exp, Width),
             Left = ifelse(logarithm, Left %>% invlogit, Left),
             Right = ifelse(logarithm, Right %>% invlogit, Right),
             # Backtransform with the reversed parameter
             Left = ifelse(reversed, 1 - Left, Left),
             Right = ifelse(reversed, 1 - Right, Right),
             p_left = ifelse(reversed, 1 - p_left, p_left),
             p_right = ifelse(reversed, 1 - p_right, p_right)) %>% 
      # We remove columns we do not need in the output
      select(-c(p_left_shelt, p_right_expos, contains("prior"), matches("max|min"))) %>%
      # Then, we choose which model is the best between the three estimated ones 
      # using the AICs
      mutate(Delta_AIC_stable = AIC_clinal - AIC_stable, Delta_AIC_linear = AIC_clinal - AIC_linear,
             Best_model = case_when(
               (Delta_AIC_linear < 0 & Delta_AIC_stable < 0) ~ "Clinal",
               (Delta_AIC_stable >= 0) & (Delta_AIC_linear < Delta_AIC_stable) ~ "Stable",
               (Delta_AIC_linear >= 0) & (Delta_AIC_stable < Delta_AIC_linear) ~ "Linear",
               TRUE ~ NA
             ),
             Best_model_AIC = ifelse(Best_model == "Stable", AIC_stable, ifelse(Best_model == "Linear", AIC_linear, AIC_clinal)),
             Second_best_model = ifelse(Best_model == "Stable" & Delta_AIC_linear > 0, "Linear",
                                        ifelse(Best_model == "Stable" & Delta_AIC_linear < 0, "Clinal",
                                               ifelse(Best_model == "Linear" & Delta_AIC_stable < 0, "Clinal",
                                                      ifelse(Best_model == "Linear" & Delta_AIC_stable > 0, "Stable",
                                                             ifelse(Best_model == "Clinal" & AIC_linear > AIC_stable, "Stable", "Linear"))))),
             Delta_AIC_second_best_model = ifelse(Best_model == "Clinal" & Second_best_model == "Linear", AIC_clinal - AIC_linear,
                                                  ifelse(Best_model == "Clinal" & Second_best_model == "Stable", AIC_clinal - AIC_stable,
                                                         ifelse(Best_model == "Linear" & Second_best_model == "Clinal", AIC_linear - AIC_clinal,
                                                                ifelse(Best_model == "Linear" & Second_best_model == "Stable", AIC_linear - AIC_stable,
                                                                       ifelse(Best_model == "Stable" & Second_best_model == "Clinal", AIC_stable - AIC_clinal, AIC_stable - AIC_linear))))),
             Significant = case_when(
               (abs(Delta_AIC_second_best_model) < 4) & (Second_best_model %in% neutral_models) ~ Second_best_model,
               (abs(Delta_AIC_second_best_model) < 4) & (Best_model %in% neutral_models) ~ Best_model,
               TRUE ~ Best_model
             )) %>%  
      # We remove the AICs
      select(-starts_with("AIC_")) %>% 
      # Rename the population and position columns
      rename(Position = Pos,
             Population = Pop) %>% 
      ungroup
    if (!is.null(write_output)){
      if (!file.exists(write_output)){
        Comp_table_i %>%
          write.table(write_output, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      }else{
        Comp_table_i %>%
          write.table(write_output, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
    }else{
      # Then, we add this to the general comp table containing all the outputs
      Comp_table <- add_table_to_df_in_iteration(Comp_table, Comp_table_i)
    }
    
    
    
  }
  
  # And return the final table. As the number of iterations is not always a round
  # number, several lines may be duplicated, so we keep only the unique ones
  return(Comp_table %>% 
           unique)
  
  print("Done !")
  
}

plot_clines <- function(cline_params, real_distance = FALSE, goodness_of_fit = FALSE, ...){
  #' Makes information for the cline plotting
  #' 
  #' This function gives a data frame containing all the information to plot the
  #' best distribution of allelic frequencies along the transect using three 
  #' models: stable, linear and clinal. 
  #' 
  #'@param cline_params (data.frame).
  #'    This data frame contains the parameters that came out of the "optimise_clines"
  #'    function containing the parameters of the optimised model parameters.
  #'@param real_distance (boolean). (default = FALSE)
  #'    This arguemnt is used to know if the graphic representation should use the
  #'    real position of the individuals along the transect rather than smoothing
  #'    the positions. This argument should be TRUE when the "goodness_of_fit"
  #'    argument is TRUE.
  #'@param goodness_of_fit (boolean). (default = FALSE)
  #'    This argument indicates if it should return the goodness of fit of the 
  #'    plotted model as well as the plotting data for the distributions of 
  #'    allelic frequencies along the transect. This argument should be used with
  #'    the "real_distance" argument set to TRUE.
  #'@param ...
  #'    Supplementary arguments to pass to the "get_gentoype_transect" function.
  #'    See this function to know which arguments to use.
  #' 
  #'@returns (data.frame).
  #'    This data frame contains the parameters to plot the best distribution of 
  #'    allelic frequencies along the transect between three models. If the 
  #'    "goodness_of_fit" argument is also given, then it should return a list
  #'    of two data frames. The first one containing the goodness of fit of the 
  #'    plotted models and then the plotting information.
  #' @export
  
  # First, we need to get the genotypes of the individuals along the transect to
  # optimise and test the models. It is recommended to use a subset of SNPs as 
  # this function runs for a long time
  genotype_transect <- get_genotype_transect(...)
  
  # Here, we initialise the data frames that we are going to fill in when
  # iterating over each SNP in each population.
  Cline_along_transect <- tibble()
  # We also initialise the goodness of fit data frame to add it for each SNP, so
  # that we can do one iteration and do both the plotting and the goodness of fit
  if (goodness_of_fit){
    if (!real_distance){
      # We add an error message if we are missing some arguments
      stop("You need to use the real distance to do a goodness of fit test")
    }
    Goodness_of_fit_reg <- tibble()
  }
  # Then, we iterate over each SNP to get the plotting information and the 
  # goodness of fit if needed.
  # We also make a counter to see the advancement
  i <- 1
  nb_positions <- cline_params$Position %>% unique %>% length
  for (position_i in cline_params$Position %>% unique){
    # Print the progress bar to see where we are in the loop
    progress_bar(iteration_number = i, nb_iterations_to_do = nb_positions)
    # Here, we initialise a temporary data frame to put the population information
    # inside and then bind it to the general data frame.
    temp_data <- tibble()
    for (population in cline_params$Population %>% unique){
      # We isolate the parameters we want to use for the SNP in the selected 
      # population.
      Cline_param_pop_pos_i <- cline_params %>% 
        filter(Population == population,
               Position == position_i)
        
      if (nrow(Cline_param_pop_pos_i) == 0){
        next
      }
      
      # For each SNP in each population, we select the genotype of the SNP for the
      # selected population.
      geno_pop_pos_i <- genotype_transect %>% 
        filter(Population == population) %>% 
        select(Sample_Name, LCmeanDist, position_i %>% matches)
      
      # In the same way as in the optimisation of the models, the cline plotting
      # works best when the allelic frequencies are p_left < p_right, so if it is
      # in the other way, we want to reverse them (1 - p), so we extract for each
      # SNP if it is reversed or not.
      is_reversed <- cline_params %>%
        ungroup %>% 
        filter(Position == position_i,
               Population == population) %>% 
        select(reversed)
      
      Cline_param_pop_pos_i <- Cline_param_pop_pos_i %>% 
        # And if the SNP is reversed (p_right > p_left), then, we reverse it
        # (1 - p)
        mutate(Left = ifelse(is_reversed, 1 - Left, Left),
               Right = ifelse(is_reversed, 1 - Right, Right))
      
      # For each SNP, we extract the model to use (clinal, linear or stable) from
      # the input data frame
      model_to_use <- (cline_params %>%
        ungroup %>% 
        filter(Position == position_i,
               Population == population) %>% 
        select(Best_model))$Best_model
      
      # Select which distance to use (real_distance or a smoothed distance along
      # the transect)
      if (real_distance){
        Distance <- geno_pop_pos_i %>% select(LCmeanDist) %>% as.vector %>% unlist %>% unname
      }else{
        Distance <- seq(from = geno_pop_pos_i %>% select(LCmeanDist) %>% min,
                        to = geno_pop_pos_i %>% select(LCmeanDist) %>% max,
                        by = 1)
      }
      
      # Then, we separate the three different models to get the plotting 
      # information for the model we want to use.
      # First, we check if the model to use is the stable one
      if (model_to_use == "Stable"){
        cline_transect_pop_pos_i <- stable(x = Distance,
                                           p_all = Cline_param_pop_pos_i$Stable_model_fit,
                                           optimisation = FALSE) %>% 
          # We rename the output of the stable function to match the name of the
          # SNP, so we can bind the same SNPs into the same column by separating
          # them by population
          rename(!!quo_name(position_i) := phen_cline,
                 LCmeanDist = position) %>% 
          mutate(Population = population)
      # Then, we check if it is the linear model
      }else if (model_to_use == "Linear"){
        cline_transect_pop_pos_i <- linear(x = Distance,
                                           p_left = Cline_param_pop_pos_i$p_left,
                                           p_right = Cline_param_pop_pos_i$p_right,
                                           optimisation = FALSE) %>% 
          # In the same way, we rename the output of the function so it matches
          # the name of the SNP.
          rename(!!quo_name(position_i) := phen_cline,
                 LCmeanDist = position) %>% 
          mutate(Population = population)
      }else{
        # The last option is that the SNP has a clinal distribution.
        cline_transect_pop_pos_i <- clinef(x = Distance,
                                           centre = Cline_param_pop_pos_i$Centre,
                                           width = Cline_param_pop_pos_i$Width,
                                           left = Cline_param_pop_pos_i$Left,
                                           right = Cline_param_pop_pos_i$Right,
                                           optimisation = FALSE) %>% 
          # In the same way, we rename the output of the function so it matches
          # the name of the SNP.
          rename(!!quo_name(position_i) := phen_cline,
                 LCmeanDist = position) %>% 
          # And if needed, we backtransform the allelic frequency.
          mutate(!!quo_name(position_i) := ifelse(rep(is_reversed, nrow(.)), 1 - !!as.name(position_i), !!as.name(position_i)),
                 Population = population)
      }
      # If we are working with the real distance along the transect, we add the 
      # names of the individuals to identify them
      if (real_distance){
        cline_transect_pop_pos_i <- cline_transect_pop_pos_i %>% 
          arrange(LCmeanDist) %>% 
          cbind(geno_pop_pos_i %>%
                  arrange(LCmeanDist) %>%
                  select(Sample_Name))
      }
      
      # Then, we store the results in a temporary data frame that we fill while
      # iterating over the populations.
      if ((temp_data %>% nrow) == 0){
        # We separate the case where the temporary data frame is empty form the 
        # case where it already has some information.
        temp_data <- cline_transect_pop_pos_i
      }else{
        temp_data <- temp_data %>% 
          rbind(cline_transect_pop_pos_i)
      }
      
      # Then ,if asked, we get the goodness of fit for each regression
      if (goodness_of_fit){
        # First, we regroup the calculated distribution of allelic frequency along
        # the transect and the genotype of the considered SNP
        cline_and_geno_along_transect_pop_pos_i <- cline_transect_pop_pos_i %>% 
          left_join(geno_pop_pos_i, by = c("Sample_Name", "LCmeanDist"), suffix = c("", "_geno")) %>% 
          drop_na %>% 
          # We rename the columns so as to make it easier to identify them
          rename(Genotype := !!quo_name(paste(position_i, "geno", sep="_")),
                 Freq := !!quo_name(position_i))
        # Then, we fit a glm to the predicted frequency to see if it can explain 
        # the observed genotype
        good_fit_pop_pos_i <- cline_and_geno_along_transect_pop_pos_i %>% 
          glm(data = ., (Genotype/2) ~ Freq, family = binomial(link = "logit"))
        
        # Once this is done, we add the name of the population, the position of
        # the SNP, the deviance of the glm, the explained deviance and the
        # difference in allelic frequency between the left and right parts of
        # the transect
        list_to_add_to_table_goodness_of_fit <- c(population, position_i, good_fit_pop_pos_i$deviance,
                                                  100 * (good_fit_pop_pos_i$null.deviance - good_fit_pop_pos_i$deviance) / good_fit_pop_pos_i$null.deviance,
                                                  Cline_param_pop_pos_i$Right - Cline_param_pop_pos_i$Left)
        # Then, we add this to the data frame containing all the goodness of fit 
        # information that we want to return.
        if((Goodness_of_fit_reg %>% nrow) == 0){
          # We separate the case where the data frame is empty from the case where
          # there is already something inside.
          Goodness_of_fit_reg <- list_to_add_to_table_goodness_of_fit %>% 
            matrix(nrow = 1, byrow = TRUE) %>% 
            as_tibble %>% 
            rename(Population = V1,
                   Position = V2,
                   Deviance = V3,
                   Explained_deviance = V4,
                   Delta_freq = V5)
        }else{
          Goodness_of_fit_reg <- Goodness_of_fit_reg %>% 
            rbind(list_to_add_to_table_goodness_of_fit)
        }
      }
    }
    # Once we iterated over the populations, we add this to the overall data 
    # frame containing the allelic frequency distributions along the transect of
    # all the SNPs.
    if ((Cline_along_transect %>% nrow) == 0){
      Cline_along_transect <- temp_data
    }else{
      # If there is already some information in the data frame, we separate the 
      # case where we are dealing with the real distance and the case where we 
      # have a smoothed distance. In these cases, the merging columns are not the
      # same.
      if (real_distance){
        Cline_along_transect <- Cline_along_transect %>% 
          left_join(temp_data, by = c("Population", "Sample_Name", "LCmeanDist"))
      }else{
        Cline_along_transect <- Cline_along_transect %>% 
          left_join(temp_data, by = c("Population", "LCmeanDist"))
      }
    }
    # Increase the position counter
    i <- i + 1
  }
  
  # Once the calculation of all the plotting parameters is done (we finished 
  # iterating over all the SNPs), we make some final adjustments to the 
  # obtained information
  Cline_return <- Cline_along_transect %>% 
    # We reorder the levels of the populations
    mutate(Population = Population %>% 
             factor(levels = c("Sweden", "France"))) %>% 
    group_by(Population) %>% 
    # And we pivot the data frame to only have a Frequency column containing the
    # allelic frequencies and a Position column containing all the names of the 
    # SNPs and their position with the adegenet format.
    pivot_longer(cols = starts_with("LG"),
                 names_to = "Position",
                 values_to = "Frequency") 
  
  
  # Finally, if we want to return the goodness of fit, we do the same modifications
  # to the obtained data for the goodness of fit and we return the necessary outputs.
  if (goodness_of_fit){
    Goodness_return <- Goodness_of_fit_reg %>% 
      mutate(Deviance = Deviance %>% as.numeric,
             Explained_deviance = Explained_deviance %>% as.numeric,
             Delta_freq = Delta_freq %>% as.numeric,
             Population = Population %>% factor(levels = c("Sweden", "France")))
    cat("\n")
    return(list("Goodness_of_fit" = Goodness_return, "Cline" = Cline_return))
  }else{
    cat("\n")
    return(Cline_return)
  }
}


################################################################################
################### 6. Run local PCA on portion of genome   ####################
################################################################################
# Order the clustering to always have the first cluster on the left and the last cluster on the right
reorder_clustering <- function(df){
  # Duplicate the column we are going to work on
  df <- df %>% 
    mutate(gr = Group)
  
  # Get the number of positions that are positive or negative in the reference group
  ref_group <- df %>%
    mutate(Group = Group %>% as.factor) %>% 
    group_by(Group) %>% 
    summarize(Mean_PC1_value = mean(Axis1)) %>% 
    arrange(Mean_PC1_value)
  
  df <- df %>% 
    mutate(gr = ifelse(
      # Smallest value
      Group == ref_group$Group[1], 1,
      ifelse(
        # Medium value
        Group == ref_group$Group[2], 2,
        # Take the last possible value
        3
      )
    ),
    Group = gr) %>%
    dplyr::select(-gr) %>% 
    return
}

# Calibrate the pca outputs to always have the exposed individuals on the left and the sheltered on the right
Calibrate_pca_outputs <- function(df, colName){
  for (col in colName){
    if ((df[[col]][which(df$Exposition == "SHELT")] %>% mean) < 0){
      df[[col]] <- - df[[col]] 
    }
  }
  df %>% 
    return
}

# Run the local pcas for all the candidate inversions
run_loca_pca_inversions <- function(df, population = NULL, genetic_data = data){
  # Initialise the loop by creating the necessary data frames
  pcas_candidate_inversions <- data.frame()
  clusters_candidate_inversions <- data.frame()
  pca_contribs_candidate_inversions <- data.frame()
  for (inversion in df$Inversion %>% unique){
    cat("\n", inversion)
    df_inv <- df %>% 
      filter(Inversion == inversion)
    # Get a sequence of all the possible positions in the candidate inversion
    if (dim(df_inv)[1] > 1){
      positions_in_inversion <- c()
      for (portion in 1:dim(df_inv[1])){
        positions_in_inversion <- c(
          positions_in_inversion,
          seq(
          from = df_inv[portion, ]$Start,
          to = df_inv[portion, ]$End,
          by = 1
          )
        )
      }
    }else{
      positions_in_inversion <- seq(
        from = (df %>% 
                  filter(Inversion == inversion))$Start,
        to = (df %>% 
                filter(Inversion == inversion))$End,
        by = 1
      )
    }
    
    # Isolate the chromosome name
    chromosome <- (df %>% 
                     filter(Inversion == inversion))$Chromosome %>% 
      unique
    # Isolate all the positions in the candidate inversion
    name_positions_in_inversion <- data@tab %>%
      colnames %>% 
      as_tibble %>% 
      rename(Position = value) %>% 
      filter(grepl(chromosome, Position)) %>% 
      transform_position_ade2tidy() %>% 
      filter(Position %in% positions_in_inversion) %>% 
      unite(Position, Chromosome, Position) %>% 
      as.vector %>% unname %>% unlist
    
    # Extract the positions and individuals from the genind object and run a pca on it
    if (!is_null(population)){
      pca_inv <- data[which(data@pop == population)][loc = name_positions_in_inversion] %>% 
        scaleGen(NA.method = "mean", scale = FALSE, center = TRUE) %>% 
        dudi.pca(scale = TRUE, nf = 5, scannf = FALSE)
    }else{
      pca_inv <- data[loc = name_positions_in_inversion] %>% 
        scaleGen(NA.method = "mean", scale = FALSE, center = TRUE) %>% 
        dudi.pca(scale = TRUE, nf = 5, scannf = FALSE)
    }
    
    # Changing the sign of the calculations depending on the population and the transect
    pca_inv_li <- pca_inv$li %>% 
      rownames_to_column("Sample_Name") %>% 
      mutate(indiv = Sample_Name) %>%
      separate(indiv, c("Species", "Cov", "Population", "Sex", "Exposition", "Size", "Sample_ID"), sep="_") %>% 
      mutate(Exposition = ifelse(Exposition == "TRANSI", "TRANS", Exposition)) %>% 
      mutate(Inversion = inversion,
             Chromosome = chromosome) %>% 
      # Calibrate the pca outputs to make it comparable between inversions
      Calibrate_pca_outputs(c("Axis1", "Axis2", "Axis3", "Axis4", "Axis5")) %>%  
      select(-c(Species, Cov, Population, Sex, Exposition, Size, Sample_ID))
    
    cat("\n    Ran PCA")
    
    # Find clusters
    if (!is_null(population)){
      clustering_pca <- data[which(data@pop == population)][loc = name_positions_in_inversion] %>%
        find.clusters(
          n.pca = 3,
          n.clust = 3
        )
    }else{
      clustering_pca <- data[loc = name_positions_in_inversion] %>%
        find.clusters(
          n.pca = 3,
          n.clust = 3
        )
    }
    
    # Reorder the clustering to be comparable between inversions
    pca_inv_li <- clustering_pca$grp %>% 
      as.data.frame %>% 
      dplyr::rename(Group = ".") %>% 
      rownames_to_column("Sample_Name") %>% 
      mutate(Inversion = inversion,
             Chromosome = chromosome,
             Group = Group %>% as.numeric) %>% 
      left_join(pca_inv_li, by = c("Sample_Name", "Chromosome", "Inversion")) %>% 
      reorder_clustering() %>%
      relocate(Group, .after = Axis5)
    
    
    cat("\n   Found clusters")
    
    # Add this to a data frame to save it
    pcas_candidate_inversions <- add_table_to_df_in_iteration(pcas_candidate_inversions,
                                                              pca_inv_li)
    pca_contribs_candidate_inversions <- add_table_to_df_in_iteration(pca_contribs_candidate_inversions,
                                                                      pca_inv$co %>% 
                                                                        as.data.frame %>% 
                                                                        rownames_to_column("Position") %>% 
                                                                        mutate(Inversion = inversion,
                                                                               Chromosome = chromosome))
  }
  list(
    "PC_scores" = pcas_candidate_inversions,
    "PC_contribs" = pca_contribs_candidate_inversions
  ) %>% 
    return
}

# Remove duplicates in the expand.grid output
remove_duplicates_expand.grid <- function(df){
  df %>% 
    rowwise() %>% 
    mutate(Inv_comb = paste0(sort(c(Inv_1, Inv_2)), collapse = "-")) %>%
    distinct() %>%
    separate(Inv_comb, c("Inv_1", "Inv_2"), "-") %>% 
    return
}

Compute_LD_inversions <- function(df, LD_version = "intra_chromosomic"){
  # Initialise the loop by defining a new data frame
  LD_between_inversions <- data.frame()
  if (LD_version == "intra_chromosomic"){
    for (chromosome in df$Chromosome %>% unique){
      combination_inversions <- (df %>% 
                                   filter(Chromosome == chromosome))$Inversion %>% unique %>% 
        expand.grid(., .) %>% 
        remove_duplicates_expand.grid() %>% 
        select(starts_with("Inv_"))
      
      
      if (nrow(combination_inversions) > 0){
        for (i in 1:nrow(combination_inversions)){
          combination_inversions_i <- c(
            combination_inversions$Inv_1[i],
            combination_inversions$Inv_2[i]
          )
          pca_out_i <- df %>% 
            filter(Inversion %in% combination_inversions_i)
          
          # Get the percentage of individuals that have the same attribution (same cluster)
          nb_indivs_same_cluster <- ((
            (pca_out_i %>% filter(Inversion == combination_inversions$Inv_1[i]))$Group == (pca_out_i %>% filter(Inversion == combination_inversions$Inv_2[i]))$Group
          ) %>% 
            sum) / (pca_out_i %>% select(Sample_Name) %>% unique %>% nrow) * 100
          
          LD_between_inversions <- add_table_to_df_in_iteration(LD_between_inversions,
                                                                combination_inversions[i, ] %>% 
                                                                  mutate(Chromosome = chromosome,
                                                                         Perc_indivs_same_cluster = nb_indivs_same_cluster))
          
        }
      }
    }
  }
  else{
    combination_inversions <- df$Inversion %>%
      unique %>% 
      expand.grid(., .) %>% 
      rename(Inv_1 = Var1,
             Inv_2 = Var2) %>% 
      remove_duplicates_expand.grid()
    
    
    if (nrow(combination_inversions) > 0){
      for (i in 1:nrow(combination_inversions)){
        combination_inversions_i <- c(
          combination_inversions$Inv_1[i],
          combination_inversions$Inv_2[i]
        )
        pca_out_i <- df %>% 
          filter(Inversion %in% combination_inversions_i)
        
        # Get the percentage of individuals that have the same attribution (same cluster)
        nb_indivs_same_cluster <- ((
          (pca_out_i %>% filter(Inversion == combination_inversions$Inv_1[i]))$Group == (pca_out_i %>% filter(Inversion == combination_inversions$Inv_2[i]))$Group
        ) %>% 
          sum) / (pca_out_i %>% select(Sample_Name) %>% unique %>% nrow) * 100
        
        LD_between_inversions <- add_table_to_df_in_iteration(LD_between_inversions,
                                                              combination_inversions[i, ] %>% 
                                                                mutate(Perc_indivs_same_cluster = nb_indivs_same_cluster))
        
      }
    }
  }
  LD_between_inversions %>% 
    return
}


################################################################################
######### 7. Subset genetic data in inversions and trace phylogeny   ###########
################################################################################


# Function to get the exposition of the individuals depending on the Group they are given in the local PCA
find_exposition_of_individuals <- function(df){
  exposed_indivs <- (df %>% 
                       filter(Group == 1))$Sample_Name %>% 
    unique
  trans_indivs <- (df %>% 
                     filter(Group == 2))$Sample_Name %>% 
    unique
  sheltered_indivs <- (df %>% 
                         filter(Group == 3))$Sample_Name %>% 
    unique
  
  return(list(
    "Exposed_indivs" = exposed_indivs,
    "Transition" = trans_indivs,
    "Sheltered" = sheltered_indivs
  ))
}

# Function to subset the vcf file to keep homokaryotype individuals and the positions located inside the inversions
Subset_genetic_data <- function(inversion, delim_invs, groups_pca, vcf_file, output_path, .is_inversion = is_inversion, .vcftools = "/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools", .bcftools = "/shared/software/miniconda/envs/bcftools-1.9/bin/bcftools",
                                force = FALSE){
  # First, create the directories to use in the function
  Create_dir_not_exist(c(paste0(output_path, "List_pos_per_inversion/"),
                         paste0(output_path, "List_indivs_inversion/"),
                         paste0(output_path, "VCF_inversion/")))
  
  
  if (force | !file.exists(paste0(output_path, "List_pos_per_inversion/", inversion, ".tsv"))){
    # Verify if we have a split inversion or not
    is_split <- delim_invs %>% 
      filter(Inversion == inversion) %>% 
      pull(Is_split) %>% 
      unique
    
    # Get the name of the chromosome we are working with
    chromosome <- (delim_invs %>% 
                     filter(Inversion == inversion))$Chromosome %>% 
      unique
    
    # Select the delimitations of the inversion to consider
    if (is_split){
      
      delims_inv_i <- delim_invs %>% 
        filter(Inversion == inversion) 
      
      Number_delims <- delims_inv_i %>% 
        group_by(Population) %>% 
        summarize(count = n()) %>% 
        pull(count) %>% 
        unique
      
      Starts_split_inversion <- c()
      Ends_split_inversion <- c()
      for (i in 1:Number_delims){
        Starts_split_inversion <- c(
          Starts_split_inversion,
          delims_inv_i %>% 
            arrange(Start) %>% 
            dplyr::slice(Number_delims * (i - 1) + 1) %>% 
            pull(Start)
        )
        Ends_split_inversion <- c(
          Ends_split_inversion,
          delims_inv_i %>% 
            arrange(End) %>% 
            dplyr::slice(Number_delims * i) %>% 
            pull(End)
        )
      }
      
      # Make a list of all the positions to keep
      pos2keep <- mapply(function(x, y) seq(x, y, 1), Starts_split_inversion, Ends_split_inversion) %>% 
        unlist
      
    }else{
      delims_inv_i <- delim_invs %>% 
        filter(Inversion == inversion) %>% 
        group_by(Inversion) %>% 
        summarize(Start = min(Start),
                  End = max(End)) %>% 
        ungroup
      
      # Make alist of all the positions to keep
      pos2keep <- seq(
        from = delims_inv_i$Start,
        to = delims_inv_i$End,
        by = 1
      )
    }
    
    pos2keep %>%
      as.data.frame %>% 
      rename(Pos = ".") %>% 
      mutate(Chr = chromosome) %>% 
      relocate(Pos, .after = Chr) %>% 
      write.table(paste0(output_path, "List_pos_per_inversion/", inversion, ".tsv"),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  # Make a list of the invdividuals to keep
  if (.is_inversion){
    list_expos_indivs <- groups_pca %>% 
      filter(Inversion == inversion) %>% 
      find_exposition_of_individuals()
  }else{
    list_expos_indivs <-  list(
      "Exposed_indivs" = (ref_individuals %>%
                            filter(Exposition == "Exposed"))$Sample_Name,
      "Sheltered" = (ref_individuals %>%
                       filter(Exposition == "Sheltered"))$Sample_Name
    )
  }
  
  if (force | !file.exists(paste0(output_path, "List_indivs_inversion/", inversion, ".txt"))){
    # Save the names of the individiuals that are homozygous for the two inversion genotypes
    list_expos_indivs$Exposed_indivs %>% 
      as.data.frame %>% 
      rbind(list_expos_indivs$Sheltered %>% 
              as.data.frame) %>% 
      write.table(paste0(output_path, "List_indivs_inversion/", inversion, ".txt"),
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  if (force | !file.exists(paste0(output_path, "VCF_inversion/", inversion, ".vcf.gz"))){
    system2(.vcftools,
            args = paste0(" --gzvcf ", vcf_file,
                          " --positions ", output_path, "List_pos_per_inversion/", inversion, ".tsv",
                          " --keep ", output_path, "List_indivs_inversion/", inversion, ".txt",
                          " --recode --stdout | ", .bcftools,
                          " view -c 2 -Oz",
                          " > ", output_path, "VCF_inversion/", inversion, ".vcf.gz"))
  }
  return(list_expos_indivs)
}

# Sample genotypes effectively
Sample_genotypes <- function(df){
  
  # First, we need all the unique combinations of locus and sample name and their count in the dataframe
  unique_comb_count <- df %>% 
    group_by(locus, Sample_Name) %>% 
    summarize(count = n()) %>% 
    ungroup
  
  # Here we initialise the output data frame and iterate over the unique counts that we have
  Sampled_genotypes <- data.frame()
  for (unique_count in unique_comb_count$count %>% unique){
    
    # Filter the unique combinations that contain the right count
    unique_comb_i <- unique_comb_count %>% 
      filter(count == unique_count) %>% 
      arrange(locus, Sample_Name)
    
    # Select these unique combinations in the raw dataframe
    df_i <- df %>% 
      filter(locus %in% unique_comb_i$locus,
             Sample_Name %in% unique_comb_i$Sample_Name) %>% 
      arrange(locus, Sample_Name)
    
    # Sample a number between 1 and the number of time the combination is used in the original dataframe
    Random_positions_to_sample <- sample(1:unique_count, size = dim(unique_comb_i)[1], replace = TRUE)
    # Here we transform the random number into the actual position of the row in the dataframe
    Random_rows_to_sample_in_df <- Random_positions_to_sample + c(0, rep(unique_count, length(Random_positions_to_sample) - 1)) %>% cumsum
    # And sample these rows in the dataframe
    Positions_sampled_count <- df_i[Random_rows_to_sample_in_df, ]
    
    # Then, add these rows to the data frame to return
    if (dim(Sampled_genotypes)[1] == 0){
      Sampled_genotypes <- Positions_sampled_count
    }else{
      Sampled_genotypes <- Sampled_genotypes %>% 
        rbind(Positions_sampled_count)
    }
  }
  
  return(Sampled_genotypes)
}

# Function to prepare a haploid genome to trace the trees
Prepare_haploid_genome <- function(inversion, output_path, genetic_data = data, .vcftools = "/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools", force = FALSE){
  # First create the directories we need to use
  Create_dir_not_exist(c(paste0(output_path, "Fastas_inversions/"),
                         paste0(output_path, "Ped_inversions/")))
  
  if (force | !file.exists(paste0(output_path, "Ped_inversions/", inversion, ".ped"))){
    # Run this command to create a ped file (genotype of the individuals for each position)
    system2(.vcftools,args =c(paste0("--gzvcf ", output_path, "VCF_inversion/", inversion, ".vcf.gz",
                                     " --plink --out ", output_path, "Ped_inversions/", inversion)))
  }
  
  if (force | !file.exists(paste0(output_path, "Fastas_inversions/", inversion, "sequences_random.txt"))){
    # Import the newly created ped file
    pedfile <- read.table(paste0(output_path, "Ped_inversions/", inversion, ".ped"),
                          header = FALSE)
    
    # Transform the table to get a table with the genotype of each individual at each position
    data_fasta <- genetic_data@loc.fac %>% 
      rbind(pedfile[, 7:dim(pedfile)[2]]) %>% 
      t %>% 
      as.data.frame
    Sample_Name <- genetic_data@tab %>% rownames
    colnames(data_fasta) <- c("locus", Sample_Name)
    
    # Make a table with the genotype of every individual at every locus
    reshape_fasta <- data_fasta %>% 
      pivot_longer(!locus,
                   names_to = "Sample_Name",
                   values_to = "Genotype") %>% 
      mutate(Genotype = str_trim(Genotype))
    
    
    # select one random SNP per individuals and marker
    data_fasta_sub <- reshape_fasta %>%
      Sample_genotypes()
    
    # Replace missing data or weird data in the subset fasta
    data_fasta_sub$Genotype[which(data_fasta_sub$Genotype %!in% c("A","T","C","G","N"))] <- "N"
    
    # Transform 3 colomn dataframe into a classic haploid genotype matrix
    data_fasta_random <- data_fasta_sub %>% 
      pivot_wider(names_from = locus, values_from = Genotype)
    
    data_fasta_random %>% 
      select(-Sample_Name) %>% 
      write.table(paste0(output_path, "Fastas_inversions/", inversion, "sequences_random.txt"),
                  sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  if (force | !file.exists(paste0(output_path, "Fastas_inversions/", inversion, "sequences_random.fasta"))){
    # Create fasta with the sampled random allele
    data_fasta_random %>% 
      select(Sample_Name) %>% 
      write.table(paste0(output_path, "Ped_inversions/", inversion, "Sample_Names.txt"), sep = "\t",
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    system(command = paste0("sed 's/$/_sequenceR/' ", output_path, "Ped_inversions/", inversion, "Sample_Names.txt",
                            " | sed 's/^/>/' > ", output_path, "Ped_inversions/", inversion, "Corrected.txt"))
    
    system(command = paste0("paste ", output_path, "Ped_inversions/", inversion, "Corrected.txt ", output_path, "Fastas_inversions/", inversion, "sequences_random.txt",
                            " | tr '\t' '\n' > ", output_path, "Fastas_inversions/", inversion, "sequences_random.fasta"))
    
    # Clean up the environment to only keep files that are useful for the next steps
    system(command = paste0("rm ",
                            output_path, "Ped_inversions/", inversion, "Sample_Names.txt ",
                            output_path, "Ped_inversions/", inversion, "Corrected.txt ",
                            output_path, "Fastas_inversions/", inversion, "sequences_random.txt ",
                            output_path, "VCF_inversion/", inversion, ".vcf.gz"))
  }
}

Get_Karyotype_name <- function(df, .list_expos_indivs = list_expos_indivs){
  df %>%
    mutate(Karyotype = ifelse((Sample_Name %in% .list_expos_indivs$Exposed) & (grepl("LOK", Sample_Name)), "SW_EE",
                              ifelse((Sample_Name %in% .list_expos_indivs$Sheltered) & (grepl("LOK", Sample_Name)), "SW_SS",
                                     ifelse((Sample_Name %in% .list_expos_indivs$Exposed) & (grepl("LAM", Sample_Name)), "FR_EE",
                                            ifelse((Sample_Name %in% .list_expos_indivs$Sheltered) & (grepl("LAM", Sample_Name)), "FR_SS", NA)))) %>%
             factor(levels = c("SW_EE", "SW_SS", "FR_EE", "FR_SS"))) %>%
    return()
}

# Function that loads a fasta file, computes the genetic distance between the sequences and draws a NJ tree
Run_and_trace_phylogeny <- function(inversion, .output_path = output_path, .list_expos_indivs = list_expos_indivs){
  
  # Load the fasta file
  data_phylo <- read.FASTA(paste0(.output_path, "Fastas_inversions/", inversion, "sequences_random.fasta"),
                           type = "DNA")
  
  # transform data
  phyDat <- phyDat(data_phylo, type = "DNA")
  fabalis <- as.phyDat(phyDat)
  
  # compute genetic distance based on polymorphic sites
  dna_dist <- dist.ml(fabalis, model = "F81")
  
  # draw NJ tree
  fabalis_NJ  <- NJ(dna_dist)
  fabalis_NJ$tip.label <- fabalis_NJ$tip.label %>% str_remove_all("_sequenceR")
  
  # Get the names of the individuals and classify them by karyotype
  Sample_Name <- fabalis_NJ$tip.label %>%
    as.data.frame %>% 
    rename(Sample_Name = ".") %>% 
    Get_Karyotype_name() %>%
    mutate(Population = ifelse(grepl("LOK", Sample_Name), "Sweden", "France"),
           Color_ID = case_when(
             Karyotype == "SW_EE" ~ met.brewer("Peru1", 7)[6],
             Karyotype == "SW_SS" ~ met.brewer("Peru1", 7)[5],
             Karyotype == "FR_EE" ~ met.brewer("Peru1", 7)[2],
             Karyotype == "FR_SS" ~ met.brewer("Peru1", 7)[3],
             TRUE ~ met.brewer("Peru1", 7)[4]
           ))
  
  # Regroup the individuals in the tree depending on their karyotype
  group_Info <- list(
    "SW_EE" = Sample_Name %>% filter(Karyotype == "SW_EE") %>% pull(Sample_Name),
    "SW_SS" = Sample_Name %>% filter(Karyotype == "SW_SS") %>% pull(Sample_Name),
    "FR_EE" = Sample_Name %>% filter(Karyotype == "FR_EE") %>% pull(Sample_Name),
    "FR_SS" = Sample_Name %>% filter(Karyotype == "FR_SS") %>% pull(Sample_Name)
  )
  
  fabalis_NJ_ <- groupOTU(fabalis_NJ, group_Info)
  
  # Plot the tree
  tree <- fabalis_NJ_ %>% 
    ggtree(aes(colour = group), layout = "unrooted", lwd = 1.05) +
    scale_color_manual(name = "Genotype",
                       values = c("FR_EE" = met.brewer("Peru1", 6)[1],
                                  "FR_SS" = met.brewer("Peru1", 6)[2],
                                  "SW_EE" = met.brewer("Peru1", 6)[6],
                                  "SW_SS" = met.brewer("Peru1", 6)[5])) +
    theme(text = element_text(size = 20))
  
  # Compute the distances between the ecotypes in the two countries
  #extract phylogenetic distance  
  phylo_dist <- cophenetic.phylo(fabalis_NJ)
  karyo_info_1 <- data.frame(cbind(fabalis_NJ$tip.label, Sample_Name$Karyotype %>% as.character))
  colnames(karyo_info_1) <- c("Var1", "Karyotype_ID1")
  karyo_info_2 <- data.frame(cbind(fabalis_NJ$tip.label, Sample_Name$Karyotype %>% as.character))
  colnames(karyo_info_2) <- c("Var2", "Karyotype_ID2")
  
  
  # #reshape phylogenetic distance matrix
  test <- as.matrix(phylo_dist)
  test[upper.tri(test, diag = TRUE)] <- NA
  test2 <- melt(test) %>% 
    filter(!is.na(value)) %>% 
    mutate(Population = ifelse((grepl("LOK", Var1) & grepl("LOK", Var2)), "Sweden", ifelse(
      (grepl("LAM", Var1) & grepl("LAM", Var2)), "France", NA)
    )) %>% 
    filter(!is.na(Population))
  test4 <- test2 %>% 
    left_join(karyo_info_2, by = "Var2") %>% 
    left_join(karyo_info_1, by = "Var1") %>% 
    rowwise() %>% 
    mutate(comp = paste0(sort(c(Karyotype_ID1, Karyotype_ID2)), collapse = "_"))
  
  #calculate mean net distance between group
  mean_divergence_karyotype <- aggregate(value ~ comp, data = test4, FUN = "mean") %>% 
    rename(Genotype = comp,
           Divergence = value) %>% 
    filter(Genotype %in% c("FR_EE_FR_SS", "SW_EE_SW_SS")) %>% 
    mutate(Divergence = Divergence %>% round(digits = 4))
  
  Divergence_karyotypes_countries <- data.frame(
    "Population" = c("Sweden", "France"),
    "Divergence" = c(
      (mean_divergence_karyotype %>% filter(grepl("SW", Genotype)))$Divergence,
      (mean_divergence_karyotype %>% filter(grepl("FR", Genotype)))$Divergence
    )
  )
  
  
  return(list(
    "Tree" = tree,
    "Divergence" = Divergence_karyotypes_countries
  ))
}
print("Finished function importation")



