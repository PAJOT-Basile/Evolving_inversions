# Import libraries
#install.packages("anyLib")
require("anyLib")
anyLib(c("tidyverse", "ggforce", "viridis", "ggnewscale", "LaplacesDemon", "tie", "bbmle", "zeallot"))

# Import functions from the `Cline_functions.R` script (the path may need to be modified)
source("./Cline_functions.R")

################################################################################
######################### 1. Useful functions  #################################
################################################################################

is_numeric_in_character <- function(x){
  #' This function checks if there is a numeric value contained in a character
  #' 
  #' Arguments:
  #'    x: string, factor, vector.
  #'      object to check
  #' Returns:
  #'    : bool
  #'      If there is a numeric value in the input, it returns TRUE, it returns
  #'      FALSE otherwise
  
  return(!is.na(x %>% as.numeric) %>% suppressWarnings)
}

is.continuous <- function(x){
  #'This function checks if a variable is discrete or continuous. 
  #' Warning: here a variable is considered discrete if it contains less than 10
  #' distinct levels. This approximation is done to simplify the distinction 
  #' between continuous or discrete variables, but it is not a real classification
  #' 
  #' Arguments:
  #'    x: numeric vector
  #'      A vector containing values to see if it is continuous or discrete
  #' Returns:
  #'    : bool.
  #'      This function returns TRUE if the variable is continuous (more than 10
  #'      levels) or FALSE if the variable is discrete (less than 10 levels)
  
  return(length(unique(x)) >= 10)
}

"%!in%" <- function(x, y){
  #' This function is the opposite of the "%in%" function. It checks if the 
  #' values contained in input x are NOT contained in y.
  #' 
  #' Arguments:
  #'  x: vector
  #'    This input contains any type of object, as long as they are the same type
  #'    as the ones in the y vector
  #'  y: vector
  #'    This input contains any type of object, as long as they are the same type
  #'    as the ones in the x vector
  #'    
  #' Returns:
  #'  : bool
  #'    This function returns TRUE for each value of x that is NOT in y and FALSE
  #'    for every value of x that is in y.
  
  return(!(x %in% y))
  }

################################################################################
############### 2. Functions for genomics with adegenet  #######################
################################################################################

transform_position_ade2tidy <- function(df){
  #' The position of each SNP in the adegenet format is as follows: 
  #' "CHROMOSOME_SNPposition.allele". 
  #' This function transforms the position column from the previous format into
  #' Two columns containing the chromosome in one column and the numeric 
  #' position in the other column.
  #' 
  #' Arguments:
  #'  df: data.frame
  #'    The dataframe containing the position column to transform
  #'  
  #' Returns:
  #'  df: data.frame
  #'    This function returns the input dataframe with two separate chromosome and
  #'    position columns

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
  #' This function selects the allele version of SNPs that have a given value 
  #' greater or equal to 0 and if this condition is not sufficient to select one
  #' allele, it picks the first one in the table.
  #' 
  #' Arguments:
  #'  df: data.frame
  #'    Dataframe containing a column called "Position" that can be duplicated 
  #'    and a column that is given to be filtered.
  #'  colName: symbol or string
  #'    It is the name of a column that is given to be filtered on. All the values
  #'    of this column that are greater or equal to 0 will be kept.
  #' 
  #' Results:
  #'  : data.frame
  #'    This dataframe contains only alleles of the first dataframe that have
  #'    values greater or equal to 0 in the selected column
  
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
  #' This function selects only clinal SNPs in at least one location. It requires
  #' that the input to have a "Position" and a "Model_to_select" column
  #' 
  #' Arguments:
  #'  df: data.frame
  #'    This dataframe must contain two columns at least with the "Position" and
  #'    the "Model_to_select" column. The "Position" column contains the position
  #'    of the SNP on the Genome. The "Model_to_select" column must contain the
  #'    name of the model to select: one of "Clinal", "Stable" or "Linear".
  #' 
  #' Results:
  #'  : data.frame
  #'    This dataframe contains all the SNPs that have at least one clinal model.
  #'    If the SNP is present several times in the dataframe, it will select all
  #'    the lines containing said SNP if at least one model is clinal.

  # We iterate over each unique SNP position in the dataframe
  for (position_i in df$Position %>% unique){
    # We count the number of rows for each SNP that use the "Clinal" model 
    is_clinal <- df %>% 
      filter(Position == position_i,
             Model_to_select == "Clinal") %>% 
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

################################################################################
########## 3. Functions for allelic frequencies from adegent objects  ##########
################################################################################

get_genotype_transect <- function(genetic_data, SNP_subset = NULL, meta_data = metadata){
  #' This function uses an adegenet object and returns the genotype of all or
  #' part of the SNPs.
  #' 
  #' Arguments:
  #'  genetic_data: Genind object 
  #'    It containing the SNPs from which to get the genotype
  #'  SNP_subset: data.frame (default = NULL)
  #'    This data frame contains columns "Position" and/or "SNP_name". If it is
  #'    given, it will only return the genotype from these SNPs. If this data
  #'    frame is used, it is required to have a column called "SNP_name".
  #'  meta_data: data.frame
  #'    This dataframe contains metadata on the individuals (for example the
  #'    position on the transect, the color, ...)
  #' 
  #' Returns:
  #'  : data.frame
  #'    This data frame contains the genotypes for both alleles of the selected
  #'    SNPs as well as some metadata on the individuals (position along the 
  #'    transect)
  
  # If the SNP subset is not selected, do this
  if (SNP_subset %>% is.null){
    # We select the table of genotypes for all SNPs along the whole genome
    genetic_data@tab %>%
      # We transform this into a data frame to keep the name of the columns
      as.data.frame %>% 
      # We add a column to the table to get the name of the individuals
      rownames_to_column("Sample_Name") %>% 
      # We merge this with the metadata using the name of the individuals
      inner_join(meta_data, by = "Sample_Name") %>% 
      # We select the columns we need to keep (position along the transect, 
      # name of the individuals, the population they are from and all the 
      # genotypes)
      select(Sample_Name, LCmeanDist, Population, starts_with("LG")) %>%
      return()
  
    # If the SNP subset is given, we use it 
  }else{
    # If there is no column called "Position" in the SNP_subset, we simply select
    # all the alleles of the selected SNPs
    if ("Position" %!in% names(SNP_subset)){
      # We select the SNPs we want in the genind object
      genetic_data[loc = SNP_subset$SNP_name]@tab %>% 
        # We transform this into a data frame to keep the row names 
        as.data.frame %>% 
        # We add a column containing the names of the individuals
        rownames_to_column("Sample_Name") %>% 
        # We merge this to the metadata using the name of the individuals
        inner_join(meta_data, by = "Sample_Name") %>% 
        # We select the required columns
        select(Sample_Name, LCmeanDist, Population, starts_with("LG")) %>% 
        return()
      
    # If there is a column called "Position" in the dataframe, we use it to 
    # select only the alleles of the SNPs we want to use
    }else{
      # The process is the same as in the previous cases
      genetic_data[loc = SNP_subset$SNP_name]@tab %>% 
        as.data.frame %>% 
        rownames_to_column("Sample_Name") %>% 
        inner_join(meta_data, by = "Sample_Name") %>% 
        # The only difference is in the following line, we select the columns
        # that contain the genotypes of only the alleles of the SNPs we want to
        # keep
        select(Sample_Name, LCmeanDist, Population, SNP_subset$Position %>% matches) %>% 
        return()
    }
  }
}

get_extreme_genotypes <- function(genetic_data, SNP_subset = NULL, Extreme_values = NULL, var = NULL, nb_extreme_indivs = 30, nb_indivs_to_keep = 20, meta_data = metadata){
  #' This function is used to get the genotypes of the individuals located on the
  #' extreme ends of the transect. If the selected individuals do not have the 
  #' right genotype, we can also use a supplementary dataset to re-filter the 
  #' individuals to select ones with more differences.
  #' 
  #' Arguments:
  #'  genetic_data: Genind object 
  #'    It containing the SNPs from which to get the genotype
  #'  SNP_subset: data.frame (default = NULL)
  #'    This data frame contains columns "Position" and/or "SNP_name". If it is
  #'    given, it will only return the genotype from these SNPs. If this data
  #'    frame is used, it is required to have a column called "SNP_name".
  #'  Extreme_values: data.frame (default = NULL)
  #'    This data frame contains the values to use to re-filter the selected
  #'    individuals. This variable works in combination with the next one.
  #'  var: string (default = NULL)
  #'    This string is the name of the column from Extreme_values to use to 
  #'    refilter the individuals.
  #'  nb_extreme_indivs: integer (default = 30)
  #'    This number is the number of individuals to sample in the first samping
  #'    process (position on the transect).
  #'  nb_indivs_to_keep: integer (default = 20)
  #'    This number is the number of individuals to re-select in the
  #'    nb_extreme_indivs using the Extreme_values dataset. It must therefore be
  #'    smaller than the number indicated by nb_extreme_indivs
  #'  meta_data: data.frame
  #'    This dataframe contains metadata on the individuals (for example the
  #'    position on the transect, the color, ...)
  #'    
  #' Results:
  #'  :list of 4 data.frames
  #'    This list contains 4 dataframes corresponding to the 4 groups of individuals
  #'    (the two extremes of each transect). Each data frame contains the 
  #'    genotypes of the selected individuals.

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
  #' This function is used to calculate the allelic frequencies of the obtained
  #' genotypes for each extreme of each transect
  #' 
  #' Arguments:
  #'  genetic_data: Genind object 
  #'    It containing the SNPs from which to get the genotype
  #'  SNP_subset: data.frame (default = NULL)
  #'    This data frame contains columns "Position" and/or "SNP_name". If it is
  #'    given, it will only return the genotype from these SNPs. If this data
  #'    frame is used, it is required to have a column called "SNP_name".
  #'  Extreme_values: data.frame (default = NULL)
  #'    This data frame contains the values to use to re-filter the selected
  #'    individuals. This variable works in combination with the next one.
  #'  var: string (default = NULL)
  #'    This string is the name of the column from Extreme_values to use to 
  #'    refilter the individuals.
  #'  nb_extreme_indivs: integer (default = 30)
  #'    This number is the number of individuals to sample in the first samping
  #'    process (position on the transect).
  #'  nb_indivs_to_keep: integer (default = 20)
  #'    This number is the number of individuals to re-select in the
  #'    nb_extreme_indivs using the Extreme_values dataset. It must therefore be
  #'    smaller than the number indicated by nb_extreme_indivs
  #'  meta_data: data.frame
  #'    This dataframe contains metadata on the individuals (for example the
  #'    position on the transect, the color, ...)
  #'    
  #' Results:
  #'  : data.frame
  #'      This data frame contains the allelic frequencies of the exposed and
  #'      sheltered parts of the transect for both populations. Each row is a SNP
  #'      and for each one, is given the population (here Sweden or France),
  #'      the allelic frequency in the sheltered and in the exposed part of the 
  #'      transect.
  
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
  #' This function calculates the allelic frequency of one allele in a population
  #' using the genotype of this allele in the population.
  #' 
  #' Arguments:
  #'  genotype: vector of numbers
  #'    This vector contains genotype values for all the individuals in the
  #'    population. It can contain missing values
  #'  
  #' Returns:
  #'  : numeric value
  #'    This value is the allelic frequency of the considered allele in the 
  #'    population.
  (sum(genotype, na.rm = TRUE) / ((length(genotype) - (is.na(genotype) %>% sum(na.rm = TRUE))) * 2)) %>% 
    return()
}

get_delta_freqs_and_F4 <- function(genetic_data, SNP_subset = NULL, Extreme_values = NULL, var = NULL, nb_extreme_indivs = 30, nb_indivs_to_keep = 20, meta_data = metadata){
  #' This function is used to calculate the differences in allelic frequencies
  #' between both ends of a transect. It also calculates the F4 statistic between
  #' two populations (the multiplication of differences in allelic frequencies).
  #' 
  #' Arguments:
  #'  genetic_data: Genind object 
  #'    It containing the SNPs from which to get the genotype
  #'  SNP_subset: data.frame (default = NULL)
  #'    This data frame contains columns "Position" and/or "SNP_name". If it is
  #'    given, it will only return the genotype from these SNPs. If this data
  #'    frame is used, it is required to have a column called "SNP_name".
  #'  Extreme_values: data.frame (default = NULL)
  #'    This data frame contains the values to use to re-filter the selected
  #'    individuals. This variable works in combination with the next one.
  #'  var: string (default = NULL)
  #'    This string is the name of the column from Extreme_values to use to 
  #'    refilter the individuals.
  #'  nb_extreme_indivs: integer (default = 30)
  #'    This number is the number of individuals to sample in the first samping
  #'    process (position on the transect).
  #'  nb_indivs_to_keep: integer (default = 20)
  #'    This number is the number of individuals to re-select in the
  #'    nb_extreme_indivs using the Extreme_values dataset. It must therefore be
  #'    smaller than the number indicated by nb_extreme_indivs
  #'  meta_data: data.frame
  #'    This dataframe contains metadata on the individuals (for example the
  #'    position on the transect, the color, ...)
  #'    
  #' Results:
  #'  : data.frame
  #'      This data frame contains the allelic frequencies of the exposed and
  #'      sheltered parts of the transect for both populations. Each row is a SNP
  #'      and for each one is given the difference in allelic frequency between 
  #'      the exposed and the sheltered part of the transect in both populations,
  #'      and the F4 statistic.
  
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
  #' This function is used to get the extreme values of a data frame with a
  #' specific condition and return a table containing the "Position" and "SNP_name"
  #' columns that can be used in the previous allelic frequency functions.
  #' 
  #' Arguments:
  #'  df: data.frame
  #'    This data frame is the data frame to filter, from which to extract
  #'    extreme values
  #'  condition: string
  #'    This string contains the condition to get the extreme values. For example
  #'    you could use '"Comp2 > 0"'. This would filter the values of column Comp2
  #'    to be greater than 0.
  #'  slice: boolean (default = FALSE)
  #'    This indicates if there is a need to get only the first n values of the
  #'    filtered data frame. If TRUE, the threshold argument is required and it
  #'    will keep only these values. If FALSE, it will simply filter.
  #'  threshold: integer (default = NULL)
  #'    This argument works in combination with the slice argument. It indicates
  #'    how many maximum values you need to keep in the final data frame.
  #' Returns:
  #'  : data.frame
  #'    This data frame contains the values that validate the given condition 
  #'    (and the most extreme values of the input data frame.)
  #'    
  #' Warning:
  #'  This function is not very replicable. It only works with some cases. It
  #'  could be perfected with more time and more rigor.

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
  #' This function creates a data frame that contains the polygon delimitations
  #' if it is required to separate chromosomes by using boxes in the visual
  #' representations
  #' 
  #' Arguments:
  #'  To_plot_manhat: data.frame
  #'    This data frame contains the data on the plotting of the values in a 
  #'    manhattan plot. 
  #'  colName: string or name
  #'    This argument is the name of the column in "To_plot_manhat" that we want
  #'    to plot.
  #'  chromosomes_centers: data.frame
  #'    This data frame contains information on the chromosomes: their name and
  #'    the position of their center is the minimum required for this function to
  #'    work.
  #'    
  #' Returns:
  #'  boxes: data.frame
  #'    This data frame contains the information on delimitations of the polygons
  #'    (position of the angles on the x and y axis). 
  
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

geom_manhattan <- function(df, mapping, thresholding = FALSE, absolute = TRUE, palette = c("grey71","orange2"), ...){
  #' This function traces a manhattan plot relying on the ggplot2 aesthetics.
  #' Some additional arguments have been added to be able to personalise the plots
  #' to the best.
  #' 
  #' Arguments:
  #'  df: data.frame
  #'    This data frame contains the information to plot. It requires at least
  #'    two columns: the position in the genome and the column to plot along the
  #'    genome.
  #'  mapping: mapping object
  #'    This is the mapping of the graph just like in the ggplot2 library. The
  #'    mapping in this function does not require an x column because the position
  #'    along the genome is taken by default.
  #'  thresholding: boolean (default = FALSE)
  #'    This is an argument that changes the output of the function. It outputs
  #'    the graph and a maximum cumulative position along the whole genome. It is
  #'    to be used with the "thresholds_manhattan" function.
  #'  absolute: boolean (default = TRUE)
  #'    This arguments simply says if you want to plot the absolute values in the
  #'    selected column or if you plot the real value
  #'  palette: vector (default = c("grey71", "orange2"))
  #'    This argument is the color palette to use to distinguish between the 
  #'    successive chromosomes. Warning: this function has been made to deal only 
  #'    with two-colored palettes.
  #'  ...
  #'    In these arguments, you can add any arguments that you would give a 
  #'    ggplot2 graph outside of the aesthetics.
  #'    
  #' Returns:
  #'  p: ggplot2 object
  #'    This returns the manhattan plot of the required column values along the
  #'    genome
  #'  : numeric value
  #'    If the "thresholding" argument is TRUE, then this function also returns
  #'    the maximum cumulative position in the whole genome.
  
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
                        arrange(as.numeric(gsub("\\D+", "", Chromosome))) %>% 
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
                        arrange(as.numeric(gsub("\\D+", "", Chromosome))) %>% 
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
    mutate(!!as.symbol(colName) := ifelse(absolute_list == TRUE, !!as.symbol(colName) %>% abs, !!as.symbol(colName))) %>% 
    # We also filter some values to lighten the plot a little bit
    filter((!!as.symbol(colName)) %>% abs > 0.05)
  
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

thresholds_manhattan <- function(df, mapping, percentages = NULL, values = NULL, ...){
  #' This function re-uses the geom_manhattan function to get a manhattan plot
  #' and it adds threshold values to it (usually represented by horizontal lines)
  #' 
  #' Arguments:
  #'  df: data.frame
  #'    This data frame contains the information to plot. It requires at least
  #'    two columns: the position in the genome and the column to plot along the
  #'    genome.
  #'  mapping: mapping object
  #'    This is the mapping of the graph just like in the ggplot2 library. The
  #'    mapping in this function does not require an x column because the position
  #'    along the genome is taken by default.
  #'  percentages: numeric value or vector of numeric values (default = NULL)
  #'    This contains values between 0 and 1. It will show the given percentages
  #'    of extreme values on the manhattan plot: if we want 20%, a line will 
  #'    appear on the plot showing where the cut of the 20% of highest values 
  #'    are. If this argument is used, do not use the values argument.
  #'  values: numeric value or vector of numeric values (default = NULL)
  #'    This contains numeric values in the range of values in the plot. It will
  #'    show the given values on the manhattan plot: if we want a cutoff value at
  #'    5, a vertical line will appear on the y axis at the value of 5.
  #'    If this argument is used, do not use the percentages argument.
  #'  ...:
  #'    supplementary arguments to pass to the geom_manhattan function
  #'    
  #' Returns:
  #'  : ggplot2 object
  #'    It will return the same plot as the geom_manhattan function, except some
  #'    vertical lines will have appeared at the given values or percentages.

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
      y_value_text <- threshold * 1.2
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
        select(Absolute) %>% min * 1.2
    }
    # For each value, we append the plot with the new line and text
    p <- p +
      # We add the horizontal line
      geom_hline(yintercept = y_value_line, color = colors_thresholds[i]) +
      # We add the text above
      annotate("text", x = max_value * 0.975, y = y_value_text,
               label = ifelse(is_null(percentages), value %>% as.character,
                              paste0((threshold * 100) %>% as.character, "%")))
    
  }
  # Once we iterated over all the possible cutoff values, we return the plot
  return(p)
}

################################################################################
#################### 5. Cline fitting and visualisation   ######################
################################################################################

optimise_clines <- function(Priors, logarithm = FALSE, ...){
  #' This function is made to optimise the best of three models (stable, linear
  #' and clinal) along a transect. It gives the best optimised parameters for
  #' all three models as well as the AIC of each model and which is better to use
  #' 
  #' Arguments:
  #'  Priors: data.frame
  #'    This data frame contains the prior values of the parameters to optimise.
  #'    This data frame has to contain at least a "Population" and "Position" 
  #'    column as well as some columns with the allelic frequencies on both ends
  #'    of the transect, the center of a cline and the width of the cline.
  #'  logarithm: boolean (default = FALSE)
  #'    This arguemnt is used to know if the optimisation process should be done
  #'    using logarithm transformations for the cline model fitting (TRUE) or 
  #'    not (FALSE)
  #'  ...:
  #'    Supplementary arguments to pass to the get_gentoype_transect function.
  #'    See this function to know which arguments to use.
  #' 
  #' Returns:
  #'  : data.frame
  #'    This data frame contains the optimised values of the three models for 
  #'    each position in the input data frame. It also tests the models on their
  #'    AIC and selects the best model to use to represent it.

  # First, we need to get the genotypes of the individuals along the transect to
  # optimise and test the models. It is recommended to use a subset of SNPs as 
  # this function runs for a long time
  genotype_transect <- get_genotype_transect(...)
  
  # This determines which function to use for the optimisation of the cline model
  # if we want to do a logarithm transformation, we use the clineflog function. 
  # Otherwise, we use the clinef function
  function_to_use <- ifelse(logarithm, clineflog, clinef)
  # Here, we make a vector of booleans to use in the tranformation of the priors
  # if we are using a logarithm transformation or not.
  logarithm_list <- rep(logarithm, (Priors %>% nrow)/2)
  
  # Then, we slighly modify the priors to be used in the optimisation process
  Priors_func <- Priors %>% 
    # We rename columns so we do not confuse them with other arguments
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
      p_left_shelt = ifelse(logarithm, p_left_shelt %>% logit, p_left_shelt),
      p_left_max = ifelse(logarithm, min(p_left_shelt + 0.1, 0.999999) %>% logit,
                          min(p_left_shelt + 0.1, 0.999999)),
      p_left_min = ifelse(logarithm, max(p_left_shelt - 0.1, 0.000001) %>% logit,
                          max(p_left_shelt - 0.1, 0.000001)),
      p_right_expos = ifelse(p_right_expos == 1, 0.999999,
                             ifelse(p_right_expos == 0, 0.000001, p_right_expos)),
      p_right_expos = ifelse(logarithm, p_right_expos %>% logit, p_right_expos),
      p_right_max = ifelse(logarithm, min(p_right_expos + 0.1, 0.999999) %>% logit,
                           min(p_right_expos + 0.1, 0.999999)),
      p_right_min = ifelse(logarithm, max(p_right_expos - 0.1, 0.000001) %>% logit,
                           max(p_right_expos - 0.1, 0.000001))
    )
  
  
  # Now that the priors are ready, we can start the optimisation process
  Clinal_model <- Priors_func %>% 
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
    left_join(Priors_func, by = c("Pop", "Pos")) %>% 
    # We select the columns we want to keep
    select(Pop, Pos, reversed, Centre, Width, Left, Right) %>% 
    # Backtransform the optimised parameters
    mutate(Width = ifelse(logarithm_list, Width %>% exp, Width),
           Left = ifelse(logarithm_list, Left %>% invlogit, Left),
           Right = ifelse(logarithm_list, Right %>% invlogit, Right),
           # Backtransform with the reversed parameter
           Left = ifelse(reversed, 1 - Left, Left),
           Right = ifelse(reversed, 1 - Right, Right)) %>% 
    rename(Position = Pos,
           Population = Pop)
  
  # Once we estimated the clinal parameters, we are now going to estimate the 
  # parameter for the stable model
  Stable_model <- Priors %>% 
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
                coef() %>% round(digits = 3)) %>% 
    rename(Position = Pos,
           Population = Pop)
  
  # Once we estimated the clinal parameters, we are now going to estimate the 
  # parameters for the linear model
  Linear_model <- Priors %>% 
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
           Position = Pos)
  
  # The last step is to merger all the estimated parameters together and calculate
  # the AICs of the models using the estimated parameters.
  Comp_table <- Clinal_model %>% 
    # First, we merge the table of the estimated parameters together
    left_join(Stable_model, by = c("Population", "Position")) %>% 
    left_join(Linear_model, by = c("Population", "Position")) %>% 
    left_join(Priors_func %>% 
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
      # Then, the logarithm transformation of the priors to do the log optimisation
      # requires that the allelic frequencies be transformed using the logit function.
      # So, in the same way as for the width, we take the logit of the allelic
      # frequencies if we want to do the logarithm transformation.
      # Wa also change fixed allele frequencies to be able to transform them with
      # the logit function
      Left = ifelse(Left == 1, 0.999999,
                            ifelse(Left == 0, 0.000001, Left)),
      Left = ifelse(logarithm, Left %>% logit, Left),
      p_left_max = ifelse(logarithm, min(Left + 0.1, 0.999999) %>% logit,
                          min(Left + 0.1, 0.999999)),
      p_left_min = ifelse(logarithm, max(Left - 0.1, 0.000001) %>% logit,
                          max(Left - 0.1, 0.000001)),
      Right = ifelse(Right == 1, 0.999999,
                             ifelse(Right == 0, 0.000001, Right)),
      Right = ifelse(logarithm, Right %>% logit, Right),
      p_right_max = ifelse(logarithm, min(Right + 0.1, 0.999999) %>% logit,
                           min(Right + 0.1, 0.999999)),
      p_right_min = ifelse(logarithm, max(Right - 0.1, 0.000001) %>% logit,
                           max(Right - 0.1, 0.000001))
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
    mutate(Width = ifelse(logarithm_list, Width %>% exp, Width),
           Left = ifelse(logarithm_list, Left %>% invlogit, Left),
           Right = ifelse(logarithm_list, Right %>% invlogit, Right),
           # Backtransform with the reversed parameter
           Left = ifelse(reversed, 1 - Left, Left),
           Right = ifelse(reversed, 1 - Right, Right)) %>% 
    # We remove columns we do not need in the output
    select(-c(p_left_shelt, p_right_expos, contains("prior"), matches("max|min"))) %>%
    # Then, we choose which model is the best between the three estimated ones 
    # using the AICs
    mutate(Delta_AIC_stable = AIC_clinal - AIC_stable, Delta_AIC_linear = AIC_clinal - AIC_linear,
           Model_to_select = ifelse((Delta_AIC_stable > Delta_AIC_linear) & (Delta_AIC_stable > 0), "Stable",
                                    ifelse((Delta_AIC_stable > Delta_AIC_linear) & (Delta_AIC_stable < 0), "Clinal",
                                           ifelse((Delta_AIC_stable < Delta_AIC_linear) & (Delta_AIC_linear > 0), "Linear",
                                                  ifelse((Delta_AIC_stable < Delta_AIC_linear) & (Delta_AIC_linear < 0), "Clinal"))))) %>% 
    # We remove the AICs
    select(-starts_with("AIC_")) %>% 
    # Rename the population and position columns
    rename(Position = Pos,
           Population = Pop) %>% 
    ungroup
  
  # And return the final table
  return(Comp_table)
  
}

plot_clines <- function(cline_params, real_distance = FALSE, goodness_of_fit = FALSE, ...){
  #' This function gives a data frame containing all the information to plot the
  #' best distribution of allelic frequencies along the transect using three 
  #' models: stable, linear and clinal. 
  #' 
  #' Arguments:
  #'  cline_params: data.frame
  #'    This data frame contains the parameters that came out of the "optimise_clines"
  #'    function containing the parameters of the optimised model parameters.
  #'  real_distance: boolean (default = FALSE)
  #'    This arguemnt is used to know if the graphic representation should use the
  #'    real position of the individuals along the transect rather than smoothing
  #'    the positions. This argument should be TRUE when the "goodness_of_fit"
  #'    argument is TRUE.
  #'  goodness_of_fit: boolean (default = FALSE)
  #'    This argument indicates if it should return the goodness of fit of the 
  #'    plotted model as well as the plotting data for the distributions of 
  #'    allelic frequencies along the transect. This argument should be used with
  #'    the "real_distance" argument set to TRUE.
  #'  ...:
  #'    Supplementary arguments to pass to the "get_gentoype_transect" function.
  #'    See this function to know which arguments to use.
  #' 
  #' Returns:
  #'  : data.frame
  #'    This data frame contains the parameters to plot the best distribution of 
  #'    allelic frequencies along the transect between three models. If the 
  #'    "goodness_of_fit" argument is also given, then it should return a list
  #'    of two data frames. The first one containing the goodness of fit of the 
  #'    plotted models and then the plotting information.
  
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
  for (position_i in cline_params$Position %>% unique){
    # Here, we initialise a temporary data frame to put the population information
    # inside and then bind it to the general data frame.
    temp_data <- tibble()
    for (population in cline_params$Population %>% unique){
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
      
      # For each SNP, we extract the model to use (clinal, linear or stable) from
      # the input data frame
      model_to_use <- cline_params %>%
        ungroup %>% 
        filter(Position == position_i,
               Population == population) %>% 
        select(Model_to_select)
      
      # Select which distance to use (real_distance or a smoothed distance along
      # the transect)
      if (real_distance){
        Distance <- geno_pop_pos_i %>% select(LCmeanDist) %>% as.vector %>% unlist %>% unname
      }else{
        Distance <- seq(from = geno_pop_pos_i %>% select(LCmeanDist) %>% min,
                        to = geno_pop_pos_i %>% select(LCmeanDist) %>% max,
                        by = 1)
      }
      # We isolate the parameters we want to use for the SNP in the selected 
      # population.
      Cline_param_pop_pos_i <- cline_params %>% 
        filter(Population == population,
               Position == position_i) %>% 
        # And if the SNP is reversed (p_right > p_left), then, we reverse it
        # (1 - p)
        mutate(Left = ifelse(is_reversed, 1 - Left, Left),
               Right = ifelse(is_reversed, 1 - Right, Right))
      
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
                                           optimisation = FALSE,
                                           plotting=TRUE) %>% 
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
    return(list("Goodness_of_fit" = Goodness_return, "Cline" = Cline_return))
  }else{
    return(Cline_return)
  }
}

print("finished importation")



