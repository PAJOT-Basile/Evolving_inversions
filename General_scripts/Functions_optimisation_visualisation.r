# Import libraries
library(tidyverse); library(ggforce); library(viridis); library(ggnewscale); library(LaplacesDemon, lib.loc = "/shared/home/bpajot/R/x86_64-conda-linux-gnu-library/4.2/"); library(tie, lib.loc = "/shared/home/bpajot/R/x86_64-conda-linux-gnu-library/4.2/"); library(bbmle); library(zeallot)

# Import functions from another script
source("/shared/projects/pacobar/finalresult/bpajot/scripts/Phenotypic_analysis/Cline_functions.R")

# Make a function that recognises if a string is a number. It returns TRUE if the character is a number and FALSE otherwise
is_numeric_in_character <- function(x){
  return(!is.na(x %>% as.numeric) %>% suppressWarnings)
}

# Function to know if a variable is discrete or continuous
is.continuous <- function(x){
  return(length(unique(x)) >= 10)
}

# Function to create the boxes and the color palette for the chromosomes if needed
geom_box_background <- function(To_plot_manhat, colName, chromosome_centers){

  colName <- as.character(colName)
  colName <- as.name(substitute(colName))
  # The height of the boxes has to be bigger than the maximal value and minimal value of the considered column
  max_col <- (To_plot_manhat %>% 
    select(as.symbol(colName)) %>% max(na.rm = TRUE)) + 0.2
  min_col <- To_plot_manhat %>% 
    select(as.symbol(colName)) %>% min(na.rm = TRUE) - 0.2
  # Create the boxes in the background
  boxes <- To_plot_manhat  %>% 
    right_join(chromosome_centers, by="Chromosome") %>% 
    filter((!!as.symbol(colName)) %>% abs > 0.1) %>% 
    group_by(Chromosome) %>% 
    summarise(Min = bp_cum %>% min,
              Max = bp_cum %>% max) %>% 
    pivot_longer(cols = c(Max, Min), values_to = "x_value", names_to = "operation") %>% 
    select(-operation) %>% 
    rbind(., .) %>% 
    arrange(Chromosome, x_value) %>% 
    cbind("y_value" = rep(c(min(c(0, min_col)), max(c(1.02, max_col)), max(c(1.02, max_col)), min(c(0, min_col))), nrow(.)/4))
  
  # Return the created variables
  return(boxes)
}

# Make a function to transform the position from a adagenet object to a usable position in tydiverse
transform_position_ade2tidy <- function(df){
  df %>%
    separate(Position, c("Position", "Allele"), sep="\\.") %>% 
    select(-Allele) %>% 
    separate(Position, c("LG", "SUPER", "SUPER_frag", "Position"), sep="_") %>% 
    unite("Chromosome", LG, SUPER, SUPER_frag) %>% 
    return()
}

#rm(df, mapping, thresholding, absolute, palette, mapping_ori, colName, to_select, i, To_plot, data_cum, To_plot_manhat, chromosome_centers, absolute_list, color_chromosome, p, function_to_call, list_default_parameters, default_param, param, color, boxes, color_polygon_chrom, list_defaults)
# Make a function that draw a manhattan plot for a specified column name in a df
geom_manhattan <- function(df, mapping, thresholding = FALSE, absolute = TRUE, palette = c("grey71","orange2"), ...){
  
  # Keeps the function from running for nothing
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
  mapping_ori <- mapping
  
  colName <- mapping$y[2] %>% as.character
  colName <- as.name(substitute(colName))
  
  # Now that all the argument importation has been done, we are going to select the columns in the input dataframe that
  # we want to keep in the analysis. To do this, we make a list of all the parameters that were called to use in the 
  # "matches" function is select
  # We initialise the string "to_select" to an empty character that will be complementeds
  to_select <- ""
  # We iterate over the names of the arguments passed to the function (except for the first one which is the name of the datafame)
  for (i in 1:length(mapping)){
    if (to_select == ""){
      to_select <- quo_name(mapping[[i]])
    }else{
      # Add the name of the parameters from the mapping to keep
      to_select <- paste0(to_select, "|", quo_name(mapping[[i]]))
    }
  }
  
  # Then, we separate the position argument in the dataframe to use it to place each SNP in the genome
  # If the positions are already numeric, we simply continue
  if (is.numeric(df$Position[1])){
    # Now, we keep only the columns of interest in the dataframe
    To_plot <- df %>% 
      select(Position, Chromosome, matches(to_select)) %>%
      mutate(Chromosome = Chromosome %>% 
               factor(levels = df %>% select(Chromosome) %>% 
                        unique %>% 
                        arrange(as.numeric(gsub("\\D+", "", Chromosome))) %>% 
                        as.vector %>% unname %>% unlist)) %>% 
      drop_na()
  }else{
    To_plot <- df %>%
      transform_position_ade2tidy() %>% 
      mutate(Chromosome = Chromosome %>% 
               factor(levels = df %>%
                        transform_position_ade2tidy() %>%
                        select(Chromosome) %>% 
                        unique %>% 
                        arrange(as.numeric(gsub("\\D+", "", Chromosome))) %>% 
                        as.vector %>% unname %>% unlist),
             Position = Position %>% as.numeric) %>% 
      select(Position, Chromosome, matches(to_select)) %>% 
      # Finally, we drop the missing values
      drop_na()
  }
  
  # Then, we make a cumulative dataframe with the positions of each end and beginning of chromosome
  data_cum <- To_plot %>% 
    group_by(Chromosome) %>% 
    summarise(max_bp = Position %>% max) %>% 
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(Chromosome, bp_add)
  
  # We join it to the filtered dataframe
  To_plot_manhat <- To_plot %>% 
    inner_join(data_cum, by="Chromosome") %>%
    mutate(Position = Position %>% as.numeric,
           bp_cum = Position + bp_add)
  
  # Then, we find the center of the chromosomes
  chromosome_centers <- To_plot_manhat %>% 
    group_by(Chromosome) %>% 
    summarise(center = bp_cum %>% mean)
  
  # We select two colors that will be used to distinguish chromosomes
  color_chromosome <- rep(palette, To_plot$Chromosome %>% unique %>% length)
  
  # Once all this is done, we create a vector of boolean to use as indicator to take or not the absolute value of the 
  # column to plot
  absolute_list <- rep(absolute, nrow(To_plot_manhat))
  
  # Finally, we mutate, if needed, the column to plot with the absolute value
  To_plot_manhat <- To_plot_manhat  %>% 
    right_join(chromosome_centers, by="Chromosome") %>% 
    mutate(!!as.symbol(colName) := ifelse(absolute_list == TRUE, !!as.symbol(colName) %>% abs, !!as.symbol(colName))) %>% 
    # We also filter some values to lighten the plot a little bit
    filter((!!as.symbol(colName)) %>% abs > 0.05)
  
  # Once this is done, we create the architecture of the plot (x axis, name of the axis and the theme to use)
  p <- ggplot(data = To_plot_manhat, aes(x = bp_cum)) +
    scale_x_continuous(labels = ((chromosome_centers$Chromosome %>% str_split_fixed(., "_", 2))[, 1] %>% 
                                   str_split_fixed(., "G", 2))[, 2],
                       breaks = chromosome_centers$center) +
    labs(x = "Linkage groups",
         y = colName) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          text = element_text(size = 20))
  
  # And we add the mapping to the plot
  list_default_parameters <- list("y" = ".", "colour" = "Chromosome")
  function_to_call <- "geom_point(aes(mapping)"
  for (param in names(list_default_parameters)){
    default_param <- list_default_parameters[[param]]
    if (param %!in% names(mapping)){
      if (param != "colour"){
        function_to_call <- paste0(function_to_call, ", ", param, " = ", default_param)
      }else{
        mapping <- c(mapping, aes(color = Chromosome %>% as.factor))
        class(mapping) <- "uneval"
      }
    }else{
      p$labels[[param]] <- mapping[[param]] %>% quo_name
    }
  }
  
  if (length(list(...))){
    function_to_call <- paste0(function_to_call, ", ...")
  }
  function_to_call <- paste0(function_to_call, ", inherit.aes = TRUE)")
  if ("colour" %in% names(mapping_ori)){
    color <- mapping$colour[2] %>% as.character
    color <- as.name(substitute(color))
    # If there is a color aesthetic to take into account, we have to find another way to show the chromosomes
    # So, we create some polygons in the background and the color palette
    boxes <- geom_box_background(To_plot_manhat, colName, chromosome_centers)
    
    
    function_to_call <- function_to_call %>% str_replace("geom_point\\(", "geom_point\\(data = To_plot_manhat, ")
    # We add this to the graph
    p <- p +
      # First, the polygons so they are in the background
      geom_polygon(data = boxes, aes(x = x_value, y = y_value, fill = Chromosome), color = NA) +
      # Then, we chose how to color the chromosomes
      scale_fill_manual(values = color_chromosome, guide="none")
    
    if (is.continuous(To_plot_manhat[[as.character(color)]])){
      p <- p +
        eval(parse(text = function_to_call)) +
        scale_color_gradientn(name = as.character(color),
                              colors = c("darkorchid4", "darkorchid", "mediumorchid1", "magenta"))
    }else{
      if (color == colName){
        To_plot_manhat <- To_plot_manhat %>% 
          mutate(color_column = !!as.symbol(color) %>% factor(levels = df %>% 
                                                                select(as.symbol(color)) %>% 
                                                                unique() %>% arrange(!!as.symbol(color)) %>% 
                                                                as.vector %>% unname %>% unlist))
        mapping$colour <- quo(color_column)
        
        color_polygon_chrom <- To_plot_manhat %>% 
          select(color_column) %>%
          n_distinct() %>% 
          plasma()
      }else{
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

      p <- p +
        eval(parse(text = function_to_call)) +
        scale_color_manual(name = as.character(color),
                           values = color_polygon_chrom,
                           drop = FALSE)
    }
    p$layers[[2]]$computed_mapping <- NULL
    p$layers[[2]]$mapping <- mapping
    
      
  }else{
    p <- p +
      eval(parse(text = function_to_call)) +
      scale_color_manual(values = color_chromosome, guide="none")
    p$layers[[1]]$computed_mapping <- NULL
    p$layers[[1]]$mapping <- mapping
    
  }
  # Finally, if there is a thresholding (i.e. if the function is embedded in the thresholding function), then, we simply return a list of parameters
  if (thresholding) {
    return(list("plot" = p, "max_value" = To_plot_manhat %>% select(bp_cum) %>% max))
  }else{
    return(p)
  }
}

# Function that plots the contribution of said variable (PC2/PC3) along the genome as well as 
# Lines allowing to see how many percent of outliers can be sampled
thresholds_manhattan <- function(df, mapping, percentages = NULL, values = NULL, ...){
  if (is_null(percentages) & is_null(values)){
    stop("At least one percentage or value threshold should be given")
  }
  
  if (is_null(percentages)){
    thresholds <- values
  }else{
    thresholds <- percentages
  }
  # First, we extract the contributions of the considered dataframe
  manhattan_plotting <- geom_manhattan(df, mapping, thresholding = TRUE, ...)
  colName <- mapping$y[2] %>% as.character
  colName <- as.name(substitute(colName))
  p <- manhattan_plotting$p
  max_value <- manhattan_plotting$max_value
  
  colors_thresholds <- plasma(n = length(thresholds))
  
  # Once this is done, we plot the contribution of the SNPs along the genome
  for (threshold in thresholds){
    i <- which(thresholds == threshold)
    # Select the y position of the line and the text on the graph
    if (is_null(percentages)){
      y_value_line <- threshold
      y_value_text <- threshold * 1.2
    }else{
      y_value_line <- df %>% 
        select(colName) %>% 
        rownames_to_column("Position") %>% 
        mutate(Absolute = (!!as.symbol(colName)) %>% abs) %>% 
        arrange(desc(Absolute)) %>% 
        dplyr::slice(1:(threshold * nrow(.))) %>% 
        select(Absolute) %>% min
      y_value_text <- df %>% 
        select(colName) %>% 
        rownames_to_column("Position") %>% 
        mutate(Absolute = (!!as.symbol(colName)) %>% abs) %>% 
        arrange(desc(Absolute)) %>% 
        dplyr::slice(1:(threshold * nrow(.))) %>% 
        select(Absolute) %>% min * 1.2
    }
    p <- p +
      geom_hline(yintercept = y_value_line, color = colors_thresholds[i]) +
      annotate("text", x = max_value * 0.975, y = y_value_text,
               label = ifelse(is_null(percentages), value %>% as.character, paste0((threshold * 100) %>% as.character, "%")))
    
  }
  return(p)
}

# Function to get the genotype along the transect of the selected SNPs if needed
#rm(genetic_data, SNP_subset)
get_genotype_transect <- function(genetic_data, SNP_subset = NULL, metadata = metadata){
  if (SNP_subset %>% is.null){
    genetic_data@tab %>%
      as.data.frame %>% 
      rownames_to_column("Sample_Name") %>% 
      inner_join(metadata, by = "Sample_Name") %>% 
      select(Sample_Name, LCmeanDist, Population, starts_with("LG")) %>%
      return()
    
  }else{
    if ("Position" %!in% names(SNP_subset)){
      genetic_data[loc = SNP_subset$SNP_name]@tab %>% 
        as.data.frame %>% 
        rownames_to_column("Sample_Name") %>% 
        inner_join(metadata, by = "Sample_Name") %>% 
        select(Sample_Name, LCmeanDist, Population, starts_with("LG")) %>% 
        return()
    }else{
      genetic_data[loc = SNP_subset$SNP_name]@tab %>% 
        as.data.frame %>% 
        rownames_to_column("Sample_Name") %>% 
        inner_join(metadata, by = "Sample_Name") %>% 
        select(Sample_Name, LCmeanDist, Population, SNP_subset$Position %>% matches) %>% 
        return()
    }
  }
}

# Make a function to calculate the allelic frequency depending on the genotype 
Calculate_frequency_transect <- function(genotype){
  (sum(genotype, na.rm = TRUE) / ((length(genotype) - (is.na(genotype) %>% sum(na.rm = TRUE))) * 2)) %>% 
    return()
}

# Make a function to get the genotypes of the individuals in the extreme locations of the transect
get_extreme_genotypes <- function(genetic_data, SNP_subset = NULL, Extreme_values = NULL, var = NULL, nb_extreme_indivs = 30, nb_indivs_to_keep = 20, metadata = metadata){
  # First, we have to get the genotype of the considered SNPs on the whole transect
  genotype_transect <- get_genotype_transect(genetic_data = genetic_data,
                                             SNP_subset = SNP_subset,
                                             metadata = metadata)
  # Then, we get the genotype of the individuals we are interested in (the
  # 30 most extreme individuals in the two environments and then we remove
  # those that have the less extreme PC scores to look for parallelism (PC2)
  #or reverse association (PC3))
  # Select the most extreme individuals and their genotype
  if (!is.null(Extreme_values) & !is.null(var)){
    LAM_EXPOS <- genotype_transect %>% 
      filter(Population == "France") %>% 
      arrange(LCmeanDist %>% desc) %>% 
      dplyr::slice(1:nb_extreme_indivs) %>% 
      left_join(Extreme_values %>% 
                  select(var %>% matches) %>% 
                  rownames_to_column("Sample_Name"),
                by = "Sample_Name") %>% 
      arrange(!!as.symbol(var) %>% desc) %>% 
      dplyr::slice(1:nb_indivs_to_keep)
    
    
    LAM_SHELT <- genotype_transect %>% 
      filter(Population == "France") %>% 
      arrange(LCmeanDist) %>% 
      dplyr::slice(1:nb_extreme_indivs) %>% 
      left_join(Extreme_values %>% 
                  select(var %>% matches) %>% 
                  rownames_to_column("Sample_Name"),
                by = "Sample_Name") %>% 
      arrange(!!as.symbol(var)) %>% 
      dplyr::slice(1:nb_indivs_to_keep)
    
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
    nb_extreme_indivs <- 20
    LAM_EXPOS <- genotype_transect %>% 
      filter(Population == "France") %>% 
      arrange(LCmeanDist %>% desc) %>% 
      dplyr::slice(1:nb_extreme_indivs)
    
    
    LAM_SHELT <- genotype_transect %>% 
      filter(Population == "France") %>% 
      arrange(LCmeanDist) %>% 
      dplyr::slice(1:nb_extreme_indivs)
    
    LOK_EXPOS <- genotype_transect %>% 
      filter(Population == "Sweden") %>% 
      arrange(LCmeanDist %>% desc) %>% 
      dplyr::slice(1:nb_extreme_indivs)
    
    LOK_SHELT <- genotype_transect %>% 
      filter(Population == "Sweden") %>% 
      arrange(LCmeanDist) %>% 
      dplyr::slice(1:nb_extreme_indivs)
  }
  
  return(list("LAM_EXPOS" = LAM_EXPOS, "LAM_SHELT" = LAM_SHELT, "LOK_EXPOS" = LOK_EXPOS, "LOK_SHELT" = LOK_SHELT))
}
# Make a function that calculates the allelic frequencies of selected SNPs between extreme individuals
#rm(list = c(genetic_data, SNP_subset, var, nb_extreme_indivs, nb_indivs_to_keep, genotypes, LAM_EXPOS, LAM_SHELT, LOK_EXPOS, LOK_SHELT, apropos("p_LAM_"), apropos("p_LOK_"), Allele_freqs))
get_allelic_frequencies <- function(genetic_data, SNP_subset = NULL, Extreme_values = NULL, var = NULL, nb_extreme_indivs = 30, nb_indivs_to_keep = 20, metadata = metadata){
  # First, we get the genotypes of extreme individuals in the transect
  genotypes <- get_extreme_genotypes(genetic_data = genetic_data, 
                                     SNP_subset = SNP_subset, 
                                     Extreme_values = Extreme_values, 
                                     var = var, 
                                     nb_extreme_indivs = nb_extreme_indivs, 
                                     nb_indivs_to_keep = nb_indivs_to_keep,
                                     metadata = metadata)
  LAM_EXPOS <- genotypes$LAM_EXPOS
  LAM_SHELT <- genotypes$LAM_SHELT
  LOK_EXPOS <- genotypes$LOK_EXPOS
  LOK_SHELT <- genotypes$LOK_SHELT
  # Calculate the allelic frequencies on the selected individuals
  p_LAM_EXPOS <- (LAM_EXPOS %>% 
                    select(starts_with("LG")) %>% 
                    colSums(na.rm = TRUE)) / ((nrow(LAM_EXPOS) - (LAM_EXPOS %>% 
                                                                    select(starts_with("LG")) %>% 
                                                                    is.na %>% 
                                                                    colSums(na.rm = TRUE))) * 2)
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
  Allele_freqs <- rbind(p_LAM_SHELT, p_LAM_EXPOS) %>% 
    matrix(nrow = p_LAM_SHELT %>% length, byrow=TRUE) %>% 
    as_tibble %>% 
    rename(p_left_shelt = V1,
           p_right_expos = V2) %>% 
    mutate(Population = "France",
           Position = p_LAM_SHELT %>% names) %>% 
    rbind(rbind(p_LOK_SHELT, p_LOK_EXPOS) %>% 
            matrix(nrow = p_LOK_SHELT %>% length, byrow=TRUE) %>% 
            as_tibble %>% 
            rename(p_left_shelt = V1,
                   p_right_expos = V2) %>% 
            mutate(Population = "Sweden",
                   Position = p_LOK_SHELT %>% names))
  
  return(Allele_freqs)
}

get_delta_freqs_and_F4 <- function(genetic_data, SNP_subset = NULL, Extreme_values = NULL, var = NULL, nb_extreme_indivs = 30, nb_indivs_to_keep = 20, metadata = metadata){
  get_allelic_frequencies(genetic_data = genetic_data,
                          SNP_subset = SNP_subset,
                          Extreme_values = Extreme_values,
                          var = var,
                          nb_extreme_indivs = nb_extreme_indivs,
                          nb_indivs_to_keep = nb_indivs_to_keep,
                          metadata = metadata) %>% 
    mutate(Delta_freq = p_right_expos - p_left_shelt)%>% 
    pivot_wider(names_from = Population, values_from = c(p_right_expos, p_left_shelt, Delta_freq)) %>% 
    mutate(F4_stat = Delta_freq_France * Delta_freq_Sweden) %>% 
    return()
}


# Isolate the extreme values of pc2 contributions for each SNP
get_extreme_values <- function(df, condition, slice = FALSE, threshold = NULL){
  var <- (condition %>% str_split_fixed(., " ", 3))[, 1]
  order_arranging <- grepl("> 0", condition)
  x <- df %>%  
    filter(!! rlang::parse_expr(condition))
  if ("Position" %!in% (x %>% names)){
    x <- x %>% 
      rownames_to_column("Position")
  }
  x <- x %>% 
    mutate(SNP_name = str_split_fixed(Position, "\\.", 2)[, 1])
  
  if (!(threshold %>% is.null) & (slice)){
    if (grepl(" ", (condition %>% str_split_fixed(., " ", 3))[, 3])){
      warning("The condition may be too complex for the slice to be effective. Please use a simple argument like 'Comp2 > 0'.")
    }
    x %>%
      select(var %>% matches) %>%
      arrange(ifelse(rep(order_arranging, nrow(.)), desc(!!as.symbol(var)), !!as.symbol(var))) %>%
      dplyr::slice(1:(threshold * nrow(.))) %>% 
      return()
  }else{
    x %>% return()
  }
}
# Make a function to optimise the clines taking into account the priors for said clines
#rm(Priors, logarithm, genetic_data, SNP_subset, genotype_transect, function_to_use, logarithm_list, Priors_func, Clinal_model, Stable_model, Linear_model, Comp_table)
optimise_clines <- function(Priors, logarithm = FALSE, ...){
  # Get the genotype along the transect
  genotype_transect <- get_genotype_transect(...)
  
  # Which function to use in the optimisation
  function_to_use <- ifelse(logarithm, clineflog, clinef)
  logarithm_list <- rep(logarithm, (Priors %>% nrow)/2)
  # Do the optimisation
  Priors_func <- Priors %>% 
    rename(Pop = Population,
           Pos = Position) %>% 
    group_by(Pop, Pos) %>% 
    mutate(# First, we prepare the width
      Width_prior = ifelse(logarithm, Width_prior %>% log, Width_prior),
      Width_max = ifelse(logarithm, Width_max %>% log, Width_max),
      Width_min = ifelse(logarithm, Width_min %>% log, Width_min),
      # Then, we do check if the association with the transect is reversed or not
      reversed = p_left_shelt > p_right_expos,
      # If the frequencies are reversed, we use 1 - freq as a prior
      p_left_shelt = ifelse(reversed, 1 - p_left_shelt, p_left_shelt),
      p_right_expos = ifelse(reversed, 1 - p_right_expos, p_right_expos),
      # Then, we look at the maximum and minimum values of the left frequency 
      p_left_shelt = ifelse(p_left_shelt == 1, 0.999999, ifelse(p_left_shelt == 0, 0.000001, p_left_shelt)),
      p_left_max = ifelse(logarithm, min(p_left_shelt + 0.1, 0.999999) %>% logit, min(p_left_shelt + 0.1, 0.999999)),
      p_left_min = ifelse(logarithm, max(p_left_shelt - 0.1, 0.000001) %>% logit, max(p_left_shelt - 0.1, 0.000001)),
      # Then we look at the maximum and minimum values of the right frequency
      p_right_expos = ifelse(p_right_expos == 1, 0.999999, ifelse(p_right_expos == 0, 0.000001, p_right_expos)),
      p_right_max = ifelse(logarithm, min(p_right_expos + 0.1, 0.999999) %>% logit, min(p_right_expos + 0.1, 0.999999)),
      p_right_min = ifelse(logarithm, max(p_right_expos - 0.1, 0.000001) %>% logit, max(p_right_expos - 0.1, 0.000001)),
      # We convert the allelic frequencies using the logit function if needed
      p_left_shelt = ifelse(logarithm, p_left_shelt %>% logit, p_left_shelt),
      p_right_expos = ifelse(logarithm, p_right_expos %>% logit, p_right_expos),
    )
  
  # Estimate the clinal model
  Clinal_model <- Priors_func %>% 
    group_by(Pop, Pos) %>% 
    bow(tie(Centre, Width, Left, Right) := mle2(function_to_use,
                                                # Priors
                                                list(centre = Centre_prior,
                                                     width = Width_prior ,
                                                     left = p_left_shelt,
                                                     right = p_right_expos),
                                                # Length along the transect and genotype
                                                data = list(x = genotype_transect %>% 
                                                              filter(Population == Pop) %>% 
                                                              select(LCmeanDist) %>% 
                                                              drop_na %>% 
                                                              as.vector %>% unlist %>% unname,
                                                            g = genotype_transect %>% 
                                                              filter(Population == Pop) %>%
                                                              select(Pos %>% matches) %>% 
                                                              rename(Geno = starts_with("LG")) %>% 
                                                              drop_na %>% 
                                                              mutate(Geno = ifelse(rep(reversed, nrow(.)), 2 - Geno, Geno)) %>% 
                                                              select(Geno) %>% 
                                                              as.vector %>% unlist %>% unname,
                                                            n = 2),
                                                # Method to use
                                                "L-BFGS-B",
                                                # Upper limits of the optimisation
                                                upper = list(centre = Centre_max,
                                                             width = Width_max,
                                                             left = p_left_max,
                                                             right = p_right_max),
                                                # Lower limits of the optimisation
                                                lower = list(centre = Centre_min,
                                                             width = Width_min,
                                                             left = p_left_min,
                                                             right = p_right_min)) %>% 
          coef() %>% 
          round(digits = 3)) %>% 
    left_join(Priors_func, by = c("Pop", "Pos")) %>% 
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
  
  # Estimate the stable model
  Stable_model <- Priors %>% 
    mutate(Mean_freq = (p_left_shelt + p_right_expos)/2) %>% 
    rename(Pop = Population,
           Pos = Position) %>% 
    group_by(Pop, Pos) %>% 
    summarise(Stable_model_fit = mle2(stable,
                                      list(p_all = Mean_freq),
                                      data = list(x = genotype_transect %>% 
                                                    filter(Population == Pop) %>% 
                                                    select(LCmeanDist) %>% 
                                                    drop_na %>% 
                                                    as.vector %>% unlist %>% unname,
                                                  g = genotype_transect %>% 
                                                    filter(Population == Pop) %>%
                                                    select(Pos %>% matches) %>% 
                                                    drop_na %>% 
                                                    as.vector %>% unlist %>% unname,
                                                  n = 2),
                                      method = "L-BFGS-B",
                                      upper = list(p_all = 0.999999),
                                      lower = list(p_all = 0.000001)) %>% 
                coef() %>% round(digits = 3)) %>% 
    rename(Position = Pos,
           Population = Pop)
  
  # Estimate the linear model
  Linear_model <- Priors %>% 
    rename(Pos = Position,
           Pop = Population) %>% 
    group_by(Pop, Pos) %>% 
    bow(tie(p_left, p_right) := mle2(linear,
                                     list(p_left = p_left_shelt,
                                          p_right = p_right_expos),
                                     data = list(x = genotype_transect %>% 
                                                   filter(Population == Pop) %>% 
                                                   select(LCmeanDist) %>% 
                                                   drop_na %>% 
                                                   as.vector %>% unlist %>% unname,
                                                 g = genotype_transect %>% 
                                                   filter(Population == Pop) %>%
                                                   select(Pos %>% matches) %>% 
                                                   drop_na %>% 
                                                   as.vector %>% unlist %>% unname,
                                                 n = 2),
                                     upper = list(p_left = 0.999999,
                                                  p_right = 0.999999),
                                     lower = list(p_left = 0.000001,
                                                  p_right = 0.000001)) %>% 
          coef() %>% round(digits = 3)) %>% 
    rename(Population = Pop,
           Position = Pos)
  
  # Merge all the models together and calculate the AIC of the models
  Comp_table <- Clinal_model %>% 
    left_join(Stable_model, by = c("Population", "Position")) %>% 
    left_join(Linear_model, by = c("Population", "Position")) %>% 
    left_join(Priors_func %>% 
                rename(Population = Pop,
                       Position = Pos), by = c("Population", "Position", "reversed")) %>% 
    rename(Pop = Population,
           Pos = Position) %>% 
    group_by(Pop, Pos) %>% 
    mutate(AIC_stable = mle2(stable, list(p_all = Stable_model_fit),
                             data = list(x = genotype_transect %>% 
                                           filter(Population == Pop) %>% 
                                           select(LCmeanDist) %>% 
                                           drop_na %>% 
                                           as.vector %>% unlist %>% unname,
                                         g = genotype_transect %>% 
                                           filter(Population == Pop) %>%
                                           select(Pos %>% matches) %>% 
                                           drop_na %>% 
                                           as.vector %>% unlist %>% unname,
                                         n = 2),
                             method = "L-BFGS-B",
                             upper = list(p_all = 0.999999),
                             lower = list(p_all = 0.000001)) %>% 
             AIC,
           AIC_linear = mle2(linear, list(p_left = p_left,
                                          p_right = p_right),
                             data = list(x = genotype_transect %>% 
                                           filter(Population == Pop) %>% 
                                           select(LCmeanDist) %>% 
                                           drop_na %>% 
                                           as.vector %>% unlist %>% unname,
                                         g = genotype_transect %>% 
                                           filter(Population == Pop) %>%
                                           select(Pos %>% matches) %>% 
                                           drop_na %>% 
                                           as.vector %>% unlist %>% unname,
                                         n = 2),
                             method = "L-BFGS-B",
                             upper = list(p_left = 0.999999,
                                          p_right = 0.999999),
                             lower = list(p_left = 0.000001,
                                          p_right = 0.000001)) %>% 
             AIC,
           AIC_clinal = mle2(function_to_use,
                             # Priors
                             list(centre = Centre_prior,
                                  width = Width_prior ,
                                  left = p_left_shelt,
                                  right = p_right_expos),
                             # Length along the transect and genotype
                             data = list(x = genotype_transect %>% 
                                           filter(Population == Pop) %>% 
                                           select(LCmeanDist) %>% 
                                           drop_na %>% 
                                           as.vector %>% unlist %>% unname,
                                         g = genotype_transect %>% 
                                           filter(Population == Pop) %>%
                                           select(Pos %>% matches) %>% 
                                           rename(Geno = starts_with("LG")) %>% 
                                           drop_na %>% 
                                           mutate(Geno = ifelse(rep(reversed, nrow(.)), 2 - Geno, Geno)) %>% 
                                           select(Geno) %>% 
                                           as.vector %>% unlist %>% unname,
                                         n = 2),
                             # Method to use
                             "L-BFGS-B",
                             # Upper limits of the optimisation
                             upper = list(centre = Centre_max,
                                          width = Width_max,
                                          left = p_left_max,
                                          right = p_right_max),
                             # Lower limits of the optimisation
                             lower = list(centre = Centre_min,
                                          width = Width_min,
                                          left = p_left_min,
                                          right = p_right_min)) %>%
             AIC) %>% 
    select(-c(p_left_shelt, p_right_expos, contains("prior"), matches("max|min"))) %>%
    mutate(Delta_AIC_stable = AIC_clinal - AIC_stable, Delta_AIC_linear = AIC_clinal - AIC_linear,
           Model_to_select = ifelse((Delta_AIC_stable > Delta_AIC_linear) & (Delta_AIC_stable > 0), "Stable",
                                    ifelse((Delta_AIC_stable > Delta_AIC_linear) & (Delta_AIC_stable < 0), "Clinal",
                                           ifelse((Delta_AIC_stable < Delta_AIC_linear) & (Delta_AIC_linear > 0), "Linear",
                                                  ifelse((Delta_AIC_stable < Delta_AIC_linear) & (Delta_AIC_linear < 0), "Clinal"))))) %>% 
    select(-starts_with("AIC_")) %>% 
    rename(Position = Pos,
           Population = Pop)
  
  return(Comp_table)
  
}


# Make a function to plot the graphs of the clines that have been optimised
#rm(cline_params, genetic_data, SNP_subset, genotype_transect, Cline_along_transect, goodness_of_fit, position_i, temp_data, population, geno_pop_pos_i, is_reversed, real_distance,Distance, Cline_param_pop_pos_i, cline_transect_pop_pos_i, model_to_use)
plot_clines <- function(cline_params, extreme_values = NULL, real_distance = FALSE, goodness_of_fit = FALSE, ...){
  # Get the genotype along the transect
  genotype_transect <- get_genotype_transect(...)
  
  # Iterate over each SNP in each position to get the things we want
  Cline_along_transect <- tibble()
  if (goodness_of_fit){
    if (!real_distance){
      stop("You need to use the real distance to do a goodness of fit test")
    }
    Goodness_of_fit_reg <- tibble()
  }
  for (position_i in cline_params$Position %>% unique){
    temp_data <- tibble()
    for (population in cline_params$Population %>% unique){
      # Select the genotype and position of each individual along the transect for each SNP in each population
      geno_pop_pos_i <- genotype_transect %>% 
        filter(Population == population) %>% 
        select(Sample_Name, LCmeanDist, position_i %>% matches)
      
      is_reversed <- cline_params %>%
        ungroup %>% 
        filter(Position == position_i,
               Population == population) %>% 
        select(reversed)
      
      # Select the model to use to plot the data
      model_to_use <- cline_params %>%
        ungroup %>% 
        filter(Position == position_i,
               Population == population) %>% 
        select(Model_to_select)

      # Select which distance to use
      if (real_distance){
        Distance <- geno_pop_pos_i %>% select(LCmeanDist) %>% as.vector %>% unlist %>% unname
      }else{
        Distance <- seq(from = geno_pop_pos_i %>% select(LCmeanDist) %>% min,
                        to = geno_pop_pos_i %>% select(LCmeanDist) %>% max,
                        by = 1)
      }
      # Isolate the parameters we need for each SNP in each population
      Cline_param_pop_pos_i <- cline_params %>% 
        filter(Population == population,
               Position == position_i) %>% 
        mutate(Left = ifelse(is_reversed, 1 - Left, Left),
               Right = ifelse(is_reversed, 1 - Right, Right))
      if (model_to_use == "Stable"){
        cline_transect_pop_pos_i <- stable(x = Distance,
                                           p_all = Cline_param_pop_pos_i$Stable_model_fit,
                                           optimisation = FALSE,
                                           plotting = TRUE) %>% 
          rename(!!quo_name(position_i) := phen_cline,
                 LCmeanDist = position) %>% 
          mutate(Population = population)
      }else if (model_to_use == "Linear"){
        cline_transect_pop_pos_i <- linear(x = Distance,
                                           p_left = Cline_param_pop_pos_i$p_left,
                                           p_right = Cline_param_pop_pos_i$p_right,
                                           optimisation = FALSE,
                                           plotting = TRUE) %>% 
        rename(!!quo_name(position_i) := phen_cline,
               LCmeanDist = position) %>% 
          mutate(Population = population)
      }else{
        cline_transect_pop_pos_i <- clinef(x = Distance,
                                         centre = Cline_param_pop_pos_i$Centre,
                                         width = Cline_param_pop_pos_i$Width,
                                         left = Cline_param_pop_pos_i$Left,
                                         right = Cline_param_pop_pos_i$Right,
                                         optimisation = FALSE,
                                         plotting=TRUE) %>% 
        rename(!!quo_name(position_i) := phen_cline,
               LCmeanDist = position) %>% 
        mutate(!!quo_name(position_i) := ifelse(rep(is_reversed, nrow(.)), 1 - !!as.name(position_i), !!as.name(position_i)),
               Population = population)
        }
      if (real_distance){
        cline_transect_pop_pos_i <- cline_transect_pop_pos_i %>% 
          arrange(LCmeanDist) %>% 
          cbind(geno_pop_pos_i %>%
                  arrange(LCmeanDist) %>%
                  select(Sample_Name))
      }
      
      # Store the results in a temporary variable
      if ((temp_data %>% nrow) == 0){
        temp_data <- cline_transect_pop_pos_i
      }else{
        temp_data <- temp_data %>% 
          rbind(cline_transect_pop_pos_i)
      }
      # Get the goodness of fit for each regression
      if (goodness_of_fit){
        if (!real_distance){
          stop("You need to use the real distance to do a goodness of fit test")
        }
        # First, we get the regression and the genotype
        cline_and_geno_along_transect_pop_pos_i <- cline_transect_pop_pos_i %>% 
          left_join(geno_pop_pos_i, by = c("Sample_Name", "LCmeanDist"), suffix = c("", "_geno")) %>% 
          drop_na %>% 
          rename(Genotype := !!quo_name(paste(position_i, "geno", sep="_")),
                 Freq := !!quo_name(position_i))
        # Then, we fit a glm to the predicted frequency
        good_fit_pop_pos_i <- cline_and_geno_along_transect_pop_pos_i %>% 
          glm(data = ., (Genotype/2) ~ Freq, family = binomial(link = "logit"))
        
        list_to_add_to_table_goodness_of_fit <- c(population, position_i, good_fit_pop_pos_i$deviance,
                                                  100 * (good_fit_pop_pos_i$null.deviance - good_fit_pop_pos_i$deviance) / good_fit_pop_pos_i$null.deviance,
                                                  Cline_param_pop_pos_i$Right - Cline_param_pop_pos_i$Left)
        if((Goodness_of_fit_reg %>% nrow) == 0){
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
    # Add this to the overall dataframe to plot this later
    if ((Cline_along_transect %>% nrow) == 0){
      Cline_along_transect <- temp_data
    }else{
      if (real_distance){
        Cline_along_transect <- Cline_along_transect %>% 
          left_join(temp_data, by = c("Population", "Sample_Name", "LCmeanDist"))
      }else{
        Cline_along_transect <- Cline_along_transect %>% 
          left_join(temp_data, by = c("Population", "LCmeanDist"))
      }
    }
  }
  
  # Once the calculation of all the plotting parameters is done, plot the calculated data
  Cline_return <- Cline_along_transect %>% 
    mutate(Population = Population %>% 
             factor(levels = c("Sweden", "France"))) %>% 
    group_by(Population) %>% 
    pivot_longer(cols = starts_with("LG"),
                 names_to = "Position",
                 values_to = "Frequency") 
  
  if (!is.null(extreme_values)){
    Cline_return <- Cline_return %>% 
    left_join(extreme_values %>% 
                select(-SNP_name),
              by = "Position")
  }
  
  # Plot the correlation between goodness of fit and delta freq
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


select_good_SNPs <- function(df, colName){
  
  colName <- as.name(substitute(colName))
  df %>% 
    arrange(Position) %>% 
    filter(!!as.symbol(colName) >= 0) %>% 
    mutate(pos_bis = lag(Position, default = "default.position")) %>% 
    filter(str_split_fixed(Position, "\\.", 2)[, 1] != str_split_fixed(pos_bis, "\\.", 2)[, 1]) %>%
    select(-pos_bis) %>% 
    return()
}


select_clinal_SNPs <- function(df){
  for (position_i in df$Position %>% unique){
    is_clinal <- df %>% 
      filter(Position == position_i,
             Model_to_select == "Clinal") %>% 
      nrow
    if (is_clinal == 0){
      df <- df %>% 
        filter(Position != position_i)
    }
  }
  return(df)
}

print("Finished function importation")

