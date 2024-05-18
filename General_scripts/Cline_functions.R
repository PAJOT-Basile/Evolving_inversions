# This file is used to store accessible functions that are used to 
# fit some data (phenotypic or genetic) to some cline models using maximum of 
# likelihood optimization.
# It is separated into two parts with in each, some functions to be used


################################################################################
################### 1. Genetic frequency variations  ###########################
################################################################################
stable <- function(x, p_all, g=0, n=0, optimisation = TRUE){
  #' This function estimates a stable frequency (constant value of allelic
  #' frequency) along the transect
  #' 
  #' Arguments:
  #'  x: vector of numeric positions of the individuals along the transect
  #'  p_all: numeric value of the allelic frequency across the transect
  #'  g: vector of genotype values for each individual along the transect
  #'  n: numeric value telling which maximal genetic value can be reached (1 or 2)
  #'  optimisation: boolean allowing the function to know if it is used in an optimisation
  #'    phase (default) or if it is used to plot the results. To plot the results of this
  #'    model, simply change the value of the optimisation to FALSE
  #' 
  #' Returns
  #'  minusll: (default return) This is opposite of the likelihood of the used parameters
  #'           to explain the observed genotypes
  #'  phen_cline: dataframe. If we choose to plot the output of the model, it returns
  #'           a dataframe containing the position of the individuals on the transect
  #'           as well as the allelic frequency along the transect
  
  
  fx <- p_all 
  # If we are optimising the function, we return the likelihood
  if (isTRUE(optimisation)){
    # The likelihood calculation relies on the hypothesis that g is binomial  
    # from a sample of n alleles with frequency fx
    # This effectively assumes HWE.
    minusll <- -sum(dbinom(g, n, fx, log=TRUE))
    return(minusll)
  }
  
  # If we are not optimising, then, we are plotting the results
  if (!isTRUE(optimisation)){
    phen_cline <- data.frame(phen_cline = fx, position = x)
    return(phen_cline)
  }
}

linear <- function(x, p_left, p_right, g=0, n=0, optimisation=TRUE){
  #' This function estimates a linear frequency variation along the transect
  #' Arguments:
  #'  x: vector of numeric positions of the individuals along the transect
  #'  p_left: numeric value of the allelic frequency on the left part of the transect
  #'  p_right: numeric value of the allelic frequency on the right part of the transect
  #'           (the optimisation of parameters works best if p_left < p_right)
  #'  g: vector of genotype values for each individual along the transect
  #'  n: numeric value telling which maximal genetic value can be reached (1 or 2)
  #'  optimisation: boolean allowing the function to know if it is used in an optimisation
  #'    phase (default) or if it is used to plot the results. To plot the results of this
  #'    model, simply change the value of the optimisation to FALSE
  #' 
  #' Returns
  #'  minusll: (default return) This is opposite of the likelihood of the used parameters
  #'           to explain the observed genotypes
  #'  phen_cline: dataframe. If we choose to plot the output of the model, it returns
  #'           a dataframe containing the position of the individuals on the transect
  #'           as well as the allelic frequency along the transect
  
  fx <- p_left + (p_right - p_left) * (x - min(x)) / (max(x) - min(x))
  # If the optimization argument is set to true (default), the calculation of 
  # maximum of likelihood is returned
  if (isTRUE(optimisation)){
    # The likelihood calculation relies on the hypothesis that g is binomial  
    # from a sample of n alleles with frequency z_x
    # This effectively assumes HWE.
    minusll <- -sum(dbinom(g, n, fx, log=TRUE))
    return(minusll)
  }
  # If the plotting argument is true, we return the data to visualize it
  if (!isTRUE(optimisation)){
    phen_cline <- data.frame(phen_cline = fx, position = x)
    return(phen_cline)
  }
}

clinef <- function(x, g=0, n=0, centre, width, left, right, optimisation=TRUE){
  #' This function estimates a sigmoid cline for allelic frequencies with no transformation. 
  #'
  #' Arguments:
  #'  x: vector of numeric positions of the individuals along the transect
  #'  centre: numeric value indicating the center of the cline
  #'  width: numeric value indicating the width of the cline
  #'  left: numeric value of the allelic frequency on the left part of the transect
  #'  right: numeric value of the allelic frequency on the right part of the transect
  #'           (the optimisation of parameters works best if p_left < p_right)
  #'  g: vector of genotype values for each individual along the transect
  #'  n: numeric value telling which maximal genetic value can be reached (1 or 2)
  #'  optimisation: boolean allowing the function to know if it is used in an optimisation
  #'    phase (default) or if it is used to plot the results. To plot the results of this
  #'    model, simply change the value of the optimisation to FALSE
  #' 
  #' Returns
  #'  minusll: (default return) This is opposite of the likelihood of the used parameters
  #'           to explain the observed genotypes
  #'  phen_cline: dataframe. If we choose to plot the output of the model, it returns
  #'           a dataframe containing the position of the individuals on the transect
  #'           as well as the allelic frequency along the transect
  
  # The first step is to get the end frequencies on a proportional scale for the
  # studied transect
  # To do this, we calculate the natural width of the cline: d
  d <- x - centre
  # p_x is the frequency cline between 0 and 1
  p_x <- 1 / (1 + exp(0 - 4 * (d) / width))
  # z_x is the classic cline function adapted to the observed data frequencies. 
  z_x <- left + (right - left) * p_x
  # If the optimization argument is set to true (default), the calculation of 
  # maximum of likelihood is returned
  if (isTRUE(optimisation)){  
    # The likelihood calculation relies on the hypothesis that g is binomial  
    # from a sample of 2 alleles with frequency z_x
    # This effectively assumes HWE. 
    minusll <- -sum(dbinom(g, n, z_x, log=TRUE))
    return(minusll)
  }
  # If the plotting argument is true, we return the data to visualize it
  if (!isTRUE(optimisation)){
    # Here, we directly add the calculated data into a data frame, so we can 
    # plot it
    phen_cline <- data.frame(phen_cline = z_x, position = x)
    return(phen_cline)
  }

}

clineflog <-  function(x, g, n, centre, width, left, right){
  #' This function estimates a sigmoid cline for allelic frequencies with a
  #' logarithmic transformation. Warning, by comparison with the other 
  #' optimisation functions, this one can not be used to plot the allelic frequency
  #' variations along the transect.
  #'
  #' Arguments:
  #'  x: vector of numeric positions of the individuals along the transect
  #'  centre: numeric value indicating the center of the cline
  #'  width: numeric value indicating the logarithm of the cline width
  #'  left: numeric value of the logit allelic frequency on the left part of the transect
  #'  right: numeric value of the logit allelic frequency on the right part of the transect
  #'           (the optimisation of parameters works best if p_left < p_right)
  #'  g: vector of genotype values for each individual along the transect
  #'  n: numeric value telling which maximal genetic value can be reached (1 or 2)
  #' 
  #' Returns
  #'  minusll: This is opposite of the likelihood of the used parameters
  #'           to explain the observed genotypes
  
  # First, we transform the values of the width and the allelic frequencies
  wi <- exp(width)
  r <- exp(right) / (1 + exp(right))
  l <- exp(left) / (1 + exp(left))
  
  # We then calculate the theoretical allelic frequency expected under a clinal
  # model using the input parameters
  d <- x - centre
  p_x <- 1/(1 + exp(0 - 4 * (d) / wi))
  z_x <- l + (r - l) * p_x  
  # likelihood calculation, assuming g is a binomial from sample of 2 with frequency z_x
  # this effectively assumes HWE
  minusll <- -sum(dbinom(g, n, z_x, log=TRUE))
  return(minusll) 
}

################################################################################
####################### 2. Phenotypic cline  ###################################
################################################################################
cline_phen <- function(x, phen, centre, width, left, right, sl, sc, sr, optimisation=TRUE){
  #' This function estimates a sigmoid cline for phenotypic variations
  #'
  #' Arguments:
  #'  x: vector of numeric positions of the individuals along the transect
  #'  centre: numeric value indicating the center of the cline
  #'  width: numeric value indicating the width of the cline
  #'  left: numeric value of the phenotypic trait on the left part of the transect
  #'  right: numeric value of the phenotypic trait on the right part of the transect
  #'           (the optimisation of parameters works best if p_left < p_right)
  #'  sl: numeric value representing the standard deviation on the left phenotypic trait value
  #'  sc: numeric value representing the standard deviation on the center of the cline
  #'  sr: numeric value representing the standard deviation on the right phenotypic trait value
  #'  n: numeric value telling which maximal genetic value can be reached (1 or 2)
  #'  optimisation: boolean allowing the function to know if it is used in an optimisation
  #'    phase (default) or if it is used to plot the results. To plot the results of this
  #'    model, simply change the value of the optimisation to FALSE
  #' 
  #' Returns
  #'  minusll: (default return) This is opposite of the likelihood of the used parameters
  #'           to explain the observed genotypes
  #'  phen_cline: dataframe. If we choose to plot the output of the model, it returns
  #'           a dataframe containing the position of the individuals on the transect,
  #'           the allelic frequency along the transect and the standard deviation of
  #'           the phenotypic cline
  
  # The first step is to calculte the natural width of the cline: d
  d <- position - centre
  # p_x is the frequency cline between 0 and 1
  p_x <- 1 / (1 + exp(0 - 4 * (d) / width))
  # z_x is the classic cline function adapted to the observed data frequencies. 
  z_x <- left + (right-left)*p_x  
  # s_x is the variance of the model (assumes variances are additive) 
  s_x <- sqrt(sl^2 + 4 * p_x * (1 - p_x) * sc^2 + (p_x^2) * (sr^2 - sl^2))
  # If the optimization argument is set to true (default), the calculation of 
  # maximum of likelihood is returned
  if (isTRUE(optimisation)){
    # The likelihood calculation relies on the hypothesis that g is binomial  
    # from a sample of 2 alleles with frequency z_x
    # This effectively assumes HWE. 
    minusll <- -sum(dnorm(phen, z_x, s_x, log=TRUE))
    return(minusll)
  }
  # If the plotting argument is true, we return the data to visualize it
  if (!isTRUE(optimisation)){
    # Here, we directly add the calculated data into a data frame, so we can 
    # plot it
    phen_cline <- data.frame(phen_cline = z_x, sd_cline = s_x, position = position)
    return(phen_cline)
  }
}




