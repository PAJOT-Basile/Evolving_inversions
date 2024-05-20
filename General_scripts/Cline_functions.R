# This file is used to store accessible functions that are used to 
# fit some data (phenotypic or genetic) to some cline models using maximum of 
# likelihood optimization.
# It is separated into two parts with in each, some functions to be used

# Libraries
#install.packages("anyLib")
require("anyLib")
anyLib("docstring")
################################################################################
################### 1. Genetic frequency variations  ###########################
################################################################################
stable <- function(x, p_all, g=0, n=0, optimisation = TRUE){
  #' Represent stable distribution
  #' 
  #' This function estimates a stable frequency (constant value of allelic
  #' frequency) along the transect
  #' 
  #'@param x (numeric vector).
  #'    Positions of the individuals along the transect
  #'@param p_all (numeric value).
  #'    Allelic frequency across the transect
  #'@param g (numeric vector)
  #'    Genotype values for each individual along the transect
  #'@param n (numeric value)
  #'    Maximal genetic value can be reached (1 or 2)
  #'@param optimisation (boolean).
  #'    This argument allows the function to know if it is used in an optimisation
  #'    phase (default) or if it is used to plot the results. To plot the
  #'    results of this model, simply change the value of the optimisation to FALSE
  #' 
  #'@returns (numeric value).(default return)
  #'    This is opposite of the likelihood of the used parameters to explain the
  #'    observed genotypes
  #'@returns (dataframe.)
  #'    If we choose to plot the output of the model, it returns a dataframe
  #'    containing the position of the individuals on the transect as well as
  #'    the allelic frequency along the transect
  #' @export
  
  
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
  #' Represent linear distribution
  #' 
  #' This function estimates a linear frequency variation along the transect.
  #' 
  #'@param x (numeric vector).
  #'    Positions of the individuals along the transect
  #'@param p_left (numeric value). 
  #'    Allelic frequency on the left part of the transect.
  #'@param p_right (numeric value). 
  #'    Allelic frequency on the right part of the transect.
  #'@param g (numeric vector)
  #'    Genotype values for each individual along the transect
  #'@param n (numeric value)
  #'    Maximal genetic value can be reached (1 or 2)
  #'@param optimisation (boolean).
  #'    This argument allows the function to know if it is used in an optimisation
  #'    phase (default) or if it is used to plot the results. To plot the
  #'    results of this model, simply change the value of the optimisation to FALSE
  #' 
  #'@returns (numeric value).(default return)
  #'    This is opposite of the likelihood of the used parameters to explain the
  #'    observed genotypes. 
  #'@returns (dataframe.)
  #'    If we choose to plot the output of the model, it returns a dataframe
  #'    containing the position of the individuals on the transect as well as
  #'    the allelic frequency along the transect
  #' @section Warning:
  #'   The optimisation of parameters for this function works best if
  #'    p_left < p_right
  #' @export
  
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
  #' Represent clinal distribution
  #' 
  #' This function estimates a sigmoid cline for allelic frequencies with no
  #' transformation. 
  #'
  #'@param x (numeric vector).
  #'    Positions of the individuals along the transect
  #'@param centre (numeric value).
  #'    Center of the cline
  #'@param width (numeric value).
  #'    Width of the cline
  #'@param left (numeric value). 
  #'    Allelic frequency on the left part of the transect.
  #'@param right (numeric value). 
  #'    Allelic frequency on the right part of the transect.
  #'@param g (numeric vector)
  #'    Genotype values for each individual along the transect
  #'@param n (numeric value)
  #'    Maximal genetic value can be reached (1 or 2)
  #'@param optimisation (boolean).
  #'    This argument allows the function to know if it is used in an optimisation
  #'    phase (default) or if it is used to plot the results. To plot the
  #'    results of this model, simply change the value of the optimisation to FALSE
  #' 
  #'@returns (numeric value).(default return)
  #'    This is opposite of the likelihood of the used parameters to explain the
  #'    observed genotypes.
  #'@returns (dataframe.)
  #'    If we choose to plot the output of the model, it returns a dataframe
  #'    containing the position of the individuals on the transect as well as
  #'    the allelic frequency along the transect
  #' @section Warning:
  #'   The optimisation of parameters for this function works best if
  #'    p_left < p_right
  #' @export
  
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
  #' Represent clinal distribution with log transformation
  #' 
  #' This function estimates a sigmoid cline for allelic frequencies with a
  #' logarithmic transformation. Warning, by comparison with the other 
  #' optimisation functions, this one can not be used to plot the allelic frequency
  #' variations along the transect.
  #'
  #'@param x (numeric vector).
  #'    Positions of the individuals along the transect
  #'@param centre (numeric value).
  #'    Center of the cline
  #'@param width (numeric value).
  #'    Logarithm of the width of the cline
  #'@param left (numeric value). 
  #'    Logit allelic frequency on the left part of the transect.
  #'@param right (numeric value). 
  #'    Logit allelic frequency on the right part of the transect.
  #'@param g (numeric vector)
  #'    Genotype values for each individual along the transect
  #'@param n (numeric value)
  #'    Maximal genetic value can be reached (1 or 2)
  #'    
  #'@returns (numeric value)
  #'    This is opposite of the likelihood of the used parameters to explain the
  #'    observed genotypes. 
  #'@section Warning:
  #'    The optimisation of parameters for this function works best if 
  #'    p_left < p_right
  #'@export
  
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
  #' Represent phenotypic cline
  #' 
  #' This function estimates a sigmoid cline for phenotypic variations
  #'
  #'@param x (numeric vector).
  #'    Positions of the individuals along the transect
  #'@param centre (numeric value).
  #'    Center of the cline
  #'@param width (numeric value).
  #'    Width of the cline
  #'@param left (numeric value)
  #'    Phenotypic trait value on the left part of the transect
  #'@param right (numeric value)
  #'    Phenotypic trait value on the right part of the transect
  #'@param sl (numeric value).
  #'    Standard deviation on the left phenotypic trait value
  #'@param sc (numeric value).
  #'    Standard deviation on the center of the cline
  #'@param sr (numeric value).
  #'    Standard deviation on the right phenotypic trait value
  #'@param n (numeric value).
  #'    Maximal genetic value can be reached (1 or 2)
  #'@param optimisation (boolean).
  #'    This argument allows the function to know if it is used in an optimisation
  #'    phase (default) or if it is used to plot the results. To plot the
  #'    results of this model, simply change the value of the optimisation to FALSE
  #' 
  #'@returns minusll (numeric value). (default return)
  #'    This is opposite of the likelihood of the used parameters to explain the
  #'    observed genotypes
  #'@returns phen_cline (dataframe).
  #'    If we choose to plot the output of the model, it returns a dataframe
  #'    containing the position of the individuals on the transect, the allelic
  #'    frequency along the transect and the standard deviation of the 
  #'    phenotypic cline
  #'@section Warning:
  #'    This function is best optimised if left < right
  #'@export
  
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

