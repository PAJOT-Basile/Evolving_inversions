# This file is used to store accessible functions that are used to 
# fit some data (phenotypic or genetic) to some cline models using maximum of 
# likelihood optimization.
# It is separated into two parts with in each, some functions to be used


################################################################################
################### 1. Genetic frequency clines  ###############################
################################################################################
stable <- function(x, p_all, g=0, n=0, optimisation = TRUE, plotting = FALSE){
  #' This function estimates a stable frequency variation along the transect
  fx <- p_all 
  # If we are optimising the function, we return the likelihood
  if (isTRUE(optimisation)){
    minusll <- -sum(dbinom(g, n, fx, log=TRUE))
    return(minusll)
  }
  if (isTRUE(plotting)){
    phen.cline <- data.frame(phen_cline = fx, position = x)
    return(phen.cline)
  }
}

linear <- function(x, p_left, p_right, g=0, n=0, optimisation=TRUE, plotting=FALSE){
  #' This function estimates a linear frequency variation along the transect
  fx <- p_left + (p_right - p_left) * (x - min(x)) / (max(x) - min(x))
  # If the optimization argument is set to true (default), the calculation of 
  # maximum of likelihood is returned
  if (isTRUE(optimisation)){
    # The likelihood calculation relies on the hypothesis that g is binomial  
    # from a sample of 2 alleles with frequency z_x
    # This effectively assumes HWE.
    minusll <- -sum(dbinom(g, n, fx, log=TRUE))
    return(minusll)
  }
  # If the plotting argument is true, we return the data to visualize it
  if (isTRUE(plotting)){
    phen_cline <- data.frame(phen_cline = fx, position = x)
    return(phen_cline)
  }
}

clinef <- function(x, g=0, n=0, centre, width, left, right, optimisation=TRUE, plotting=FALSE){
  #' This function estimates a sigmoid cline with no transformation. It can also
  #' return the data that can be used to visualize the optimization
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
  if (isTRUE(plotting)){
    # Here, we directly add the calculated data into a data frame, so we can 
    # plot it
    phen_cline <- data.frame(phen_cline = z_x, position = x)
    return(phen_cline)
  }

}
clineflog <-  function(x, g, n, centre, width, left, right){
  wi <- exp(width)
  c <- exp(right) / (1 + exp(right))
  w <- exp(left) / (1 + exp(left))
  d <- x - centre
  p_x <- 1/(1 + exp(0 - 4 * (d) / wi))
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right
  z_x <- w + (c - w) * p_x  
  # z_x is expected frequency
  minusll <- -sum(dbinom(g, n, z_x, log=TRUE))
  # likelihood calculation, assuming g is a binomial from sample of 2 with frequency z_x
  # this effectively assumes HWE
  return(minusll) 
}

################################################################################
####################### 2. Phenotypic clines  ##################################
################################################################################
cline_phen <- function(phen, position, centre, w, left, right, sl, sc, sr, optimisation=TRUE, plotting=FALSE){
  #' This function allows to optimise a cline model using phenotypic values
  # The first step is to calculte the natural width of the cline: d
  d <- position - centre
  # p_x is the frequency cline between 0 and 1
  p_x <- 1 / (1 + exp(0 - 4 * (d) / w))
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
  if (isTRUE(plotting)){
    phen_cline <- data.frame(phen_cline = z_x, sd_cline = s_x, position = position)
    # Here, we directly add the calculated data into a data frame, so we can 
    # plot it
    return(phen_cline)
  }
}
