#! /usr/bin/Rscript

#########################################
## Detection of genomic islands        ##
## from Fst values (log10-transformed) ##
## using a Hidden Markov model with    ##
## 2 normally distributed states       ##
#########################################

# based on Hofer et al. 2012 BMC Genomics and Soria-Carrasco et al. 2014 Science
# adapted by David Alexander Marques
# please cite both Hofer et al. 2012 BMC Genomics, Soria-Carrasco et al. 2014 Science
# and Marques et al. 2016 Molecular Ecology

# Changelog
# 2018-01-31: fixed typos

# USAGE
# HMM_log_delta_freqs_modif_2norms.R filename.suffix ncores [nstart]
# filename.suffix is a file with FST values with a header 'x'
# ncores is the number of cores allowed for parallel computation
# nstart (optional) is the number of random starting parameters in the parameter estimation [default=1000]

# REQUIRED
# R-libraries HiddenMarkov, foreach, doParallel

# OUTPUT
# - filename_inputdata_histogram.pdf : Histogram of the z-transformed q-values (input data), should look approximately normally distributed and continuous, otherwise consider data filtering (minor allele frequency filter) or using an other appraoch.
# - filename_HMMparest.txt : Text file with Baum-Welch algorithm optimized parameters from nstart [default=1000] random starting values.
# - filename_2state_parest_plot.pdf : Scatterplot of the maximized parameter estimates against their likelihood
# - filename_zscore_vs_statedistributions.pdf : Histogram of the z-transformed q-values with the distribution and parameters of the three states.
# - filename_bestparameters.txt : Text file with the best HMM parameters used for the state reconstruction
# - filename_2state_HMMstates.txt : Text file with the HMM states from the state reconstruction using the best parameters

##################
# Load libraries #
##################
# Add a miror to download the libraries if needed
utils::setRepositories(ind = 0, addURLs = c(CRAN = "https://cran.irsn.fr/"))
# Install libraries if needed and load them
libraries <- c("HiddenMarkov", "foreach", "doParallel")

if (!require("pacman")) install.packages("pacman")
for (lib in libraries){
  pacman::p_load(lib, character.only = TRUE)
}
#library("HiddenMarkov")
#library("foreach")
#library("doParallel")

##################################
# Read input parameters and data #
##################################

# Reads input parameters
args <- commandArgs(trailingOnly = T)

# Checks if enough input parameters were given
if (length(args) < 2){
  stop("Usage: HMM_log_delta_freqs_modif_2norms.R filename.suffix ncores [nstart]\n\t\t- filename.suffix is a file with FST values with a header 'x'\n\t\t- ncores is the number of cores allowed for parallel computation\n\t\t- nstart (optional) is the number of random starting parameters in the parameter estimation [default=1000]") 
}

# Defines base name for output files
elements <- strsplit(args[1], split = "\\.")
base <- paste(elements[[1]][-length(elements[[1]])], collapse = ".")
rm(elements)

# Loads input data
data <- read.table(file = args[1], header = T, quote = "\"")

##############################
# LOG10-transform FST values #
##############################

# use FST+1 instead of FST to avoid negative numbers in transformation
lfst <- log10(abs(data$WEIR_AND_COCKERHAM_FST) + 1)

# Prints distribution of log10-transformed FST+1 values
pdf(paste(base, "_inputdata_histogram.pdf", sep = ""))
hist(lfst, breaks = 50, col = 1, xlab = "log10(FST+1)", main = "HMM input data distribution", freq = F)
dev.off()

#############################
## PARAMETER ESTIMATION BY ##
## BAUM-WELCH ALGORITHM    ##
#############################

# Parallelizes parameter estimation from args[3] random starting parameters

# Defines the number of cores to be used
registerDoParallel(cores = args[2])

# Default number of random starting parameters
if(is.na(args[3])){args[3] = 1000}

# Runs parameter estimation in parallel on args[2] cores
parout <- foreach(i = 1:as.integer(args[3]), .packages = "HiddenMarkov", .combine = 'c') %dopar% {
  print(i); set.seed(i)
  # Samples random initial parameters for the transition matrix (trans),
  #   marginal/initial probabilities (init), for the state means (means)
  #   and state standard deviations (sdevs)
  #prob <- runif(4, min = 0, max = 1)
  #prob2 <- runif(1, min = 0, max = (1 - prob[2]))
  #prob4 <- runif(1, min = 0, max = (1 - prob[4]))
  #mean2 <- runif(1, min = quantile(lfst, c(0.4)), max = quantile(lfst, c(0.6)))
  #mean1 <- runif(1, min = quantile(lfst, c(0.2)), max = mean2)
  #mean3 <- runif(1, min = mean2, max = quantile(lfst, c(0.8)))
  #sdevs <- runif(3, min = 0.5, max = 2)
  #trans <- matrix(c(prob[1], 1 - prob[1], 0,
  #                prob[2], prob2,(1 - (prob[2] + prob2)),
  #                0, 1 - prob[3], prob[3]), byrow = T, nrow = 3)
  #init <- c(prob[4], prob4, (1 - (prob[4] + prob4)))
  #means <- c(mean1, mean2, mean3)
  prob <- runif(3, min = 0, max = 1)
  mean <- runif(1, min = quantile(lfst, c(0.1)), max = quantile(lfst, c(0.5)))
  mean2 <- runif(1, min = mean, max = quantile(lfst, c(0.8)))
  sdevs <- runif(2, min = 0.5, max = 2)
  trans <- matrix(c(prob[1], 1 - prob[1],
                    prob[2], 1 - prob[2]),
                  byrow = TRUE, nrow = 2)
  init <- c(prob[3], 1 - prob[3])
  means <- c(mean, mean2)
  
  # Builds Hidden Markov Model with random initial parameters
  myhmm <- dthmm(lfst, trans, init, "norm", list(mean = means, sd = sdevs), discrete = F)
  
  # Optimizes parameters with Baum-Welch algorithm, with 3 additional runs to find maximal estimates
  # Baum-Welch configuration
  a <- bwcontrol(maxiter = 1000, tol = 1e-07, prt = F, posdiff = F)
  bwhmm2 <- try(BaumWelch(myhmm, control = a), silent = T)
  bwhmm1 <- try(BaumWelch(bwhmm2, control = a), silent = T)
  bwhmm <- try(BaumWelch(bwhmm1, control = a), silent = T)
  
  # Output parameters
  if(length(bwhmm) > 1){
    c(bwhmm$Pi, bwhmm$delta, bwhmm$pm$mean, bwhmm$pm$sd, bwhmm$LL, bwhmm$iter)
  }else{
    rep(NA, 20)
  }
}

hmmpar <- as.data.frame(matrix(parout, ncol = 20,byrow = T))
names(hmmpar) <- c("a11", "a21", "a12", "a22", "i1", "i2", "m1", "m2",
                   "sd1", "sd2", "LL", "iter")

# Writes all parameter estimates to a single file
write.table(hmmpar, file = paste(base, "_HMMparest.txt", sep = ""))
hmmpar <- read.table(paste(base, "_HMMparest.txt", sep = ""), stringsAsFactors = F)

########################################
## EVALUATION OF PARAMETER ESTIMATION ##
########################################

# Plots parameter estimates against their likelihood
pdf(file = paste(base, "_2state_parest_plot.pdf", sep=""), height = 16, width = 10)
par(mfrow = c(6, 2))
attach(hmmpar)
for(k in colnames(hmmpar)[-(length(colnames(hmmpar)) - c(0, 1))]){
  plot(hmmpar$LL, get(k), xlab = "Likelihood", ylab = k, pch = 20, col = "#88888888", xlim = c(min(hmmpar$LL, na.rm = T),
       max(hmmpar$LL, na.rm = T)))
  abline(h = with(hmmpar, get(k)[LL == max(LL, na.rm = T)]), col = "red")
}
dev.off()

# Store the best initial parameters (the one with the highest likelihood)
bestpar <- head(hmmpar[with(hmmpar, order(-hmmpar$LL)), ], 1)

# Run another parameter optimization round to ensure the maximum likelihood is reached
{
  # Builds Hidden Markov Model with best parameters
  myhmm <- dthmm(lfst, matrix(c(bestpar$a11, bestpar$a12,
                             bestpar$a21, bestpar$a22), byrow=T, nrow=2),
               c(bestpar$i1, bestpar$i2), "norm",
               list(mean = c(bestpar$m1, bestpar$m2), sd = c(bestpar$sd1, bestpar$sd2)), discrete = F)
  
  # Run Baum-Welch algorithm 3 times in a row to maximize parameter estimates
  # Baum-Welch configuration
  a <- bwcontrol(maxiter = 1000, tol = 1e-07, prt = F, posdiff = F)
  bwhmm2 <- try(BaumWelch(myhmm, control = a), silent = T)
  bwhmm1 <- try(BaumWelch(bwhmm2, control = a), silent = T)
  bwhmm <- try(BaumWelch(bwhmm1, control = a), silent = T)  
}

# Saves new best parameters to file
bestpar <- c(bwhmm$Pi, bwhmm$delta, bwhmm$pm$mean, bwhmm$pm$sd, bwhmm$LL, bwhmm$iter)
names(bestpar) <- c("a11", "a21", "a12", "a22",
                    "i1", "i2", "m1", "m2", "sd1", "sd2", "LL", "iter")
write.table(bestpar, file = paste(base, "_bestparameters.txt", sep = ""))

# Plots the data (log10-transformed FST values) against the inferred state distributions
pdf(file = paste(base, "_log10FST+1_vs_statedistributions.pdf", sep = ""))
par(mfrow = c(1, 1))
hist(lfst,breaks = 50, col = "grey", main = "Data vs. Emission & Transition Probabilities",
     border = F, xlab = "log10(FST+1)", freq = F)
for(i in 1:length(bwhmm$pm$mean)){
  points(seq(par("usr")[1], par("usr")[2], length.out = 1000),
         dnorm(seq(par("usr")[1], par("usr")[2], length.out = 1000),
         mean = bwhmm$pm$mean[i], sd = bwhmm$pm$sd[i]), type = "l", col = c(4, 1, 2)[i], lwd = 2)
}
legend("topright", legend = c(paste("state1 mean =", round(bwhmm$pm$mean[1],2), " sd =", round(bwhmm$pm$sd[1], 2)),
                           paste("state2 mean =", round(bwhmm$pm$mean[2], 2), " sd =", round(bwhmm$pm$sd[2], 2)),
                           paste("1 > 1: ", round(bwhmm$Pi[1, 1], 2), "; 1 > 2: ", round(bwhmm$Pi[1, 2], 2)),
                           paste("2 > 1: ", round(bwhmm$Pi[2, 1], 2), "; 2 > 2: ", round(bwhmm$Pi[2, 2], 2))),
       fill = c(4, 1, 2, NA, NA, NA), border = c(1, 1, 1, 0, 0, 0), bty = "n")
dev.off()

##############################
## RECONSTRUCTION OF STATES ##
##############################

# Reconstruction of states with the Viterbi algorithm
states <- Viterbi(bwhmm)

# Output state reconstruction results
write.table(states, file = paste(base, "_2state_HMMstates.txt", sep = ""))
