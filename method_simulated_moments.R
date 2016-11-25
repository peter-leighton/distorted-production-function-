# Name:     master_peter.R
# Date:     22/11/2016 (SC)
# Note:     This solves the model for optimal I*


########################################################
# PACKAGES
########################################################
#install.packages("rootSolve")
#install.packages("nleqslv")
library("rootSolve")
library("nleqslv")
library(haven)

########################################################
# DIRECTORIES 
########################################################
dir     <- c("P:/ECD Colombia/DistortedPF")
setwd(dir)

########################################################
# LOAD DATA 
########################################################
data <- read_dta("distortedPF.dta")
n       <- nrow(data)
options(scipen=10)


########################################################
# DEFINE PARAMETER VALUES
########################################################

#theta      <- 0.5 
#phi        <- 0.5 
#alpha      <- 0.5 
#beta       <- 0.5 
#var_xi     <- 1 
param_set  <- c(theta, phi, alpha, beta, var_xi)

########################################################
# CREATE SOME DATA TO TEST THE PROGRAM 
########################################################

#n         <- 100
#H0        <- runif(n, min = 0, max = 10)
#Y         <- runif(n, min = 10, max = 100)
#P         <- 1
#xi        <- rlnorm(n, mean = 0, sd = var_xi)


########################################################
# DEFINING VARIABLES FROM DATA  
########################################################

## Take data set columns and put them into vectors ###

Y <-  data[, 1]
Y <- unlist(Y)
Y <- as.vector(Y, mode='numeric')

P <-  data[, 2]
P <- unlist(P)
P <- as.vector(P, mode='numeric')

H0 <-  data[, 3]
H0 <- unlist(H0)
H0 <- as.vector(H0, mode='numeric')




true_moments <- c(0.24168653, 2.96648602, 0.04080262, 30.16635837, 505914.42609943, 0.01764332, 0.66210755, 505914.42609943) #set population moments in data# 

## This calculates simulated moments (as before) and then returns distance to specified true moments ##

ma.msm.objective <- function(param_set, Y, P, HO) {
  
  set.seed(24112016)
  xi <- rlnorm(n, mean = 0, sd = var_xi)
  
  foc.inv <- function(I) {
    func <-  1/(beta*(1 - theta) - 1) * (log(P[i]) - theta * log(Y[i] - P[i]*I) - log(beta) - alpha*(1 - phi) * log(H0[i]) - alpha*(1 - phi)*log(xi[i])) - log(I)
    func2 <- Vectorize(func) # Need to turn into a vector to get to work with uniroot.all # 
    return(func2)
  }
  
  optinv <- numeric(length = n)
  
  # Solving for everyone in the data # 
  for (i in 1:n) {
    optinv[i]   <- uniroot.all(foc.inv, lower =- 1, upper = 49)
  }
  
  H <- numeric(length = n)
  for (i in 1:n) { 
    H[i] <- xi[i] * H0[i]^alpha * optinv[i]^beta
  } 
  
  mean_I <- mean(optinv, na.rm=TRUE)
  mean_H <- mean(H, na.rm=TRUE)
  
  
  var_I <- var(optinv, na.rm=TRUE)
  var_H <- var(H, na.rm=TRUE)
  
  covI_Y <- cov(optinv, Y)
  covI_H0 <- cov(optinv, H0)
  covI_H <- cov(optinv, H)
  
  simulated_moments <- c(mean_I, mean_H, var_I, var_H, covI_Y, covI_H0, covI_H, covI_Y)
  
  diff <- simulated_moments - true_moments 
  length <- length(diff)
  weight <- diag(length) #weighting function - for now identity matrix # 

  objective <- t(diff)%*%weight%*%diff #distance between simulated and true moments#
  
  return(objective)
}


## This should choose parameters to minimise the above objective function 

ma.msm.est <- function(param_set, Y, P, HO, xi) {
  
  # Taking initial guess at parameters # 

  theta.0      <- 0.6 
  phi.0        <- 0.6 
  alpha.0      <- 0.6 
  beta.0       <- 0.6
  var_xi.0     <- 1.2
  
  fit <- nlm(ma.msm.objective, c(theta.0, phi.0, alpha.0, beta.0, var_xi.0), Y=Y, P=P)
  return(fit)
}

