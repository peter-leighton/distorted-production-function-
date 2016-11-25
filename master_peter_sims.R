# Name:     master_peter_sims.R
# Date:     22/11/2016 (SC)
# Note:     This solves the model for optimal I*, but does so using 'n2' simulations for each individual.
            # I is then taken as the average over all these simulations of xi. 


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
#rm(list = ls())
dir     <- c("P:/ECD Colombia/DistortedPF")
setwd(dir)

########################################################
# LOAD DATA 
########################################################
data <- read_dta("distortedPF.dta")



########################################################
# DEFINE PARAMETER VALUES
########################################################
var_xi     <- 1 
theta      <- 0.5 
phi        <- 0.5 
alpha      <- 0.5 
beta       <- 0.5 
param_set  <- c(theta, phi, alpha, beta, var_xi)

n1 <- nrow(data)  # number of individuals#
n2 <- 100  #number of simulations per person #
N <- n1*n2  #total dimensionality of matrix # 

########################################################
# CREATE SOME DATA TO TEST THE PROGRAM 
########################################################
#n1         <- 100 
#n2        <- 100
#N <- n1*n2
#H0        <- runif(n1, min = 0, max = 10)
#Y         <- runif(n1, min = 10, max = 100)
#P         <- 1
#xi        <- rlnorm(n1, mean = 0, sd = var_xi)






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

## Creating vector of n2 realisations of xi for each person (h0, Y and P stay constant) ##

H02 <- rep(H0, each = n2) #creating n2 identical values for each individual in vector #
Y2 <- rep(Y, each = n2)
P2 <- rep(P, each = n2)
xi <- rlnorm(N, mean =0, sd = var_xi)


########################################################
# WRITE FUNCTION TO SOLVE THE MODEL 
########################################################

foc.inv <- function(I) {
  func <-  1/(beta*(1 - theta) - 1) * (log(P2[i]) - theta * log(Y2[i] - P2[i]*I) - log(beta) - alpha*(1 - phi) * log(H02[i]) - alpha*(1 - phi)*log(xi[i])) - log(I)
  func2 <- Vectorize(func) # Need to turn into a vector to get to work with uniroot.all # 
  return(func2)
}




########################################################
# SOLVE INVESTMENT FOR EVERYONE IN THE DATA 
########################################################

### for all realizations ###
optinv_full <- numeric(length = N)
for (i in 1:N) {
  optinv_full[i]   <- uniroot.all(foc.inv, lower =- 1, upper = 49)
}



## calculating an average over all realizations for each individual ##
optinv_i <- numeric(length = n1)

for (i in 1:n1) { 
  low <- i*n2 - (n2-1)  ## each individual is 1-100, 101-200 etc - so taking these rows#
  high <- i*n2
  vec <- optinv_full[low:high]
  optinv_i[i] <- mean(vec) # creating mean choice of I over 100 realizations of xi for each individual ##
} 

optinv_i

##### Calculating human capital at each individuals average xi ###


xi_mean <- numeric(length = n1)
for (i in 1:n1) { 
  low <- i*n2 - (n2-1)  ## each individual is 1-100, 101-200 etc - so taking these rows#
  high <- i*n2
  vec <- xi[low:high]
  xi_mean[i] <- mean(vec) # creating mean choice of I over 100 realizations of xi for each individual ##
} 

H <- numeric(length = n1)
for (i in 1:n1) { 
  H[i] <- xi_mean[i] * H0[i]^alpha * optinv_i[i]^beta
} 


#### Calculating human capital for all possible iterations in the data #### 
#H_full <- numeric(length = N)
#for (i in 1:N) { 
  #H_full[i] <- xi[i]*H02[i]^alpha * optinv_full[i]^beta # what to do here - what xi to include? #,
#} 


### ## calculating an average over all realizations for each individual ##
#H_i <- numeric(length = n1)

#for (i in 1:n1) { 
  #low <- i*n2 - (n2-1)  ## each individual is 1-100, 101-200 etc - so taking these rows#
  #high <- i*n2
  #vec <- H_full[low:high]
  #H_i[i] <- mean(vec) # creating mean choice of I over 100 realizations of xi for each individual ##
#} 





########################################################
# SIMULATED MOMENTS 
########################################################
mean_I <- mean(optinv_i)
mean_H <- mean(H)


var_I <- var(optinv_i)
var_H <- var(H)

covI_Y <- cov(optinv_i, Y)
covI_H0 <- cov(optinv_i, H0)
covI_H <- cov(optinv_i, H)

simulated_moments2 <- c(mean_I, mean_H, var_I, var_H, covI_Y, covI_H0, covI_H, covI_Y)

## Main difference between simulating many times for individuals and not is the variance of I and it's covariance with H - get a lot smaller variance when do this ##

simulated_moments
simulated_moments2

