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

options(scipen=10) # making so format doesn't include e's # 


########################################################
# DEFINE PARAMETER VALUES
########################################################
var_xi     <- 1 
theta      <- 0.5 
phi        <- 0.5 
alpha      <- 0.5 
beta       <- 0.5 
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


set.seed(24112016)
xi <- rlnorm(n, mean = 0, sd = var_xi)

########################################################
# WRITE FUNCTION TO SOLVE THE MODEL 
########################################################
foc.inv       <- function(I){
  func <-  1/(beta*(1 - theta) - 1) * (log(P[i]) - theta * log(Y[i] - P[i]*I) - log(beta) - alpha*(1 - phi) * log(H0[i]) - alpha*(1 - phi)*log(xi[i])) - log(I)
  func2 <- Vectorize(func) # Need to turn into a vector to get to work with uniroot.all # 
  return(func2)
}



########################################################
# SOLVE INVESTMENT FOR EVERYONE IN THE DATA 
########################################################

#creating vector to be filled#
optinv <- numeric(length = n)

# Solving for everyone in the data # 
  for (i in 1:n) {
    optinv[i]   <- uniroot.all(foc.inv, lower =- 1, upper = 49)
  }
optinv 

#### Calculating human capital for all in the data #### 
H <- numeric(length = n)
for (i in 1:n) { 
  H[i] <- xi[i] * H0[i]^alpha * optinv[i]^beta
} 




########################################################
# SIMULATED MOMENTS 
########################################################

mean_I <- mean(optinv)
mean_H <- mean(H)


var_I <- var(optinv)
var_H <- var(H)

covI_Y <- cov(optinv, Y)
covI_H0 <- cov(optinv, H0)
covI_H <- cov(optinv, H)

simulated_moments <- c(mean_I, mean_H, var_I, var_H, covI_Y, covI_H0, covI_H, covI_Y)
length(simulated_moments)

simulated_moments





