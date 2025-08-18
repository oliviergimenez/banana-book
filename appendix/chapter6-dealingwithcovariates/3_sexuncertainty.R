# Uncertainty in sex

# STATES
# alive as male (M)
# alive as female (F)
# dead (D)

# OBS
# not seen (1)
# judged from copulation to be M (2)
# judged from begging food to be M (3)
# judged from coutship feeding to be M (4)
# judged from body size to be M (5)
# judged from copulation to be F (6)
# judged from begging food to be F (7)
# judged from coutship feeding to be F (8)
# judged from body size to be F (9)
# not judged (10)

library(nimble)
library(MCMCvis)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# read in data
y <- as.matrix(read.table("audouin-gull.txt"))

hmm.sex <- nimbleCode({
  
  # priors
  phiF ~ dunif(0, 1) # prior survival male
  phiM ~ dunif(0, 1) # prior survival female
  p ~ dunif(0, 1) # prior encounter prob
  pi ~ dunif(0, 1) # prob init state male
  e ~ dunif(0,1)
  for (j in 1:4){
    x[j] ~ dunif(0,1)
  }
  m[1] ~ dunif(0,1)
  m[2] <- (1 - m[1]) * alpha
  alpha ~ dunif(0,1)
  m[3] <- 1 - m[1] - m[2]
  m[4] ~ dunif(0,1)
  
  # HMM ingredients
  delta[1] <- pi
  delta[2] <- 1 - pi
  delta[3] <- 0
  
  gamma[1,1] <- phiM 
  gamma[1,2] <- 0 
  gamma[1,3] <- 1 - phiM
  gamma[2,1] <- 0    
  gamma[2,2] <- phiF        
  gamma[2,3] <- 1 - phiF   
  gamma[3,1] <- 0      
  gamma[3,2] <- 0    
  gamma[3,3] <- 1       
  
  omega[1,1] <- 1 - p                                  
  omega[1,2] <- p * e * (1 - m[4]) * m[1] * x[1]      
  omega[1,3] <- p * e * (1 - m[4]) * m[2] * x[2]      
  omega[1,4] <- p * e * (1 - m[4]) * m[3] * x[3]      
  omega[1,5] <- p * e * m[4] * x[4]      
  omega[1,6] <- p * e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega[1,7] <- p * e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega[1,8] <- p * e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega[1,9] <- p * e * m[4] * (1 - x[4])      
  omega[1,10] <- p * (1 - e)      
  
  omega[2,1] <- 1 - p                                  
  omega[2,2] <- p * e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega[2,3] <- p * e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega[2,4] <- p * e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega[2,5] <- p * e * m[4] * (1 - x[4])      
  omega[2,6] <- p * e * (1 - m[4]) * m[1] * x[1]      
  omega[2,7] <- p * e * (1 - m[4]) * m[2] * x[2]     
  omega[2,8] <- p * e * (1 - m[4]) * m[3] * x[3]   
  omega[2,9] <- p * e * m[4] * x[4]   
  omega[2,10] <- p * (1 - e)      
  
  omega[3,1] <- 1                                  
  omega[3,2] <- 0     
  omega[3,3] <- 0     
  omega[3,4] <- 0     
  omega[3,5] <- 0     
  omega[3,6] <- 0     
  omega[3,7] <- 0     
  omega[3,8] <- 0     
  omega[3,9] <- 0     
  omega[3,10] <- 0     
  
  omega.init[1,1] <- 0                                 
  omega.init[1,2] <- e * (1 - m[4]) * m[1] * x[1]      
  omega.init[1,3] <- e * (1 - m[4]) * m[2] * x[2]      
  omega.init[1,4] <- e * (1 - m[4]) * m[3] * x[3]      
  omega.init[1,5] <- e * m[4] * x[4]      
  omega.init[1,6] <- e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega.init[1,7] <- e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega.init[1,8] <- e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega.init[1,9] <- e * m[4] * (1 - x[4])      
  omega.init[1,10] <- (1 - e)      
  omega.init[2,1] <- 0                                  
  omega.init[2,2] <- e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega.init[2,3] <- e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega.init[2,4] <- e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega.init[2,5] <- e * m[4] * (1 - x[4])      
  omega.init[2,6] <- e * (1 - m[4]) * m[1] * x[1]      
  omega.init[2,7] <- e * (1 - m[4]) * m[2] * x[2]     
  omega.init[2,8] <- e * (1 - m[4]) * m[3] * x[3]   
  omega.init[2,9] <- e * m[4] * x[4]   
  omega.init[2,10] <- (1 - e)      
  omega.init[3,1] <- 1                                  
  omega.init[3,2] <- 0     
  omega.init[3,3] <- 0     
  omega.init[3,4] <- 0     
  omega.init[3,5] <- 0     
  omega.init[3,6] <- 0     
  omega.init[3,7] <- 0     
  omega.init[3,8] <- 0     
  omega.init[3,9] <- 0     
  omega.init[3,10] <- 0     
  
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:10])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3])
      y[i,j] ~ dcat(omega[z[i,j], 1:10])
    }
  }
})

first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), K = ncol(y), first = first) 
my.data <- list(y = y + 1)


# function to create initial latent states z
create_initial_z <- function(y, first) {
  N <- nrow(y)
  K <- ncol(y)
  zinit <- matrix(NA, nrow = N, ncol = K)
  
  for (i in 1:N) {
    # determine sex at first detection
    first_obs <- y[i, first[i]]  
    if (first_obs %in% 1:5) {
      initial_sex <- 1  # male
    } else if (first_obs %in% 6:8) {
      initial_sex <- 2  # female
    } else {
      initial_sex <- sample(1:2, 1)  # random sex if uncertain
    }
    
    zinit[i, first[i]] <- initial_sex
    
    # for subsequent occasion
    current_state <- initial_sex
    if (first[i] < K) {
      for (j in (first[i] + 1):K) {
        obs <- y[i, j]
        
        if (obs == 0) {
          zinit[i, j] <- current_state
        } else {
          if (current_state == 3) {
            current_state <- sample(1:2, 1)
          }
          zinit[i, j] <- current_state
        }
      }
    }
  }
  
  return(zinit)
}

# function to generate valid m parameters
generate_valid_m <- function() {
  m <- rep(NA, 4)
  
  # m[4]
  m[4] <- runif(1, 0.1, 0.4)
  
  # m[1] and m[2] w/ constraint m[1] + m[2] + m[3] = 1 and m[3] > 0
  repeat {
    m[1] <- runif(1, 0.1, 0.4)
    alpha <- runif(1, 0.1, 0.4)
    m[2] <- (1 - m[1]) * alpha
    if (m[1] + m[2] <= 0.8) break  # m[3] >= 0.2
  }
  m[3] <- 1 - m[1] - m[2]
  
  return(m)
}

# function to create valid initial parameter values
create_initial_params <- function() {
  phiM <- runif(1, 0.75, 0.95)  
  phiF <- runif(1, 0.75, 0.95)  
  p <- runif(1, 0.4, 0.8)       
  pi <- runif(1, 0.45, 0.55)    
  e <- runif(1, 0.3, 0.7)       
  x <- 1 - runif(4, 0.01, 0.3)   
  m <- generate_valid_m()
  
  return(list(
    phiM = phiM,
    phiF = phiF, 
    p = p,
    pi = pi,
    e = e,
    x = x,
    m = m
  ))
}

# Create the complete initial values function
initial.values <- function() {
  params <- create_initial_params()
  zinit <- create_initial_z(y, first)
  return(c(params, list(z = zinit)))
}

parameters.to.save <- c("phiM", 
                        "phiF", 
                        "p", 
                        "pi",
                        "e", 
                        "m", 
                        "x")

n.iter <- 2500
n.burnin <- 500
n.chains <- 2

out <- nimbleMCMC(code = hmm.sex, 
                  constants = my.constants,
                  data = my.data,              
                  inits = initial.values(),
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin, 
                  nchains = n.chains)

MCMCsummary(out, round = 2)

MCMCtrace(out, pdf = FALSE, params = c("p"))
MCMCtrace(out, pdf = FALSE, params = c("e"))
MCMCtrace(out, pdf = FALSE, params = c("m"))
MCMCtrace(out, pdf = FALSE, params = c("x"))
MCMCtrace(out, pdf = FALSE, params = c("phiM", "phiF"))
MCMCtrace(out, pdf = FALSE, params = c("pi"))


