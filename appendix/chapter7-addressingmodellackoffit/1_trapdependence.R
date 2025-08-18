# Accounting for trap-dependence

library(nimble)
library(MCMCvis)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# read in data
y <- as.matrix(read.table("calonectris_diomedea.txt"))

# gof
library(R2ucare)
size <- rep(1, nrow(y))
overall_CJS(y, size)
# test for transience
test3sr(y, size)
# test for trap-dependence:
test2ct(y, size)

# OBS
# 1 not seen
# 2 seen

# STATES
# 1 trap aware
# 2 trap unaware
# 3 dead

hmm.td <- nimbleCode({
  
  # priors
  beta[1] ~ dunif(0, 1) # phi1
  beta[2] ~ dunif(0, 1) # phi2
  phi1 <- beta[1]
  phi2 <- beta[2]
  p ~ dunif(0, 1) # prior detection
  pprime ~ dunif(0, 1) # prior detection
  #pi ~ dunif(0, 1) # prob init state 1
  
  # HMM ingredients
  for (i in 1:N){
    for (t in first[i]:(K-1)){
      phi[i,t] <- beta[age[i,t]] # beta1 = phi1, beta2 = phi1+
      gamma[1,1,i,t] <- phi[i,t] * pprime 
      gamma[1,2,i,t] <- phi[i,t] * (1 - pprime)    
      gamma[1,3,i,t] <- 1 - phi[i,t]              
      gamma[2,1,i,t] <- phi[i,t] * p           
      gamma[2,2,i,t] <- phi[i,t] * (1 - p)     
      gamma[2,3,i,t] <- 1 - phi[i,t]            
      gamma[3,1,i,t] <- 0          
      gamma[3,2,i,t] <- 0           
      gamma[3,3,i,t] <- 1                
    }
  }
  
  delta[1] <- 1                  
  delta[2] <- 0               
  delta[3] <- 0                
  
  omega[1,1] <- 0                 
  omega[1,2] <- 1                
  omega[2,1] <- 1                    
  omega[2,2] <- 0                    
  omega[3,1] <- 1                  
  omega[3,2] <- 0                    
  omega.init[1,1] <- 0               
  omega.init[1,2] <- 1             
  omega.init[2,1] <- 1              
  omega.init[2,2] <- 0            
  omega.init[3,1] <- 1             
  omega.init[3,2] <- 0            
  
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:2])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3, i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# age effect via an individual by time covariate and use nested indexing 
# to distinguish survival over the interval after first detection from survival afterwards: 
first <- apply(y, 1, function(x) min(which(x != 0)))
age <- matrix(NA, nrow = nrow(y), ncol = ncol(y) - 1)
for (i in 1:nrow(age)){
  for (j in 1:ncol(age)){
    if (j == first[i]) age[i,j] <- 1 # age = 1
    if (j > first[i]) age[i,j] <- 2  # age > 1
  }
}

my.constants <- list(N = nrow(y), K = ncol(y), first = first, age = age)
my.data <- list(y = y + 1)

# Viterbi-like initialization
create_initial_z <- function(y, first, phi = 0.9, p = 0.5, pprime = 0.5, pi = 1) {
  N <- nrow(y)
  K <- ncol(y)
  zinit <- matrix(3, N, K)  # initialize all as dead
  
  # transition matrix
  gamma <- matrix(0, nrow = 3, ncol = 3)
  gamma[1, ] <- c(phi * pprime, phi * (1 - pprime), 1 - phi)
  gamma[2, ] <- c(phi * p,      phi * (1 - p),      1 - phi)
  gamma[3, ] <- c(0,            0,                 1)
  
  # observation matrix: y is 1 = not seen, 2 = seen
  omega <- matrix(0, nrow = 3, ncol = 2)
  omega[1, ] <- c(0, 1)  # trap-aware: only seen
  omega[2, ] <- c(1, 0)  # trap-unaware: only not seen
  omega[3, ] <- c(1, 0)  # dead: only not seen
  
  for (i in 1:N) {
    if (is.na(first[i]) || first[i] > K) next
    
    T_len <- K - first[i] + 1
    obs_seq <- y[i, first[i]:K] + 1
    
    # forward algorithm (unnormalized alpha)
    alpha <- matrix(0, nrow = 3, ncol = T_len)
    
    # initial distribution
    init_probs <- c(pi, 1 - pi, 0)
    alpha[, 1] <- init_probs * omega[, obs_seq[1]]
    
    if (sum(alpha[, 1]) > 0) {
      alpha[, 1] <- alpha[, 1] / sum(alpha[, 1])
    } else {
      alpha[, 1] <- c(1, 0, 0)  # fallback: only trap-aware
    }
    
    if (T_len > 1) {
      for (t in 2:T_len) {
        for (s in 1:3) {
          trans_prob <- sum(alpha[, t - 1] * gamma[, s])
          alpha[s, t] <- trans_prob * omega[s, obs_seq[t]]
        }
        if (sum(alpha[, t]) > 0) {
          alpha[, t] <- alpha[, t] / sum(alpha[, t])
        } else {
          alpha[, t] <- c(1, 0, 0)  # fallback to avoid NA
        }
      }
    }
    
    # backtrace the most likely path
    path <- rep(NA, T_len)
    path[T_len] <- which.max(alpha[, T_len])
    
    if (T_len > 1) {
      for (t in (T_len - 1):1) {
        probs <- alpha[, t] * gamma[, path[t + 1]]
        if (sum(probs) == 0) {
          path[t] <- 1  # fallback to trap-aware
        } else {
          path[t] <- which.max(probs)
        }
      }
    }
    
    # assign path to zinit
    zinit[i, first[i]:K] <- path
  }
  
  return(zinit)
}

initial.values <- function() {
  beta <- runif(2, 0, 1)
  p <- runif(1, 0, 1)
  pprime <- runif(1, 0, 1)
  pi <- runif(1, 0, 1)
  zinit <- create_initial_z(y, first, phi = beta[1], p, pprime, pi = 1)
  list(beta = beta, p = p, pprime = pprime, z = zinit)
}

parameters.to.save <- c("phi1", "phi2", "p", "pprime")

n.iter <- 15000
n.burnin <- 2500
n.chains <- 2

out <- nimbleMCMC(code = hmm.td,
                  constants = my.constants,
                  data = my.data,
                  inits = initial.values,
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin,
                  nchains = n.chains)

MCMCsummary(out, round = 2)

MCMCtrace(out, pdf = FALSE)

#------------------------------------------------------------------------------------

# NO trap-dependence

# read in data
y <- as.matrix(read.table("calonectris_diomedea.txt"))

# OBS
# 1 not seen
# 2 seen

# STATES
# 1 alive
# 2 dead

hmm <- nimbleCode({
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    for (t in first[i]:(K-1)){
      phi[i,t] <- beta[age[i,t]] # beta1 = phi1, beta2 = phi1+
      gamma[1,1,i,t] <- phi[i,t]      # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t]  # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0           # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1           # Pr(dead t -> dead t+1)
    }
  }
  beta[1] ~ dunif(0, 1) # phi1
  beta[2] ~ dunif(0, 1) # phi1+
  phi1 <- beta[1]
  phi2 <- beta[2]
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# age effect via an individual by time covariate and use nested indexing 
# to distinguish survival over the interval after first detection from survival afterwards: 
age <- matrix(NA, nrow = nrow(y), ncol = ncol(y) - 1)
for (i in 1:nrow(age)){
  for (j in 1:ncol(age)){
    if (j == first[i]) age[i,j] <- 1 # age = 1
    if (j > first[i]) age[i,j] <- 2  # age > 1
  }
}

first <- apply(y, 1, function(x) min(which(x != 0)))
my.constants <- list(N = nrow(y), K = ncol(y), first = first, age = age)
my.data <- list(y = y + 1)

zinits <- y
zinits[zinits == 0] <- 1
initial.values <- function() list(beta = runif(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

parameters.to.save <- c("phi1", "phi2", "p")

n.iter <- 2500
n.burnin <- 500
n.chains <- 2

out <- nimbleMCMC(code = hmm,
                  constants = my.constants,
                  data = my.data,
                  inits = initial.values,
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin,
                  nchains = n.chains)

MCMCsummary(out, round = 2)

MCMCtrace(out, pdf = FALSE)



