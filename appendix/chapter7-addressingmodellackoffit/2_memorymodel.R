# Allowing your Markov models to remember

library(nimble)
library(nimbleEcology)
library(MCMCvis)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# read in data
y <- as.matrix(read.table("all-geese-2sites.txt"))

# get the occasion of first capture for each individual
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

# filter out individuals that are first captured at last occasion. 
# These individuals do not contribute to parameter estimation, and also they cause problems with nimbleEcology. 
mask <- which(first!=ncol(y)) # individuals that are not first encountered at last occasion
y <- y[mask, ]                # keep only these
first <- first[mask]

# let's write the model. Note that the likelihood is simpler due to the use of the function `dHMM`.
multisite.marginalized <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phi111: survival-mov probability from state 11 to state 11
  # phi112: survival-mov probability from state 11 to state 12
  # phi121: survival-mov probability from state 12 to state 21
  # phi122: survival-mov probability from state 12 to state 22
  # phi211: survival-mov probability from state 21 to state 11
  # phi212: survival-mov probability from state 21 to state 12
  # phi221: survival-mov probability from state 22 to state 21
  # phi222: survival-mov probability from state 22 to state 22
  # det1: detection probability site 1
  # det2: detection probability site 2
  # pi11: init stat prob 11
  # pi12: init stat prob 12
  # pi21: init stat prob 21
  # pi22: init stat prob 22
  # -------------------------------------------------
  # States (S):
  # 1 alive 11
  # 2 alive 12
  # 3 alive 21
  # 4 alive 22
  # 5 dead
  # Observations (O):  
  # 1 not seen
  # 2 seen at site 1 
  # 3 seen at site 2
  # -------------------------------------------------
  
  # priors
  det1 ~ dunif(0, 1)
  det2 ~ dunif(0, 1)
  phi11[1:3] ~ ddirch(alpha[1:3]) # phi111, phi112, 1-sum
  phi12[1:3] ~ ddirch(alpha[1:3]) # phi121, phi122, 1-sum
  phi21[1:3] ~ ddirch(alpha[1:3]) # phi211, phi212, 1-sum
  phi22[1:3] ~ ddirch(alpha[1:3]) # phi221, phi222, 1-sum
  
  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phi11[1]
  gamma[1,2] <- phi11[2]
  gamma[1,3] <- 0
  gamma[1,4] <- 0
  gamma[1,5] <- phi11[3]
  gamma[2,1] <- 0
  gamma[2,2] <- 0
  gamma[2,3] <- phi12[1]
  gamma[2,4] <- phi12[2]
  gamma[2,5] <- phi12[3]
  gamma[3,1] <- phi21[1]
  gamma[3,2] <- phi21[2]
  gamma[3,3] <- 0
  gamma[3,4] <- 0
  gamma[3,5] <- phi21[3]
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- phi22[1]
  gamma[4,4] <- phi22[2]
  gamma[4,5] <- phi22[3]
  gamma[5,1] <- 0
  gamma[5,2] <- 0
  gamma[5,3] <- 0
  gamma[5,4] <- 0
  gamma[5,5] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- 1 - det1
  omega[1,2] <- det1
  omega[1,3] <- 0
  omega[2,1] <- 1 - det2
  omega[2,2] <- 0
  omega[2,3] <- det2
  omega[3,1] <- 1 - det1
  omega[3,2] <- det1
  omega[3,3] <- 0
  omega[4,1] <- 1 - det2
  omega[4,2] <- 0
  omega[4,3] <- det2
  omega[5,1] <- 1
  omega[5,2] <- 0
  omega[5,3] <- 0
  
  # initial state probs
  for(i in 1:N) {
    init[i, 1:5] <- gamma[ y[i, first[i] ] - 1, 1:5 ] # First state propagation
  }
  
  # likelihood 
  for (i in 1:N){
    y[i,(first[i]+1):K] ~ dHMM(init = init[i,1:5],  # count data from first[i] + 1
                               probObs = omega[1:5,1:3],     # observation matrix
                               probTrans = gamma[1:5,1:5],   # transition matrix
                               len = K - first[i],           # nb of occasions
                               checkRowSums = 0)             # do not check whether elements in a row sum tp 1
  }
})

# data in a list
my.data <- list(y = y + 1, 
                alpha = c(1, 1, 1))

# constants in a list
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

# initial values - note that we do not need initial values for the latent states anymore. 
initial.values <- function(){list(phi11 = rdirch(1, rep(1,3)), 
                                  phi12 = rdirch(1, rep(1,3)), 
                                  phi21 = rdirch(1, rep(1,3)), 
                                  phi22 = rdirch(1, rep(1,3)),
                                  det1 = runif(1, 0, 1), 
                                  det2 = runif(1, 0, 1))}  

# parameters to monitor.
parameters.to.save <- c("phi11", "phi12", "phi21", "phi22", "det1", "det2")

# MCMC settings.
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# run NIMBLE
out <- nimbleMCMC(code = multisite.marginalized, 
                  constants = my.constants,
                  data = my.data,              
                  inits = initial.values(),
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin, 
                  nchains = n.chains)

MCMCsummary(out, round = 2)

MCMCplot(out)

MCMCtrace(out, params = c("det1","det2"), pdf = FALSE)

