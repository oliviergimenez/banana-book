# Accomodating individual heterogeneity (w/ finite mixture)

library(nimble)
library(MCMCvis)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# simulate some data
set.seed(1234) # for reproducibility

# generate heterogeneous detection probabilities
phi <- 0.7 # survival
prop_class1 <- 0.2 # pi
p_class1 <- 0.8 # low prop are highly detectable
p_class2 <- 0.3
nind <- 400 # nb of ind
nyear <- 10 # duration of the study
expit <- function(x){exp(x)/(1+exp(x))} # reciprocal of the logit function 
z <- data <- x <- matrix(NA, nrow = nind, ncol = nyear)
first <- rep(1,nind)
detection <- rep(NA,nind)
which_mixture <- rep(NA,nind)
# first assign ind to a class, then use corresponding detection
for (i in 1:nind){
  which_mixture[i] <- rbinom(1,1,prop_class1) # assign ind i to a class with prob pi
  if (which_mixture[i] == 1){ detection[i] <- p_class1  
  } else {
    detection[i] <- p_class2} 
}

# generate encounter histories
for(i in 1:nind){
  z[i,first[i]] <- x[i,first[i]] <- 1 
  for(j in (first[i]+1):nyear){
    z[i,j] <- rbinom(1,1,phi*z[i,j-1])
    x[i,j] <- rbinom(1,1,z[i,j]*detection[i]) }
}
y <- x 
y[is.na(y)] <- 0

# write down model
# two alive states for two classes of individuals
# detection probabilities dependent on states and constant survival
hmm.phipmix <- nimbleCode({
  
  # priors
  phi ~ dunif(0, 1) # prior survival
  one ~ dconstraint(pp1 < pp2) # to avoid label switching
  pp1 ~ dunif(0, 1) # prior detection
  pp2 ~ dunif(0, 1) # prior detection
  pi ~ dunif(0, 1) # prob init state 1
  # transition matrix
  gamma[1,1] <- phi      # A1(t)->A1(t+1)
  gamma[1,2] <- 0        # A1(t)->A2(t+1)
  gamma[1,3] <- 1 - phi  # A1(t)->D(t+1)
  gamma[2,1] <- 0        # A2(t)->A1(t+1)
  gamma[2,2] <- phi      # A2(t)->A2(t+1)
  gamma[2,3] <- 1 - phi  # A2(t)->D(t+1)
  gamma[3,1] <- 0        # D(t)->A1(t+1)
  gamma[3,2] <- 0        # D(t)->A2(t+1)
  gamma[3,3] <- 1        # D(t)->D(t+1)
  
  # vector of initial state probs
  delta[1] <- pi         # A1(first)
  delta[2] <- 1 - pi     # A2(first)
  delta[3] <- 0          # D(first)
  
  # observation matrix
  omega[1,1] <- 1 - pp1   # A1(t)->0(t)
  omega[1,2] <- pp1       # A1(t)->1(t)
  omega[2,1] <- 1 - pp2   # A2(t)->0(t)
  omega[2,2] <- pp2       # A2(t)->1(t)
  omega[3,1] <- 1         # D(t)->0(t)
  omega[3,2] <- 0         # D(t)->1(t)
  
  # first encounter t = first
  omegae[1,1] <- 0        # A1(t)->0(t)
  omegae[1,2] <- 1        # A1(t)->1(t)
  omegae[2,1] <- 0        # A2(t)->0(t)
  omegae[2,2] <- 1        # A2(t)->1(t)
  omegae[3,1] <- 1        # D(t)->0(t)
  omegae[3,2] <- 0        # D(t)->1(t)
  
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omegae[z[i,first[i]], 1:2])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# get date of first capture
first <- apply(y, 1, function(x) min(which(x != 0)))

# put constants in a list
my.constants <- list(N = nrow(y), 
                     K = ncol(y), 
                     first = first)

# put data in a list
my.data <- list(y = y + 1, one = 1)

# initial values 
zinit <- y
for (i in 1:nrow(y)) {
  for (j in first[i]:ncol(y)) {
    if (j == first[i]) zinit[i,j] <- which(rmultinom(1, 1, c(2/3,1/3))==1) # pick alive state
    if (j > first[i]) zinit[i,j] <- zinit[i,j-1]
  }     
}
zinit <- as.matrix(zinit)
initial.values <- function() list(phi = runif(1,0,1),
                                  pp1 = runif(1,0,1),
                                  pp2 = runif(1,0,1),
                                  pi = runif(1,0,1),
                                  z = as.matrix(zinit))

# parameters to be monitored
parameters.to.save <- c("phi", "pp1", "pp2", "pi")

# MCMC details
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# model fitting
mcmc.phipmix <- nimbleMCMC(code = hmm.phipmix, 
                           constants = my.constants,
                           data = my.data,              
                           inits = initial.values,
                           monitors = parameters.to.save,
                           niter = n.iter,
                           nburnin = n.burnin, 
                           nchains = n.chains)

# numerical summaries
MCMCsummary(mcmc.phipmix, round = 2)

# mean   sd 2.5%  50% 97.5% Rhat n.eff
# phi 0.70 0.01 0.67 0.70  0.73 1.00  1074
# pi  0.70 0.03 0.63 0.70  0.76 1.01   722
# pp1 0.45 0.02 0.40 0.45  0.49 1.00   492
# pp2 0.49 0.03 0.43 0.48  0.55 1.01   632

# trace and posterior distribution
MCMCtrace(mcmc.phipmix, 
          pdf = FALSE, 
          params = c("pi", "pp1", "pp2", "phi"))

# everything looks fine, but posterior means are far from values we used
# to simulate data

# let's try marginalized likelihood
# use forward algo for HMM, and pooling of encounter histories. 

# get rid of individuals for which first==K
mask <- which(first!=ncol(y)) # individuals that are not first encountered at last occasion
y <- y[mask, ]                # keep only these
first <- first[mask]

# pool encounter histories
y_weighted <- y %>% 
  as_tibble() %>% 
  group_by_all() %>% 
  summarise(size = n()) %>% 
  relocate(size) %>% 
  as.matrix()
head(y_weighted)
size <- y_weighted[,1] # nb of individuals w/ a particular encounter history
y <- y_weighted[,-1] # pooled data

# assemble data and constants
my.data <- list(y = y + 1)
my.constants <- list(N = nrow(y), 
                     K = ncol(y), 
                     first = first,
                     size = size,
                     one = 1)

# NIMBLE functions
dwolfHMM <- nimbleFunction(
  run = function(x = double(1), 
                 probInit = double(1), # vector of initial states
                 probObs = double(2), #observation matrix
                 probObse = double(2),
                 probTrans = double(2), # transition matrix
                 size = double(0),
                 len = double(0, default = 0), # number of sampling occasions
                 log = integer(0, default = 0)) {
    alpha <- probInit[1:3] * probObse[1:3,x[1]]# * probObs[1:3,x[1]] == 1 due to conditioning on first detection
    for (t in 2:len) {
      alpha[1:3] <- (alpha[1:3] %*% probTrans[1:3,1:3]) * probObs[1:3,x[t]]
    }
    logL <- size * log(sum(alpha[1:3]))
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

rwolfHMM <- nimbleFunction(
  run = function(n = integer(),
                 probInit = double(1),
                 probObs = double(2),
                 probObse = double(2),
                 probTrans = double(2),
                 size = double(0),
                 len = double(0, default = 0)) {
    returnType(double(1))
    z <- numeric(len)
    z[1] <- rcat(n = 1, prob = probInit[1:3]) # all individuals alive at t = 0
    y <- z
    y[1] <- 2 # all individuals are detected at t = 0
    for (t in 2:len){
      # state at t given state at t-1
      z[t] <- rcat(n = 1, prob = probTrans[z[t-1],1:3]) 
      # observation at t given state at t
      y[t] <- rcat(n = 1, prob = probObs[z[t],1:2]) 
    }
    return(y)
  })

assign('dwolfHMM', dwolfHMM, .GlobalEnv)
assign('rwolfHMM', rwolfHMM, .GlobalEnv)

# NIMBLE code again
hmm.phipmix <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # pi: initial state probability C1
  # phi: survival probability
  # pp1: recapture probability C1
  # pp2: recapture probability C2
  # -------------------------------------------------
  # States (S):
  # 1 alive (A1)
  # 2 alive (A2)
  # 5 dead (D)
  # Observations (O):
  # 1 neither seen nor recovered (0)
  # 2 seen alive (1)
  # -------------------------------------------------
  
  # priors
  phi ~ dunif(0, 1) # prior survival
  one ~ dconstraint(pp1 < pp2) # to avoid label switching
  pp1 ~ dunif(0, 1) # prior detection
  pp2 ~ dunif(0, 1) # prior detection
  pi ~ dunif(0, 1) # prob init state 1
  # transition matrix
  gamma[1,1] <- phi      # A1(t)->A1(t+1)
  gamma[1,2] <- 0        # A1(t)->A2(t+1)
  gamma[1,3] <- 1 - phi  # A1(t)->D(t+1)
  gamma[2,1] <- 0        # A2(t)->A1(t+1)
  gamma[2,2] <- phi      # A2(t)->A2(t+1)
  gamma[2,3] <- 1 - phi  # A2(t)->D(t+1)
  gamma[3,1] <- 0        # D(t)->A1(t+1)
  gamma[3,2] <- 0        # D(t)->A2(t+1)
  gamma[3,3] <- 1        # D(t)->D(t+1)
  
  # vector of initial state probs
  delta[1] <- pi         # A1(first)
  delta[2] <- 1 - pi     # A2(first)
  delta[3] <- 0          # D(first)
  
  # observation matrix
  omega[1,1] <- 1 - pp1   # A1(t)->0(t)
  omega[1,2] <- pp1       # A1(t)->1(t)
  omega[2,1] <- 1 - pp2   # A2(t)->0(t)
  omega[2,2] <- pp2       # A2(t)->1(t)
  omega[3,1] <- 1        # D(t)->0(t)
  omega[3,2] <- 0        # D(t)->1(t)
  
  omegae[1,1] <- 0   # A1(t)->0(t)
  omegae[1,2] <- 1   # A1(t)->1(t)
  omegae[2,1] <- 0   # A2(t)->0(t)
  omegae[2,2] <- 1   # A2(t)->1(t)
  omegae[3,1] <- 1   # D(t)->0(t)
  omegae[3,2] <- 0   # D(t)->1(t)
  
  # likelihood
  for(i in 1:N) {
    y[i,first[i]:K] ~ dwolfHMM(probInit = delta[1:3], # count data from first[i] + 1
                               probObs = omega[1:3,1:2], # observation matrix
                               probObse = omegae[1:3,1:2], # observation matrix
                               probTrans = gamma[1:3,1:3], # transition matrix
                               size = size[i],
                               len = K - first[i] + 1) # nb of occasions
  }
})

# initial values 
# cool thing is that we do not need inits for the latent states anymore
initial.values <- function() list(phi = runif(1,0,1),
                                  pp1 = 0.3,
                                  pp2 = 0.8,
                                  pi = runif(1,0,1))

# parameters to be monitored
parameters.to.save <- c("phi", "pp1", "pp2", "pi")

# MCMC details
n.iter <- 5000*2
n.burnin <- 1000*2
n.chains <- 2

# run NIMBLE
mcmc.phipmix.marginalized <- nimbleMCMC(code = hmm.phipmix, 
                                        constants = my.constants,
                                        data = my.data,              
                                        inits = initial.values,
                                        monitors = parameters.to.save,
                                        niter = n.iter,
                                        nburnin = n.burnin, 
                                        nchains = n.chains)

# numerical summaries
MCMCsummary(mcmc.phipmix.marginalized, round = 2)

# mean   sd 2.5%  50% 97.5% Rhat n.eff
# phi 0.72 0.02 0.69 0.72  0.75    1  1120
# pi  0.72 0.08 0.54 0.73  0.87    1   298
# pp1 0.27 0.05 0.18 0.27  0.35    1   343
# pp2 0.79 0.07 0.65 0.79  0.92    1   336

# we're back in business! Remember, true survival is 0.7 
# and we have prop 0.2 with detection 0.8 and 0.8 with detection 0.3
