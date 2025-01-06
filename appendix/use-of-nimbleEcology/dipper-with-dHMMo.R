################################### MODEL PHI,P ################################

#---------------------------------- with dCJSss --------------------------------

# load nimbleEcology and a few other packages
library(nimbleEcology)
library(tidyverse)
library(MCMCvis)

# data
y <- read_csv("dipper.csv") %>%
  select(year_1981:year_1987) %>%
  as.matrix()

# get occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))

# get rid of individuals for which first detection = last capture
y <- y[ first != ncol(y), ] 

# re-calculate occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))

# NIMBLE code
hmm.phip <- nimbleCode({
  phi ~ dunif(0, 1) # survival prior
  p ~ dunif(0, 1)   # detection prior
  # likelihood
  for (i in 1:N){
    y[i, first[i]:T] ~ dCJS_ss(probSurvive = phi, 
                               probCapture = p, 
                               len = T - first[i] + 1)
  }
})
# constants
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first)
# data
my.data <- list(y = y) # format is 0 for non-detected and 1 for detected
# initial values; we use the marginalized likelihood, so no latent states 
# in it, therefore no need for initial values for the latent states
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1))
initial.values <- function() list(phi = 0.6,
                                  p = 0.9)
# parameters to monitor
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

survival <- nimbleModel(code = hmm.phip,
                        data = my.data,
                        constants = my.constants,
                        inits = initial.values())
survival$calculate()
# [1] -334.7412

# run NIMBLE
mcmc.phip <- nimbleMCMC(code = hmm.phip, 
                        constants = my.constants,
                        data = my.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = n.iter,
                        nburnin = n.burnin, 
                        nchains = n.chains)
# numerical summaries
MCMCsummary(mcmc.phip, round = 2)

# mean   sd 2.5%  50% 97.5% Rhat n.eff
# p   0.90 0.03 0.83 0.90  0.95 1.01   577
# phi 0.56 0.03 0.51 0.56  0.61 1.00   619


#---------------------------------- with dHMMo --------------------------------

# load nimbleEcology and a few other packages
library(nimbleEcology)
library(tidyverse)
library(MCMCvis)

# data
y <- read_csv("dipper.csv") %>%
  select(year_1981:year_1987) %>%
  as.matrix()

# get occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))

# get rid of individuals for which first detection = last capture
y <- y[ first != ncol(y), ] 

# re-calculate occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))

# NIMBLE code (fixed by Daniel)
hmm.survival <- nimbleCode({
  phi ~ dunif(0, 1) # prior survival
  p ~ dunif(0, 1) # prior detection
  # likelihood
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
  omega[1,1,1] <- 0        # Pr(alive first -> non-detected first)
  omega[1,2,1] <- 1        # Pr(alive first -> detected first)
  omega[2,1,1] <- 1        # Pr(dead first -> non-detected first)
  omega[2,2,1] <- 0        # Pr(dead first -> detected first)
  for (t in 2:K){
    omega[1,1,t] <- 1 - p    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1        # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0        # Pr(dead t -> detected t)
  }
  for (i in 1:N){
    y[i,first[i]:K] ~ dHMMo(init = delta[1:2], # vector of initial state probabilities
                            probObs = omega[1:2,1:2,1:(K - first[i] + 1)], # observation matrix
                            probTrans = gamma[1:2,1:2], # transition matrix
                            len = K - first[i] + 1, # nb of sampling occasions
                            checkRowSums = 0) # skip validity checks
  }
})


# constants
my.constants <- list(N = nrow(y), 
                     K = ncol(y), 
                     first = first)

# data
my.data <- list(y = y + 1) # format is 1 for non-detected and 2 for detected
# initial values; we use the marginalized likelihood, so no latent states 
# in it, therefore no need for initial values for the latent states
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1))
initial.values <- function() list(phi = 0.6,
                                  p = 0.9)

# parameters to monitor
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

# create model as an R object (uncompiled model)
survival <- nimbleModel(code = hmm.survival,
                        data = my.data,
                        constants = my.constants,
                        inits = initial.values())

# Warning message:
#   In nimble::getParam(.self, node, param) :
#   Note that getParam is not available for parameters of dimension greater than two, 
#   but other calculations with models are unaffected

survival$getNodeNames()
survival$getParam(node = "phi", param = "max")
survival$getParam(node = "y[248, 6:7]", param = "probTrans")
survival$getParam(node = "y[248, 6:7]", param = "probObs") # warning is here

survival$calculate() # value different from what we get with dCJSss
# [1] -334.7412


# run NIMBLE
mcmc.phip <- nimbleMCMC(code = hmm.survival, 
                        constants = my.constants,
                        data = my.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = n.iter,
                        nburnin = n.burnin, 
                        nchains = n.chains)

# Warning message:
#   In nimble::getParam(.self, node, param) :
#   Note that getParam is not available for parameters of dimension greater than two,
#   but other calculations with models are unaffected

# numerical summaries
MCMCsummary(mcmc.phip, round = 2)

# mean   sd 2.5%  50% 97.5% Rhat n.eff
# p   0.90 0.03 0.83 0.90  0.95 1.01   620
# phi 0.56 0.03 0.51 0.56  0.61 1.00   634



################################### MODEL PHI,P(t) ################################

#---------------------------------- with NIMBLE and full likelihood -----------------

# load packages
library(nimble) # NIMBLE 
library(tidyverse) # manipulate/visualize data
library(MCMCvis) # post-process MCMC outputs
# read in data
dipper <- read_csv("dipper.csv")
head(dipper)
# format the data
y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()

# get rid of individuals for which first detection = last capture
y <- y[ first != ncol(y), ] 

# re-calculate occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))


# NIMBLE code 
hmm.phipt <- nimbleCode({
  phi ~ dunif(0, 1) # prior survival
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    p[t] ~ dunif(0, 1) # prior detection
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1        # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0        # Pr(dead t -> detected t)
  }
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
    }
  }
})

# constants
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first)
# data
my.data <- list(y = y + 1)
# initial values
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)
# parameters to monitor
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2
# run NIMBLE
mcmc.phipt <- nimbleMCMC(code = hmm.phipt, 
                         constants = my.constants,
                         data = my.data,              
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin, 
                         nchains = n.chains)
# numerical summaries
MCMCsummary(object = mcmc.phipt, round = 2)

# 
# mean   sd 2.5%  50% 97.5% Rhat n.eff
# p[1] 0.74 0.12 0.48 0.75  0.94 1.00   493
# p[2] 0.85 0.08 0.66 0.86  0.97 1.01   331
# p[3] 0.84 0.07 0.69 0.85  0.96 1.02   368
# p[4] 0.89 0.05 0.78 0.90  0.97 1.00   441
# p[5] 0.91 0.05 0.81 0.92  0.98 1.00   396
# p[6] 0.90 0.07 0.73 0.91  0.99 1.01   144
# phi  0.56 0.03 0.51 0.56  0.61 1.02   614


#---------------------------------- with nimbleEcology ---------------------------

# load nimbleEcology and a few other packages
library(nimbleEcology)
library(tidyverse)
library(MCMCvis)

# data
y <- read_csv("dipper.csv") %>%
  select(year_1981:year_1987) %>%
  as.matrix()

# get occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))

# get rid of individuals for which first detection = last capture
y <- y[ first != ncol(y), ] 

# re-calculate occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))

# NIMBLE code
hmm.survival <- nimbleCode({
  phi ~ dunif(0, 1) # prior survival
  p[1] <- 1
  for (j in 2:K){
    p[j] ~ dunif(0,1)
  }
  # likelihood
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
  omega[1,1,1] <- 1 - p[1]        # Pr(alive first -> non-detected first)
  omega[1,2,1] <- p[1]        # Pr(alive first -> detected first)
  omega[2,1,1] <- 1        # Pr(dead first -> non-detected first)
  omega[2,2,1] <- 0        # Pr(dead first -> detected first)
  for (t in 2:K){
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1        # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0        # Pr(dead t -> detected t)
  }
  for (i in 1:N){
    y[i,first[i]:K] ~ dHMMo(init = delta[1:2], # vector of initial state probabilities
                            probObs = omega[1:2,1:2,1:(K - first[i] + 1)], # observation matrix
                            probTrans = gamma[1:2,1:2], # transition matrix
                            len = K - first[i] + 1, # nb of sampling occasions
                            checkRowSums = 0) # skip validity checks
  }
})


# constants
my.constants <- list(N = nrow(y), 
                     K = ncol(y), 
                     first = first)

# data
my.data <- list(y = y + 1) # format is 1 for non-detected and 2 for detected
# initial values; we use the marginalized likelihood, so no latent states 
# in it, therefore no need for initial values for the latent states
initial.values <- function() list(phi = runif(1,0,1),
                                  p = c(NA,runif(ncol(y)-1,0,1)))

# parameters to monitor
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

# create model as an R object (uncompiled model)
survival <- nimbleModel(code = hmm.survival,
                        data = my.data,
                        constants = my.constants,
                        inits = initial.values())

survival$calculate() 

# run NIMBLE
mcmc.phipt <- nimbleMCMC(code = hmm.survival, 
                        constants = my.constants,
                        data = my.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = n.iter,
                        nburnin = n.burnin, 
                        nchains = n.chains)

# numerical summaries
MCMCsummary(mcmc.phipt, round = 2)

# mean   sd 2.5%  50% 97.5% Rhat n.eff
# p[1] 1.00 0.00 1.00 1.00  1.00  NaN     0
# p[2] 0.89 0.04 0.81 0.89  0.95 1.01   667
# p[3] 0.85 0.06 0.73 0.86  0.94 1.00   630
# p[4] 0.94 0.05 0.81 0.95  1.00 1.03   350
# p[5] 0.84 0.12 0.54 0.86  0.99 1.01   424
# p[6] 0.75 0.17 0.36 0.78  0.99 1.02   539
# p[7] 0.36 0.25 0.02 0.30  0.91 1.00   545
# phi  0.57 0.03 0.52 0.57  0.62 1.01   582

