# load nimbleEcology in case we forgot previously
library(crdata)
library(nimbleEcology)
library(tidyverse)
library(MCMCvis)

# data
y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()
# occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))
y <- y[ first != ncol(y), ] # get rid of individuals for which first detection = last capture
first <- apply(y, 1, function(x) min(which(x !=0)))

#----------------------------------------------------------------------
# NIMBLE code
hmm.phip.nimbleecology <- nimbleCode({
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

survival <- nimbleModel(code = hmm.phip.nimbleecology,
                        data = my.data,
                        constants = my.constants,
                        inits = initial.values())
survival$calculate()
# [1] -334.7412

# run NIMBLE
mcmc.phip.nimbleecology <- nimbleMCMC(code = hmm.phip.nimbleecology, 
                                      constants = my.constants,
                                      data = my.data,              
                                      inits = initial.values,
                                      monitors = parameters.to.save,
                                      niter = n.iter,
                                      nburnin = n.burnin, 
                                      nchains = n.chains)
# numerical summaries
MCMCsummary(mcmc.phip.nimbleecology, round = 2)



#------------------------------------------------------------

# load nimbleEcology in case we forgot previously
library(crdata)
library(nimbleEcology)
library(tidyverse)
library(MCMCvis)

# data
y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()
# occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))
y <- y[ first != ncol(y), ] # get rid of individuals for which first detection = last capture
first <- apply(y, 1, function(x) min(which(x !=0)))

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
                            probObs = omega[1:2,1:2,first[i]:K], # observation matrix
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

survival$getNodeNames()

survival$getParam(node = "phi", param = "max")
survival$getParam(node = "y[248, 6:7]", param = "probTrans")

survival$getParam(node = "y[248, 6:7]", param = "probObs")




survival$getVarNames()

survival$calculate()
# [1] -359.2902

#survival$plotGraph()

maps <- survival$modelDef$maps
ls(maps)

maps$elementNames      ## a name for every scalar element

maps$nodeNames         ## a name for everything in the graph

maps$types             ## labels of node types





# run NIMBLE
mcmc.phip.nimbleecology <- nimbleMCMC(code = hmm.survival, 
                                      constants = my.constants,
                                      data = my.data,              
                                      inits = initial.values,
                                      monitors = parameters.to.save,
                                      niter = n.iter,
                                      nburnin = n.burnin, 
                                      nchains = n.chains)

# numerical summaries
MCMCsummary(mcmc.phip.nimbleecology, round = 2)



#--------------------------------------------------


# load nimbleEcology in case we forgot previously
library(nimbleEcology)
library(tidyverse)
library(MCMCvis)



set.seed(2022) # for reproducibility
nocc <- 7 # nb of winters or sampling occasions
nind <- 50 # nb of animals
p <- 0.6 # detection prob
phi <- 0.8 # survival prob
# Vector of initial states probabilities
delta <- c(1,0) # all individuals are alive in first winter
# Transition matrix
Gamma <- matrix(NA, 2, 2)
Gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
Gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
Gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
Gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
# Observation matrix 
Omega <- matrix(NA, 2, 2)
Omega[1,1] <- 1 - p      # Pr(alive t -> non-detected t)
Omega[1,2] <- p          # Pr(alive t -> detected t)
Omega[2,1] <- 1          # Pr(dead t -> non-detected t)
Omega[2,2] <- 0          # Pr(dead t -> detected t)
# Matrix of states
z <- matrix(NA, nrow = nind, ncol = nocc)
y <- z
y[,1] <- 2 # all individuals are detected in first winter, as we condition on first detection
for (i in 1:nind){
  z[i,1] <- rcat(n = 1, prob = delta) # 1 for sure
  for (t in 2:nocc){
    # state at t given state at t-1
    z[i,t] <- rcat(n = 1, prob = Gamma[z[i,t-1],1:2]) 
    # observation at t given state at t
    y[i,t] <- rcat(n = 1, prob = Omega[z[i,t],1:2]) 
  }
}
y
y <- y - 1 # non-detection = 0, detection = 1

# occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))
y <- y[ first != ncol(y), ] # get rid of individuals for which first detection = last capture
first <- apply(y, 1, function(x) min(which(x !=0)))

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
    y[i,1:K] ~ dHMMo(init = delta[1:2], # vector of initial state probabilities
                            probObs = omega[1:2,1:2,1:K], # observation matrix
                            probTrans = gamma[1:2,1:2], # transition matrix
                            len = K, # nb of sampling occasions
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

survival$getNodeNames()

survival$getParam(node = "phi", param = "max")
survival$getParam(node = "y[248, 6:7]", param = "probTrans")

survival$getParam(node = "y[248, 6:7]", param = "probObs")




survival$getVarNames()

survival$calculate()
# [1] -359.2902

#survival$plotGraph()

maps <- survival$modelDef$maps
ls(maps)

maps$elementNames      ## a name for every scalar element

maps$nodeNames         ## a name for everything in the graph

maps$types             ## labels of node types





# run NIMBLE
mcmc.phip.nimbleecology <- nimbleMCMC(code = hmm.survival, 
                                      constants = my.constants,
                                      data = my.data,              
                                      inits = initial.values,
                                      monitors = parameters.to.save,
                                      niter = n.iter,
                                      nburnin = n.burnin, 
                                      nchains = n.chains)

# numerical summaries
MCMCsummary(mcmc.phip.nimbleecology, round = 2)

