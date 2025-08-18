library(nimble)
library(MCMCvis)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# Alive and dead {#survival}

## Introduction

## The Cormack-Jolly-Seber (CJS) model

## Capture-recapture data {#crdataeg}

# The data look like: 
dipper <- read_csv("dipper.csv")
dipper

dipper %>%
  ggplot() + 
  aes(x = wing_length, fill = sex) +
  geom_histogram(color = "white", binwidth = 1) +
  labs(x = "wing length", fill = "") + 
  scale_fill_manual(labels = c("females", "males"),
                    values = c("darkgreen", "orange"))

## Fitting the CJS model to the dipper data with NIMBLE

# Overall, the code looks like:
hmm.phitpt <- nimbleCode({
  # parameters
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    phi[t] ~ dunif(0, 1)        # prior survival
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
    p[t] ~ dunif(0, 1)          # prior detection
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0           # Pr(dead t -> detected t)
  }
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
    }
  }
})

# We extract the detections and non-detections from the data:
y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()

# We get the occasion of first capture for all individuals, by finding the position of detections in the encounter history (`which(x !=0)`), and keeping the first one:
first <- apply(y, 1, function(x) min(which(x !=0)))

# Now we specify the constants:
my.constants <- list(N = nrow(y),   # number of animals
                     T = ncol(y),   # number of sampling occasions
                     first = first) # first capture for all animales

# We then put the data in a list. We add 1 to the data to code 
# non-detections as 1's detections as 2's (see Section \@ref(fittinghmmnimble)). 
my.data <- list(y = y + 1)

# Let's write a function for the initial values. For the latent states, 
# we go the easy way, and say that all individuals are alive through 
# the study period:
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(my.constants$T-1,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)

# We specify the parameters we would like to monitor, survival and 
# detection probabilities here:
parameters.to.save <- c("phi", "p")

# We provide MCMC details:
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# And we run NIMBLE:
mcmc.phitpt <- nimbleMCMC(code = hmm.phitpt,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

# We may have a look to the numerical summaries: 
MCMCsummary(mcmc.phitpt, params = c("phi","p"), round = 2)

# You may have noticed the small effective sample size for the last survival 
# (`phi[6]`) and detection (`p[6]`) probabilities. Let's have a look 
# to mixing for parameter `phi[6]` for example:
priors <- runif(3000, 0, 1)
MCMCtrace(object = mcmc.phitpt,
          ISB = FALSE,
          exact = TRUE, 
          params = c("phi[6]"),
          pdf = FALSE, 
          priors = priors)

## CJS model derivatives {#cjsderivatives}

# But I realize that we did not actually fit the model with constant 
# parameters from Chapter \@ref(hmmcapturerecapture). Let's do it. 
# You should be familiar with the process by now: 

# NIMBLE code 
hmm.phip <- nimbleCode({
  phi ~ dunif(0, 1)      # prior survival
  p ~ dunif(0, 1)        # prior detection
  # likelihood
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})
# occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))
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
                                  p = runif(1,0,1),
                                  z = zinits)
# parameters to monitor
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2
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

# Let's now turn to the model with time-varying survival and constant detection. 
# We modify the CJS model NIMBLE code by no longer having the observation 
# matrix time-specific. I'm just providing the model code to save space:
hmm.phitp <- nimbleCode({
  for (t in 1:(T-1)){
    phi[t] ~ dunif(0, 1) # prior survival
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0        # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1        # Pr(dead t -> dead t+1)
  }
  p ~ dunif(0, 1) # prior detection
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

initial.values <- function() list(phi = runif(my.constants$T-1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")
parameters.to.save

n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

mcmc.phitp <- nimbleMCMC(code = hmm.phitp, 
                         constants = my.constants,
                         data = my.data,              
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin, 
                         nchains = n.chains)

# We obtain the following numerical summaries for parameters, 
# confirming high detection and temporal variation in survival:
MCMCsummary(object = mcmc.phitp, params = c("phi","p"), round = 2)

# Now the model with time-varying detection and constant survival, 
# for which the NIMBLE code has a constant over time transition matrix:
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

initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")
parameters.to.save

n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

mcmc.phipt <- nimbleMCMC(code = hmm.phipt, 
                         constants = my.constants,
                         data = my.data,              
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin, 
                         nchains = n.chains)

# Numerical summaries for the parameters are:
MCMCsummary(object = mcmc.phipt, params = c("phi","p"), round = 2)

# Before we move to the next section, you might ask how to fit these 
# models with `nimbleEcology` as in Section \@ref(nimblemarginalization). 
# ... For example, to implement the model with constant survival and 
# detection probabilities, you would use `dCJSss()`:

# load nimbleEcology in case we forgot previously
library(nimbleEcology)
# data
y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()
y <- y[ first != ncol(y), ] # get rid of individuals for which first detection = last capture
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
# occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))
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
# parameters to monitor
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2
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

## Model comparison with WAIC {#waic}

# Re-run the four models to calculate the WAIC value for each of them

hmm.phitpt <- nimbleCode({
  # parameters
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    phi[t] ~ dunif(0, 1)        # prior survival
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
    p[t] ~ dunif(0, 1)          # prior detection
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0           # Pr(dead t -> detected t)
  }
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
    }
  }
})
first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y),   # number of animals
                     T = ncol(y),   # number of sampling occasions
                     first = first) # first capture for all animales
my.data <- list(y = y + 1)
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(my.constants$T-1,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2
mcmc.phitpt <- nimbleMCMC(code = hmm.phitpt,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains,
                          WAIC = TRUE)
hmm.phip <- nimbleCode({
  phi ~ dunif(0, 1)      # prior survival
  p ~ dunif(0, 1)        # prior detection
  # likelihood
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})
first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first)
my.data <- list(y = y + 1)
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2
mcmc.phip <- nimbleMCMC(code = hmm.phip, 
                        constants = my.constants,
                        data = my.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = n.iter,
                        nburnin = n.burnin, 
                        nchains = n.chains,
                        WAIC = TRUE)
hmm.phitp <- nimbleCode({
  for (t in 1:(T-1)){
    phi[t] ~ dunif(0, 1) # prior survival
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0        # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1        # Pr(dead t -> dead t+1)
  }
  p ~ dunif(0, 1) # prior detection
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})
initial.values <- function() list(phi = runif(my.constants$T-1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")
parameters.to.save
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2
mcmc.phitp <- nimbleMCMC(code = hmm.phitp, 
                         constants = my.constants,
                         data = my.data,              
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin, 
                         nchains = n.chains,
                         WAIC = TRUE)
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
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")
parameters.to.save
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2
mcmc.phipt <- nimbleMCMC(code = hmm.phipt, 
                         constants = my.constants,
                         data = my.data,              
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin, 
                         nchains = n.chains,
                         WAIC = TRUE)

#save(mcmc.phip, mcmc.phitp, mcmc.phipt, mcmc.phitpt, file = 'dipper_waic.RData')

data.frame(model = c("both survival & detection constant",
                     "time-dependent survival, constant detection",
                     "constant survival, time-dependent detection",
                     "both survival & detection time-dependent"),
           WAIC = c(mcmc.phip$WAIC$WAIC,
                    mcmc.phitp$WAIC$WAIC,
                    mcmc.phipt$WAIC$WAIC,
                    mcmc.phitpt$WAIC$WAIC))

## Goodness of fit {#gof}

### Posterior predictive checks

### Classical tests

# To install the R2ucare package:
if(!require(devtools)) install.packages("devtools")
devtools::install_github("oliviergimenez/R2ucare")

# We load the package `R2ucare`:
library(R2ucare)

# We get the capture-recapture data:

# capture-recapture data
dipper <- read_csv("dipper.csv")
# get encounter histories
y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()
# number of birds with a particular capture-recapture history
size <- rep(1, nrow(y))

# The overall test shows that we cannot reject the hypothesis that 
# the CJS models fits the data well
overall_CJS(y, size)

# We may perform a test specifically to assess a transient effect:
test3sr(y, size)

# Or trap-dependence:
test2ct(y, size)

### Design considerations

## Covariates {#covariates}

### Temporal covariates

#### Discrete 

y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()
first <- apply(y, 1, function(x) min(which(x !=0)))
my.data <- list(y = y + 1)

# Let's use a covariate `flood` that contains 1s and 2s, 
# indicating whether we are in a flood or nonflood year for each year: 
# 1 if nonflood year, and 2 if flood year. 
flood <- c(1, # 1981-1982 (nonflood)
           2, # 1982-1983 (flood)
           2, # 1983-1984 (flood)
           1, # 1984-1985 (nonflood)
           1, # 1985-1986 (nonflood)
           1) # 1986-1987 (nonflood)

# Then we write the model code: 
hmm.phifloodp <- nimbleCode({
  delta[1] <- 1                        # Pr(alive t = 1) = 1
  delta[2] <- 0                        # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    gamma[1,1,t] <- phi[flood[t]]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[flood[t]]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0                  # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1                  # Pr(dead t -> dead t+1)
  }
  p ~ dunif(0, 1)        # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  phi[1] ~ dunif(0, 1)   # prior for survival in nonflood years
  phi[2] ~ dunif(0, 1)   # prior for survival in flood years
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# The data

# Let's provide the constants in a list: 
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     flood = flood)

# Now a function to generate initial values:
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

# The parameters to be monitored:
parameters.to.save <- c("p", "phi")

# The MCMC details:
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# We're all set, and we run NIMBLE: 
mcmc.phifloodp <- nimbleMCMC(code = hmm.phifloodp, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

# The numerical summaries are given by:
MCMCsummary(mcmc.phifloodp, round = 2)

# Having a look to the numerical summaries, we see that as expected, 
# survival in flood years (`phi[2]`) is much lower than survival in 
# non-flood years (`phi[1]`). You could formally test this difference 
# by considering the difference `phi[1] - phi[2]`. Alternatively, 
# this can be done afterwards and calculating the probability that 
# this difference is positive (or `phi[1] > phi[2]`). Using a single 
# chain for convenience, we do:
phi1 <- mcmc.phifloodp$chain1[,'phi[1]']
phi2 <- mcmc.phifloodp$chain1[,'phi[2]']
mean(phi1 - phi2 > 0)

# Another point of attention is the prior we will assign to regression 
# coefficients. We no longer assign a prior to survival directly like 
# previously, but we need to assign prior to the $\beta$'s which will 
# induce some prior on survival. And sometimes, your priors on the 
# regression coefficients are non-informative, but the prior on survival 
# is not. Consider for example the case of a single intercept with no 
# covariate. If you assign as a prior to this regression coefficient a 
# normal distribution with mean 0 and large standard deviation (left 
# figure below), which would be my first reflex, then you end up with a 
# very informative prior on survival with a bathtub shape, putting much 
# importance on low and high values (right figure below):

# 1000 random values from a N(0,10)
intercept <- rnorm(1000, mean = 0, sd = 10) 
# plogis() is the inverse-logit function in R
survival <- plogis(intercept) 
df <- data.frame(intercept = intercept, survival = survival)
plot1 <- df %>%
  ggplot(aes(x = intercept)) +
  geom_density() +
  labs(x = "prior on intercept")
plot2 <- df %>%
  ggplot(aes(x = survival)) +
  geom_density() +
  labs(x = "prior on survival")
plot1 + plot2

# Now if you go for a lower standard deviation for the intercept 
# prior (left figure below), e.g. 1.5, the prior on survival is 
# non-informative, looking like a uniform distribution between 
# 0 and 1 (right figure below):

set.seed(123)
# 1000 random values from a N(0,1.5)
intercept <- rnorm(1000, mean = 0, sd = 1.5) 
# plogis() is the inverse-logit function in R
survival <- plogis(intercept) 
df <- data.frame(intercept = intercept, survival = survival)
plot1 <- df %>%
  ggplot(aes(x = intercept)) +
  geom_density() +
  labs(x = "prior on intercept")
plot2 <- df %>%
  ggplot(aes(x = survival)) +
  geom_density() +
  labs(x = "prior on survival")
plot1 + plot2

# Now let's go back to our model. We first define our flood covariate 
# with 0 if nonflood year, and 1 if flood year:
flood <- c(0, # 1981-1982 (nonflood)
           1, # 1982-1983 (flood)
           1, # 1983-1984 (flood)
           0, # 1984-1985 (nonflood)
           0, # 1985-1986 (nonflood)
           0) # 1986-1987 (nonflood)

# Then we write the NIMBLE code:
hmm.phifloodp <- nimbleCode({
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    logit(phi[t]) <- beta[1] + beta[2] * flood[t]
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
  }
  p ~ dunif(0, 1)               # prior detection
  omega[1,1] <- 1 - p           # Pr(alive t -> non-detected t)
  omega[1,2] <- p               # Pr(alive t -> detected t)
  omega[2,1] <- 1               # Pr(dead t -> non-detected t)
  omega[2,2] <- 0               # Pr(dead t -> detected t)
  beta[1] ~ dnorm(0, sd = 1.5)  # prior intercept
  beta[2] ~ dnorm(0, sd = 1.5)  # prior slope
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# Let's put our constants in a list:
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     flood = flood)

# Then our function for generating initial values: 
initial.values <- function() list(beta = rnorm(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

# And the parameters to be monitored: 
parameters.to.save <- c("beta", "p", "phi")

# Finaly, we run NIMBLE: 
mcmc.phifloodp <- nimbleMCMC(code = hmm.phifloodp, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

# You may check that we get the same numerical summaries as above 
# for survival in nonflood years (`phi[1]`, `phi[4]`, `phi[5]` and 
# `phi[6]`) and flood years (`phi[2]` and `phi[3]`):
MCMCsummary(mcmc.phifloodp, round = 2)

# You may also check how to go from the $\beta$'s to the survival 
# probabilities $\phi$. Let's get the draws from the posterior distribution 
# of the $\beta$'s first: 
beta1 <- c(mcmc.phifloodp$chain1[,'beta[1]'], # beta1 chain 1
           mcmc.phifloodp$chain2[,'beta[1]']) # beta1 chain 2
beta2 <- c(mcmc.phifloodp$chain1[,'beta[2]'], # beta2 chain 1
           mcmc.phifloodp$chain2[,'beta[2]']) # beta2 chain 2

# Then apply the inverse-logit function to get survival in nonflood 
# years, e.g. its posterior mean and credible interval:
mean(plogis(beta1))
quantile(plogis(beta1), probs = c(2.5, 97.5)/100)

# Same thing for survival in flood years:
mean(plogis(beta1 + beta2))
quantile(plogis(beta1 + beta2), probs = c(2.5, 97.5)/100)

#### Continuous

# We build a covariate with water flow in liters per second measured 
# during the March to May period each year, starting with year 1982:

# water flow in L/s
water_flow <- c(443,  # March-May 1982
                1114, # March-May 1983
                529,  # March-May 1984
                434,  # March-May 1985
                627,  # March-May 1986
                466)  # March-May 1987

# Importantly, we standardize our covariate to improve convergence:
water_flow_st <- (water_flow - mean(water_flow))/sd(water_flow)

# Now we write the model code:
hmm.phiflowp <- nimbleCode({
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    logit(phi[t]) <- beta[1] + beta[2] * flow[t] 
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
  }
  p ~ dunif(0, 1)               # prior detection
  omega[1,1] <- 1 - p           # Pr(alive t -> non-detected t)
  omega[1,2] <- p               # Pr(alive t -> detected t)
  omega[2,1] <- 1               # Pr(dead t -> non-detected t)
  omega[2,2] <- 0               # Pr(dead t -> detected t)
  beta[1] ~ dnorm(0, 1.5)       # prior intercept
  beta[2] ~ dnorm(0, 1.5)       # prior slope
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# We put the constants in a list: 
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     flow = water_flow_st)

# Initial values as usual: 
initial.values <- function() list(beta = rnorm(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

# And parameters to be monitored:
parameters.to.save <- c("beta", "p", "phi")

# Eventually, we run NIMBLE: 
mcmc.phiflowp <- nimbleMCMC(code = hmm.phiflowp, 
                          constants = my.constants,
                          data = my.data,              
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin, 
                          nchains = n.chains)

# We can have a look to the results 
MCMCsummary(mcmc.phiflowp, round = 2)

# Through a caterpillar plot 
# of the regression parameters: 
MCMCplot(object = mcmc.phiflowp, params = "beta")

# Let's inspect the time-dependent survival probability:
MCMCplot(object = mcmc.phiflowp, params = "phi", ISB = TRUE)

### Individual covariates

#### Discrete

# We first consider a covariate `sex` that contains 1's and 2's 
# indicating the sex of each bird: 1 if male, and 2 if female. 
# We implement the model with sex effect using nested indexing, 
# similarly to the model with flood vs. nonflood years. The section 
# of the NIMBLE code that needs to be amended is:
hmm.phisexpt.ni <- nimbleCode({
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    phi[i] <- beta[sex[i]]
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0        # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1        # Pr(dead t -> dead t+1)
  }
  beta[1] ~ dunif(0,1)
  beta[2] ~ dunif(0,1)
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first,
                     sex = if_else(dipper$sex == "M", 1, 2)) # beta[1] male survival, beta[2] female survival 
my.data <- list(y = y + 1)

initial.values <- function() list(beta = runif(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
parameters.to.save <- c("beta", "p")

n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

mcmc.phisexp.ni <- nimbleMCMC(code = hmm.phisexpt.ni, 
                              constants = my.constants,
                              data = my.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin, 
                              nchains = n.chains)

# After running NIMBLE, we get:
MCMCsummary(object = mcmc.phisexp.ni, round = 2)


#### Continuous 

# We consider wing length here, and more precisely its measurement 
# at first detection. We first  standardize the covariate: 
wing.length.st <- as.vector(scale(dipper$wing_length))
head(wing.length.st)

# Now we write the model:
hmm.phiwlp <- nimbleCode({
    p ~ dunif(0, 1)             # prior detection
    omega[1,1] <- 1 - p         # Pr(alive t -> non-detected t)
    omega[1,2] <- p             # Pr(alive t -> detected t)
    omega[2,1] <- 1             # Pr(dead t -> non-detected t)
    omega[2,2] <- 0             # Pr(dead t -> detected t)
  for (i in 1:N){
    logit(phi[i]) <- beta[1] + beta[2] * winglength[i]
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
  }
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# We put the constants in a list:
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first,
                     winglength = wing.length.st)

# We write a function for generating initial values:
initial.values <- function() list(beta = rnorm(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

# And we run NIMBLE:
mcmc.phiwlp <- nimbleMCMC(code = hmm.phiwlp, 
                          constants = my.constants,
                          data = my.data,              
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin, 
                          nchains = n.chains)

# Let's inspect the numerical summaries for the regression parameters: 
MCMCsummary(mcmc.phiwlp, params = "beta", round = 2)

# Let's plot the relationship between survival and wing length. 
# First, we gather the values generated from the posterior distribution 
# of the regression parameters in the two chains:
beta1 <- c(mcmc.phiwlp$chain1[,'beta[1]'], # intercept, chain 1
           mcmc.phiwlp$chain2[,'beta[1]']) # intercept, chain 2
beta2 <- c(mcmc.phiwlp$chain1[,'beta[2]'], # slope, chain 1
           mcmc.phiwlp$chain2[,'beta[2]']) # slope, chain 2

# Then we define a grid of values for wing length, and predict survival 
# for each MCMC iteration:
predicted_survival <- matrix(NA, 
                             nrow = length(beta1), 
                             ncol = length(my.constants$winglength))
for (i in 1:length(beta1)){
  for (j in 1:length(my.constants$winglength)){
    predicted_survival[i,j] <- plogis(beta1[i] + 
                               beta2[i] * my.constants$winglength[j])
  }
}

# Now we calculate posterior mean and the credible interval 
# (note the ordering):
mean_survival <- apply(predicted_survival, 2, mean)
lci <- apply(predicted_survival, 2, quantile, prob = 2.5/100)
uci <- apply(predicted_survival, 2, quantile, prob = 97.5/100)
ord <- order(my.constants$winglength)
df <- data.frame(wing_length = my.constants$winglength[ord],
                 survival = mean_survival[ord],
                 lci = lci[ord],
                 uci = uci[ord])

# Now time to visualize: 
df %>%
  ggplot() + 
  aes(x = wing_length, y = survival) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lci, ymax = uci), 
              fill = "grey70", 
              alpha = 0.5) + 
  ylim(0,1) + 
  labs(x = "wing length", y = "estimated survival")

### Several covariates

# You may wish to have an effect of both sex and wing length in a model. 
# The NIMBLE code is:
hmm.phisexwlp <- nimbleCode({
  p ~ dunif(0, 1)               # prior detection
  omega[1,1] <- 1 - p           # Pr(alive t -> non-detected t)
  omega[1,2] <- p               # Pr(alive t -> detected t)
  omega[2,1] <- 1               # Pr(dead t -> non-detected t)
  omega[2,2] <- 0               # Pr(dead t -> detected t)
  for (i in 1:N){
    logit(phi[i]) <- beta[1] + beta[2] * sex[i] + beta[3] * winglength[i]
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
  }
  beta[1] ~ dnorm(mean = 0, sd = 1.5) # intercept male
  beta[2] ~ dnorm(mean = 0, sd = 1.5) # difference bw male and female
  beta[3] ~ dnorm(mean = 0, sd = 1.5) # slope wing length
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# We put constants and data in lists:
first <- apply(y, 1, function(x) min(which(x !=0)))
wing.length.st <- as.vector(scale(dipper$wing_length))
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first,
                     winglength = wing.length.st,
                     sex = if_else(dipper$sex == "M", 0, 1))
my.data <- list(y = y + 1)

# We write a fuction to generate initial values:
zinits <- y
zinits[zinits == 0] <- 1
initial.values <- function() list(beta = rnorm(3,0,2),
                                  p = runif(1,0,1),
                                  z = zinits)

# We specify the parameters to be monitored:
parameters.to.save <- c("beta", "p")

# MCMC details
n.iter <- 5000*4
n.burnin <- 1000
n.chains <- 2

# And now we run NIMBLE:
mcmc.phisexwlp <- nimbleMCMC(code = hmm.phisexwlp, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

# Let's display numerical summaries for all parameters:
MCMCsummary(mcmc.phisexwlp, round = 2)

# Let's visualize survival as a function of wing length for both sexes. 
beta1 <- c(mcmc.phisexwlp$chain1[,'beta[1]'], # beta1 chain 1
           mcmc.phisexwlp$chain2[,'beta[1]']) # beta1 chain 2
beta2 <- c(mcmc.phisexwlp$chain1[,'beta[2]'], # beta2 chain 1
           mcmc.phisexwlp$chain2[,'beta[2]']) # beta2 chain 2
beta3 <- c(mcmc.phisexwlp$chain1[,'beta[3]'], # beta3 chain 1
           mcmc.phisexwlp$chain2[,'beta[3]']) # beta3 chain 2
predicted_survivalM <- matrix(NA, nrow = length(beta1), 
                              ncol = length(my.constants$winglength))
predicted_survivalF <- matrix(NA, nrow = length(beta1), 
                              ncol = length(my.constants$winglength))
for (i in 1:length(beta1)){
  for (j in 1:length(my.constants$winglength)){
    predicted_survivalM[i,j] <- plogis(beta1[i] + 
                                beta3[i] * my.constants$winglength[j]) 
    predicted_survivalF[i,j] <- plogis(beta1[i] + 
                                beta2[i] + 
                                beta3[i] * my.constants$winglength[j])
  }
}
mean_survivalM <- apply(predicted_survivalM, 2, mean)
lciM <- apply(predicted_survivalM, 2, quantile, prob = 2.5/100)
uciM <- apply(predicted_survivalM, 2, quantile, prob = 97.5/100)
mean_survivalF <- apply(predicted_survivalF, 2, mean)
lciF <- apply(predicted_survivalF, 2, quantile, prob = 2.5/100)
uciF <- apply(predicted_survivalF, 2, quantile, prob = 97.5/100)
ord <- order(my.constants$winglength)
df <- data.frame(wing_length = c(my.constants$winglength[ord], 
                                 my.constants$winglength[ord]),
                 survival = c(mean_survivalM[ord], 
                              mean_survivalF[ord]),
                 lci = c(lciM[ord],lciF[ord]),
                 uci = c(uciM[ord],uciF[ord]),
                 sex = c(rep("male", length(mean_survivalM)), 
                         rep("female", length(mean_survivalF))))
df %>%
  ggplot() + 
  aes(x = wing_length, y = survival, color = sex) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = sex), alpha = 0.5) + 
  ylim(0,1) + 
  labs(x = "wing length", y = "estimated survival", color = "", fill = "")

# Note that the two curves are not exactly parallel because we 
# back-transformed the linear part of the relationship between survival 
# and wing length. You may check that parallelism occurs on the logit scale:
predicted_lsurvivalM <- matrix(NA, nrow = length(beta1), 
                              ncol = length(my.constants$winglength))
predicted_lsurvivalF <- matrix(NA, nrow = length(beta1), 
                              ncol = length(my.constants$winglength))
for (i in 1:length(beta1)){
  for (j in 1:length(my.constants$winglength)){
    predicted_lsurvivalM[i,j] <- beta1[i] + beta3[i] * my.constants$winglength[j] 
    predicted_lsurvivalF[i,j] <- beta1[i] + beta2[i] + beta3[i] * my.constants$winglength[j]
  }
}
mean_lsurvivalM <- apply(predicted_lsurvivalM, 2, mean)
mean_lsurvivalF <- apply(predicted_lsurvivalF, 2, mean)
ord <- order(my.constants$winglength)
df <- data.frame(wing_length = c(my.constants$winglength[ord], 
                                 my.constants$winglength[ord]),
                 survival = c(mean_lsurvivalM[ord], 
                              mean_lsurvivalF[ord]),
                 sex = c(rep("male", length(mean_lsurvivalM)), 
                         rep("female", length(mean_lsurvivalF))))
df %>%
  ggplot() + 
  aes(x = wing_length, y = survival, color = sex) + 
  geom_line() + 
  ylim(-2,2) + 
  labs(x = "wing length", 
       y = "estimated survival (on the logit scale)", 
       color = "")

### Random effects {#randomeffects}

# Let's go back to our dipper example with wing length as a covariate, 
# and write the NIMBLE code: 
hmm.phiwlrep <- nimbleCode({
    p ~ dunif(0, 1)             # prior detection
    omega[1,1] <- 1 - p         # Pr(alive t -> non-detected t)
    omega[1,2] <- p             # Pr(alive t -> detected t)
    omega[2,1] <- 1             # Pr(dead t -> non-detected t)
    omega[2,2] <- 0             # Pr(dead t -> detected t)
  for (i in 1:N){
    logit(phi[i]) <- beta[1] + beta[2] * winglength[i] + eps[i]
    eps[i] ~ dnorm(mean = 0, sd = sdeps)
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
  }
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  sdeps ~ dunif(0, 10)
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# We now write a function for generating initial values:
initial.values <- function() list(beta = rnorm(2,0,1.5),
                                  sdeps = runif(1,0,3),
                                  p = runif(1,0,1),
                                  z = zinits)

# We specify the parameters to be monitored:
parameters.to.save <- c("beta", "sdeps", "p")

# We increase the number of iterations and the length of the burn-in 
# period to reach better convergence:
n.iter <- 10000
n.burnin <- 5000
n.chains <- 2

# And finally, we run NIMBLE:
mcmc.phiwlrep <- nimbleMCMC(code = hmm.phiwlrep, 
                            constants = my.constants,
                            data = my.data,              
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            nburnin = n.burnin, 
                            nchains = n.chains)

# We inspect the numerical summaries: 
MCMCsummary(mcmc.phiwlrep, round = 2)

# In the *centered parameterization* we've used so far, the random effect 
# is expressed directly in terms of its variance (or standard deviation). 
# This often creates strong correlations between the random effects and their 
# variance parameter, which in turn can slow down mixing and reduce the effective 
# sample size. In the *non-centered parameterization*, we instead separate the 
# structure of the random effect from its scale. In practice, we introduce a 
# standardized random effect $\varepsilon_i \sim N(0,1)$ and then multiply it by 
# the standard deviation $\sigma$. This decouples the randomness from the scale 
# parameter and usually improves sampling efficiency. In code, this looks like:
hmm.phiwlrep <- nimbleCode({
  p ~ dunif(0, 1)             # prior detection
  omega[1,1] <- 1 - p         # Pr(alive t -> non-detected t)
  omega[1,2] <- p             # Pr(alive t -> detected t)
  omega[2,1] <- 1             # Pr(dead t -> non-detected t)
  omega[2,2] <- 0             # Pr(dead t -> detected t)
  for (i in 1:N){
    logit(phi[i]) <- beta[1] + beta[2] * winglength[i] + + sdeps * eps[i]
    eps[i] ~ dnorm(mean = 0, sd = 1)
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1           # Pr(dead t -> dead t+1)
  }
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  sdeps ~ dunif(0, 10)
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# We now write a function for generating initial values:
initial.values <- function() list(beta = rnorm(2,0,1.5),
                                  sdeps = runif(1,0,3),
                                  p = runif(1,0,1),
                                  z = zinits)

# We specify the parameters to be monitored:
parameters.to.save <- c("beta", "sdeps", "p")

# We increase the number of iterations and the length of the burn-in 
# period to reach better convergence:
n.iter <- 10000
n.burnin <- 5000
n.chains <- 2

# And finally, we run NIMBLE:
mcmc.phiwlrep <- nimbleMCMC(code = hmm.phiwlrep, 
                            constants = my.constants,
                            data = my.data,              
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            nburnin = n.burnin, 
                            nchains = n.chains)

# After running NIMBLE, we inspect the numerical summaries, and 
# see that effective sample sizes are much better:
MCMCsummary(mcmc.phiwlrep, round = 2)

sdeps <- c(mcmc.phiwlrep$chain1[,"sdeps"], mcmc.phiwlrep$chain2[,"sdeps"])
sdeps %>%
  as_tibble() %>%
  ggplot() + 
  aes(x = value) + 
  geom_histogram(color = "white", binwidth = .03, fill = "gray70") + 
  geom_density(aes(y = .03 * ..count..))

### Individual time-varying covariates {#agecov}

# The convenient thing is that age has no missing value because age at 
# $t+1$ is just age at $t$ to which we add 1. This suggests a way 
# to code the age effect in NIMBLE as follows: 
hmm.phiage.in <- nimbleCode({
  p ~ dunif(0, 1)                     # prior detection
  omega[1,1] <- 1 - p                 # Pr(alive t -> non-detected t)
  omega[1,2] <- p                     # Pr(alive t -> detected t)
  omega[2,1] <- 1                     # Pr(dead t -> non-detected t)
  omega[2,2] <- 0                     # Pr(dead t -> detected t)
  for (i in 1:N){
    for (t in first[i]:(T-1)){
    # phi1 = beta1 + beta2; phi1+ = beta1
    logit(phi[i,t]) <- beta[1] + beta[2] * equals(t, first[i]) 
    gamma[1,1,i,t] <- phi[i,t]        # Pr(alive t -> alive t+1)
    gamma[1,2,i,t] <- 1 - phi[i,t]    # Pr(alive t -> dead t+1)
    gamma[2,1,i,t] <- 0               # Pr(dead t -> alive t+1)
    gamma[2,2,i,t] <- 1               # Pr(dead t -> dead t+1)
    }
  }
  beta[1] ~ dnorm(mean = 0, sd = 1.5) # phi1+
  beta[2] ~ dnorm(mean = 0, sd = 1.5) # phi1 - phi1+
  phi1plus <- plogis(beta[1])         # phi1+
  phi1 <- plogis(beta[1] + beta[2])   # phi1
  delta[1] <- 1                       # Pr(alive t = 1) = 1
  delta[2] <- 0                       # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# We put all constants in a list:
first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first)

# And the data in a list:
my.data <- list(y = y + 1)

# We write a function to generate initial values:
zinits <- y
zinits[zinits == 0] <- 1
initial.values <- function() list(beta = rnorm(2,0,5),
                                  p = runif(1,0,1),
                                  z = zinits)

# And specify the parameters to be monitored:
parameters.to.save <- c("phi1", "phi1plus", "p")

# We now run NIMBLE:
mcmc.phi.age.in <- nimbleMCMC(code = hmm.phiage.in, 
                              constants = my.constants,
                              data = my.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin, 
                              nchains = n.chains)

# We display the results: 
MCMCsummary(mcmc.phi.age.in, round = 2)

# Another method to include an age effect is to create an individual by time 
# covariate and use nested indexing (as in the flood/nonflood example) 
# to distinguish survival over the interval after first detection from survival afterwards: 
age <- matrix(NA, nrow = nrow(y), ncol = ncol(y) - 1)
for (i in 1:nrow(age)){
  for (j in 1:ncol(age)){
    if (j == first[i]) age[i,j] <- 1 # age = 1
    if (j > first[i]) age[i,j] <- 2  # age > 1
  }
}

# Now we may write the NIMBLE code for this model. We just need to remember that 
# survival is no longer defined on the logit scale as in the previous model, 
# so we simply use uniform priors:  
hmm.phiage.out <- nimbleCode({
  p ~ dunif(0, 1)                   # prior detection
  omega[1,1] <- 1 - p               # Pr(alive t -> non-detected t)
  omega[1,2] <- p                   # Pr(alive t -> detected t)
  omega[2,1] <- 1                   # Pr(dead t -> non-detected t)
  omega[2,2] <- 0                   # Pr(dead t -> detected t)
  for (i in 1:N){
    for (t in first[i]:(T-1)){
    phi[i,t] <- beta[age[i,t]]      # beta1 = phi1, beta2 = phi1+
    gamma[1,1,i,t] <- phi[i,t]      # Pr(alive t -> alive t+1)
    gamma[1,2,i,t] <- 1 - phi[i,t]  # Pr(alive t -> dead t+1)
    gamma[2,1,i,t] <- 0             # Pr(dead t -> alive t+1)
    gamma[2,2,i,t] <- 1             # Pr(dead t -> dead t+1)
    }
  }
  beta[1] ~ dunif(0, 1) # phi1
  beta[2] ~ dunif(0, 1) # phi1+
  phi1 <- beta[1]
  phi1plus <- beta[2]
  delta[1] <- 1                     # Pr(alive t = 1) = 1
  delta[2] <- 0                     # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# We put all constants in a list, including the age covariate:
first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first,
                     age = age)

# We re-write a function to generate initial values: 
zinits <- y
zinits[zinits == 0] <- 1
initial.values <- function() list(beta = runif(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

# And we run NIMBLE:
mcmc.phi.age.out <- nimbleMCMC(code = hmm.phiage.out, 
                               constants = my.constants,
                               data = my.data,              
                               inits = initial.values,
                               monitors = parameters.to.save,
                               niter = n.iter,
                               nburnin = n.burnin, 
                               nchains = n.chains)

# We display numerical summaries for the model parameters, and acknowledge 
# that we obtain similar results to the other parameterization:
MCMCsummary(mcmc.phi.age.out, round = 2)

## Why Bayes? Incorporate prior information {#elicitprior}

### Prior elicitation

# Now let's assume that we had only the three first years of data, what would 
# have happened? We fit the model with constant parameters with both the non-informative 
# and informative priors to the dataset from which we delete the final 4 years 
# of data. Now the benefit of using the prior information becomes clear as the 
# credible interval when prior information is ignored has a width of 0.53, which 
# is more than twice as much as when prior information is used (0.24), illustrating 
# the increased precision provided by the prior. We may assess visually this gain 
# in precision by comparing the survival posterior distributions with and without informative prior:
  
hmm.phipriorp <- nimbleCode({
  phi ~ dnorm(mean = 0.57, sd = 0.073) # prior survival
  p ~ dunif(0, 1) # prior detection
  # likelihood
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  #  omega_init[1,1] <- 0    
  #  omega_init[1,2] <- 1        
  #  omega_init[2,1] <- 1        
  #  omega_init[2,2] <- 0        
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    #    y[i,first[i]] ~ dcat(omega_init[z[i,first[i]], 1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

hmm.phip <- nimbleCode({
  phi ~ dunif(0, 1) # prior survival
  p ~ dunif(0, 1) # prior detection
  # likelihood
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  #  omega_init[1,1] <- 0    
  #  omega_init[1,2] <- 1        
  #  omega_init[2,1] <- 1        
  #  omega_init[2,2] <- 0        
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    #    y[i,first[i]] ~ dcat(omega_init[z[i,first[i]], 1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

yshort <- y[,1:3]
mask <- apply(yshort, 1, sum)
yshort <- yshort[mask!=0,]
first <- apply(yshort, 1, function(x) min(which(x !=0)))

my.constants <- list(N = nrow(yshort), T = ncol(yshort), first = first)
my.constants
my.data <- list(y = yshort + 1)

zinits <- yshort + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")
parameters.to.save

n.iter <- 5000*4
n.burnin <- 1000
n.chains <- 2

mcmc.phip <- nimbleMCMC(code = hmm.phip, 
                        constants = my.constants,
                        data = my.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = n.iter,
                        nburnin = n.burnin, 
                        nchains = n.chains)

MCMCsummary(mcmc.phip, round = 2)
phinoprior <- c(mcmc.phip$chain1[,"phi"], mcmc.phip$chain2[,"phi"])


first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), T = ncol(y), first = first)
my.constants
my.data <- list(y = y + 1)

zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")

mcmc.phip <- nimbleMCMC(code = hmm.phip, 
                        constants = my.constants,
                        data = my.data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = n.iter,
                        nburnin = n.burnin, 
                        nchains = n.chains)

MCMCsummary(mcmc.phip, round = 2)
phiprior <- c(mcmc.phip$chain1[,"phi"], mcmc.phip$chain2[,"phi"])

df <- data.frame(posterior = c(phinoprior, phiprior),
                 type = c(rep("w/ vague prior", length(phinoprior)),
                          rep("w/ informative prior", length(phiprior))))
df %>%
  ggplot() +
  aes(x = posterior, fill = type) +
  geom_density(aes(y = ..density..),
               bins = 40,
               color = "white",
               alpha = 0.6) +
  labs(x = "survival", fill = "") +
  scale_fill_manual(values = wesanderson::wes_palette("Royal1")[2:1])

### Moment matching

# For our model, that means:
(alpha <- ( (1 - 0.57)/(0.073*0.073) - (1/0.57) )*0.57^2)
(beta <- alpha * ( (1/0.57) - 1))

# Now we simply have to use `phi ~ dbeta(25.6,19.3)` as a prior instead of 
# our `phi ~ dnorm(0.57, sd = 0.073)`.

## Summary

## Suggested reading
