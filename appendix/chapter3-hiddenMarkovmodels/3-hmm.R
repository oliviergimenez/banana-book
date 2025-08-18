library(nimble)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# Hidden Markov models {#hmmcapturerecapture}

## Introduction

## Longitudinal data

# Let us denote $z_{i,t} = 1$ if individual $i$ alive at winter $t$, and 
# $z_{i,t} = 2$ if dead. Then longitudinal data look like in the table below. 
# Each row is an individual $i$, and columns are for winters $t$, or sampling 
# occasions. Variable $z$ is indexed by both $i$ and $t$, and takes value 1 if 
# individual $i$ is alive in winter $t$, and 2 otherwise.

# 1 = alive, 2 = dead
nind <- 57
nocc <- 5
phi <- 0.8 # survival probability
delta <- c(1,0) # (Pr(alive at t = 1), Pr(dead at t = 1))
Gamma <- matrix(NA, 2, 2) # transition matrix
Gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
Gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
Gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
Gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
z <- matrix(NA, nrow = nind, ncol = nocc)
set.seed(2022)
for (i in 1:nind){
  z[i,1] <- nimble::rcat(n = 1, prob = delta) # 1 for sure
  for (t in 2:nocc){
    z[i,t] <- nimble::rcat(n = 1, prob = Gamma[z[i,t-1],1:2]) 
  }
}
colnames(z) <- paste0("winter ", 1:nocc)
z %>%
  as_tibble() %>%
  add_column(id = 1:nind, .before = "winter 1") #%>%
#  kableExtra::kable() %>%
#  kableExtra::scroll_box(width = "100%", height = "400px")
#  kableExtra::kable_styling(font_size = 8,
#                            latex_options = "scale_down")

## A Markov model for longitudinal data

### Assumptions

### Transition matrix

### Initial states

### Likelihood

### Example

## Bayesian formulation

rcat(n = 1, prob = c(0.1, 0.3, 0.6))

rcat(n = 20, prob = c(0.1, 0.1, 0.4, 0.2, 0.2))

## NIMBLE implementation

# How to implement in NIMBLE the Markov model we just built? We need to put 
# in place a few bricks before running our model. 

# Puting everything together, the NIMBLE code for our Markov model is:
markov.survival <- nimbleCode({
  phi ~ dunif(0, 1) # prior
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
  # likelihood
  for (i in 1:N){ # loop over individual i
    z[i,1] ~ dcat(delta[1:2]) # t = 1
    for (j in 2:T){ # loop over time t
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2]) # t = 2,...,T
    } # t
  } # i
})

# Now we're ready to resume our NIMBLE workflow. First we read in data. 
# The code I used to simulate the $z$ with survival $\phi = 0.8$ is as follows: 

# 1 = alive, 2 = dead
nind <- 57
nocc <- 5
phi <- 0.8 # survival probability
delta <- c(1,0) # (Pr(alive at t = 1), Pr(dead at t = 1))
Gamma <- matrix(NA, 2, 2) # transition matrix
Gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
Gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
Gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
Gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
z <- matrix(NA, nrow = nind, ncol = nocc)
set.seed(2022)
for (i in 1:nind){
  z[i,1] <- rcat(n = 1, prob = delta) # 1 for sure
  for (t in 2:nocc){
    z[i,t] <- rcat(n = 1, prob = Gamma[z[i,t-1],1:2]) 
  }
}
head(z) 

# Note the resemblance with the model NIMBLE code above. Because we have loops 
# and indices that do not change, we use constants as explained in 
# Section \@ref(start-nimble):
my.data <- list(z = z)
my.constants <- list(N = 57, T = 5)

# We also specify initial values for survival with a function:
initial.values <- function() list(phi = runif(1,0,1))
initial.values()

# There is a single parameter to monitor:
parameters.to.save <- c("phi")
parameters.to.save

# We run 2 chains with 5000 iterations including 1000 iterations as burnin:
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# Let's run NIMBLE:
mcmc.output <- nimbleMCMC(code = markov.survival,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

# Let's calculate the usual posterior numerical summaries for survival:
MCMCsummary(mcmc.output, round = 2)

## Hidden Markov models

### Capture-recapture data {#capturerecapturedata}

# Let's get back to the data we had in the previous section.
# The truth is in $z$ which contains the fate of all individuals with 
# $z = 1$ for alive, and $z = 2$ for dead:

colnames(z) <- paste0("winter ", 1:nocc)
z %>%
  as_tibble() %>%
  add_column(id = 1:nind, .before = "winter 1")

# The easiest connection is with dead animals which go undetected for sure. 
# Therefore when an animal is dead i.e. $z = 2$, it cannot be detected, 
# therefore $y = 0$:
z %>%
  as_tibble() %>%
  replace(. == 2, 0) %>%
  add_column(id = 1:nind, .before = "winter 1")

# Now alive animals may be detected or not. If an animal is alive $z = 1$, 
# it is detected $y = 1$ with probability $p$ or not $y = 0$ with probability 
# $1-p$. In our example, first detection coincides with first winter for 
# all individuals. 
set.seed(2022)
nocc <- 5
nind <- 57
p <- 0.6
phi <- 0.8
delta <- c(1,0)
Gamma <- matrix(NA, 2, 2)
Omega <- matrix(NA, 2, 2)
Gamma[1,1] <- phi # Pr(alive t -> alive t+1)
Gamma[1,2] <- 1 - phi # Pr(alive t -> dead t+1)
Gamma[2,1] <- 0 # Pr(dead t -> alive t+1)
Gamma[2,2] <- 1 # Pr(dead t -> dead t+1)
Omega[1,1] <- 1 - p # Pr(alive t -> non-detected t)
Omega[1,2] <- p # Pr(alive t -> detected t)
Omega[2,1] <- 1 # Pr(dead t -> non-detected t)
Omega[2,2] <- 0 # Pr(dead t -> detected t)
z <- matrix(NA, nrow = nind, ncol = nocc)
y <- z
y[,1] <- 2 # all animals detected in first winter
for (i in 1:nind){
 z[i,1] <- nimble::rcat(n = 1, prob = delta) # 1 for sure
 for (t in 2:nocc){
 z[i,t] <- nimble::rcat(n = 1, prob = Gamma[z[i,t-1],1:2]) 
 y[i,t] <- nimble::rcat(n = 1, prob = Omega[z[i,t],1:2]) 
 }
}
y <- y - 1 # non-detection = 0, detection = 1
colnames(y) <- paste0("winter ", 1:nocc)
nobs <- sum(apply(y,1,sum) != 0)
y <- y[apply(y,1,sum) !=0, ] # remove rows w/ non-detections only
y %>%
 as_tibble() %>%
 add_column(id = 1:nobs, .before = "winter 1")

### Observation matrix

### Hidden Markov model

### Likelihood {#likelihoodhmm}

## Fitting HMM with NIMBLE {#fittinghmmnimble}

# If we denote *first* the time of first detection, then our model so far is 
# written as follows:
hmm.survival <- nimbleCode({
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
  for (i in 1:N){
    z[i,1] ~ dcat(delta[1:2])
    for (j in 2:T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# Now we specify the constants:
my.constants <- list(N = nrow(y), T = 5)
my.constants

# The data are made of 0's for non-detections and 1's for detections. 
# To simulate the $y$, here is the code I used: 
set.seed(2022) # for reproducibility
nocc <- 5 # nb of winters or sampling occasions
nind <- 57 # nb of animals
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

# To use the categorical distribution, we need to code 1's and 2's. 
# We simply add 1 to get the correct format, that is $y = 1$ for non-detection 
# and $y = 2$ for detection: 
my.data <- list(y = y + 1)

# Now let's write a function for the initial values:
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

# We specify the parameters we'd like to monitor:
parameters.to.save <- c("phi", "p")
parameters.to.save

# We provide MCMC details:
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# At last, we're ready to run NIMBLE:
start_time <- Sys.time()
mcmc.output <- nimbleMCMC(code = hmm.survival,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
end_time <- Sys.time()
end_time - start_time

# We can have a look to numerical summaries:
MCMCsummary(mcmc.output, round = 2)

## Marginalization {#marginalization}

### Brute-force approach

### Forward algorithm {#forward-algorithm}

### NIMBLE implementation {#nimblemarginalization}

#### Do it yourself {#diymarginalisation}

# In NIMBLE, we use functions to implement the forward algorithm. The only 
# differences with the theory above is that i) we work on the log scale for 
# numerical stability and ii) we use a matrix formulation of the recurrence. 

# First we write the density function:
dHMMhomemade <- nimbleFunction(
  run = function(x = double(1), 
                 probInit = double(1), # vector of initial states
                 probObs = double(2), # observation matrix
                 probTrans = double(2), # transition matrix
                 len = double(0, default = 0), # number of sampling occasions
                 log = integer(0, default = 0)) {
    alpha <- probInit[1:2] # * probObs[1:2,x[1]] == 1 due to conditioning on first detection
    for (t in 2:len) {
      alpha[1:2] <- (alpha[1:2] %*% probTrans[1:2,1:2]) * probObs[1:2,x[t]]
    }
    logL <- log(sum(alpha[1:2]))
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

# Then we write a function to simulate values from a HMM:
rHMMhomemade <- nimbleFunction(
  run = function(n = integer(),
                 probInit = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = double(0, default = 0)) {
    returnType(double(1))
    z <- numeric(len)
    z[1] <- rcat(n = 1, prob = probInit[1:2]) # all individuals alive at t = 0
    y <- z
    y[1] <- 2 # all individuals are detected at t = 0
    for (t in 2:len){
      # state at t given state at t-1
      z[t] <- rcat(n = 1, prob = probTrans[z[t-1],1:2]) 
      # observation at t given state at t
      y[t] <- rcat(n = 1, prob = probObs[z[t],1:2]) 
    }
    return(y)
  })

# We assign these functions to the global R environment:
assign('dHMMhomemade', dHMMhomemade, .GlobalEnv)
assign('rHMMhomemade', rHMMhomemade, .GlobalEnv)

# Now we resume our workflow: 

# code
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
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    y[i,1:T] ~ dHMMhomemade(probInit = delta[1:2], 
                            probObs = omega[1:2,1:2], # observation matrix
                            probTrans = gamma[1:2,1:2], # transition matrix
                            len = T) # nb of sampling occasions
  }
})
# constants
my.constants <- list(N = nrow(y), T = 5)
# data
my.data <- list(y = y + 1)
# initial values - no need to specify values for z anymore
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1))
# parameters to save
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# And run NIMBLE:
start_time <- Sys.time()
mcmc.output <- nimbleMCMC(code = hmm.survival,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
end_time <- Sys.time()
end_time - start_time

# The numerical summaries are similar to those we obtained with the complete 
# likelihood, and effective samples sizes are larger denoting better mixing:
MCMCsummary(mcmc.output, round = 2)

#### Do it with `nimbleEcology` {#nimbleecologyintro}

# We load the package:
library(nimbleEcology)

# The NIMBLE code is:
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
  for (t in 2:5){
    omega[1,1,t] <- 1 - p    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1        # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0        # Pr(dead t -> detected t)
  }
  for (i in 1:N){
    y[i,1:5] ~ dHMMo(init = delta[1:2], # vector of initial state probabilities
                     probObs = omega[1:2,1:2,1:5], # observation matrix
                     probTrans = gamma[1:2,1:2], # transition matrix
                     len = 5, # nb of sampling occasions
                     checkRowSums = 0) # skip validity checks
  }
})

# constants
my.constants <- list(N = nrow(y))
# data
my.data <- list(y = y + 1)
# initial values - no need to specify values for z anymore
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1))
# parameters to save
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# Now we run NIMBLE:
start_time <- Sys.time()
mcmc.output <- nimbleMCMC(code = hmm.survival,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
end_time <- Sys.time()
end_time - start_time

# Now we display the numerical summaries of the posterior distributions: 
MCMCsummary(mcmc.output, round = 2)

## Pooled encounter histories {#pooled-likelihood}

# In this section, we amend the NIMBLE functions we wrote for marginalizing latent 
# states in Section \@ref(marginalization) to express the likelihood using 
# pooled encounter histories. We use a vector `size` that contains the number of 
# individuals with the same encounter history. 

# The density function is the function `dHMMhomemade` to which we add a 
# `size` argument, and raise the individual likelihood to the power `size`, 
# or multiply by `size` as we work on the log scale `log(sum(alpha[1:2])) * size`:
dHMMpooled <- nimbleFunction(
  run = function(x = double(1), 
                 probInit = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = double(0),
                 size = double(0),
                 log = integer(0, default = 0)) {
    alpha <- probInit[1:2]
    for (t in 2:len) {
      alpha[1:2] <- (alpha[1:2] %*% probTrans[1:2,1:2]) * probObs[1:2,x[t]]
    }
    logL <- log(sum(alpha[1:2])) * size
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

# The `rHMMhomemade` function is renamed `rHMMpooled` for compatibility but 
# remains unchanged:
rHMMpooled <- nimbleFunction(
  run = function(n = integer(),
                 probInit = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = double(0),
                 size = double(0)) {
    returnType(double(1))
    z <- numeric(len)
    z[1] <- rcat(n = 1, prob = probInit[1:2]) # all individuals alive at t = 0
    y <- z
    y[1] <- 2 # all individuals are detected at t = 0
    for (t in 2:len){
      # state at t given state at t-1
      z[t] <- rcat(n = 1, prob = probTrans[z[t-1],1:2]) 
      # observation at t given state at t
      y[t] <- rcat(n = 1, prob = probObs[z[t],1:2]) 
    }
    return(y)
  })

# We assign these two function to the global R environment so that we can use them:
assign('dHMMpooled', dHMMpooled, .GlobalEnv)
assign('rHMMpooled', rHMMpooled, .GlobalEnv)

# You can now plug your pooled HMM density function in your NIMBLE code:
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
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    y[i,1:T] ~ dHMMpooled(probInit = delta[1:2], 
                          probObs = omega[1:2,1:2], # observation matrix
                          probTrans = gamma[1:2,1:2], # transition matrix
                          len = T, # nb of sampling occasions
                          size = size[i]) # number of individuals with encounter history i
  }
})

# Before running NIMBLE, we need to actually pool individuals with the same 
# encounter history together:
y_pooled <- y %>% 
  as_tibble() %>% 
  group_by_all() %>% # group
  summarise(size = n()) %>% # count
  relocate(size) %>% # put size in front
  arrange(-size) %>% # sort along size
  as.matrix()
y_pooled

# For example, we have `r y_pooled[1,1]` individuals with encounter history (`r y_pooled[1,2:6]`).
y_pooled[1,1]

# Now you can resume the NIMBLE workflow:
my.constants <- list(N = nrow(y_pooled), T = 5, size = y_pooled[,'size'])
my.data <- list(y = y_pooled[,-1] + 1) # delete size from dataset
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1))
parameters.to.save <- c("phi", "p")
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2
start_time <- Sys.time()
mcmc.output <- nimbleMCMC(code = hmm.survival,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
end_time <- Sys.time()
end_time - start_time
MCMCsummary(mcmc.output, round = 2)

## Decoding after marginalization {#decoding}

### Theory {#viterbi-theory}

set.seed(2022) # for reproducibility
nocc <- 5 # nb of winters or sampling occasions
nind <- 57 # nb of animals
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
y[,1] <- 2 # all individuals are detected in first winter
for (i in 1:nind){
  z[i,1] <- rcat(n = 1, prob = delta) # 1 for sure
  for (t in 2:nocc){
    # state at t given state at t-1
    z[i,t] <- rcat(n = 1, prob = Gamma[z[i,t-1],1:2]) 
    # observation at t given state at t
    y[i,t] <- rcat(n = 1, prob = Omega[z[i,t],1:2]) 
  }
}
y <- y - 1

### Implementation

# Let's write a R function to implement the Viterbi algorithm. As parameters, 
# our function will take the transition and observation matrices, the vector 
# of initial state probabilities and the observed sequence of detections and 
# non-detections for which you aim to compute the sequence of states from which 
# it was most likely generated:

# getViterbi() returns sequence of states that most likely generated sequence of observations
# adapted from https://github.com/vbehnam/viterbi
getViterbi <- function(Omega, Gamma, delta, y) {
# Omega: transition matrix
# Gamma: observation matrix
# delta: vector of initial state probabilities
# y: observed sequence of detections and non-detections
  
# get number of states and sampling occasions
N <- nrow(Gamma)
T <- length(y)
  
# nu is the corresponding likelihood
nu <- matrix(0, nrow = N, ncol = T)
# zz contains the most likely states up until this point
zz <- matrix(0, nrow = N, ncol = T)
firstObs <- y[1]
  
# fill in first columns of both matrices
#nu[,1] <- initial * emission[,firstObs]
#zz[,1] <- 0
nu[,1] <- c(1,0) # initial = (1, 0) * emission[,firstObs] = (1, 0)
zz[,1] <- 1 # alive at first occasion

for (i in 2:T) {
    for (j in 1:N) {
      obs <- y[i]
      # initialize to -1, then overwritten by for loop coz all possible values are >= 0
      nu[j,i] <- -1
      # loop to find max and argmax for k
      for (k in 1:N) {
        value <- nu[k,i-1] * Gamma[k,j] * Omega[j,obs]
        if (value > nu[j,i]) {
          # maximizing for k
          nu[j,i] <- value
          # argmaximizing for k
          zz[j,i] <- k
        }
      }
    }
  }
  # mlp = most likely path
  mlp <- numeric(T)
  # argmax for stateSeq[,T]
  am <- which.max(nu[,T])
  mlp[T] <- zz[am,T]
  
  # backtrace using backpointers
  for (i in T:2) {
    zm <- which.max(nu[,i])
    mlp[i-1] <- zz[zm,i]
  }
  return(mlp)
}

# Let's test our `getViterbi()` function with our previous example. 
delta # Vector of initial states probabilities
Gamma # Transition matrix
Omega # Observation matrix
getViterbi(Omega = Omega, 
           Gamma = Gamma, 
           delta = delta, 
           y = y[15,] + 1)

# For both options, we will need the values from the posterior distributions 
# of survival and detection probabilities:
phi <- c(mcmc.output$chain1[,'phi'], mcmc.output$chain2[,'phi'])
p <- c(mcmc.output$chain1[,'p'], mcmc.output$chain2[,'p'])

### Compute first, average after {#compute-average}

# First option is to apply Viterbi to each MCMC sample, then to compute median 
# of the MCMC Viterbi paths for each observed sequence:
niter <- length(p)
T <- 5
res <- matrix(NA, nrow = nrow(y), ncol = T)
for (i in 1:nrow(y)){
  res_mcmc <- matrix(NA, nrow = niter, ncol = T)
  for (j in 1:niter){
    # Initial states
    delta <- c(1, 0)
    # Transition matrix
    transition <- matrix(NA, 2, 2)
    transition[1,1] <- phi[j]      # Pr(alive t -> alive t+1)
    transition[1,2] <- 1 - phi[j]  # Pr(alive t -> dead t+1)
    transition[2,1] <- 0        # Pr(dead t -> alive t+1)
    transition[2,2] <- 1        # Pr(dead t -> dead t+1)
    # Observation matrix 
    emission <- matrix(NA, 2, 2)
    emission[1,1] <- 1 - p[j]      # Pr(alive t -> non-detected t)
    emission[1,2] <- p[j]          # Pr(alive t -> detected t)
    emission[2,1] <- 1          # Pr(dead t -> non-detected t)
    emission[2,2] <- 0          # Pr(dead t -> detected t)
    res_mcmc[j,1:T] <- getViterbi(emission, transition, delta, y[i,] + 1)
  }
  res[i, 1:length(y[1,])] <- apply(res_mcmc, 2, median)
}

# You can compare the Viterbi decoding to the actual states $z$:
df <- expand.grid(X = 1:nrow(z), Y = 1:ncol(z))
df$Z <- as_factor(c(z - res))
myPallette <- RColorBrewer::brewer.pal(name = "RdBu", n = 4)
df %>% 
  ggplot() + 
  aes(x = X, y = Y, fill = Z) + 
  geom_tile() +
  scale_fill_manual(values = c("0" = "#9597f0", 
                               "1" = "#d4b4f6", 
                               "-1" = "#ff688c"),
                    labels = c("0" = "Decoding is correct", 
                               "1" = "Dead is decoded alive", 
                               "-1" = "Alive is decoded dead"))+
  labs(x = "individuals", y = "winters", fill = NULL) +
  theme_light()

### Average first, compute after

# Second option is to compute the posterior mean of the observation and 
# transition matrices, then to apply Viterbi:

# Initial states
delta <- c(1, 0)
# Transition matrix
transition <- matrix(NA, 2, 2)
transition[1,1] <- mean(phi)      # Pr(alive t -> alive t+1)
transition[1,2] <- 1 - mean(phi)  # Pr(alive t -> dead t+1)
transition[2,1] <- 0              # Pr(dead t -> alive t+1)
transition[2,2] <- 1              # Pr(dead t -> dead t+1)
# Observation matrix 
emission <- matrix(NA, 2, 2)
emission[1,1] <- 1 - mean(p)      # Pr(alive t -> non-detected t)
emission[1,2] <- mean(p)          # Pr(alive t -> detected t)
emission[2,1] <- 1                # Pr(dead t -> non-detected t)
emission[2,2] <- 0                # Pr(dead t -> detected t)
res <- matrix(NA, nrow = nrow(y), ncol = T)
for (i in 1:nrow(y)){
  res[i, 1:length(y[1,]) ] <- getViterbi(emission, transition, delta, y[i,] + 1)
}

# Again, you can compare the result of the Viterbi decoding to the actual 
# states we simulated and used to generate the observations:
df <- expand.grid(X = 1:nrow(z), Y = 1:ncol(z))
df$Z <- as_factor(c(z - res))
myPallette <- RColorBrewer::brewer.pal(name = "RdBu", n = 4)
df %>% 
  ggplot() + 
  aes(x = X, y = Y, fill = Z) + 
  geom_tile() +
  scale_fill_manual(values = c("0" = "#9597f0", 
                               "1" = "#d4b4f6", 
                               "-1" = "#ff688c"),
                    labels = c("0" = "Decoding is correct", 
                               "1" = "Dead is decoded alive", 
                               "-1" = "Alive is decoded dead"))+
  labs(x = "individuals", y = "winters", fill = NULL) +
  theme_light()

## Summary

## Suggested reading
