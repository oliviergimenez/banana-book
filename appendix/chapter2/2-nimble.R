library(nimble)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# NIMBLE tutorial {#intronimble}

## Introduction

## What is NIMBLE?

## Getting started {#start-nimble}

# To run NIMBLE, you will need to:  
# 1. Build a model consisting of a likelihood and priors.   
# 2. Read in some data.   
# 3. Specify parameters you want to make inference about.   
# 4. Pick initial values for parameters to be estimated (for each chain).   
# 5. Provide MCMC details namely the number of chains, the length of the burn-in period and the number of iterations following burn-in.

library(nimble)

# Now let's go back to our example on animal survival from the previous chapter. 
# First step is to build our model by specifying the binomial likelihood and a 
# uniform prior on survival probability `theta`. We use the `nimbleCode()` 
# function and wrap code within curly brackets:
model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
  # derived quantity
  lifespan <- -1/log(theta)
})

# You can check that the `model` R object contains your code:
model

# You can think of models in NIMBLE as graphs... Such graphs are 
# called directed acyclic graph or DAG.
mc <- nimbleModel(model, data = list(released = 57, survived = 19))
#mc$getVarNames()
#mc$getNodeNames()
#mc$getNodeNames(determOnly = TRUE)
#mc$getNodeNames(stochOnly = TRUE)
#mc$getNodeNames(dataOnly = TRUE)
#mc$getDependencies("theta")
mc$plotGraph()

# Second step in our workflow is to read in some data. 
# We use a list in which each component corresponds to a known quantity 
# in the model:
my.data <- list(released = 57, survived = 19)

# Third step is to tell NIMBLE which nodes in your model you would like to 
# keep track of, in other words the quantities you'd like to do inference about. 
# In our model we want survival `theta` and `lifespan`:
parameters.to.save <- c("theta", "lifespan")

# Fourth step is to specify initial values for all model parameters. 
init1 <- list(theta = 0.1)
init2 <- list(theta = 0.5)
init3 <- list(theta = 0.9)
initial.values <- list(init1, init2, init3)
initial.values

# Alternatively, you can write an R function that generates random initial values:
initial.values <- function() list(theta = runif(1,0,1))
initial.values()

# If you are using a function to generate random initial values, it's always a 
# good idea to set the seed in your code before you draw the initial values. 
my.seed <- 666
set.seed(my.seed)

# Fifth and last step, you need to tell NIMBLE the number of chains to run, 
# say `n.chain`, how long the burn-in period should be, say `n.burnin`, 
# and the number of iterations following the burn-in period to be used for
# posterior inference:
n.iter <- 5000
n.burnin <- 1000
n.chains <- 3

# We now have all the ingredients to run our model, that is to sample from the 
# posterior distribution of model parameters using MCMC simulations. 
# This is accomplished using function `nimbleMCMC()`: 
mcmc.output <- nimbleMCMC(code = model,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

# Now let's inspect what we have in `mcmc.output`: 
str(mcmc.output)

# The R object `mcmc.output` is a list with three components, 
# one for each MCMC chain. Let's have a look to `chain1` for example:
dim(mcmc.output$chain1)
head(mcmc.output$chain1)

# From there, you can compute the posterior mean of `theta`:
mean(mcmc.output$chain1[,'theta'])

# You can also obtain the 95% credible interval for `theta`:
quantile(mcmc.output$chain1[,'theta'], probs = c(2.5, 97.5)/100)

# Let's visualise the posterior distribution of `theta` with a histogram: 
mcmc.output$chain1[,"theta"] %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = value), color = "white") +
  labs(x = "survival probability")

# Let's load the package `MCMCvis`:
library(MCMCvis)

# To get the most common numerical summaries, the function `MCMCsummary()` does the job:
MCMCsummary(object = mcmc.output, round = 2)

# You can use a caterpillar plot to visualise the posterior distributions of 
# `theta` with `MCMCplot()`:
MCMCplot(object = mcmc.output, 
         params = 'theta')

# Visualization of a MCMC chain itself, i.e. the values of posterior samples 
# plotted against iteration number, is called a trace. The trace and posterior 
# density of theta can be obtained with `MCMCtrace()`:
MCMCtrace(object = mcmc.output,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "theta")

# You can also add the diagnostics of convergence we discussed in the previous chapter:
MCMCtrace(object = mcmc.output,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE, # add Rhat
          n.eff = TRUE, # add eff sample size
          params = "theta")

# In our example, all you need is samples  from the posterior distribution of 
# `theta`, which we pool between the three chains with:
theta_samples <- c(mcmc.output$chain1[,'theta'], 
                   mcmc.output$chain2[,'theta'],
                   mcmc.output$chain3[,'theta'])

# To get samples from the posterior distribution of lifespan, we apply the 
# function to calculate lifespan to the samples from the posterior distribution 
# of survival:
lifespan <- -1/log(theta_samples)

# As usual then, you can calculate the posterior mean and 95% credible interval:
mean(lifespan)
quantile(lifespan, probs = c(2.5, 97.5)/100)

# You can also visualise the posterior distribution of lifespan:
lifespan %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = value), color = "white") +
  labs(x = "lifespan")

# For convenience I have summarized the steps above in the box below. 
# The NIMBLE workflow provided with `nimbleMCMC()` allows you to build 
# models and make inference. This is what you can achieve with other 
# software like WinBUGS or JAGS. 
# **NIMBLE workflow:**

# model building
model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
  # derived quantity
  lifespan <- -1/log(theta)
})
# read in data
my.data <- list(released = 57, survived = 19)
# specify parameters to monitor
parameters.to.save <- c("theta", "lifespan")
# pick initial values
initial.values <- function() list(theta = runif(1,0,1))
# specify MCMC details
n.iter <- 5000
n.burnin <- 1000
n.chains <- 3
# run NIMBLE
mcmc.output <- nimbleMCMC(code = model,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
# calculate numerical summaries
MCMCsummary(object = mcmc.output, round = 2)
# visualize parameter posterior distribution
MCMCplot(object = mcmc.output, 
         params = 'theta')
# check convergence
MCMCtrace(object = mcmc.output,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "theta")

## Programming {#functions-in-nimble}

### NIMBLE functions

# NIMBLE provides `nimbleFunctions` for programming. A `nimbleFunction` is like 
# an R function, plus it can be compiled for faster computation. Going back to 
# our animal survival example, we can write a `nimbleFunction` to compute lifespan:
computeLifespan <- nimbleFunction(
    run = function(theta = double(0)) { # type declarations
        ans <- -1/log(theta)
        return(ans)
        returnType(double(0))  # return type declaration
    } )

# You can use your `nimbleFunction` in R:
computeLifespan(0.8)

# You can compile it and use the C++ code for faster computation: 
CcomputeLifespan <- compileNimble(computeLifespan)
CcomputeLifespan(0.8)

# You can also use your `nimbleFunction` in a model:
model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
  # derived quantity
  lifespan <- computeLifespan(theta)
})

# The rest of the workflow remains the same: 
my.data <- list(survived = 19, released = 57)
parameters.to.save <- c("theta", "lifespan")
initial.values <- function() list(theta = runif(1,0,1))
n.iter <- 5000
n.burnin <- 1000
n.chains <- 3
mcmc.output <- nimbleMCMC(code = model,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
MCMCsummary(object = mcmc.output, round = 2)

### Calling R/C++ functions {#callrfninnimble}

# As an example, imagine you'd like to use an R function `myfunction()`, 
# either a function you wrote yourself, or a function available 
# in your favorite R package:
myfunction <- function(x) {
  -1/log(x)
}

# Now wrap this function using `nimbleRcall()` or `nimbleExternalCall()` for a 
# C or C++ function:
Rmyfunction <- nimbleRcall(prototype = function(x = double(0)){}, 
                           Rfun = 'myfunction',
                           returnType = double(0))

# Now you can call your R function from a model (or any `nimbleFunctions`):
model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
  lifespan <- Rmyfunction(theta)
})

# The rest of the workflow remains the same: 
my.data <- list(survived = 19, released = 57)
parameters.to.save <- c("theta", "lifespan")
initial.values <- function() list(theta = runif(1,0,1))
n.iter <- 5000
n.burnin <- 1000
n.chains <- 3
mcmc.output <- nimbleMCMC(code = model,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
MCMCsummary(object = mcmc.output, round = 2)

### User-defined distributions

# With `nimbleFunctions` you can provide user-defined distributions to NIMBLE. 
# You need to write functions for density (`d`) and simulation (`r`) for your 
# distribution. As an example, we write our own binomial distribution:

# density
dmybinom <- nimbleFunction(
  run = function(x = double(0), 
                 size = double(0), 
                 prob = double(0), 
                 log = integer(0, default = 1)) {
    returnType(double(0))
    # compute binomial coefficient = size! / [x! (n-x)!] and take log
    lchoose <- lfactorial(size) - lfactorial(x) - lfactorial(size - x)
    # binomial density function = size! / [x! (n-x)!] * prob^x * (1-prob)^(size-x) and take log
    logProb <- lchoose + x * log(prob) + (size - x) * log(1 - prob)
    if(log) return(logProb)
    else return(exp(logProb)) 
  })
# simulation using the coin flip method (p. 524 in Devroye 1986)
# note: the n argument is required by NIMBLE but is not used, default is 1
rmybinom <- nimbleFunction(
  run = function(n = integer(0, default = 1),
                 size = double(0),
                 prob = double(0)) {
      x <- 0
      y <- runif(n = size, min = 0, max = 1)
      for (j in 1:size){
        if (y[j] < prob){
          x <- x + 1
        }else{
          x <- x
        }
      }
    returnType(double(0))
    return(x)    
  })

# You need to define the `nimbleFunctions` in R's global environment for them 
# to be accessed: 
assign('dmybinom', dmybinom, .GlobalEnv)
assign('rmybinom', rmybinom, .GlobalEnv)

# You can try out your function and simulate a single random value (n = 1 by default) from a binomial distribution with size 5 and probability 0.1: 
rmybinom(size = 5, prob = 0.1)

# All set. You can run your workflow:
model <- nimbleCode({
 # likelihood
 survived ~ dmybinom(prob = theta, size = released)
 # prior
 theta ~ dunif(0, 1)
})
my.data <- list(released = 57, survived = 19)
initial.values <- function() list(theta = runif(1,0,1))
n.iter <- 5000
n.burnin <- 1000
n.chains <- 3
mcmc.output <- nimbleMCMC(code = model,
                          data = my.data,
                          inits = initial.values,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
MCMCsummary(mcmc.output)

## Under the hood {#under-the-hood}

# We write the model code, read in data and pick initial values as before:
model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
  # derived quantity
  lifespan <- -1/log(theta)
})
my.data <- list(survived = 19, released = 57)
initial.values <- list(theta = 0.5)

# First step is to create the model as an R object (uncompiled model) with `nimbleModel()`:
survival <- nimbleModel(code = model,
                        data = my.data,
                        inits = initial.values)

# You can look at its nodes:
survival$getNodeNames()

# You can look at the values stored at each node:
survival$theta
survival$survived
survival$lifespan 
# this is -1/log(0.5)

# We can also calculate the log-likelihood at the initial value for `theta`:
survival$calculate()
# this is dbinom(x = 19, size = 57, prob = 0.5, log = TRUE)

# The ability in NIMBLE to access the nodes of your model and to evaluate 
# the model likelihood can help you in identifying bugs in your code. 
# For example, if we provide a negative initial value for `theta`, 
# `survival$calculate()` returns NA:
survival <- nimbleModel(code = model,
                        data = my.data,
                        inits = list(theta = -0.5))
survival$calculate()

# As another example, if we convey in the data the information that more 
# animals survived than were released, we'll get an infinity value for the 
# log-likelihood:
my.data <- list(survived = 61, released = 57)
initial.values <- list(theta = 0.5)
survival <- nimbleModel(code = model,
                        data = my.data,
                        inits = initial.values)
survival$calculate()

# As a check that the model is correctly initialized and that your code is 
# without bugs, the call to `model$calculate()` should return a number and 
# not NA or -Inf:
my.data <- list(survived = 19, released = 57)
initial.values <- list(theta = 0.5)
survival <- nimbleModel(code = model,
                        data = my.data,
                        inits = initial.values)
survival$calculate()

# You can obtain the graph of the model as in Figure \@ref(fig:dag-survival) with:
survival$plotGraph()

# Second we compile the model with `compileNimble()`:
Csurvival <- compileNimble(survival)

# With `compileNimble()`, the C++ code is generated, compiled and loaded 
# back into R so that it can be used in R (compiled model):
Csurvival$theta

# Now you have two versions of the model, `survival` is in R and `Csurvival` 
# in C++. Being able to separate the steps of model building and parameter 
# estimation is a strength of NIMBLE. This gives you a lot of flexibility at 
# both steps. For example, imagine you would like to fit your model with maximum 
# likelihood, then you can do it by wrapping your model in an R function that 
# gets the likelihood and maximise this function. Using the C version of the 
# model, you can write:

# function for negative log-likelihood to minimize
f <- function(par) {
    Csurvival[['theta']] <- par # assign par to theta 
    ll <- Csurvival$calculate() # update log-likelihood with par value
    return(-ll) # return negative log-likelihood
}
# evaluate function at 0.5 and 0.9
f(0.5)
f(0.9)
# minimize function
out <- optimize(f, interval = c(0,1))
round(out$minimum, 2)

# By maximising the likelihood (or minimising the negative log-likelihood), 
# you obtain the maximum likelihood estimate of animal survival, which is 
# exactly 19 surviving animals over 57 released animals or `r round(19/57, 2)`.

# Third we create a MCMC configuration for our model with `configureMCMC()`:
survivalConf <- configureMCMC(survival)

# To monitor `lifespan` in addition to `theta`, you write:
survivalConf$addMonitors(c("lifespan"))
survivalConf

# Third, we create a MCMC function with `buildMCMC()` and compile it 
# with `compileNimble()`:
survivalMCMC <- buildMCMC(survivalConf)
CsurvivalMCMC <- compileNimble(survivalMCMC, project = survival)

# Fourth, we run NIMBLE with `runMCMC()`:
n.iter <- 5000
n.burnin <- 1000
samples <- runMCMC(mcmc = CsurvivalMCMC, 
                   niter = n.iter,
                   nburnin = n.burnin)

# You can look into `samples` which contains values simulated from the posterior 
# distribution of the parameters we monitor:
head(samples)

# From here, you can obtain numerical summaries with `samplesSummary()` 
# (or `MCMCvis::MCMCsummary()`):
samplesSummary(samples)

# I have summarized the steps above in the box below. 
# **Detailed NIMBLE workflow:**

# model building
model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
  # derived quantity
  lifespan <- -1/log(theta)
})
# read in data
my.data <- list(released = 57, survived = 19)
# pick initial values
initial.values <- function() list(theta = runif(1,0,1))
# create model as an R object (uncompiled model)
survival <- nimbleModel(code = model,
                        data = my.data,
                        inits = initial.values())
# compile model
Csurvival <- compileNimble(survival)
# create a MCMC configuration
survivalConf <- configureMCMC(survival)
# add lifespan to list of parameters to monitor
survivalConf$addMonitors(c("lifespan"))
# create a MCMC function and compile it
survivalMCMC <- buildMCMC(survivalConf)
CsurvivalMCMC <- compileNimble(survivalMCMC, project = survival)
# specify MCMC details
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2
# run NIMBLE
samples <- runMCMC(mcmc = CsurvivalMCMC, 
                   niter = n.iter,
                   nburnin = n.burnin,
                   nchain = n.chains)
# calculate numerical summaries
MCMCsummary(object = samples, round = 2)
# visualize parameter posterior distribution
MCMCplot(object = samples, 
         params = 'theta')
# check convergence
MCMCtrace(object = samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "theta")

## MCMC samplers

### Default samplers {#change-sampler}

# What is the default sampler used by NIMBLE in our example? 
survivalConf$printSamplers()

# Now that we have control on the MCMC configuration, let's mess it up. 
# We start by removing the default sampler:
survivalConf$removeSamplers(c('theta'))
survivalConf$printSamplers()

# And we change it for a slice sampler:
survivalConf$addSampler(target = c('theta'),
                        type = 'slice')
survivalConf$printSamplers()

# Now you can resume the workflow:

# create a new MCMC function and compile it:
survivalMCMC2 <- buildMCMC(survivalConf)
CsurvivalMCMC2 <- compileNimble(survivalMCMC2, 
                                project = survival,
                                resetFunctions = TRUE) # to compile new functions 
                                                       # into existing project, 
                                                       # need to reset nimbleFunctions
# run NIMBLE:
samples2 <- runMCMC(mcmc = CsurvivalMCMC2, 
                    niter = n.iter,
                    nburnin = n.burnin)
# obtain numerical summaries:
samplesSummary(samples2)

# NIMBLE implements many samplers, and a list is available with 
?samplers

### User-defined samplers

# Allowing you to code your own sampler is another topic on which NIMBLE thrives. 
# As an example, we focus on the Metropolis algorithm of Section 
# \@ref(metropolis-algorithm) which we coded in R. In this section, we make it 
# a `nimbleFunction` so that we can use it within our model: 
my_metropolis <- nimbleFunction(
  name = 'my_metropolis', # fancy name for our MCMC sampler
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # i) get dependencies for 'target' in 'model'
    calcNodes <- model$getDependencies(target) 
    # ii) get sd of proposal distribution
    scale <- control$scale 
  },
  run = function() {
    # (1) log-lik at current value
    initialLP <- model$getLogProb(calcNodes) 
    # (2) current parameter value
    current <- model[[target]] 
    # (3) logit transform
    lcurrent <- log(current / (1 - current))
    # (4) propose candidate value
    lproposal <- lcurrent  + rnorm(1, mean = 0, scale) 
    # (5) back-transform
    proposal <- plogis(lproposal)
    # (6) plug candidate value in model 
    model[[target]] <<- proposal 
    # (7) log-lik at candidate value
    proposalLP <- model$calculate(calcNodes)
    # (8) compute lik ratio on log scale
    lMHR <- proposalLP - initialLP 
    # (9) spin continuous spinner and compare to ratio
    if(runif(1,0,1) < exp(lMHR)) { 
      # (10) if candidate value is accepted, update current value
      copy(from = model, to = mvSaved, nodes = calcNodes, logProb = TRUE, row = 1)
    } else {
      ## (11) if candidate value is accepted, keep current value
      copy(from = mvSaved, to = model, nodes = calcNodes, logProb = TRUE, row = 1)
    }
  },
  methods = list(
    reset = function() {}
  )
)

# Now that we have our user-defined MCMC algorithm, we can change the default 
# sampler for our new sampler as in Section \@ref(change-sampler). We start 
# from scratch:
model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
})
my.data <- list(survived = 19, released = 57)
initial.values <- function() list(theta = runif(1,0,1))
survival <- nimbleModel(code = model, 
                        data = my.data, 
                        inits = initial.values())
Csurvival <- compileNimble(survival)
survivalConf <- configureMCMC(survival)

# We print the samplers used by default, remove the default sampler for 
# `theta`, replace it with our `my_metropolis()` sampler with the standard 
# deviation of the proposal distribution set to 0.1, and print again to make 
# sure NIMBLE now uses our new sampler:
survivalConf$printSamplers()
survivalConf$removeSamplers(c('theta'))
survivalConf$addSampler(target = 'theta', 
                        type = 'my_metropolis', 
                        control = list(scale = 0.1)) # standard deviation
                                                     # of proposal distribution
survivalConf$printSamplers()

# The rest of the workflow is unchanged:
survivalMCMC <- buildMCMC(survivalConf)
CsurvivalMCMC <- compileNimble(survivalMCMC, 
                               project = survival)
samples <- runMCMC(mcmc = CsurvivalMCMC, 
                   niter = 5000, 
                   nburnin = 1000)
samplesSummary(samples)

# You can re-run the analysis by setting the standard deviation of the proposal 
# to different values, say 1 and 10, and compare the results to traceplots we 
# obtained with our R implementation of the Metropolis algorithm in the previous 
# chapter:

# standard deviation of proposal is 0.1
scale <- 0.1
Rmodel <- nimbleModel(code = model, data = my.data, inits = initial.values())
conf <- configureMCMC(Rmodel, monitors = c('theta'), print = FALSE)
conf$removeSamplers(c('theta'))
conf$addSampler(target = 'theta', type = 'my_metropolis', control = list(scale = scale))
Rmcmc <- buildMCMC(conf)
out <- compileNimble(list(model = Rmodel, mcmc = Rmcmc))
Cmcmc <- out$mcmc
samples_sd01 <- runMCMC(Cmcmc, niter = 10000, nburnin = 9000, progressBar = FALSE)
# standard deviation of proposal is 1
scale <- 1
Rmodel <- nimbleModel(code = model, data = my.data, inits = initial.values())
conf <- configureMCMC(Rmodel, monitors = c('theta'), print = FALSE)
conf$removeSamplers(c('theta'))
conf$addSampler(target = 'theta', type = 'my_metropolis', control = list(scale = scale))
Rmcmc <- buildMCMC(conf)
out <- compileNimble(list(model = Rmodel, mcmc = Rmcmc))
Cmcmc <- out$mcmc
samples_sd1 <- runMCMC(Cmcmc, niter = 10000, nburnin = 9000, progressBar = FALSE)
# standard deviation of proposal is 10
scale <- 10
Rmodel <- nimbleModel(code = model, data = my.data, inits = initial.values())
conf <- configureMCMC(Rmodel, monitors = c('theta'), print = FALSE)
conf$removeSamplers(c('theta'))
conf$addSampler(target = 'theta', type = 'my_metropolis', control = list(scale = scale))
Rmcmc <- buildMCMC(conf)
out <- compileNimble(list(model = Rmodel, mcmc = Rmcmc))
Cmcmc <- out$mcmc
samples_sd10 <- runMCMC(Cmcmc, niter = 10000, nburnin = 9000, progressBar = FALSE)
# trace plot for scenario with standard deviation 0.1
plot01 <- samples_sd01 %>%
  as_tibble() %>%
  ggplot() + 
  aes(x = 9001:10000, y = theta) +
  geom_line() + 
  labs(x = "iterations", title = "scale = 0.1")
# trace plot for scenario with standard deviation 1
plot1 <- samples_sd1 %>%
  as_tibble() %>%
  ggplot() + 
  aes(x = 9001:10000, y = theta) +
  geom_line() + 
  labs(x = "iterations", title = "scale = 1")
# trace plot for scenario with standard deviation 10
plot10 <- samples_sd10 %>%
  as_tibble() %>%
  ggplot() + 
  aes(x = 9001:10000, y = theta) +
  geom_line() + 
  labs(x = "iterations", title = "scale = 10")
# Assemble all three trace plots
library(patchwork)
plot01 + plot1 + plot10

## Tips and tricks

### Precision vs standard deviation

### Indexing

### Faster compilation

model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
})
my.data <- list(survived = 19, released = 57)
initial.values <- function() list(theta = runif(1,0,1))
survival <- nimbleModel(code = model, 
                        data = my.data, 
                        inits = initial.values(),
                        calculate = FALSE) # first tip
Csurvival <- compileNimble(survival)
survivalConf <- configureMCMC(survival)
survivalMCMC <- buildMCMC(survivalConf, useConjugacy = FALSE) # second tip
CsurvivalMCMC <- compileNimble(survivalMCMC, 
                               project = survival)
samples <- runMCMC(mcmc = CsurvivalMCMC, 
                   niter = 5000, 
                   nburnin = 1000)
samplesSummary(samples)

### Updating MCMC chains

# Sometimes it is useful to run your MCMC chains a little bit longer to 
# improve convergence. Re-starting from the run in previous section, you can use:
niter_ad <- 6000
CsurvivalMCMC$run(niter_ad, reset = FALSE)

# Then you can extract the matrix of previous MCMC samples augmented with new 
# ones and obtain numerical summaries:
more_samples <- as.matrix(CsurvivalMCMC$mvSamples)
samplesSummary(more_samples)

### Reproducibility {#tipreproducibility}

model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
})
my.data <- list(survived = 19, released = 57)
initial.values <- list(list(theta = 0.1),
                       list(theta = 0.5),
                       list(theta = 0.9))

# If you want your results to be reproducible, you can control the state of 
# R the random number generator with the `setSeed` argument in functions 
# `nimbleMCMC()` and `runMCMC()`. Going back to the animal survival example, 
# you can check that two calls to `nimbleMCMC()` give the same results when 
# `setSeed` is set to the same value. Note that we need to specify a seed 
# for each chain, hence a vector of three components here: 

# first call to nimbleMCMC()
mcmc.output1 <- nimbleMCMC(code = model,
                           data = my.data,
                           inits = initial.values,
                           niter = 5000,
                           nburnin = 1000,
                           nchains = 3,
                           summary = TRUE,
                           setSeed = c(2024, 2025, 2026))

# second call to nimbleMCMC()
mcmc.output2 <- nimbleMCMC(code = model,
                           data = my.data,
                           inits = initial.values,
                           niter = 5000,
                           nburnin = 1000,
                           nchains = 3,
                           summary = TRUE,
                           setSeed = c(2024, 2025, 2026))

# outputs from both calls are the same
mcmc.output1$summary$all.chains
mcmc.output2$summary$all.chains

# Note that to make your workflow reproducible, you need to set the seed 
# not only within the `nimbleMCMC()` call, but also before setting your 
# initial values if you are using a randomized function for that. 

### Parallelization

library(parallel)

# First you create a cluster using the total amount of cores you have but 
# one to make sure your computer can go on working:

nbcores <- detectCores() - 1
my_cluster <- makeCluster(nbcores)

# Then you wrap your workflow in a function to be run in parallel:

workflow <- function(seed, data) {
  
  library(nimble)
  
  model <- nimbleCode({
    # likelihood
    survived ~ dbinom(theta, released)
    # prior
    theta ~ dunif(0, 1)
  })
  
  set.seed(123) # for reproducibility
  initial.values <- function() list(theta = runif(1,0,1))
  
  survival <- nimbleModel(code = model, 
                          data = data, 
                          inits = initial.values())
  Csurvival <- compileNimble(survival)
  survivalMCMC <- buildMCMC(Csurvival)
  CsurvivalMCMC <- compileNimble(survivalMCMC)
  
  samples <- runMCMC(mcmc = CsurvivalMCMC, 
                     niter = 5000, 
                     nburnin = 1000,
                     setSeed = seed)
  
  return(samples)
}

# Now we run the code using `parLapply()`, which uses cluster nodes to execute our workflow:
output <- parLapply(cl = my_cluster, 
                    X = c(2022, 666),
                    fun = workflow, 
                    data = list(survived = 19, released = 57))

stopCluster(my_cluster)

str(output)

MCMCsummary(output)

### Incomplete and incorrect initialization

model <- nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
})
#initial.values <- list(theta = runif(1,0,1))
survival <- nimbleModel(code = model, 
                        data = list(survived = 19, released = 57))

survival$calculate() # gives NA
survival$initializeInfo()

survival$theta <- 0.5 # assign initial value to theta
survival$calculate() 

Csurvival <- compileNimble(survival)
survivalMCMC <- buildMCMC(Csurvival)
CsurvivalMCMC <- compileNimble(survivalMCMC)

samples <- runMCMC(mcmc = CsurvivalMCMC, 
                   niter = 5000, 
                   nburnin = 1000)

samplesSummary(samples)

### Vectorization

## Summary

## Suggested reading
