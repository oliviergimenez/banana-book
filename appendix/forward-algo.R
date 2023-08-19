
#-- model with constant survival and detection probabilities

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

# density
dHMM <- nimbleFunction(
  run = function(x = double(1), 
                 firstEncounter = double(0),
                 probInit = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = double(0, default = 0),
                 log = integer(0, default = 0)) {
    alpha <- probInit[1:2]
    for (t in (firstEncounter + 1):len) {
      if (firstEncounter == len) next
      alpha[1:2] <- (alpha[1:2] %*% probTrans[1:2,1:2]) * probObs[1:2,x[t]]
    }
    logL <- log(sum(alpha[1:2]))
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)


# random generator
rHMM <- nimbleFunction(
  run = function(n = integer(),
                 firstEncounter = double(0),
                 probInit = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = double(0, default = 0)) {
    returnType(double(1))
    z <- numeric(len)
    z[firstEncounter] <- rcat(n = 1, prob = probInit[1:2]) # all individuals alive at t = 0
    y <- z
    y[firstEncounter] <- 2 # all individuals are detected at t = 0
    for (t in (firstEncounter + 1):len){
      if (firstEncounter == len) next
      # state at t given state at t-1
      z[t] <- rcat(n = 1, prob = probTrans[z[t-1],1:2]) 
      # observation at t given state at t
      y[t] <- rcat(n = 1, prob = probObs[z[t],1:2]) 
    }
    return(y)
  })

# global env
assign('dHMM', dHMM, .GlobalEnv)
assign('rHMM', rHMM, .GlobalEnv)

# seems to work
dHMM(x = c(1,1,1,1,2),
     firstEncounter = 5,
     probInit = delta,
     probObs = Gamma,
     probTrans = Omega,
     len = 5)
rHMM(firstEncounter = 5,
     probInit = delta,
     probObs = Gamma,
     probTrans = Omega,
     len = 5)

# NIMBLE code 
hmm.phip <- nimbleCode({
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
    y[i,1:T] ~ dHMM(firstEncounter = first[i],
                    probInit = delta[1:2], 
                    probObs = omega[1:2,1:2], # observation matrix
                    probTrans = gamma[1:2,1:2], # transition matrix
                    len = T) # nb of sampling occasions
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
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1))
# parameters to monitor
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 2500
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
MCMCsummary(object = mcmc.phip, round = 2)
# caterpillar plot 
MCMCtrace(object = mcmc.phip, pdf = FALSE)



#-- model with time-varying survival and detection probabilities

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

#' We filter out individuals that are first captured at last occasion. 
#' These individuals do not contribute to parameter estimation, and also they cause problems with nimbleEcology.
#' 
#' # occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))

mask <- which(first!=ncol(y)) # individuals that are not first encountered at last occasion
y <- y[mask, ]                # keep only these
first <- first[mask]


# density
dHMM <- nimbleFunction(
  run = function(x = double(1), 
                 firstEncounter = double(0),
                 probInit = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 len = double(0, default = 0),
                 log = integer(0, default = 0)) {
    alpha <- probInit[1:2]
    for (t in (firstEncounter + 1):len) {
      if (firstEncounter == len) next
      alpha[1:2] <- (alpha[1:2] %*% probTrans[1:2,1:2,t-1]) * probObs[1:2,x[t],t-1]
    }
    logL <- log(sum(alpha[1:2]))
    returnType(double(0))
    if (log) return(logL)
    return(exp(logL))
  }
)

# random generator
rHMM <- nimbleFunction(
  run = function(n = integer(),
                 firstEncounter = double(0),
                 probInit = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 len = double(0, default = 0)) {
    returnType(double(1))
    z <- numeric(len)
    z[firstEncounter] <- rcat(n = 1, prob = probInit[1:2]) # all individuals alive at t = 0
    y <- z
    y[firstEncounter] <- 2 # all individuals are detected at t = 0
    for (t in (firstEncounter + 1):len){
      if (firstEncounter == len) next
      # state at t given state at t-1
      z[t] <- rcat(n = 1, prob = probTrans[z[t-1],1:2,t-1]) 
      # observation at t given state at t
      y[t] <- rcat(n = 1, prob = probObs[z[t],1:2,t-1]) 
    }
    return(y)
  })

# global env
assign('dHMM', dHMM, .GlobalEnv)
assign('rHMM', rHMM, .GlobalEnv)

# seems to work
dHMM(x = c(1,1,1,1,2),
     firstEncounter = 5,
     probInit = delta,
     probObs = array(Gamma, dim = c(2,2,5)),
     probTrans = array(Omega, dim = c(2,2,5)),
     len = 5)
rHMM(firstEncounter = 5,
     probInit = delta,
     probObs = array(Gamma, dim = c(2,2,5)),
     probTrans = array(Omega, dim = c(2,2,5)),
     len = 5)

# NIMBLE code 
hmm.phip <- nimbleCode({
  # likelihood
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0        # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1        # Pr(dead t -> dead t+1)
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1        # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0        # Pr(dead t -> detected t)
    phi[t] ~ dunif(0, 1) # prior survival
    p[t] ~ dunif(0, 1) # prior detection
  }
  
  for (i in 1:N){
    y[i,1:T] ~ dHMM(firstEncounter = first[i],
                    probInit = delta[1:2], 
                    probObs = omega[1:2,1:2,1:(T-1)], # observation matrix
                    probTrans = gamma[1:2,1:2,1:(T-1)], # transition matrix
                    len = T) # nb of sampling occasions
  }
})


# constants
my.constants <- list(N = nrow(y), 
                     T = ncol(y), 
                     first = first)
# data
my.data <- list(y = y + 1)
# initial values 
initial.values <- function() list(phi = runif(my.constants$T-1,0,1),
                                  p = runif(my.constants$T-1,0,1))


model <- nimbleModel(code = hmm.phip, 
            constants = my.constants,
            data = my.data,              
            inits = initial.values())


# parameters to monitor
parameters.to.save <- c("phi", "p")
# MCMC details
n.iter <- 2500
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
MCMCsummary(object = mcmc.phip, round = 2)
# caterpillar plot 
MCMCtrace(object = mcmc.phip, pdf = FALSE)

