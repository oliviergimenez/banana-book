library(nimble)

# simulate data

set.seed(2022) # for reproducibility

nocc <- 5 # nb of sampling occasions
nind <- 57 # nb of animals
p <- 0.6 # detection prob
phi <- 0.8 # survival prob

# Vector of initial states probabilities
delta <- c(1,0) # all individuals are alive in first occasion

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

# Matrix of states and observations
# z = 1 for alive, z = 2 for dead
# y = 1 for non-detected, y = 2 for detected
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
y
y <- y - 1 # non-detection = 0, detection = 1

# density
dHMM <- nimbleFunction(
  run = function(x = double(1), 
                 probInit = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = double(0, default = 0),
                 log = integer(0, default = 0)) {
    alpha <- probInit[1:2]
    for (t in 2:len) {
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

# global env
assign('dHMM', dHMM, .GlobalEnv)
assign('rHMM', rHMM, .GlobalEnv)

# seems to work
dHMM(x = c(2,2,1,1,1),
     probInit = delta,
     probObs = Gamma,
     probTrans = Omega,
     len = 5)
rHMM(probInit = delta,
     probObs = Gamma,
     probTrans = Omega,
     len = 5)

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
    y[i,1:T] ~ dHMM(probInit = delta[1:2], 
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

# run NIMBLE
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

# numerical summaries
library(MCMCvis)
MCMCsummary(mcmc.output)
