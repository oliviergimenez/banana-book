
library(tidyverse)
library(nimble)

######---------- MEMORY
geese <- read_csv(here::here("dat", "geese.csv"), col_names = TRUE)
#geese <- read_csv2(here::here("slides", "dat", "allgeese.csv"), col_names = TRUE) %>%
#  uncount(weights = `158`)

y <- geese %>%
  as.matrix()
y2 <- y
y2[y2==3] <- 0
mask <- apply(y2, 1, sum)
y <- y2[mask!=0,]
dim(y)
first <- apply(y, 1, function(x) min(which(x != 0)))

# memory model

# To fit a model in which the survival-transition probabilities 
# depend of the two sites previously occupied we have to consider 
# a specific set of states composed of all the couples of sites 
# (sites occupied at times t-1 and t) plus the state “dead”. 


# Obs = non-detection, detected in site 1, detected in site 2 = 1, 2, 3
# States = couples of sites occupied: { 11, 12, 21, 22, dead } = {1 2 3 4 5}

memory <- nimbleCode({
  
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
  pi[1:4] ~ ddirch(beta[1:4]) # pi11, pi12, pi21, pi22
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
  
  delta[1] <- pi[1]
  delta[2] <- pi[2]
  delta[3] <- pi[3]
  delta[4] <- pi[4]
  delta[5] <- 0
  
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
  
  omega.init[1,1] <- 0
  omega.init[1,2] <- 1
  omega.init[1,3] <- 0
  omega.init[2,1] <- 0
  omega.init[2,2] <- 0
  omega.init[2,3] <- 1
  omega.init[3,1] <- 0
  omega.init[3,2] <- 1
  omega.init[3,3] <- 0
  omega.init[4,1] <- 0
  omega.init[4,2] <- 0
  omega.init[4,3] <- 1
  omega.init[5,1] <- 1
  omega.init[5,2] <- 0
  omega.init[5,3] <- 0
  
  # Likelihood 
  for (i in 1:N){
    # Define latent state at first capture
    z[i,first[i]] ~ dcat(delta[1:5])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]],1:3])
    for (t in (first[i]+1):K){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:5])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(omega[z[i,t],1:3])
    }
  }
})


# initial values for unknown z 
# not optimal since it ignores impossible transitions in the surv-mov matrix
zinit <- y
for (i in 1:nrow(y)) {
  for (j in 1:ncol(y)) {
    if (j > first[i] & y[i,j]==0) {zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j > first[i] & y[i,j]==1) {zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j > first[i] & y[i,j]==2) {zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
    #    if (j > first[i] & y[i,j]==1) {zinit[i,j] <- which(rmultinom(1, 1, c(1/2,0,1/2,0))==1)}
    #    if (j > first[i] & y[i,j]==2) {zinit[i,j] <- which(rmultinom(1, 1, c(0,1/2,0,1/2))==1)}
    if (j < first[i]) {zinit[i,j] <- 0}
  }
}
zinit <- as.matrix(zinit)

# constants
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

# Initial values
initial.values <- function(){list(phi11 = rdirch(1, rep(1,3)), 
                                  phi12 = rdirch(1, rep(1,3)), 
                                  phi21 = rdirch(1, rep(1,3)), 
                                  phi22 = rdirch(1, rep(1,3)),
                                  det1 = runif(1, 0, 1), 
                                  det2 = runif(1, 0, 1), 
                                  pi = rdirch(1, rep(1,4)),
                                  z = zinit)}  

# data
my.data <- list(y = y + 1,
                alpha = rep(1,3),
                beta = rep(1,4))

# MCMC settings
parameters.to.save <- c("phi11", "phi12","phi21", "phi22", "det1", "det2", "pi")
n.iter <- 2500
n.burnin <- 500
n.chains <- 2

mcmc.memory <- nimbleMCMC(code = memory,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

save(mcmc.memory, file = "memory.RData")
library(MCMCvis)
MCMCsummary(mcmc.memory, round = 2)
MCMCplot(mcmc.memory)
MCMCtrace(mcmc.memory, params = c("det1","det2"), pdf = FALSE)

# NE MARCHE PLUS
# [Note] Infinite values were detected in model variables: logProb_z, logProb_y.
# SUREMENT ECNORE DES PBS DANS INITIALISATION DE Z

## PROBLEM: p1 and p2 seem to be updated, but Rhat  = Inf and neff = 0; inspecting their posterior distribution
## it feels like they are not updated, two modes corresp to the two chains. 
# Usually, it means that these parameters are not used in the model. 



#' 
#' ## Marginalization with nimbleEcology
#' 
#' Let's get back to the analysis of the Canada geese data, with 2 sites. 
## ------------------------------------------------------------------------

library(nimbleEcology)

geese <- read_csv(here::here("dat", "geese.csv"), col_names = TRUE)
#geese <- read_csv2(here::here("slides", "dat", "allgeese.csv"), col_names = TRUE) %>%
#  uncount(weights = `158`)

y <- geese %>%
  as.matrix()
y2 <- y
y2[y2==3] <- 0
mask <- apply(y2, 1, sum)
y <- y2[mask!=0,]

# Get the occasion of first capture for each individual.
## ------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

#' 
#' We filter out individuals that are first captured at last occasion. These individuals do not contribute to parameter estimation, and also they cause problems with nimbleEcology. 
## ------------------------------------------------------------------------
mask <- which(first!=ncol(y)) # individuals that are not first encountered at last occasion
y <- y[mask, ]                # keep only these
first <- first[mask]

#' 
#' Let's write the model. Note that the likelihood is simpler due to the use of the function `dHMM`.
## ------------------------------------------------------------------------
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

#' 
#' Data in a list. 
## ------------------------------------------------------------------------
my.data <- list(y = y + 1, 
                alpha = c(1, 1, 1))

#' 
#' Constants in a list. 
## ------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

#' 
#' Initial values. Note that we do not need initial values for the latent states anymore. 
## ------------------------------------------------------------------------
initial.values <- function(){list(phi11 = rdirch(1, rep(1,3)), 
                                  phi12 = rdirch(1, rep(1,3)), 
                                  phi21 = rdirch(1, rep(1,3)), 
                                  phi22 = rdirch(1, rep(1,3)),
                                  det1 = runif(1, 0, 1), 
                                  det2 = runif(1, 0, 1))}  

#' 
#' Parameters to monitor.
## ------------------------------------------------------------------------
parameters.to.save <- c("phi11", "phi12", "phi21", "phi22", "det1", "det2")

#' 
#' MCMC settings.
## ------------------------------------------------------------------------
n.iter <- 2500
n.burnin <- 500
n.chains <- 2

#' 
#' Run nimble.
## ----eval = FALSE--------------------------------------------------------
out <- nimbleMCMC(code = multisite.marginalized, 
                  constants = my.constants,
                  data = my.data,              
                  inits = initial.values(),
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin, 
                  nchains = n.chains)

save(out, file = "memory.RData")

library(MCMCvis)
MCMCsummary(out, round = 2)
MCMCplot(out)
MCMCtrace(out, params = c("det1","det2"), pdf = FALSE)


#################### full dataset?

# geese <- read_csv2(here::here("dat", "allgeese.csv"), col_names = FALSE)
# 
# y <- geese %>%
#   as.matrix()
# 
# # Separate sequence and count
# sequence <- y[, 1:(ncol(y)-1)]
# counts <- y[, ncol(y)]
# 
# # delete site 3
# y2 <- sequence
# y2[y2==3] <- 0
# mask <- apply(y2, 1, sum)
# y <- y2[mask!=0,]
# n <- counts[mask!=0]
# 
# dim(y)
# length(n)
# 
# # Expand rows according to counts
# expanded <- y[rep(1:nrow(y), n), ]
# 
# # Convert to data.frame (optional)
# expanded_df <- as.data.frame(expanded)
# 
# # Write to file (optional)
# write.table(expanded_df, file = "all-geese-2sites.txt", row.names = FALSE, col.names = FALSE)
# 



library(nimbleEcology)

y <- as.matrix(read.table(here::here("dat", "all-geese-2sites.txt")))

# Get the occasion of first capture for each individual.
## ------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

#' 
#' We filter out individuals that are first captured at last occasion. These individuals do not contribute to parameter estimation, and also they cause problems with nimbleEcology. 
## ------------------------------------------------------------------------
mask <- which(first!=ncol(y)) # individuals that are not first encountered at last occasion
y <- y[mask, ]                # keep only these
first <- first[mask]

#' 
#' Let's write the model. Note that the likelihood is simpler due to the use of the function `dHMM`.
## ------------------------------------------------------------------------
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

#' 
#' Data in a list. 
## ------------------------------------------------------------------------
my.data <- list(y = y + 1, 
                alpha = c(1, 1, 1))

#' 
#' Constants in a list. 
## ------------------------------------------------------------------------
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

#' 
#' Initial values. Note that we do not need initial values for the latent states anymore. 
## ------------------------------------------------------------------------
initial.values <- function(){list(phi11 = rdirch(1, rep(1,3)), 
                                  phi12 = rdirch(1, rep(1,3)), 
                                  phi21 = rdirch(1, rep(1,3)), 
                                  phi22 = rdirch(1, rep(1,3)),
                                  det1 = runif(1, 0, 1), 
                                  det2 = runif(1, 0, 1))}  

#' 
#' Parameters to monitor.
## ------------------------------------------------------------------------
parameters.to.save <- c("phi11", "phi12", "phi21", "phi22", "det1", "det2")

#' 
#' MCMC settings.
## ------------------------------------------------------------------------
n.iter <- 2500
n.burnin <- 500
n.chains <- 2

#' 
#' Run nimble.
## ----eval = FALSE--------------------------------------------------------
out <- nimbleMCMC(code = multisite.marginalized, 
                  constants = my.constants,
                  data = my.data,              
                  inits = initial.values(),
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin, 
                  nchains = n.chains)

save(out, file = "memory.RData")

library(MCMCvis)
MCMCsummary(out, round = 2)
MCMCplot(out)
MCMCtrace(out, params = c("det1","det2"), pdf = FALSE)

