library(tidyverse)
theme_set(theme_light(base_size = 14))




###### titis example with no uncertainty
# We use (historical) data collected between 1940 and 1957 by Lance Richdale on Sooty shearwaters Puffinus griseus (titis hereafter). 
# These data were reanalyzed using multistate capture-recapture models by Scofield et al. (2001) who kindly provided us with the data. 
# Following the way the data were collected, four states were originally considered:
# (1) Breeder (SB);
# (2) Accompanied by another bird in a burrow;
# (3) Alone in a burrow;
# (4) On the surface.
# Because of numerical issues, we pooled 2-3-4 together in a Non-Breeder state (NSB) that includes 
# failed breeders (birds that had bred previously – skip reproduction or divorce) and pre-breeders (birds that had yet to breed). 
# Note that because burrows were not checked before hatching, some birds in the category NSB might have already failed. 
# We therefore regard those birds in the SB state as successful breeders, and those in the NSB state as nonbreeders + prebreeders and failed breeders.
# We artificially generated uncertainty on the assignment of the two states Non-Successful-Breeder and Successful-Breeder as follows. 
# For each individual at each detection occasion: 
# 1) a Successful-Breeder was ascertained as Successful-Breeder with probability 0.2 
# and assigned to obs 1, and not ascertained with the complementary probability 0.8 and assigned to obs 3, while 
# 2) a Non-Successful-Breeder was ascertained as Non-Successful-Breeder with probability 0.7 
# and assigned to obs 2, and not ascertained with the complementary probability 0.3 and assigned to obs 3.
# This procedure was implemented in R, see the script make-uncertain.R in the codes/ directory
# Try and fit a model with uncertainty in the state assignment. Consider 2 probabilities to ascertain the breeding status of an individual 
# encountered as Non-Successful-Breeder (deltaNB) or Successful-Breeder (deltaB).

# To artificially generate uncertainty on both the states non-breeder and breeder, 
# we used the R script below to alter the raw capture-recapture data from file titis.inp; the resulting file is titis2.inp.

# read in data
titi <- read.table('dat/titis.txt')
# 1 seen as breeder
# 2 seen as non-breeder
# 0 not seen

# nb of capture occasions
ny <- ncol(titi)
# nb of individuals
nind <- nrow(titi)

titi2 <- titi
for (i in 1:nind)
{
  for (j in 1:ny){
    # 1 seen and ascertained Breeder (with probability .2)
    # 2 seen and ascertained Non-Breeder (with probability .7)
    # 3 seen but not ascertained (Non-Breeders with probability .8 + Breeders with probability .3)
    # 0 not seen
    
    # Breeders are ascertained with probability .2
    if (titi[i,j] == 1)
    {
      temp <- rbinom(1,size=1,prob=.2) 
      if (temp == 1) titi2[i,j] <- 1 # if ascertained B, event = 1
      if (temp == 0) titi2[i,j] <- 3 # if not ascertained, event = 3
    }
    
    # Non-Breeders are ascertained with probability .7 (event = 1), 
    # or not ascertained with probability .3 (event = 2)
    if (titi[i,j] == 2) 
    {
      temp <- rbinom(1,size=1,prob=.7)
      if (temp == 1) titi2[i,j] <- 2 # if ascertained NB, event = 2
      if (temp == 0) titi2[i,j] <- 3 # if not ascertained, event = 3
    }
    
  }
}

# ->>> titis_with_uncertainty.csv

# 1 seen as breeder
# 2 seen as non-breeder
# 0 not seen

titis <- read_csv2(here::here("slides", "dat", "titis_with_uncertainty.csv"), col_names = FALSE)
y <- titis %>%
  as.matrix()


library(nimble)
multievent <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phiB: survival probability state B
  # phiNB: survival probability state NB
  # psiBNB: transition probability from B to NB
  # psiNBB: transition probability from NB to B
  # pB: recapture probability B
  # pNB: recapture probability NB
  # piB prob. of being in initial state breeder
  # betaNB prob to ascertain the breeding status of an individual encountered as non-breeder
  # betaB prob to ascertain the breeding status of an individual encountered as breeder
  # -------------------------------------------------
  # States (S):
  # 1 alive B
  # 2 alive NB
  # 3 dead
  # Observations (O):  
  # 1 = non-detected
  # 2 = seen and ascertained as breeder
  # 3 = seen and ascertained as non-breeder
  # 4 = not ascertained
  # -------------------------------------------------
  
  # Priors
  phiB ~ dunif(0, 1)
  phiNB ~ dunif(0, 1)
  psiBNB ~ dunif(0, 1)
  psiNBB ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  pNB ~ dunif(0, 1)
  piB ~ dunif(0, 1)
  betaNB ~ dunif(0, 1)
  betaB ~ dunif(0, 1)
  
  # vector of initial stats probs
  delta[1] <- piB # prob. of being in initial state B
  delta[2] <- 1 - piB # prob. of being in initial state NB
  delta[3] <- 0 # prob. of being in initial state dead
  
  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phiB * (1 - psiBNB)
  gamma[1,2] <- phiB * psiBNB
  gamma[1,3] <- 1 - phiB
  gamma[2,1] <- phiNB * psiNBB
  gamma[2,2] <- phiNB * (1 - psiNBB)
  gamma[2,3] <- 1 - phiNB
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- 1 - pB             # Pr(alive B t -> non-detected t)
  omega[1,2] <- pB * betaB         # Pr(alive B t -> detected B t)
  omega[1,3] <- 0                  # Pr(alive B t -> detected NB t)
  omega[1,4] <- pB * (1 - betaB)   # Pr(alive B t -> detected U t)
  omega[2,1] <- 1 - pNB            # Pr(alive NB t -> non-detected t)
  omega[2,2] <- 0                  # Pr(alive NB t -> detected B t)
  omega[2,3] <- pNB * betaNB       # Pr(alive NB t -> detected NB t)
  omega[2,4] <- pNB * (1 - betaNB) # Pr(alive NB t -> detected U t)
  omega[3,1] <- 1                  # Pr(dead t -> non-detected t)
  omega[3,2] <- 0                  # Pr(dead t -> detected N t)
  omega[3,3] <- 0                  # Pr(dead t -> detected NB t)
  omega[3,4] <- 0                  # Pr(dead t -> detected U t)
  
  omega.init[1,1] <- 0                  # Pr(alive B t = 1 -> non-detected t = 1)
  omega.init[1,2] <- betaB              # Pr(alive B t = 1 -> detected B t = 1)
  omega.init[1,3] <- 0                  # Pr(alive B t = 1 -> detected NB t = 1)
  omega.init[1,4] <- 1 - betaB          # Pr(alive B t = 1 -> detected U t = 1)
  omega.init[2,1] <- 0                  # Pr(alive NB t = 1 -> non-detected t = 1)
  omega.init[2,2] <- 0                  # Pr(alive NB t = 1 -> detected B t = 1)
  omega.init[2,3] <- betaNB             # Pr(alive NB t = 1 -> detected NB t = 1)
  omega.init[2,4] <- 1 - betaNB         # Pr(alive NB t = 1 -> detected U t = 1)
  omega.init[3,1] <- 1                  # Pr(dead t = 1 -> non-detected t = 1)
  omega.init[3,2] <- 0                  # Pr(dead t = 1 -> detected N t = 1)
  omega.init[3,3] <- 0                  # Pr(dead t = 1 -> detected NB t = 1)
  omega.init[3,4] <- 0                  # Pr(dead t = 1 -> detected U t = 1)
  
  # Likelihood 
  for (i in 1:N){
    # Define latent state at first capture
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]],1:4])
    for (t in (first[i]+1):K){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:3])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
})

# constants
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

# data
my.data <- list(y = y + 1)

# initial values for unknown z
zinit <- y
zinit[zinit==3] <- sample(c(1,2), sum(zinit==3), replace = TRUE)
for (i in 1:nrow(y)) {
  for (j in 1:ncol(y)) {
    if (j > first[i] & y[i,j]==0) {zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j < first[i]) {zinit[i,j] <- 0}
  }
}
zinit <- as.matrix(zinit)


# Initial values
initial.values <- function(){list(phiNB = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  psiNBB = runif(1, 0, 1), 
                                  psiBNB = runif(1, 0, 1), 
                                  pNB = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1),
                                  piB = runif(1, 0, 1),
                                  betaB = runif(1, 0, 1),
                                  betaNB = runif(1, 0, 1),
                                  z = zinit)}  

# MCMC settings
parameters.to.save <- c("phiB", 
                        "phiNB", 
                        "psiNBB", 
                        "psiBNB", 
                        "piB", 
                        "pB", 
                        "pNB", 
                        "betaNB", 
                        "betaB")

n.iter <- 5000
n.burnin <- 2500
n.chains <- 2

mcmc.multievent <- nimbleMCMC(code = multievent, 
                              constants = my.constants,
                              data = my.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin, 
                              nchains = n.chains)

save(mcmc.multievent, file = "titisuncertain.RData")
library(MCMCvis)
MCMCsummary(mcmc.multievent, round = 2)
MCMCplot(mcmc.multievent)






############# DISEASE

# OBS
# 0 = not seen
# 1 = captured healthy
# 2 = captured ill
# 3 = health status is unknown (i.e. seen at distance

# STATES
# 1 alive healthy
# 2 alive ill
# 3 dead

y <- as.matrix(read.table(here::here("dat", "hofi_2000.txt")))

library(nimble)

hmm.disease <- nimbleCode({
  
  # priors
  phiH ~ dunif(0, 1) # prior survival healthy
  phiI ~ dunif(0, 1) # prior survival ill
  psiIH ~ dunif(0, 1) # prior transition ill -> healthy
  psiHI ~ dunif(0, 1) # prior transition healthy -> ill
  pH ~ dunif(0, 1) # prior detection healthy
  pI ~ dunif(0, 1) # prior detection ill
  pi ~ dunif(0, 1) # prob init state healthy
  betaH ~ dunif(0,1)
  betaI ~ dunif(0,1)

  # HMM ingredients
  delta[1] <- pi         # Pr(healthy t = 1) = pi
  delta[2] <- 1 - pi     # Pr(ill t = 1) = 1 - pi
  delta[3] <- 0          # Pr(dead t = 1) = 0
  
  gamma[1,1] <- phiH * (1 - psiHI)      # Pr(H t -> H t+1)
  gamma[1,2] <- phiH * psiHI            # Pr(H t -> I t+1)
  gamma[1,3] <- 1 - phiH                # Pr(alive t -> dead t+1)
  gamma[2,1] <- phiI * psiIH            # Pr(I t -> H t+1)
  gamma[2,2] <- phiI * (1 - psiIH)      # Pr(I t -> I t+1)
  gamma[2,3] <- 1 - phiI                # Pr(alive t -> dead t+1)
  gamma[3,1] <- 0                       # Pr(dead t -> alive t+1)
  gamma[3,2] <- 0                       # Pr(dead t -> alive t+1)
  gamma[3,3] <- 1                       # Pr(dead t -> dead t+1)
  
  omega[1,1] <- 1 - pH           # Pr(H t -> non-detected t)
  omega[1,2] <- pH * betaH       # Pr(H t -> detected H t)
  omega[1,3] <- 0                # Pr(H t -> detected I t)
  omega[1,4] <- pH * (1 - betaH) # Pr(H t -> detected U t)
  omega[2,1] <- 1 - pI           # Pr(I t -> non-detected t)
  omega[2,2] <- 0                # Pr(I t -> detected H t)
  omega[2,3] <- pI * betaI       # Pr(I t -> detected I t)
  omega[2,4] <- pI * (1 - betaI) # Pr(I t -> detected U t)
  omega[3,1] <- 1                # Pr(dead t -> non-detected t)
  omega[3,2] <- 0                # Pr(dead t -> detected H t)
  omega[3,3] <- 0                # Pr(dead t -> detected I t)
  omega[3,4] <- 0                # Pr(dead t -> detected U t)
  
  omega.init[1,1] <- 0                # Pr(H t = 1 -> non-detected t = 1)
  omega.init[1,2] <- betaH            # Pr(H t = 1 -> detected H t = 1)
  omega.init[1,3] <- 0                # Pr(H t = 1 -> detected I t = 1)
  omega.init[1,4] <- 1 - betaH        # Pr(H t = 1 -> detected U t = 1)
  omega.init[2,1] <- 0                # Pr(I t = 1 -> non-detected t = 1)
  omega.init[2,2] <- 0                # Pr(I t = 1 -> detected H t = 1)
  omega.init[2,3] <- betaI            # Pr(I t = 1 -> detected I t = 1)
  omega.init[2,4] <- 1 - betaI        # Pr(I t = 1 -> detected U t = 1)
  omega.init[3,1] <- 1                # Pr(dead t = 1 -> non-detected t = 1)
  omega.init[3,2] <- 0                # Pr(dead t = 1 -> detected H t = 1)
  omega.init[3,3] <- 0                # Pr(dead t = 1 -> detected I t = 1)
  omega.init[3,4] <- 0                # Pr(dead t = 1 -> detected U t = 1)

  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:4])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3])
      y[i,j] ~ dcat(omega[z[i,j], 1:4])
    }
  }
})

first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), K = ncol(y), first = first)
my.data <- list(y = y + 1)

# initial values for unknown z
zinit <- y
zinit[zinit==3] <- sample(c(1,2), sum(zinit==3), replace = TRUE)
for (i in 1:nrow(y)) {
  for (j in 1:ncol(y)) {
    if (j > first[i] & y[i,j]==0) {zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j < first[i]) {zinit[i,j] <- 0}
  }
}
zinit <- as.matrix(zinit)

initial.values <- function() list(phiH = runif(1, 0, 1),
                                  phiI = runif(1, 0, 1),
                                  pH = runif(1, 0, 1),
                                  pI = runif(1, 0, 1),
                                  pi = runif(1, 0, 1),
                                  betaH = runif(1, 0, 1),
                                  betaI = runif(1, 0, 1),
                                  psiHI = runif(1, 0, 1),
                                  psiIH = runif(1, 0, 1),
                                  z = zinit)

parameters.to.save <- c("phiH", 
                        "phiI", 
                        "pH", 
                        "pI", 
                        "pi",
                        "betaH", 
                        "betaI", 
                        "psiHI", 
                        "psiIH")

n.iter <- 40000
n.burnin <- 25000
n.chains <- 2

out <- nimbleMCMC(code = hmm.disease, 
                           constants = my.constants,
                           data = my.data,              
                           inits = initial.values,
                           monitors = parameters.to.save,
                           niter = n.iter,
                           nburnin = n.burnin, 
                           nchains = n.chains)

library(MCMCvis)

MCMCsummary(out, round = 2)
MCMCtrace(out, pdf = FALSE, params = c("pH", "pI"))
MCMCtrace(out, pdf = FALSE, params = c("betaH", "betaI"))
MCMCtrace(out, pdf = FALSE, params = c("phiH", "phiI"))
MCMCtrace(out, pdf = FALSE, params = c("psiHI", "psiIH"))
MCMCtrace(out, pdf = FALSE, params = c("pi"))

save(out, file = "disease.RData")





############# WOLF

# read in data
y <- as.matrix(read.table(here::here("slides", "dat", "wolf.txt")))
# 1 seen
# 0 not seen

library(nimble)

hmm.phipmix <- nimbleCode({
  
  # priors
  phi ~ dunif(0, 1) # prior survival
  p1 ~ dunif(0, 1) # prior detection
  p2 ~ dunif(0, 1) # prior detection
  pi ~ dunif(0, 1) # prob init state 1

  # HMM ingredients
  gamma[1,1] <- phi      # Pr(alive 1 t -> alive 1 t+1)
  gamma[1,2] <- 0        # Pr(alive 1 t -> alive 2 t+1) // no transition
  gamma[1,3] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(alive 1 t -> alive 1 t+1) // no transition
  gamma[2,2] <- phi      # Pr(alive 1 t -> alive 2 t+1) 
  gamma[2,3] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[3,1] <- 0        # Pr(dead t -> alive 1 t+1)
  gamma[3,2] <- 0        # Pr(dead t -> alive 1 t+1)
  gamma[3,3] <- 1        # Pr(dead t -> dead t+1)
  delta[1] <- pi         # Pr(alive t = 1) = pi
  delta[2] <- 1 - pi     # Pr(alive t = 1) = 1 - pi
  delta[3] <- 0          # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p1   # Pr(alive state 1 t -> non-detected t)
  omega[1,2] <- p1       # Pr(alive state 1 t -> detected t)
  omega[2,1] <- 1 - p2   # Pr(alive state 2 t -> non-detected t)
  omega[2,2] <- p2       # Pr(alive state 2 t -> detected t)
  omega[3,1] <- 1        # Pr(dead t -> non-detected t)
  omega[3,2] <- 0        # Pr(dead t -> detected t)
  omega.init[1,1] <- 0   # Pr(alive state 1 t -> non-detected t)
  omega.init[1,2] <- 1       # Pr(alive state 1 t -> detected t)
  omega.init[2,1] <- 0   # Pr(alive state 2 t -> non-detected t)
  omega.init[2,2] <- 1       # Pr(alive state 2 t -> detected t)
  omega.init[3,1] <- 1        # Pr(dead t -> non-detected t)
  omega.init[3,2] <- 0        # Pr(dead t -> detected t)
  
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:2])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

first <- apply(y, 1, function(x) min(which(x != 0)))
my.constants <- list(N = nrow(y), K = ncol(y), first = first)
my.data <- list(y = y + 1)

# initial values for unknown z
zinit <- y
for (i in 1:nrow(y)) {
  for (j in first[i]:ncol(y)) {
    if (j == first[i]) zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1) # pick alive state
    if (j > first[i]) zinit[i,j] <- zinit[i,j-1] # because no transition, state remains the same
  }
}
zinit <- as.matrix(zinit)

initial.values <- function() list(phi = runif(1,0,1),
                                  p1 = runif(1,0,1),
                                  p2 = runif(1,0,1),
                                  pi = runif(1,0,1),
                                  z = zinit)

parameters.to.save <- c("phi", "p1", "p2", "pi")

n.iter <- 10000
n.burnin <- 5000
n.chains <- 2

mcmc.phipmix <- nimbleMCMC(code = hmm.phipmix, 
                           constants = my.constants,
                           data = my.data,              
                           inits = initial.values,
                           monitors = parameters.to.save,
                           niter = n.iter,
                           nburnin = n.burnin, 
                           nchains = n.chains)

library(MCMCvis)
MCMCsummary(mcmc.phipmix, round = 2)
MCMCtrace(mcmc.phipmix, pdf = FALSE, params = c("pi", "p1", "p2"))
MCMCtrace(mcmc.phipmix, pdf = FALSE, params = c("phi"))
save(mcmc.phipmix, file = "wolf_het.RData")


# read in data
y <- as.matrix(read.table(here::here("slides", "dat", "wolf.txt")))
# 1 seen
# 0 not seen

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

  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), T = ncol(y), first = first)
my.data <- list(y = y + 1)

zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")

n.iter <- 10000
n.burnin <- 5000
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
MCMCplot(mcmc.phip)
save(mcmc.phip, file = "wolf_homg.RData")



# ADD TRANSITION BETWEEN HETEROGENEIT STATES?


######---------- MEMORY
library(tidyverse)
library(nimble)
geese <- read_csv(here::here("slides", "dat", "geese.csv"), col_names = TRUE)
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
    if (j > first[i] & y[i,j]==0) {zinit[i,j] <- which(rmultinom(1, 1, c(1/4,1/4,1/4,1/4))==1)}
    if (j > first[i] & y[i,j]==1) {zinit[i,j] <- which(rmultinom(1, 1, c(1/4,1/4,1/4,1/4))==1)}
    if (j > first[i] & y[i,j]==2) {zinit[i,j] <- which(rmultinom(1, 1, c(1/4,1/4,1/4,1/4))==1)}
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
n.iter <- 5000
n.burnin <- 2500
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



## PROBLEM: p1 and p2 seem to be updated, but Rhat  = Inf and neff = 0; inspecting their posterior distribution
## it feels like they are not updated, two modes corresp to the two chains. 
# Usually, it means that these parameters are not used in the model. 





hmm.phip <- nimbleModel(code = memory,
                        constants = my.constants,
                        data = my.data,
                        inits = initial.values())
phip.mcmc <- buildMCMC(hmm.phip)
phip.model <- compileNimble(hmm.phip) 
c.phip.mcmc <- compileNimble(phip.mcmc, project = phip.model)
c.phip.mcmc$run(10000)
samples <- as.matrix(c.phip.mcmc$mvSamples)
summary(samples[,"det1"])




####################################################  age uncertainty

# OBS
# not detected (1)
# detected as cub by both criteria (2)
# detected as cub by criterion 1 (3)
# detected as adult by both criteria (4)

# STATES
# alive as cub (1)
# alive as adult (2)
# dead (3)

# read in data
y_raw <- as.matrix(read.table(here::here("dat", "apennine-bear.txt")))
# encounter histories
y <- apply(y_raw[, 1:7], 2, as.numeric)

# sex and age
sex <- as.numeric(as.factor(y_raw[, 8])) # 1 is female, 2 is male
age <- as.numeric(factor(as.factor(y_raw[, 9]), levels = c("Unknown", "Cub", "Adult"))) - 1

library(nimble)

hmm.age <- nimbleCode({
  
  # priors
  phiC ~ dunif(0, 1) # prior survival cub
  phiA[1] ~ dunif(0, 1) # prior survival adult female
  phiA[2] ~ dunif(0, 1) # prior survival adult male
  p ~ dunif(0, 1) # prior detection cub
  pi ~ dunif(0, 1) # prob init state cub
  betaCCC ~ dunif(0,1)
  betaCCA ~ dunif(0,1)
  betaACC ~ dunif(0,1)
  betaACA ~ dunif(0,1)
  
  # HMM ingredients
  
  # As we knew the age-class at first detection for some bears, we fixed 
  # their initial state in the model, thus providing the model with a subset of individuals of known age.
  for (i in 1:N){
    # age[i] == 0 → unknown → use probabilistic delta: pi, 1 - pi
    # age[i] == 1 → known cub → force delta to [1, 0, 0]
    # age[i] == 2 → known adult → force delta to [0, 1, 0]
    delta[1, i] <- equals(age[i], 0) * pi + equals(age[i], 1)
    delta[2, i] <- equals(age[i], 0) * (1 - pi) + equals(age[i], 2)
    delta[3, i] <- 0
  }
  
  for (i in 1:N){
    gamma[1,1,i] <- 0                       # Pr(C t -> C t+1)
    gamma[1,2,i] <- phiC                    # Pr(C t -> A t+1)
    gamma[1,3,i] <- 1 - phiC                # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0                       # Pr(A t -> C t+1)
    gamma[2,2,i] <- phiA[sex[i]]            # Pr(A t -> A t+1)
    gamma[2,3,i] <- 1 - phiA[sex[i]]        # Pr(alive t -> dead t+1)
    gamma[3,1,i] <- 0                       # Pr(dead t -> alive t+1)
    gamma[3,2,i] <- 0                       # Pr(dead t -> alive t+1)
    gamma[3,3,i] <- 1                       # Pr(dead t -> dead t+1)
  }
  
  omega[1,1] <- 1 - p                                  # Pr(H t -> non-detected t)
  omega[1,2] <- p * betaCCC                            # Pr(H t -> detected H t)
  omega[1,3] <- p * betaCCA                            # Pr(H t -> detected I t)
  omega[1,4] <- p * (1 - betaCCC - betaCCA)            # Pr(H t -> detected U t)
  omega[2,1] <- 1 - p                                  # Pr(I t -> non-detected t)
  omega[2,2] <- p * betaACC                            # Pr(I t -> detected H t)
  omega[2,3] <- p * betaACA                            # Pr(I t -> detected I t)
  omega[2,4] <- p * (1 - betaACC - betaACA)            # Pr(I t -> detected U t)
  omega[3,1] <- 1                                      # Pr(dead t -> non-detected t)
  omega[3,2] <- 0                                      # Pr(dead t -> detected H t)
  omega[3,3] <- 0                                      # Pr(dead t -> detected I t)
  omega[3,4] <- 0                                      # Pr(dead t -> detected U t)
  
  omega.init[1,1] <- 0                      # Pr(H t = 1 -> non-detected t = 1)
  omega.init[1,2] <- betaCCC                # Pr(H t = 1 -> detected H t = 1)
  omega.init[1,3] <- betaCCA                # Pr(H t = 1 -> detected I t = 1)
  omega.init[1,4] <- 1 - betaCCC - betaCCA  # Pr(H t = 1 -> detected U t = 1)
  omega.init[2,1] <- 0                      # Pr(I t = 1 -> non-detected t = 1)
  omega.init[2,2] <- betaACC                # Pr(I t = 1 -> detected H t = 1)
  omega.init[2,3] <- betaACA                # Pr(I t = 1 -> detected I t = 1)
  omega.init[2,4] <- 1 - betaACC - betaACA  # Pr(I t = 1 -> detected U t = 1)
  omega.init[3,1] <- 1                      # Pr(dead t = 1 -> non-detected t = 1)
  omega.init[3,2] <- 0                      # Pr(dead t = 1 -> detected H t = 1)
  omega.init[3,3] <- 0                      # Pr(dead t = 1 -> detected I t = 1)
  omega.init[3,4] <- 0                      # Pr(dead t = 1 -> detected U t = 1)
  
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3,i])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:4])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:4])
    }
  }
})

first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), K = ncol(y), first = first, sex = sex, age = age) 
# sex = 1 is female, 2 is male
# age = 0 if unknown, 1 if cub and 2 if adult
my.data <- list(y = y + 1)

# initial values for unknown z
zinit <- y
zinit[zinit==2] <- sample(c(1,3), sum(zinit==2), replace = TRUE) # if unsure, pick cub or adult at random
zinit[zinit==3] <- 2 # 3s become 2s for alive as adults
# fill in 0s w/ cub/adult alive states
for (i in 1:nrow(y)) { 
  for (j in 1:ncol(y)) {
    if (j > first[i] & y[i,j]==0) {zinit[i,j] <- 2}
    if (j < first[i]) {zinit[i,j] <- 0}
  }
}
zinit <- as.matrix(zinit)

# for those individuals for which state is known at first detection, use it
for (i in 1:nrow(y)) {
  if (age[i] !=0) zinit[i,first[i]] <- age[i]
}

# Function to generate valid initial beta values
generate_valid_betas <- function() {
  # Ensure beta parameters allow valid probabilities
  repeat {
    betaCCC <- runif(1, 0.1, 0.8)
    betaCCA <- runif(1, 0.1, 0.8)
    if (betaCCC + betaCCA <= 0.9) break
  }
  
  repeat {
    betaACC <- runif(1, 0.1, 0.8)
    betaACA <- runif(1, 0.1, 0.8)
    if (betaACC + betaACA <= 0.9) break
  }
  
  return(list(betaCCC = betaCCC, betaCCA = betaCCA, 
              betaACC = betaACC, betaACA = betaACA))
}

initial.values <- function() {
  betas <- generate_valid_betas()
  list(phiC = runif(1, 0.5, 0.9),
       phiA = runif(2, 0.7, 0.95),
       p = runif(1, 0.3, 0.8),
       pi = runif(1, 0.2, 0.6),
       betaCCC = betas$betaCCC,
       betaCCA = betas$betaCCA,
       betaACC = betas$betaACC,
       betaACA = betas$betaACA,
       z = zinit)
}

parameters.to.save <- c("phiC", 
                        "phiA", 
                        "p", 
                        "pi",
                        "betaCCC", 
                        "betaCCA", 
                        "betaACC", 
                        "betaACA")

n.iter <- 40000
n.burnin <- 25000
n.chains <- 2

out <- nimbleMCMC(code = hmm.age, 
                  constants = my.constants,
                  data = my.data,              
                  inits = initial.values,
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin, 
                  nchains = n.chains)

save(out, file = "age.RData")

library(MCMCvis)

MCMCsummary(out, round = 2)
MCMCtrace(out, pdf = FALSE, params = c("p"))
MCMCtrace(out, pdf = FALSE, params = c("betaCCC", "betaCCA"))
MCMCtrace(out, pdf = FALSE, params = c("betaACC", "betaACA"))
MCMCtrace(out, pdf = FALSE, params = c("phiC", "phiA"))
MCMCtrace(out, pdf = FALSE, params = c("pi"))

save(out, file = "age.RData")





####################################################  sex uncertainty

# STATES
# alive as male (M)
# alive as female (F)
# dead (D)

# OBS
# not seen (1)
# judged from copulation to be M (2)
# judged from begging food to be M (3)
# judged from coutship feeding to be M (4)
# judged from body size to be M (5)
# judged from copulation to be F (6)
# judged from begging food to be F (7)
# judged from coutship feeding to be F (8)
# judged from body size to be F (9)
# not judged (10)

# # Read the data (as a matrix, assuming space-separated values in a file)
# dat <- as.matrix(read.table("/Users/oliviergimenez/Dropbox/OG/GITHUB/book/matos/case-studies/covariates/sex-uncertainty/sexo_audouines9/Sexaudouin9.text"))  # or use here::here("path", "to", "file.txt")
# 
# # Separate sequence and count
# sequence <- dat[, 1:(ncol(dat)-1)]
# counts <- dat[, ncol(dat)]
# 
# ind <- c(1, 176, 991, 177, 178, 179, 468, 667, 668, 669, 180, 670, 469, 671, 181, 182, 830, 183, 917, 184, 918, 992, 672, 993) 
# length(ind)
# 
# 
# 
# geneticsex <- c(rep('F', 11), rep('M', 13))
# 
# # Expand rows according to counts
# expanded <- sequence[rep(1:nrow(sequence), counts), ]
# 
# # Convert to data.frame (optional)
# expanded_df <- as.data.frame(expanded)
# 
# # Write to file (optional)
# write.table(expanded_df, file = "audouin-gull.txt", row.names = FALSE, col.names = FALSE)
# 

# read in data
y <- as.matrix(read.table(here::here("dat", "audouin-gull.txt")))

library(nimble)

hmm.sex <- nimbleCode({
  
  # priors
  phiF ~ dunif(0, 1) # prior survival male
  phiM ~ dunif(0, 1) # prior survival female
  p ~ dunif(0, 1) # prior encounter prob
  pi ~ dunif(0, 1) # prob init state male
  e ~ dunif(0,1)
  for (j in 1:4){
    x[j] ~ dunif(0,1)
  }
  m[1] ~ dunif(0,1)
  m[2] <- (1 - m[1]) * alpha
  alpha ~ dunif(0,1)
  m[3] <- 1 - m[1] - m[2]
  m[4] ~ dunif(0,1)
  
  # HMM ingredients
  delta[1] <- pi
  delta[2] <- 1 - pi
  delta[3] <- 0
  
  gamma[1,1] <- phiM            # Pr(C t -> C t+1)
  gamma[1,2] <- 0               # Pr(C t -> A t+1)
  gamma[1,3] <- 1 - phiM        # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0               # Pr(A t -> C t+1)
  gamma[2,2] <- phiF            # Pr(A t -> A t+1)
  gamma[2,3] <- 1 - phiF        # Pr(alive t -> dead t+1)
  gamma[3,1] <- 0               # Pr(dead t -> alive t+1)
  gamma[3,2] <- 0               # Pr(dead t -> alive t+1)
  gamma[3,3] <- 1               # Pr(dead t -> dead t+1)
  
  omega[1,1] <- 1 - p                                  
  omega[1,2] <- p * e * (1 - m[4]) * m[1] * x[1]      
  omega[1,3] <- p * e * (1 - m[4]) * m[2] * x[2]      
  omega[1,4] <- p * e * (1 - m[4]) * m[3] * x[3]      
  omega[1,5] <- p * e * m[4] * x[4]      
  omega[1,6] <- p * e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega[1,7] <- p * e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega[1,8] <- p * e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega[1,9] <- p * e * m[4] * (1 - x[4])      
  omega[1,10] <- p * (1 - e)      
  omega[2,1] <- 1 - p                                  
  omega[2,2] <- p * e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega[2,3] <- p * e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega[2,4] <- p * e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega[2,5] <- p * e * m[4] * (1 - x[4])      
  omega[2,6] <- p * e * (1 - m[4]) * m[1] * x[1]      
  omega[2,7] <- p * e * (1 - m[4]) * m[2] * x[2]     
  omega[2,8] <- p * e * (1 - m[4]) * m[3] * x[3]   
  omega[2,9] <- p * e * m[4] * x[4]   
  omega[2,10] <- p * (1 - e)      
  omega[3,1] <- 1                                  
  omega[3,2] <- 0     
  omega[3,3] <- 0     
  omega[3,4] <- 0     
  omega[3,5] <- 0     
  omega[3,6] <- 0     
  omega[3,7] <- 0     
  omega[3,8] <- 0     
  omega[3,9] <- 0     
  omega[3,10] <- 0     

  omega.init[1,1] <- 0                                 
  omega.init[1,2] <- e * (1 - m[4]) * m[1] * x[1]      
  omega.init[1,3] <- e * (1 - m[4]) * m[2] * x[2]      
  omega.init[1,4] <- e * (1 - m[4]) * m[3] * x[3]      
  omega.init[1,5] <- e * m[4] * x[4]      
  omega.init[1,6] <- e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega.init[1,7] <- e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega.init[1,8] <- e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega.init[1,9] <- e * m[4] * (1 - x[4])      
  omega.init[1,10] <- (1 - e)      
  omega.init[2,1] <- 0                                  
  omega.init[2,2] <- e * (1 - m[4]) * m[1] * (1 - x[1])      
  omega.init[2,3] <- e * (1 - m[4]) * m[2] * (1 - x[2])      
  omega.init[2,4] <- e * (1 - m[4]) * m[3] * (1 - x[3])      
  omega.init[2,5] <- e * m[4] * (1 - x[4])      
  omega.init[2,6] <- e * (1 - m[4]) * m[1] * x[1]      
  omega.init[2,7] <- e * (1 - m[4]) * m[2] * x[2]     
  omega.init[2,8] <- e * (1 - m[4]) * m[3] * x[3]   
  omega.init[2,9] <- e * m[4] * x[4]   
  omega.init[2,10] <- (1 - e)      
  omega.init[3,1] <- 1                                  
  omega.init[3,2] <- 0     
  omega.init[3,3] <- 0     
  omega.init[3,4] <- 0     
  omega.init[3,5] <- 0     
  omega.init[3,6] <- 0     
  omega.init[3,7] <- 0     
  omega.init[3,8] <- 0     
  omega.init[3,9] <- 0     
  omega.init[3,10] <- 0     
  
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:10])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3])
      y[i,j] ~ dcat(omega[z[i,j], 1:10])
    }
  }
})

first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), K = ncol(y), first = first) 
my.data <- list(y = y + 1)


# # Function to create initial latent states z
# #create_initial_z <- function(y, first) {
#   N <- nrow(y)
#   K <- ncol(y)
#   zinit <- matrix(NA, nrow = N, ncol = K)
#   
#   for (i in 1:N) {
#     # Détermination du sexe à la première détection
#     first_obs <- y[i, first[i]]  # valeurs déjà de 0-9, pas besoin de +1
#     
#     # Classification basée sur les codes d'observation
#     # 0 = non observé, 1-5 = mâle, 6-8 = femelle, 9 = incertain
#     if (first_obs %in% 1:5) {
#       initial_sex <- 1  # mâle
#     } else if (first_obs %in% 6:8) {
#       initial_sex <- 2  # femelle
#     } else {
#       initial_sex <- sample(1:2, 1)  # sexe aléatoire si incertain (9)
#     }
#     
#     zinit[i, first[i]] <- initial_sex
#     
#     # Pour les occasions suivantes
#     current_state <- initial_sex
#     if (first[i] < K) {
#       for (j in (first[i] + 1):K) {
#         obs <- y[i, j]
#         
#         if (obs == 0) {  # non observé
#           # Probabilité de survie élevée si pas d'évidence de mort
#           if (runif(1) < 0.85) {
#             zinit[i, j] <- current_state  # survit
#           } else {
#             zinit[i, j] <- 3  # meurt
#             current_state <- 3
#           }
#         } else {  # observé vivant (1-9)
#           # Maintient l'état vivant
#           if (current_state == 3) {
#             # Si était mort, ne peut pas revenir vivant - problème dans les données
#             current_state <- sample(1:2, 1)
#           }
#           zinit[i, j] <- current_state
#         }
#       }
#     }
#   }
#   
#   return(zinit)
# }

# # Function to create initial latent states z
create_initial_z <- function(y, first) {
  N <- nrow(y)
  K <- ncol(y)
  zinit <- matrix(NA, nrow = N, ncol = K)
  
  for (i in 1:N) {
    # Détermination du sexe à la première détection
    first_obs <- y[i, first[i]]  # valeurs déjà de 0-9, pas besoin de +1
    
    # Classification basée sur les codes d'observation
    # 0 = non observé, 1-5 = mâle, 6-8 = femelle, 9 = incertain
    if (first_obs %in% 1:5) {
      initial_sex <- 1  # mâle
    } else if (first_obs %in% 6:8) {
      initial_sex <- 2  # femelle
    } else {
      initial_sex <- sample(1:2, 1)  # sexe aléatoire si incertain (9)
    }
    
    zinit[i, first[i]] <- initial_sex
    
    # Pour les occasions suivantes
    current_state <- initial_sex
    if (first[i] < K) {
      for (j in (first[i] + 1):K) {
        obs <- y[i, j]
        
        if (obs == 0) {  # non observé
            zinit[i, j] <- current_state  # survit
        } else {  # observé vivant (1-9)
          # Maintient l'état vivant
          if (current_state == 3) {
            # Si était mort, ne peut pas revenir vivant - problème dans les données
            current_state <- sample(1:2, 1)
          }
          zinit[i, j] <- current_state
        }
      }
    }
  }
  
  return(zinit)
}

# Function to generate valid m parameters
generate_valid_m <- function() {
  m <- rep(NA, 4)
  
  # m[4] : proportion de la 4ème catégorie
  m[4] <- runif(1, 0.1, 0.4)
  
  # m[1] et m[2] avec contrainte que m[1] + m[2] + m[3] = 1 et m[3] > 0
  repeat {
    m[1] <- runif(1, 0.1, 0.4)
    alpha <- runif(1, 0.1, 0.4)
    m[2] <- (1 - m[1]) * alpha
    if (m[1] + m[2] <= 0.8) break  # assure que m[3] >= 0.2
  }
  m[3] <- 1 - m[1] - m[2]
  
  return(m)
}

# Function to create valid initial parameter values - VERSION CORRIGEE
create_initial_params <- function() {
  # Paramètres avec des valeurs plus réalistes
  phiM <- runif(1, 0.75, 0.95)  # survie mâle
  phiF <- runif(1, 0.75, 0.95)  # survie femelle
  p <- runif(1, 0.4, 0.8)       # probabilité de détection
  pi <- runif(1, 0.45, 0.55)    # ratio des sexes
  e <- runif(1, 0.3, 0.7)       # probabilité d'encounter given detected
  x <- 1 - runif(4, 0.01, 0.3)       # probabilités de sexage correct
  m <- generate_valid_m()
  
  return(list(
    phiM = phiM,
    phiF = phiF, 
    p = p,
    pi = pi,
    e = e,
    x = x,
    m = m
  ))
}

# Create the complete initial values function
initial.values <- function() {
  params <- create_initial_params()
  zinit <- create_initial_z(y, first)
  return(c(params, list(z = zinit)))
}

# # Fonction de diagnostic pour vérifier les valeurs initiales
# check_initial_values <- function() {
#   init_vals <- initial.values()
#   
#   cat("Vérification des paramètres initiaux:\n")
#   cat("phiM:", init_vals$phiM, "\n")
#   cat("phiF:", init_vals$phiF, "\n")
#   cat("p:", init_vals$p, "\n")
#   cat("pi:", init_vals$pi, "\n")
#   cat("e:", init_vals$e, "\n")
#   cat("x:", init_vals$x, "\n")
#   cat("m:", init_vals$m, "\n")
#   cat("Somme m[1:3]:", sum(init_vals$m[1:3]), "\n")
#   
#   # Vérifier quelques probabilités omega
#   e_val <- init_vals$e
#   m_val <- init_vals$m
#   x_val <- init_vals$x
#   p_val <- init_vals$p
#   
#   # Probabilité pour état 1, observation 2
#   prob_example <- p_val * e_val * (1 - m_val[4]) * m_val[1] * x_val[1]
#   cat("Exemple probabilité omega[1,2]:", prob_example, "\n")
#   
#   return(init_vals)
# }

parameters.to.save <- c("phiM", 
                        "phiF", 
                        "p", 
                        "pi",
                        "e", 
                        "m", 
                        "x")

n.iter <- 2500
n.burnin <- 500
n.chains <- 2

# # Diagnostic avant création du modèle
# cat("=== DIAGNOSTIC DES VALEURS INITIALES ===\n")
# test_init <- check_initial_values()
# 
# # Création du modèle avec gestion d'erreur
# cat("\n=== CRÉATION DU MODÈLE ===\n")
# set.seed(123)
# 

# survival <- nimbleModel(code = hmm.sex,
#                         data = my.data,
#                         constants = my.constants,
#                         inits = initial.values())
#   
# survival$calculate()
# 
# for (i in 1:N) {
#   for (j in first[i]:K) {
#     z_ij <- survival$z[i,j]
#     y_ij <- survival$y[i,j]
#     if (z_ij %in% 1:3 && y_ij %in% 1:10) {
#       prob <- survival$omega[z_ij, y_ij]
#       if (prob == 0) {
#         cat("⚠️ omega[", z_ij, ",", y_ij, "] == 0 at y[", i, ",", j, "]\n")
#       }
#     }
#   }
# }
# 
# survival$getNodeNames()
# 
# survival$e
# survival$p
# survival$m
# survival$x
# 
# apply(survival$omega,1,sum)
# apply(survival$omega.init,1,sum)
# 
# survival$gamma
# apply(survival$gamma, 1, sum)
# 
# survival$delta

out <- nimbleMCMC(code = hmm.sex, 
                  constants = my.constants,
                  data = my.data,              
                  inits = initial.values(),
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin, 
                  nchains = n.chains)

save(out, file = "sex.RData")

library(MCMCvis)

MCMCsummary(out, round = 2)
MCMCtrace(out, pdf = FALSE, params = c("p"))
MCMCtrace(out, pdf = FALSE, params = c("e"))
MCMCtrace(out, pdf = FALSE, params = c("m"))
MCMCtrace(out, pdf = FALSE, params = c("x"))
MCMCtrace(out, pdf = FALSE, params = c("phiM", "phiF"))
MCMCtrace(out, pdf = FALSE, params = c("pi"))



########################################## RJMCM storks

# white stork RJMCMC
# rainfall at several meteorological stations in Africa

# dakar	diourbel (2)	gao (3)	kandi	kayes (5)	kita (6)	koutiala	maradi (8)	mopti (9)	natitingou	ouahigouya (11)	sikasso	ségou (13)	tahoua (14)	tombouctou (15) 
rainfall_raw <- as.matrix(read.table(here::here("dat", "pluies.txt")))
# select Diourbel, Gao, Kayes, Kita, Maradi, Mopti, Ouahigouya, Ségou, Tahoua, and Tombouctou
filter <- c(2,3,5,6,8,9,11,13,14,15)
rainfall_raw <- rainfall_raw[,filter]
mean_rainfall <- apply(rainfall_raw, 2, mean)
sd_rainfall <- apply(rainfall_raw, 2, sd)
rainfall <- (rainfall_raw - matrix(rep(mean_rainfall, 16), ncol = 10, byrow = T)) / matrix(rep(sd_rainfall, 16), ncol = 10, byrow = T)

# dat <- as.matrix(read.table(here::here("dat", "stork_grouped.txt")))
# 
# # Separate sequence and count
# sequence <- dat[, 1:(ncol(dat)-1)]
# counts <- dat[, ncol(dat)]
# 
# # Expand rows according to counts
# expanded <- sequence[rep(1:nrow(sequence), counts), ]
# 
# # Convert to data.frame (optional)
# expanded_df <- as.data.frame(expanded)
# 
# # Write to file (optional)
# write.table(expanded_df, file = "stork.txt", row.names = FALSE, col.names = FALSE)
# 
# R2ucare::marray(sequence, counts)

library(nimble)

y <- as.matrix(read.table(here::here("dat", "stork.txt")))

first <- apply(y, 1, function(x) min(which(x !=0)))

# Then we write the NIMBLE code:
hmm.phixp <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    logit(phi[t]) <- intercept + inprod(beta[1:10], x[t, 1:10])
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
  }
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (j in 1:10){
    beta[j] ~ dnorm(0, sd = 1.5) # prior slopes
  }
  intercept ~ dnorm(0, sd = 1.5) # prior intercept
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

#Let's put our constants in a list:
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     x = rainfall)

#Then our function for generating initial values: 
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(intercept = rnorm(1,0,1),
                                  beta = rnorm(10,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

#And the parameters to be monitored: 
parameters.to.save <- c("beta", "p", "intercept")

my.data <- list(y = y + 1)

#Finaly, we run NIMBLE: 
n.iter <- 15000
n.burnin <- 5000
n.chains <- 2
# out <- nimbleMCMC(code = hmm.phixp, 
#                   constants = my.constants,
#                   data = my.data,              
#                   inits = initial.values,
#                   monitors = parameters.to.save,
#                   niter = n.iter,
#                   nburnin = n.burnin, 
#                   nchains = n.chains)
# 
# library(MCMCvis)
# MCMCsummary(out, round = 2)


#### with RJMCMC, cf manual section 7.10.1 


# Then we write the NIMBLE code:
hmm.phirjmcmcp <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  for (j in 1:10){
    betaksi[j] <- beta[j] * ksi[j]
    ksi[j] ~ dbern(psi) ## indicator variable associated with betaj
  }
  psi ~ dbeta(1, 1) ## hyperprior on inclusion probability
  for (t in 1:(T-1)){
    logit(phi[t]) <- intercept + inprod(betaksi[1:10], x[t, 1:10])
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
  }
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (j in 1:10){
    beta[j] ~ dnorm(0, sd = 1.5) # prior slopes
  }
  intercept ~ dnorm(0, sd = 1.5) # prior intercept
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

#Let's put our constants in a list:
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     x = rainfall)

#Then our function for generating initial values: 
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(intercept = rnorm(1,0,1),
                                  beta = rnorm(10,0,1),
                                  p = runif(1,0,1),
                                  z = zinits,
                                  eta = sample(c(0,1), 10, replace = T),
                                  psi = 0.5)

## build the model and configure default MCMC
RJsurvival <- nimbleModel(code = hmm.phirjmcmcp, 
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values())

Csurvival <- compileNimble(RJsurvival)

RJsurvivalConf <- configureMCMC(RJsurvival)
RJsurvivalConf$addMonitors('ksi')

# modify the MCMC configuration object, RJsurvivalConf, to use reversible
# jump samplers for selection on betas
configureRJ(conf = RJsurvivalConf,
            targetNodes = 'beta',
            indicatorNodes = 'ksi',
            control = list(mean = 0, scale = .2))

survivalMCMC <- buildMCMC(RJsurvivalConf)
CsurvivalMCMC <- compileNimble(survivalMCMC, project = RJsurvival)

# run NIMBLE:
samples <- runMCMC(mcmc = CsurvivalMCMC, 
                   niter = n.iter,
                   nburnin = n.burnin,
                   nchains = n.chains)
## obtain numerical summaries:
#samplesSummary(samples)

library(tidyverse)

# individual inclusion proportion
out <- rbind(samples[[1]], samples[[2]])

grepl('ksi', colnames(out)) %>% 
  which() %>%
  out[, .] %>%
  colMeans()

#ksi[1]  ksi[2]  ksi[3]  ksi[4]  ksi[5]  ksi[6]  ksi[7]  ksi[8]  ksi[9] ksi[10] 
#0.02575 0.01720 0.01560 0.96120 0.04430 0.05350 0.02025 0.03360 0.02630 0.02965

grepl('beta', colnames(out)) %>% 
  which() %>%
  out[, .] %>%
  colMeans()

#beta[1]       beta[2]       beta[3]       beta[4]       beta[5]       beta[6]       beta[7]       beta[8]       beta[9]      beta[10] 
#2.748795e-03 -9.999704e-06  2.316042e-04  3.482938e-01 -6.244761e-03  7.479767e-03 -1.501896e-03  4.541015e-03  2.183717e-03 -3.139205e-03 
 
save(out, file = "rjmcmc.RData")

## check BAPE book for content

