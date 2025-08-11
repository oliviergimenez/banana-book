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



