library(nimble)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# Sites and states {#dispersal}

## Introduction

## The Arnason-Schwarz (AS) model {#ASmodel}

### Theory

## Fitting the AS model {#ASmodelfitting}

### Geese data

# You may see the data below: 
y <- read_csv("geese.csv") %>% as.matrix()
head(y)

### NIMBLE implementation

# We replace all 3's by 0's in the dataset:
y[y==3] <- 0
mask <- apply(y, 1, sum)
y <- y[mask!=0,]


# NIMBLE code
multisite <- nimbleCode({
  # -------------------------------------------------
  # Parameters:
  # phi1: survival probability site 1
  # phi2: survival probability site 2
  # psi12: movement probability from site 1 to site 2
  # psi21: movement probability from site 2 to site 1
  # p1: detection probability site 1
  # p2: detection probability site 2
  # -------------------------------------------------
  # States (z):
  # 1 alive at site 1
  # 2 alive at site 2
  # 3 dead
  # Observations (y):
  # 1 not seen
  # 2 seen at site 1
  # 3 seen at site 2
  # -------------------------------------------------
  
  # Priors
  phi1 ~ dunif(0, 1)
  phi2 ~ dunif(0, 1)
  psi12 ~ dunif(0, 1)
  psi21 ~ dunif(0, 1)
  p1 ~ dunif(0, 1)
  p2 ~ dunif(0, 1)
  
  # Define state-transition and observation matrices
  
  # Define probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phi1 * (1 - psi12)
  gamma[1,2] <- phi1 * psi12
  gamma[1,3] <- 1 - phi1
  gamma[2,1] <- phi2 * psi21
  gamma[2,2] <- phi2 * (1 - psi21)
  gamma[2,3] <- 1 - phi2
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 1
  
#  delta[1] <- pi1          # Pr(alive in 1 at t = first)
#  delta[2] <- 1 - pi1      # Pr(alive in 2 at t = first)
#  delta[3] <- 0            # Pr(dead at t = first) = 0
  
  # Define probabilities of y(t) given z(t)
  omega[1,1] <- 1 - p1     # Pr(alive 1 t -> non-detected t)
  omega[1,2] <- p1         # Pr(alive 1 t -> detected site 1 t)
  omega[1,3] <- 0          # Pr(alive 1 t -> detected site 2 t)
  omega[2,1] <- 1 - p2     # Pr(alive 2 t -> non-detected t)
  omega[2,2] <- 0          # Pr(alive 2 t -> detected site 1 t)
  omega[2,3] <- p2         # Pr(alive 2 t -> detected site 2 t)
  omega[3,1] <- 1          # Pr(dead t -> non-detected t)
  omega[3,2] <- 0          # Pr(dead t -> detected site 1 t)
  omega[3,3] <- 0          # Pr(dead t -> detected site 2 t)
  
  # Likelihood 
  for (i in 1:N){
    # latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1
    for (t in (first[i]+1):K){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:3])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:3])
    }
  }
})

# We need to provide NIMBLE with constants, data, initial values,
# some parameters to monitor and MCMC details:

# occasions of first capture
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)
first

# constants
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))
# data
my.data <- list(y = y + 1)

# initial values 
zinits <- y # say states are observations, detections in A or B are taken as alive in same sites 
zinits[zinits==0] <- sample(c(1,2), sum(zinits==0), replace = TRUE) # non-detections become alive in site A or B
initial.values <- function(){list(phi1 = runif(1, 0, 1), 
                                  phi2 = runif(1, 0, 1), 
                                  psi12 = runif(1, 0, 1), 
                                  psi21 = runif(1, 0, 1), 
                                  p1 = runif(1, 0, 1), 
                                  p2 = runif(1, 0, 1), 
                                  #pi1 = runif(1, 0, 1),
                                  z = zinits)}
# parameters to monitor
#parameters.to.save <- c("phi1", "phi2","psi12", "psi21", "p1", "p2", "pi1")
parameters.to.save <- c("phi1", "phi2","psi12", "psi21", "p1", "p2")

# MCMC details
n.iter <- 20000
n.burnin <- 1000
n.chains <- 2

# Now we may run NIMBLE:
mcmc.multisite <- nimbleMCMC(code = multisite, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

# We may have a look to the results via a caterpillar plot: 
MCMCplot(mcmc.multisite)

# Remember mid--Atlantic is site 1, and Chesapeake site 2. Detection in mid--Atlantic (around 0.5) is higher than in Cheasapeake (around 0.4) although it comes with more uncertainty (wider credible interval). Survival in both sites are estimated at around 0.6--0.7. Note that by going multisite, we could make these parameters site-specific and differences might reflect habitat quality for example. Now the novelty lies in our capability to estimate movements from site 1 to site 2 and from site 2 to site 1 from a winter to the next. The annual probability of remaining in the same site for two successive winters, used as a measure of site fidelity, was lower in the mid--Atlantic ($1-\psi_{12}$ around 0.8) than in the Chesapeake ($1-\psi_{21}$ around 0.9). The estimated probability of moving to the Chesapeake from the mid--Atlantic was four times as high as the probability of moving in the opposite direction.

# We may also have a look to numerical summaries, which confirm our 
# ecological interpretation of the model parameter estimates:
MCMCsummary(mcmc.multisite, round = 2)

# I will illustrate the use of `nimbleEcology` to fit the AS to the 
# geese data with 2 sites (see Section \@ref(nimblemarginalization)). 
# Using the function `dHMM()` which implements HMM with time-independent 
# observation and transition matrices, we have:

# read in data
geese <- read_csv("geese.csv", col_names = TRUE)
y <- as.matrix(geese)
# drop Carolinas
y[y==3] <- 0 # act as if there was no detection in site 3 Carolinas
mask <- apply(y, 1, sum)
y <- y[mask!=0,] # remove rows w/ 0s only
# get occasions of first encounter
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)
# filter out individuals that are first captured at last occasion. 
# These individuals do not contribute to parameter estimation, 
# and also they cause problems with nimbleEcology
mask <- which(first!=ncol(y)) 
y <- y[mask, ]                # keep only these
first <- first[mask]

# NIMBLE code 
multisite.marginalized <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phi1: survival probability site 1
  # phi2: survival probability site 2
  # psi12: movement probability from site 1 to site 2
  # psi21: movement probability from site 2 to site 1
  # p1: recapture probability site 1
  # p2: recapture probability site 2
  # pi1: prop of being in site 1 at initial capture
  # -------------------------------------------------
  # States (z):
  # 1 alive at 1
  # 2 alive at 2
  # 3 dead
  # Observations (y):  
  # 1 not seen
  # 2 seen at site 1 
  # 3 seen at site 2
  # -------------------------------------------------
  
  # priors
  phi1 ~ dunif(0, 1)
  phi2 ~ dunif(0, 1)
  psi12 ~ dunif(0, 1)
  psi21 ~ dunif(0, 1)
  p1 ~ dunif(0, 1)
  p2 ~ dunif(0, 1)
  #pi1 ~ dunif(0, 1)
  
  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phi1 * (1 - psi12)
  gamma[1,2] <- phi1 * psi12
  gamma[1,3] <- 1 - phi1
  gamma[2,1] <- phi2 * psi21
  gamma[2,2] <- phi2 * (1 - psi21)
  gamma[2,3] <- 1 - phi2
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- 1 - p1     # Pr(alive 1 t -> non-detected t)
  omega[1,2] <- p1         # Pr(alive 1 t -> detected 1 t)
  omega[1,3] <- 0          # Pr(alive 1 t -> detected 2 t)
  omega[2,1] <- 1 - p2     # Pr(alive 2 t -> non-detected t)
  omega[2,2] <- 0          # Pr(alive 2 t -> detected 1 t)
  omega[2,3] <- p2         # Pr(alive 2 t -> detected 2 t)
  omega[3,1] <- 1          # Pr(dead t -> non-detected t)
  omega[3,2] <- 0          # Pr(dead t -> detected 1 t)
  omega[3,3] <- 0          # Pr(dead t -> detected 2 t)
  
  # likelihood 
  # initial state probs
  for(i in 1:N) {
    init[i, 1:3] <- gamma[ y[i, first[i] ] - 1, 1:3 ]        # first state propagation
  }
  for (i in 1:N){
    y[i,(first[i]+1):K] ~ dHMM(init = init[i,1:3],           # count data from first[i] + 1
                               probObs = omega[1:3,1:3],     # observation matrix
                               probTrans = gamma[1:3,1:3],   # transition matrix
                               len = K - first[i],           # nb of occasions
                               checkRowSums = 0)             # do not check whether elements in a row sum up to 1
  }
})
# constants
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))
# data
my.data <- list(y = y + 1)
# initial values 
initial.values <- function(){list(phi1 = runif(1, 0, 1), 
                                  phi2 = runif(1, 0, 1), 
                                  psi12 = runif(1, 0, 1), 
                                  psi21 = runif(1, 0, 1), 
                                  p1 = runif(1, 0, 1), 
                                  p2 = runif(1, 0, 1))}
# parameters to monitor
parameters.to.save <- c("phi1", "phi2","psi12", "psi21", "p1", "p2")

# MCMC details
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# run NIMBLE
mcmc.multisite.marginalized <- nimbleMCMC(code = multisite.marginalized, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

# We obtain:
MCMCsummary(mcmc.multisite.marginalized, round = 2)

### Goodness of fit {#gofas}

# Going back to the geese data, the WBWA test is implemented as 
# follows on the whole dataset:
library(R2ucare)
geese <- read_csv2("allgeese.csv") %>% as.matrix()
y <- geese[,1:6]
size <- geese[,7]
wbwa <- test3Gwbwa(y, size)
wbwa$test3Gwbwa

## What if more than 2 sites?

### Dirichlet prior {#dirichletprior}

library(gtools) # to make the rdirichlet() function available
library(ggtern) # to visually represent multidim prob distribution 
set.seed(123) # for reproducibility
n <- 1000 # number of values drawn from Dirichlet distribution
alpha1 <- c(.1, .1, .1)
p1 <- rdirichlet(n, alpha1)
alpha2 <- c(1, 1, 1)
p2 <- rdirichlet(n, alpha2)
alpha3 <- c(10, 10, 10)
p3 <- rdirichlet(n, alpha3)
df <- cbind(rbind(p1, p2, p3), c(rep("alpha = c(0.1, 0.1, 0.1)", n),
                                     rep("alpha = c(1, 1, 1)", n),
                                     rep("alpha = c(10, 10, 10)", n))) %>%
  as_tibble() %>%
  mutate(x = as.numeric(V1),
         y = as.numeric(V2),
         z = as.numeric(V3),
         alpha = V4)

df %>%
  ggtern(aes(x = x, y = y, z = z)) +
  stat_density_tern(aes(fill=..level.., alpha=..level..),
                    geom = 'polygon',
                    bdl = 0.005) + # a 2D kernel density estimation of the distribution
  scale_fill_viridis_b() +
  geom_point(alpha = 0.3, pch = "+") +
  theme_showarrows() +
  scale_T_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = seq(0, 1, by = 0.2)) +
  scale_L_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = seq(0, 1, by = 0.2)) +
  scale_R_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = seq(0, 1, by = 0.2)) +
  labs(x = "",
       y = "",
       z = "",
       Tarrow = "psi11",
       Larrow = "psi12",
       Rarrow = "psi13") +
  guides(color = "none", fill = "none", alpha = "none") +
  facet_wrap(~alpha, ncol = 3)

# Going back to our example, in NIMBLE, we consider a Dirichlet 
# prior for each triplet of movement parameters, from site 
# 1 ($\psi^{11}$, $\psi^{12}$ and $\psi^{13}$), from site 2 
# ($\psi^{21}$, $\psi^{22}$ and $\psi^{23}$) and from site 
# 3 ($\psi^{31}$, $\psi^{32}$ and $\psi^{33}$). 

# Read in data
geese <- read_csv("geese.csv", col_names = TRUE)
y <- geese %>%
  as.matrix()

# Let us get the occasion of first capture for each individual
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)
first

# Multisite model
multisite <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phi1: survival probability site 1
  # phi2: survival probability site 2
  # phi3: survival probability site 3
  # psi11 = psi1[1]: movement probability from site 1 to site 1 (reference)
  # psi12 = psi1[2]: movement probability from site 1 to site 2
  # psi13 = psi1[3]: movement probability from site 1 to site 3 
  # psi21 = psi2[1]: movement probability from site 2 to site 1
  # psi22 = psi2[2]: movement probability from site 2 to site 2 (reference)
  # psi23 = psi2[3]: movement probability from site 2 to site 3
  # psi31 = psi3[1]: movement probability from site 3 to site 1
  # psi32 = psi3[2]: movement probability from site 3 to site 2
  # psi33 = psi3[3]: movement probability from site 3 to site 3 (reference)
  # p1: recapture probability site 1
  # p2: recapture probability site 2
  # p3: recapture probability site 3
  # -------------------------------------------------
  # States (z):
  # 1 alive at 1
  # 2 alive at 2
  # 2 alive at 3
  # 3 dead
  # Observations (y):  
  # 1 not seen
  # 2 seen at 1 
  # 3 seen at 2
  # 3 seen at 3
  
  # Priors
  phi1 ~ dunif(0, 1)
  phi2 ~ dunif(0, 1)
  phi3 ~ dunif(0, 1)
  p1 ~ dunif(0, 1)
  p2 ~ dunif(0, 1)
  p3 ~ dunif(0, 1)
  # transitions: Dirichlet priors
  psi1[1:3] ~ ddirch(alpha[1:3]) # psi11, psi12, psi13
  psi2[1:3] ~ ddirch(alpha[1:3]) # psi21, psi22, psi23
  psi3[1:3] ~ ddirch(alpha[1:3]) # psi31, psi32, psi33
  # Define state-transition and observation matrices
  # Define probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phi1 * psi1[1]
  gamma[1,2] <- phi1 * psi1[2]
  gamma[1,3] <- phi1 * psi1[3]
  gamma[1,4] <- 1 - phi1
  gamma[2,1] <- phi2 * psi2[1]
  gamma[2,2] <- phi2 * psi2[2]
  gamma[2,3] <- phi2 * psi2[3]
  gamma[2,4] <- 1 - phi2
  gamma[3,1] <- phi3 * psi3[1]
  gamma[3,2] <- phi3 * psi3[2]
  gamma[3,3] <- phi3 * psi3[3]
  gamma[3,4] <- 1 - phi3
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- 1
  
  # Define probabilities of y(t) given z(t)
  omega[1,1] <- 1 - p1     # Pr(alive A t -> non-detected t)
  omega[1,2] <- p1         # Pr(alive A t -> detected A t)
  omega[1,3] <- 0          # Pr(alive A t -> detected B t)
  omega[1,4] <- 0          # Pr(alive A t -> detected C t)
  omega[2,1] <- 1 - p2     # Pr(alive B t -> non-detected t)
  omega[2,2] <- 0          # Pr(alive B t -> detected A t)
  omega[2,3] <- p2         # Pr(alive B t -> detected B t)
  omega[2,4] <- 0          # Pr(alive B t -> detected C t)
  omega[3,1] <- 1 - p3     # Pr(alive C t -> non-detected t)
  omega[3,2] <- 0          # Pr(alive C t -> detected A t)
  omega[3,3] <- 0          # Pr(alive C t -> detected B t)
  omega[3,4] <- p3         # Pr(alive C t -> detected C t)
  omega[4,1] <- 1          # Pr(dead t -> non-detected t)
  omega[4,2] <- 0          # Pr(dead t -> detected A t)
  omega[4,3] <- 0          # Pr(dead t -> detected B t)
  omega[4,4] <- 0          # Pr(dead t -> detected C t)
  
  # Likelihood 
  for (i in 1:N){
    # Define latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1
    for (t in (first[i]+1):K){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
})

# data
my.data <- list(y = y + 1)

# constants
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y),
                     alpha = c(1, 1, 1))

# Initial values
zinits <- y
zinits[zinits==0] <- sample(c(1,2,3), sum(zinits==0), replace = TRUE)
initial.values <- function(){list(phi1 = runif(1, 0, 1), 
                                  phi2 = runif(1, 0, 1), 
                                  phi3 = runif(1, 0, 1), 
                                  psi1 = rep(1/3, 3),
                                  psi2 = rep(1/3, 3),
                                  psi3 = rep(1/3, 3),
                                  p1 = runif(1, 0, 1), 
                                  p2 = runif(1, 0, 1), 
                                  p3 = runif(1, 0, 1),
                                  z = zinits)}  

# MCMC settings
parameters.to.save <- c("phi1", "phi2", "phi3", "psi1", "psi2", "psi3","p1", "p2", "p3")
parameters.to.save

n.iter <- 10000
n.burnin <- 2500
n.chains <- 2

mcmc.multisite <- nimbleMCMC(code = multisite, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

MCMCsummary(mcmc.multisite, round = 2)

MCMCplot(mcmc.multisite)


# The distribution gamma($\alpha$,$\theta$) for different values 
# of $\alpha$ and $\theta$. 
# https://blog.revolutionanalytics.com/2015/10/parameters-and-percentiles-the-gamma-distribution.html
plotGamma <- function(shape=2, rate=0.5, to=0.99, p=c(0.1, 0.9), cex=1, ...){
  to <- qgamma(p = to, shape = shape, rate = rate)
  curve(dgamma(x, shape, rate), from = 0, to = to, n = 500, type="l", 
        main = sprintf("gamma(%1.1f, %1.1f)", shape, rate),
        bty = "n", xaxs = "i", yaxs = "i", col = "blue", xlab = "", ylab = "", 
        las = 1, lwd = 2, cex = cex, cex.axis = cex, cex.main = cex, ...)
}
par(mfrow = c(2, 3))
plotGamma(1, 0.1)
plotGamma(3, 0.1)
plotGamma(6, 0.1)
plotGamma(1, 1)
plotGamma(3, 1)
plotGamma(6, 1)

### Multinomial logit {#multinomiallogit}

geese <- read_csv("geese.csv", col_names = TRUE)
y <- geese %>%
  as.matrix()

# Let us get the occasion of first capture for each individual
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

# Multisite model
multisite <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phi1: survival probability site 1
  # phi2: survival probability site 2
  # phi3: survival probability site 3
  # psi11 = psi1[1]: movement probability from site 1 to site 1 (reference)
  # psi12 = psi1[2]: movement probability from site 1 to site 2
  # psi13 = psi1[3]: movement probability from site 1 to site 3 
  # psi21 = psi2[1]: movement probability from site 2 to site 1
  # psi22 = psi2[2]: movement probability from site 2 to site 2 (reference)
  # psi23 = psi2[3]: movement probability from site 2 to site 3
  # psi31 = psi3[1]: movement probability from site 3 to site 1
  # psi32 = psi3[2]: movement probability from site 3 to site 2
  # psi33 = psi3[3]: movement probability from site 3 to site 3 (reference)
  # p1: recapture probability site 1
  # p2: recapture probability site 2
  # p3: recapture probability site 3
  # -------------------------------------------------
  # States (z):
  # 1 alive at 1
  # 2 alive at 2
  # 2 alive at 3
  # 3 dead
  # Observations (y):  
  # 1 not seen
  # 2 seen at 1 
  # 3 seen at 2
  # 3 seen at 3
  # -------------------------------------------------
  
  # Priors
  phi1 ~ dunif(0, 1)
  phi2 ~ dunif(0, 1)
  phi3 ~ dunif(0, 1)
  p1 ~ dunif(0, 1)
  p2 ~ dunif(0, 1)
  p3 ~ dunif(0, 1)
  # transitions: multinomial logit
  # normal priors on logit of all but one transition probs
  for (i in 1:2){
    beta1[i] ~ dnorm(0, sd = 1.5)
    beta2[i] ~ dnorm(0, sd = 1.5)
    beta3[i] ~ dnorm(0, sd = 1.5)
  }
  # constrain the transitions such that their sum is < 1
  for (i in 1:2){
    psi1[i] <- exp(beta1[i]) / (1 + exp(beta1[1]) + exp(beta1[2]))
    psi2[i] <- exp(beta2[i]) / (1 + exp(beta2[1]) + exp(beta2[2]))
    psi3[i] <- exp(beta3[i]) / (1 + exp(beta3[1]) + exp(beta3[2]))
  }
  # last transition probability
  psi1[3] <- 1 - psi1[1] - psi1[2]
  psi2[3] <- 1 - psi2[1] - psi2[2]
  psi3[3] <- 1 - psi3[1] - psi3[2]
  # Define state-transition and observation matrices
  # Define probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phi1 * psi1[1]
  gamma[1,2] <- phi1 * psi1[2]
  gamma[1,3] <- phi1 * psi1[3]
  gamma[1,4] <- 1 - phi1
  gamma[2,1] <- phi2 * psi2[1]
  gamma[2,2] <- phi2 * psi2[2]
  gamma[2,3] <- phi2 * psi2[3]
  gamma[2,4] <- 1 - phi2
  gamma[3,1] <- phi3 * psi3[1]
  gamma[3,2] <- phi3 * psi3[2]
  gamma[3,3] <- phi3 * psi3[3]
  gamma[3,4] <- 1 - phi3
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- 1
  
  # Define probabilities of y(t) given z(t)
  omega[1,1] <- 1 - p1     # Pr(alive A t -> non-detected t)
  omega[1,2] <- p1         # Pr(alive A t -> detected A t)
  omega[1,3] <- 0          # Pr(alive A t -> detected B t)
  omega[1,4] <- 0          # Pr(alive A t -> detected C t)
  omega[2,1] <- 1 - p2     # Pr(alive B t -> non-detected t)
  omega[2,2] <- 0          # Pr(alive B t -> detected A t)
  omega[2,3] <- p2         # Pr(alive B t -> detected B t)
  omega[2,4] <- 0          # Pr(alive B t -> detected C t)
  omega[3,1] <- 1 - p3     # Pr(alive C t -> non-detected t)
  omega[3,2] <- 0          # Pr(alive C t -> detected A t)
  omega[3,3] <- 0          # Pr(alive C t -> detected B t)
  omega[3,4] <- p3         # Pr(alive C t -> detected C t)
  omega[4,1] <- 1          # Pr(dead t -> non-detected t)
  omega[4,2] <- 0          # Pr(dead t -> detected A t)
  omega[4,3] <- 0          # Pr(dead t -> detected B t)
  omega[4,4] <- 0          # Pr(dead t -> detected C t)
  
  # Likelihood 
  for (i in 1:N){
    # Define latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1
    for (t in (first[i]+1):K){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
})

# data
my.data <- list(y = y + 1)

# constants
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

# Initial values
zinits <- y
zinits[zinits==0] <- sample(c(1,2,3), sum(zinits==0), replace = TRUE)
initial.values <- function(){list(phi1 = runif(1, 0, 1), 
                                  phi2 = runif(1, 0, 1), 
                                  phi3 = runif(1, 0, 1), 
                                  beta1 = rnorm(2, 0, 10),
                                  beta2 = rnorm(2, 0, 10),
                                  beta3 = rnorm(2, 0, 10),
                                  p1 = runif(1, 0, 1), 
                                  p2 = runif(1, 0, 1), 
                                  p3 = runif(1, 0, 1),
                                  z = zinits)}  

# MCMC settings
parameters.to.save <- c("phi1", "phi2", "phi3", "psi1", "psi2", "psi3","p1", "p2", "p3")
parameters.to.save

n.iter <- 15000
n.burnin <- 2500
n.chains <- 2

mcmc.multisite <- nimbleMCMC(code = multisite, 
                             constants = my.constants,
                             data = my.data,              
                             inits = initial.values,
                             monitors = parameters.to.save,
                             niter = n.iter,
                             nburnin = n.burnin, 
                             nchains = n.chains)

MCMCsummary(mcmc.multisite, round = 2)

MCMCplot(mcmc.multisite)

## Sites may be states {#states}

### Titis data 

# You may see the data below: 
titis <- read_csv2("titis.csv", col_names = FALSE)
titis %>%
  rename(year_1942 = X1,
         year_1943 = X2,
         year_1944 = X3,
         year_1949 = X4,
         year_1952 = X5,
         year_1953 = X6,
         year_1956 = X7)

### The AS model for states

### NIMBLE implementation

y <- titis %>%
  as.matrix()

multistate <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phiB: survival probability state B
  # phiNB: survival probability state NB
  # psiBNB: transition probability from B to NB
  # psiNBB: transition probability from NB to B
  # pB: recapture probability B
  # pNB: recapture probability NB
  # -------------------------------------------------
  # States (S):
  # 1 alive B
  # 2 alive NB
  # 3 dead
  # Observations (O):  
  # 1 not seen
  # 2 seen as B 
  # 3 seen as NB
  # -------------------------------------------------
  
  # Priors
  phiB ~ dunif(0, 1)
  phiNB ~ dunif(0, 1)
  psiBNB ~ dunif(0, 1)
  psiNBB ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  pNB ~ dunif(0, 1)
  
  # Define state-transition and observation matrices
  
  # Define probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phiB * (1 - psiBNB)
  gamma[1,2] <- phiB * psiBNB
  gamma[1,3] <- 1 - phiB
  gamma[2,1] <- phiNB * psiNBB
  gamma[2,2] <- phiNB * (1 - psiNBB)
  gamma[2,3] <- 1 - phiNB
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 1
  
  # Define probabilities of y(t) given z(t)
  omega[1,1] <- 1 - pB    # Pr(alive B t -> non-detected t)
  omega[1,2] <- pB        # Pr(alive B t -> detected B t)
  omega[1,3] <- 0         # Pr(alive B t -> detected NB t)
  omega[2,1] <- 1 - pNB   # Pr(alive NB t -> non-detected t)
  omega[2,2] <- 0         # Pr(alive NB t -> detected B t)
  omega[2,3] <- pNB       # Pr(alive NB t -> detected NB t)
  omega[3,1] <- 1         # Pr(dead t -> non-detected t)
  omega[3,2] <- 0         # Pr(dead t -> detected N t)
  omega[3,3] <- 0         # Pr(dead t -> detected NB t)
  
  # Likelihood 
  for (i in 1:N){
    # Define latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1
    for (t in (first[i]+1):K){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:3])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(omega[z[i,t],1:3])
    }
  }
})


# data
my.data <- list(y = y + 1)

# constants
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)
first
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y))

# Initial values
zinits <- y
zinits[zinits == 0] <- sample(c(1,2), sum(zinits == 0), replace = TRUE)
initial.values <- function(){list(phiNB = runif(1, 0, 1), 
                                  phiB = runif(1, 0, 1), 
                                  psiNBB = runif(1, 0, 1), 
                                  psiBNB = runif(1, 0, 1), 
                                  pNB = runif(1, 0, 1), 
                                  pB = runif(1, 0, 1), 
                                  z = zinits)}  

# MCMC settings
parameters.to.save <- c("phiNB", "phiB","psiNBB", "psiBNB", "pNB", "pB")
parameters.to.save

n.iter <- 10000
n.burnin <- 1000
n.chains <- 2

mcmc.multistate <- nimbleMCMC(code = multistate, 
                              constants = my.constants,
                              data = my.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin, 
                              nchains = n.chains)

MCMCsummary(mcmc.multistate, round = 2)

MCMCplot(mcmc.multistate)

# First we gather the values generated for $\phi^B$ and $\phi^{NB}$ for the two chains:
phiB <- c(mcmc.multistate$chain1[,"phiB"], mcmc.multistate$chain2[,"phiB"])
phiNB <- c(mcmc.multistate$chain1[,"phiNB"], mcmc.multistate$chain2[,"phiNB"])
df <- data.frame(param = c(rep("phiB", length(phiB)), 
                           rep("phiNB", length(phiB))), 
                 value = c(phiB, phiNB))

# Then, we plot the two posterior distributions:
df %>%
  ggplot(aes(x = value, fill = param)) +
  geom_density(color = "white", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(fill = "", x = "survival")

# What about a potential trade-off on reproduction?
psiBNB <- c(mcmc.multistate$chain1[,"psiBNB"], mcmc.multistate$chain2[,"psiBNB"])
psiBB <- 1 - psiBNB
psiNBB <- c(mcmc.multistate$chain1[,"psiNBB"], mcmc.multistate$chain2[,"psiNBB"])
df <- data.frame(param = c(rep("psiBB", length(phiB)), 
                           rep("psiNBB", length(phiB))), 
                 value = c(psiBB, psiNBB))
df %>%
  ggplot(aes(x = value, fill = param)) +
  geom_density(color = "white", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(fill = "", x = "breeding probabilities")

## Issue of local minima {#localminima}

dat <- matrix(c(2, 0, 2, 1, 2, 0, 2, 
                2, 0, 2, 1, 2, 0, 2,
                2, 0, 2, 1, 2, 0, 2,
                2, 0, 2, 1, 2, 0, 2,
                1, 1, 1, 0, 1, 0, 1,
                1, 1, 1, 0, 1, 0, 1,
                1, 1, 1, 0, 1, 0, 1,
                1, 1, 1, 0, 1, 0, 1,
                2, 0, 2, 0, 2, 0, 1,
                2, 0, 2, 0, 2, 0, 1,
                2, 0, 2, 0, 2, 0, 1,
                2, 0, 2, 0, 2, 0, 1,
                1, 0, 1, 0, 1, 0, 1,
                1, 0, 1, 0, 1, 0, 1,
                1, 0, 1, 0, 1, 0, 1,
                1, 0, 1, 0, 1, 0, 1,
                2, 0, 2, 0, 2, 0, 2,
                2, 0, 2, 0, 2, 0, 2,
                2, 0, 2, 0, 2, 0, 2,
                2, 0, 2, 0, 2, 0, 2,
                1, 0, 1, 0, 1, 0, 2,
                1, 0, 1, 0, 1, 0, 2,
                1, 0, 1, 0, 1, 0, 2,
                1, 0, 1, 0, 1, 0, 2,
                2, 2, 0, 1, 0, 2, 1,
                2, 2, 0, 1, 0, 2, 1,
                2, 2, 0, 1, 0, 2, 1,
                2, 2, 0, 1, 0, 2, 1,
                2, 1, 0, 2, 0, 1, 1,
                2, 1, 0, 2, 0, 1, 1,
                2, 1, 0, 2, 0, 1, 1,
                2, 1, 0, 2, 0, 1, 1),
              byrow = T,
              ncol = 7)

# Let us get the occasion of first capture for each individual
get.first <- function(x) min(which(x != 0))
first <- apply(dat, 1, get.first)

# Multistate model
multistate <- nimbleCode({
  
  # -------------------------------------------------
  # Parameters:
  # phi: survival probability
  # psiAB: movement probability from site A to site B
  # psiBA: movement probability from site B to site A
  # p: recapture probability
  # -------------------------------------------------
  # States (S):
  # 1 alive at A
  # 2 alive at B
  # 3 dead
  # Observations (O):  
  # 1 not seen
  # 2 seen at A 
  # 3 seen at B
  # -------------------------------------------------
  
  # Priors
  phi ~ dunif(0, 1)
  psi12 ~ dunif(0, 1)
  psi21 ~ dunif(0, 1)
  p ~ dunif(0, 1)
  
  # Define state-transition and observation matrices
  
  # Define probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phi * (1 - psi12)
  gamma[1,2] <- phi * psi12
  gamma[1,3] <- 1 - phi
  gamma[2,1] <- phi * psi21
  gamma[2,2] <- phi * (1 - psi21)
  gamma[2,3] <- 1 - phi
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 1
  
  # Define probabilities of y(t) given z(t)
  omega[1,1] <- 1 - p     # Pr(alive A t -> non-detected t)
  omega[1,2] <- p         # Pr(alive A t -> detected A t)
  omega[1,3] <- 0         # Pr(alive A t -> detected B t)
  omega[2,1] <- 1 - p     # Pr(alive B t -> non-detected t)
  omega[2,2] <- 0         # Pr(alive B t -> detected A t)
  omega[2,3] <- p         # Pr(alive B t -> detected B t)
  omega[3,1] <- 1         # Pr(dead t -> non-detected t)
  omega[3,2] <- 0         # Pr(dead t -> detected A t)
  omega[3,3] <- 0         # Pr(dead t -> detected B t)
  
  # Likelihood 
  for (i in 1:N){
    # Define latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1
    for (t in (first[i]+1):K){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:3])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(omega[z[i,t],1:3])
    }
  }
})

# data
my.data <- list(y = dat + 1)

# constants
my.constants <- list(first = first, 
                     K = ncol(dat), 
                     N = nrow(dat))

# Initial values
zinits <- dat
zinits[zinits == 0] <- sample(c(1,2), sum(zinits == 0), replace = TRUE)
initial.values <- function(){list(phi = runif(1, 0, 1), 
                                  psi12 = runif(1, 0, 1), 
                                  psi21 = runif(1, 0, 1), 
                                  p = runif(1, 0, 1), 
                                  z = zinits)}  

# MCMC settings
parameters.to.save <- c("phi", "psi12", "psi21", "p")
parameters.to.save

n.iter <- 10000
n.burnin <- 2500
n.chains <- 2

mcmc.multisite <- nimbleMCMC(code = multistate, 
                              constants = my.constants,
                              data = my.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin, 
                              nchains = n.chains)

MCMCsummary(mcmc.multisite, round = 2)

psiBA <- c(mcmc.multisite$chain1[,"psi21"], mcmc.multisite$chain2[,"psi21"])
plotpsiBA <- psiBA %>%
  as_tibble() %>%
  ggplot() +
  aes(x = 1:length(value), y = value) +
  geom_line(color = "gray70") +
  labs(y = "psi21", x = "iterations") +
  geom_hline(yintercept = 0.85, lty = 2, color = "blue")

psiAB <- c(mcmc.multisite$chain1[,"psi12"],mcmc.multisite$chain2[,"psi12"])
plotpsiAB <- psiAB %>%
  as_tibble() %>%
  ggplot() +
  aes(x = 1:length(value), y = value) +
  geom_line(color = "gray70") +
  labs(y = "psi12", x = "iterations") +
  geom_hline(yintercept = 0.6, lty = 2, color = "blue")

det <- c(mcmc.multisite$chain1[,"p"], mcmc.multisite$chain2[,"p"])
plotdet <- det %>%
  as_tibble() %>%
  ggplot() +
  aes(x = 1:length(value), y = value) +
  geom_line(color = "gray70") +
  labs(y = "detection", x = "iterations") +
  geom_hline(yintercept = 0.6, lty = 2, color = "blue")

library(patchwork)
plotdet / (plotpsiAB + plotpsiBA)

plotpsiBA <- psiBA %>%
  as_tibble() %>%
  ggplot() +
  aes(x = value) +
  geom_histogram(color = "white", binwidth = .03, fill = "gray70") +
  geom_density(aes(y = .03 * ..count..)) +
  labs(x = "psi21", y = "") +
  geom_vline(xintercept = 0.85, lty = 2, color = "blue")

plotpsiAB <- psiAB %>%
  as_tibble() %>%
  ggplot() +
  aes(x = value) +
  geom_histogram(color = "white", binwidth = .03, fill = "gray70") +
  geom_density(aes(y = .03 * ..count..)) +
  labs(x = "psi12", y = "") +
  geom_vline(xintercept = 0.6, lty = 2, color = "blue")

plotdet <- det %>%
  as_tibble() %>%
  ggplot() +
  aes(x = value) +
  geom_histogram(color = "white", binwidth = .03, fill = "gray70") +
  geom_density(aes(y = .03 * ..count..)) +
  labs(x = "detection", y = "") +
  geom_vline(xintercept = 0.6, lty = 2, color = "blue")

plotdet / (plotpsiAB + plotpsiBA)

## Uncertainty {#multievent}

### Breeding states {#breedingmultievent}

titis <- read_csv2("titis_with_uncertainty.csv", col_names = FALSE)
y <- titis %>%
  as.matrix()

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

n.iter <- 15000
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

MCMCsummary(mcmc.multievent, round = 2)

MCMCplot(mcmc.multievent)


### Disease states {#diseasemultievent}


# model building
y <- as.matrix(read.table("hofi_2000.txt"))

model <- nimbleCode({
  
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
set.seed(250375)
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

n.iter <- 30000
n.burnin <- 2500
n.chains <- 2

init1 <- initial.values()
init2 <- initial.values()

# create model as an R object (uncompiled model)
survival <- nimbleModel(code = model,
                        data = my.data,
                        constants = my.constants,
                        inits = list(init1, init2))
# compile model
Csurvival <- compileNimble(survival)
# create a MCMC configuration
survivalConf <- configureMCMC(survival)
# create a MCMC function and compile it
survivalMCMC <- buildMCMC(survivalConf)
CsurvivalMCMC <- compileNimble(survivalMCMC, project = survival)

survivalConf$printSamplers()

# Now that we have control on the MCMC configuration, let's mess it up. 
# We start by removing the default sampler:
survivalConf$removeSamplers(c('phiH','phiI','psiIH','psiHI','betaH','betaI'))
survivalConf$printSamplers()

# add block samplers for coupled params
survivalConf$addSampler(target = c('phiH','phiI'), type = 'RW_block')
survivalConf$addSampler(target = c('psiIH','psiHI'), type = 'RW_block')
survivalConf$addSampler(target = c('betaH','betaI'), type = 'RW_block')

survivalConf$printSamplers()

# Now you can resume the workflow:

# create a new MCMC function and compile it:
survivalMCMC2 <- buildMCMC(survivalConf)
CsurvivalMCMC2 <- compileNimble(survivalMCMC2, 
                                project = survival,
                                resetFunctions = TRUE) 
# run NIMBLE:
samples2 <- runMCMC(mcmc = CsurvivalMCMC2, 
                    niter = n.iter,
                    nburnin = n.burnin,
                    nchains = n.chains)

# obtain numerical summaries:
MCMCvis::MCMCsummary(samples2)

## Summary

## Suggested reading

