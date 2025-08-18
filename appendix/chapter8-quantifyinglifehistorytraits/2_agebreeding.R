# Quantifying age-specific breeding

# OBS
# data from 2004 to 2016
# 1 = not seen
# 2 = seen alone
# 3 = seen with a young of the year
# 4 = seen with a calf (1 to 3 years old)

# STATES
# NB Nonbreeding adult female
# Byoy Breeding adult female with a young‐of‐the‐year
# Bc1 Breeding adult female with a 1‐year‐old calf 
# Bc3 Breeding adult female with a 3‐year‐old calf
# Bc2 Breeding adult female with a 2‐year‐old calf
# D Dead female

# INTERMEDIATE STATES
# Byoy‐D Breeding adult female that had lost her young‐of‐the‐year
# Bc1‐D Breeding adult female that had lost her 1‐year‐ old calf
# Bc2‐D Breeding adult female that had lost her 2‐year‐ old calf
# Bc3‐leave Breeding adult female that raised her calf to the age of 3

library(nimble)
library(MCMCvis)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# read in data
dolphins <- read_csv2("dolphins.csv", col_names = TRUE)
y <- dolphins %>%
  as.matrix()

# model code
hmm <- nimbleCode({
  
  # priors
  youngsurvival[1] ~ dunif(0, 1)
  youngsurvival[2] ~ dunif(0, 1)
  adultsurvival ~ dunif(0, 1)
  beta[1] ~ dunif(0, 1)
  beta[2] ~ dunif(0, 1)
  p[1] ~ dunif(0, 1) # detection prob
  p[2] ~ dunif(0, 1) # detection prob
  q[1] ~ dunif(0, 1) # obs prob
  q[2] ~ dunif(0, 1) # obs prob
  
  pi[1:5] ~ ddirch(alpha[1:5])

  # vector of initial stats probs
  delta[1] <- pi[1] 
  delta[2] <- pi[2] 
  delta[3] <- pi[3] 
  delta[4] <- pi[4] 
  delta[5] <- pi[5] 
  delta[6] <- 0 # prob. of being in initial state dead
  
  # probabilities of state z(t+1) given z(t)
  
  # PHIA is 6x6
  for (i in 1:5) {
    PHIA[i, i] <- adultsurvival
    PHIA[i, 6] <- 1 - adultsurvival
  }
  PHIA[6, 6] <- 1
  
  # PHIY is 6x9
  # Row 1: NB stays NB
  PHIY[1, 1] <- 1
  
  # Row 2: Byoy
  PHIY[2, 2] <- youngsurvival[1]        # YOY survives → becomes Bc1
  PHIY[2, 3] <- 1 - youngsurvival[1]    # YOY dies     → becomes Byoy-D
  
  # Row 3: Bc1
  PHIY[3, 4] <- youngsurvival[1]        # 1-year-old survives → Bc2
  PHIY[3, 5] <- 1 - youngsurvival[1]    # dies → Bc1-D
  
  # Row 4: Bc2
  PHIY[4, 6] <- youngsurvival[2]        # 2-year-old survives → Bc3-L (emancipation)
  PHIY[4, 7] <- 1 - youngsurvival[2]    # dies → Bc2-D
  
  # Row 5: Bc3 always emancipates
  PHIY[5, 8] <- 1
  
  # Row 6: dead female stays dead
  PHIY[6, 9] <- 1
  
  # AGING is 9x9

  # Deterministic transitions
  AGING[1, 1] <- 1   # NB → NB
  AGING[2, 2] <- 1   # Byoy → Bc1
  AGING[3, 5] <- 1   # Byoy-D → Byoy-D
  AGING[4, 3] <- 1   # Bc1 → Bc2
  AGING[5, 6] <- 1   # Bc1-D → Bc1-D
  AGING[6, 4] <- 1   # Bc2 → Bc3
  AGING[7, 7] <- 1   # Bc2-D → Bc2-D
  AGING[8, 8] <- 1   # Bc3-L → Bc3-L
  AGING[9, 9] <- 1   # D → D
  
  # BREED is 9x6

  # Row 1: NB → Byoy or NB
  BREED[1, 2] <- beta[1]
  BREED[1, 1] <- 1 - beta[1]
  
  # Row 2: Bc1 → Bc1
  BREED[2, 3] <- 1
  
  # Row 3: Bc2 → Bc2
  BREED[3, 4] <- 1
  
  # Row 4: Bc3 → Bc3
  BREED[4, 5] <- 1
  
  # Row 5: Byoy-D → Byoy or NB
  BREED[5, 2] <- beta[1]
  BREED[5, 1] <- 1 - beta[1]
  
  # Row 6: Bc1-D → Byoy or NB
  BREED[6, 2] <- beta[1]
  BREED[6, 1] <- 1 - beta[1]
  
  # Row 7: Bc2-D → Byoy or NB
  BREED[7, 2] <- beta[2]
  BREED[7, 1] <- 1 - beta[2]
  
  # Row 8: Bc3-L → Byoy or NB
  BREED[8, 2] <- beta[2]
  BREED[8, 1] <- 1 - beta[2]
  
  # Row 9: D → D
  BREED[9, 6] <- 1
  
  gamma[1:6,1:6] <- PHIA[1:6,1:6] %*% PHIY[1:6,1:9] %*% AGING[1:9,1:9] %*% BREED[1:9,1:6]
  
  # probabilities of y(t) given z(t)
  
  # DETECTION is 6x6

  # detection probabilities
  DETECTION[1, 1] <- 1 - p[1]     # NB not detected
  DETECTION[1, 2] <- p[1]         # NB detected
  
  DETECTION[2, 1] <- 1 - p[2]
  DETECTION[2, 3] <- p[2]        # Byoy detected
  
  DETECTION[3, 1] <- 1 - p[2]
  DETECTION[3, 4] <- p[2]         # Bc1 detected
  
  DETECTION[4, 1] <- 1 - p[1]
  DETECTION[4, 5] <- p[1]         # Bc2 detected
  
  DETECTION[5, 1] <- 1 - p[1]
  DETECTION[5, 6] <- p[1]         # Bc3 detected
  
  DETECTION[6, 1] <- 1           # Dead never detected
  
  # OBS is 6x4

  # Not detected → not sighted
  OBS[1, 1] <- 1
  
  # Detected NB → alone
  OBS[2, 2] <- 1
  
  # Detected Byoy → YOY visible or not
  OBS[3, 2] <- 1 - q[1]    # alone
  OBS[3, 3] <- q[1]        # with YOY
  
  # Detected Bc1 → calf visible or not
  OBS[4, 2] <- 1 - q[2]
  OBS[4, 4] <- q[2]
  
  # Detected Bc2
  OBS[5, 2] <- 1 - q[2]
  OBS[5, 4] <- q[2]
  
  # Detected Bc3
  OBS[6, 2] <- 1 - q[2]
  OBS[6, 4] <- q[2]
  
  omega[1:6,1:4] <- DETECTION[1:6,1:6] %*% OBS[1:6,1:4]
  
  # DETECTION is 6x6
  
  # detection probabilities
  DETECTION.init[1, 1] <- 0     # NB not detected
  DETECTION.init[1, 2] <- 1         # NB detected
  
  DETECTION.init[2, 1] <- 0
  DETECTION.init[2, 3] <- 1        # Byoy detected
  
  DETECTION.init[3, 1] <- 0
  DETECTION.init[3, 4] <- 1         # Bc1 detected
  
  DETECTION.init[4, 1] <- 0
  DETECTION.init[4, 5] <- 1         # Bc2 detected
  
  DETECTION.init[5, 1] <- 0
  DETECTION.init[5, 6] <- 1         # Bc3 detected
  
  DETECTION.init[6, 1] <- 1           # Dead never detected

  omega.init[1:6,1:4] <- DETECTION.init[1:6,1:6] %*% OBS[1:6,1:4]
  
  # Likelihood 
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:6])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]],1:4])
    for (t in (first[i]+1):K){
      z[i,t] ~ dcat(gamma[z[i,t-1],1:6])
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
})

# constants
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

PHIA <- matrix(0, 6, 6)
PHIY <- matrix(0, 6, 9)
AGING <- matrix(0, 9, 9)
BREED <- matrix(0, 9, 6)
DETECTION <- matrix(0, 6, 6)
OBS <- matrix(0, 6, 4)
DETECTION.init <- matrix(0, 6, 6)

my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y),
                     alpha = rep(1,5),
                     PHIA = PHIA,
                     PHIY = PHIY,
                     AGING = AGING,
                     BREED = BREED,
                     DETECTION = DETECTION,
                     DETECTION.init = DETECTION.init,
                     OBS = OBS)

# data
my.data <- list(y = y + 1)

# MCMC settings
parameters.to.save <- c("adultsurvival", 
                        "youngsurvival", 
                        "beta", 
                        "pi", 
                        "p", 
                        "q")

n.iter <- 15000
n.burnin <- 5000
n.chains <- 2

# initial values
start_params <- list(
  adultsurvival = 0.90,
  youngsurvival = c(0.70, 0.80),
  beta = c(0.40, 0.30),
  p = c(0.70, 0.80),
  q = c(0.70, 0.80),
  # initial state probabilities over 5 live states (NB,Byoy,Bc1,Bc2,Bc3)
  pi = c(0.50, 0.15, 0.10, 0.10, 0.15)
)

build_matrices <- function(par) {
  # --- PHIA 6x6 ---
  PHIA <- matrix(0, 6, 6)
  for(i in 1:5){ PHIA[i,i] <- par$adultsurvival; PHIA[i,6] <- 1 - par$adultsurvival }
  PHIA[6,6] <- 1
  
  # --- PHIY 6x9 ---
  PHIY <- matrix(0, 6, 9)
  PHIY[1,1] <- 1
  PHIY[2,2] <- par$youngsurvival[1]; PHIY[2,3] <- 1 - par$youngsurvival[1]
  PHIY[3,4] <- par$youngsurvival[1]; PHIY[3,5] <- 1 - par$youngsurvival[1]
  PHIY[4,6] <- par$youngsurvival[2]; PHIY[4,7] <- 1 - par$youngsurvival[2]
  PHIY[5,8] <- 1
  PHIY[6,9] <- 1
  
  # --- AGING 9x9 (deterministic, correct wiring) ---
  AGING <- matrix(0, 9, 9)
  AGING[1,1] <- 1   # NB -> NB
  AGING[2,2] <- 1   # Byoy -> Bc1
  AGING[3,5] <- 1   # Byoy-D -> Byoy-D
  AGING[4,3] <- 1   # Bc1 -> Bc2
  AGING[5,6] <- 1   # Bc1-D -> Bc1-D
  AGING[6,4] <- 1   # Bc2 -> Bc3-L
  AGING[7,7] <- 1   # Bc2-D -> Bc2-D
  AGING[8,8] <- 1   # Bc3-L -> Bc3-L
  AGING[9,9] <- 1   # D -> D
  
  # --- BREED 9x6 ---
  BREED <- matrix(0, 9, 6)
  BREED[1,2] <- par$beta[1]; BREED[1,1] <- 1 - par$beta[1]
  BREED[2,3] <- 1
  BREED[3,4] <- 1
  BREED[4,5] <- 1
  BREED[5,2] <- par$beta[1]; BREED[5,1] <- 1 - par$beta[1]
  BREED[6,2] <- par$beta[1]; BREED[6,1] <- 1 - par$beta[1]
  BREED[7,2] <- par$beta[2]; BREED[7,1] <- 1 - par$beta[2]
  BREED[8,2] <- par$beta[2]; BREED[8,1] <- 1 - par$beta[2]
  BREED[9,6] <- 1
  
  # --- γ 6x6 ---
  gamma <- PHIA %*% PHIY %*% AGING %*% BREED
  gamma <- sweep(gamma, 1, pmax(rowSums(gamma), .Machine$double.eps), "/")
  
  # --- DETECTION 6x6 ---
  DET <- matrix(0, 6, 6)
  DET[1,1] <- 1 - par$p[1]; DET[1,2] <- par$p[1]
  DET[2,1] <- 1 - par$p[2]; DET[2,3] <- par$p[2]
  DET[3,1] <- 1 - par$p[2]; DET[3,4] <- par$p[2]
  DET[4,1] <- 1 - par$p[1]; DET[4,5] <- par$p[1]
  DET[5,1] <- 1 - par$p[1]; DET[5,6] <- par$p[1]
  DET[6,1] <- 1
  
  # --- OBS 6x4 ---
  OBS <- matrix(0, 6, 4)
  OBS[1,1] <- 1
  OBS[2,2] <- 1
  OBS[3,2] <- 1 - par$q[1]; OBS[3,3] <- par$q[1]
  OBS[4,2] <- 1 - par$q[2]; OBS[4,4] <- par$q[2]
  OBS[5,2] <- 1 - par$q[2]; OBS[5,4] <- par$q[2]
  OBS[6,2] <- 1 - par$q[2]; OBS[6,4] <- par$q[2]
  
  # --- ω 6x4 ---
  omega <- DET %*% OBS
  omega <- pmax(omega, 0)  # guard numerical noise
  omega <- sweep(omega, 1, pmax(rowSums(omega), .Machine$double.eps), "/")
  
  list(gamma = gamma, omega = omega)
}

viterbi_path <- function(obs, delta, gamma, omega) {
  Tlen <- length(obs)
  S <- nrow(gamma)
  # tiny floor to avoid -Inf in logs for init only
  eps <- 1e-12
  gL <- log(pmax(gamma, eps))
  oL <- log(pmax(omega, eps))
  dL <- log(pmax(delta, eps))
  
  V <- matrix(-Inf, nrow = S, ncol = Tlen)
  ptr <- matrix(NA_integer_, nrow = S, ncol = Tlen)
  
  # init
  V[,1] <- dL + oL[, obs[1]]
  ptr[,1] <- 0L
  
  # forward
  for(t in 2:Tlen){
    for(s in 1:S){
      cand <- V[, t-1] + gL[, s]
      j <- which.max(cand)
      V[s, t] <- cand[j] + oL[s, obs[t]]
      ptr[s, t] <- j
    }
  }
  # backtrack
  z <- integer(Tlen)
  z[Tlen] <- which.max(V[, Tlen])
  for(t in (Tlen-1):1){
    z[t] <- ptr[z[t+1], t+1]
  }
  z
}

make_delta <- function(pi5){
  c(pi5[1:5], 0)
}

# generate init values with Viterbi
generate_inits_viterbi <- function(Y, first, par){
  mats <- build_matrices(par)
  gamma <- mats$gamma
  omega <- mats$omega
  
  N <- nrow(Y); K <- ncol(Y)
  z <- matrix(NA_integer_, N, K)
  delta <- make_delta(par$pi)
  
  for(i in 1:N){
    t0 <- first[i]
    obs_i <- Y[i, t0:K]
    z_i <- viterbi_path(obs_i, delta, gamma, omega)
    z[i, t0:K] <- z_i
  }
  
  list(
    adultsurvival = par$adultsurvival,
    youngsurvival = par$youngsurvival,
    beta = par$beta,
    pi = par$pi,
    p = par$p,
    q = par$q,
    z = z
  )
}

inits <- generate_inits_viterbi(my.data$y, my.constants$first, start_params)

# run NIMBLE
out <- nimbleMCMC(code = hmm, 
                  constants = my.constants,
                  inits = inits,
                  data = my.data,
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin, 
                  nchains = n.chains)

# inspect results
MCMCsummary(out, round = 2)

