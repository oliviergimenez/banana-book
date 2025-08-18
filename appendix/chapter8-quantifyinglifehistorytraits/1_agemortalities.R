# Assessing age- and cause-specific mortalities

library(nimble)
library(MCMCvis)
library(tidyverse)
theme_set(theme_light(base_size = 14))

# read lines from file
lines <- readLines("roedeer.txt")

# remove trailing semicolon
lines <- gsub(";", "", lines)

# split into components
split_lines <- strsplit(lines, " ")

# number of rows (individuals)
n <- length(split_lines)

# initialize ind vector
ind <- character(n)

# get length of obs (assume all same length)
obs_length <- nchar(split_lines[[1]][1])

# initialize matrix for obs
obs_matrix <- matrix(NA, nrow = n, ncol = obs_length)

# initialize sex vector
# fill matrix and sex vector
for (i in seq_along(split_lines)) {
  obs_str <- split_lines[[i]][1]
  
  # fill observation matrix with digits
  obs_matrix[i, ] <- as.integer(strsplit(obs_str, "")[[1]])
  ind[i] <- as.integer(split_lines[[i]][6])
  
  # determine sex
  if (ind[i] == 1) {
    ind[i] <- "uncensored"
  } else {
    ind[i] <- "censored"
  }
}

head(obs_matrix)
nrow(obs_matrix) # number of female roe deer (of known age, marked as fawn)
ncol(obs_matrix) # number of capture occasions, 1985 to 2013
2013-1985+1

sum(obs_matrix==2) # number of human-related mortalities
# deadly injured during handling, 
# victim of car collisions, or 
# recovered and reported by hunters

ind
sum(ind == "censored")
sum(ind == "uncensored")

# OBS
# 0, not seen; 
# 1, captured for the first time or recaptured; 
# 2, killed by human activities and reported. 

# STATES
# A, alive; 
# H, individual just died from human causes; 
# NH, individual just died from a natural (non-human) cause; 
# D, individual already dead.

hmm <- nimbleCode({
  
  # ----- priors -----
  # A -> H (human-caused): two constants by age class
  alpha_H_01  ~ dnorm(0, sd = 2)   # age < 1
  alpha_H_ge1 ~ dnorm(0, sd = 2)   # age >= 1
  
  # A -> NH (natural): piecewise (0-1), (1-2), then linear >=2
  alpha_NH_01 ~ dnorm(0, sd = 2)   # 0 <= age < 1
  alpha_NH_12 ~ dnorm(0, sd = 2)   # 1 <= age < 2
  beta0_NH    ~ dnorm(0, sd = 2)   # intercept for age >= 2
  beta1_NH    ~ dnorm(0, sd = 2)   # slope for age >= 2  (on logit scale)
  
  # recovery constant
  ll ~ dunif(0, 1)
  
  # detection p: time-dependent, 2 age classes (<=1 vs >1), on logit scale
  for(t in 1:K){
    lp_ageLE1[t] ~ dnorm(0, sd = 2)
    lp_ageGT1[t] ~ dnorm(0, sd = 2)
    p_ageLE1[t] <- ilogit(lp_ageLE1[t])
    p_ageGT1[t] <- ilogit(lp_ageGT1[t])
  }
  
  # initial state probabilities (start alive)
  delta[1] <- 1; delta[2] <- 0; delta[3] <- 0; delta[4] <- 0
  
  # first-time obs 
  omega_init[1,1] <- 0; omega_init[1,2] <- 1; omega_init[1,3] <- 0   # A: captured at first time
  omega_init[2,1] <- 0; omega_init[2,2] <- 0; omega_init[2,3] <- 1   # H: recovered
  omega_init[3,1] <- 1; omega_init[3,2] <- 0; omega_init[3,3] <- 0   # NH: not seen
  omega_init[4,1] <- 1; omega_init[4,2] <- 0; omega_init[4,3] <- 0   # D:  not seen
  
  # ----- Build time- and individual-specific transitions & observations -----

  for(i in 1:N){
    for(t in 1:(K-1)) {
      
      # age effects
      etaH[i,t] <- (alpha_H_01 * age_lt1[i,t] + alpha_H_ge1 * age_ge1[i,t])
      
      etaNH[i,t] <- (alpha_NH_01 * age_lt1[i,t] +
                       alpha_NH_12 * age_1to2[i,t] +
                       (beta0_NH + beta1_NH * age_std_mu2[i,t]) * age_ge2[i,t])
      
      # denominator 
      den[i,t] <- 1 + exp(etaH[i,t]) + exp(etaNH[i,t])
      
      # transitions from A
      gamma[i,t,1,1] <- 1 / den[i,t]  # default: stay alive
      gamma[i,t,1,2] <- exp(etaH[i,t]) / den[i,t]   
      gamma[i,t,1,3] <- exp(etaNH[i,t]) / den[i,t]  
      gamma[i,t,1,4] <- 0  
      
      # from H, NH, D (always go to D)
      gamma[i,t,2,1] <- 0; gamma[i,t,2,2] <- 0; gamma[i,t,2,3] <- 0; gamma[i,t,2,4] <- 1
      gamma[i,t,3,1] <- 0; gamma[i,t,3,2] <- 0; gamma[i,t,3,3] <- 0; gamma[i,t,3,4] <- 1
      gamma[i,t,4,1] <- 0; gamma[i,t,4,2] <- 0; gamma[i,t,4,3] <- 0; gamma[i,t,4,4] <- 1
      
      # detection probability 
      pA[i,t] <- p_ageLE1[t] * age_lt1[i,t] + p_ageGT1[t] * age_ge1[i,t]
      
      # observation for A
      omega[i,t,1,1] <- (1 - pA[i,t])  # Default: not seen
      omega[i,t,1,2] <- pA[i,t]  # Default: 0
      omega[i,t,1,3] <- 0  # Always 0
      
      # observation for H
      omega[i,t,2,1] <- 1 - ll
      omega[i,t,2,2] <- 0
      omega[i,t,2,3] <- ll
      
      # observation for NH and D
      omega[i,t,3,1] <- 1; omega[i,t,3,2] <- 0; omega[i,t,3,3] <- 0
      omega[i,t,4,1] <- 1; omega[i,t,4,2] <- 0; omega[i,t,4,3] <- 0
    }
  }
  
  # ----- Likelihood -----
  for (i in 1:N) {
    z[i, first[i]] ~ dcat(delta[1:4])
    y[i, first[i]] ~ dcat(omega_init[z[i, first[i]], 1:3])
    
    if (first[i] < K) {
      for (t in (first[i] + 1):K) {
        z[i, t] ~ dcat(gamma[i, t-1, z[i, t-1], 1:4])
        y[i, t] ~ dcat(omega[i, t-1, z[i, t], 1:3])
      }
    }
  }
})

# constants
get.first <- function(x) min(which(x != 0))
first <- apply(obs_matrix, 1, get.first)

# years of the study 
years <- 1985:2013
K <- length(years)

# build an N x K age matrix
age <- matrix(0, nrow = nrow(obs_matrix), ncol = K)
for (i in 1:nrow(obs_matrix)) {
  for (t in first[i]:K) {
    age[i, t] <- t - first[i]   # age 0 at first capture
  }
}

# standardize age (center and scale by SD) over non-NA values
age_mean <- mean(age, na.rm = TRUE)
age_sd   <- sd(age, na.rm = TRUE)
age_std  <- (age - age_mean) / age_sd

# for mu[2] effect only from age >= 2
age_std_mu2 <- age_std
age_std_mu2[age < 3] <- 0  # slope kicks in only from age >= 2

age_ge2 <- age_1to2 <- matrix(0, nrow(obs_matrix), K)
age_lt1 <- age_ge1 <- matrix(0, nrow(obs_matrix), K)

for (i in 1:nrow(obs_matrix)) {
  for (t in first[i]:K) {
    age_val <- age[i, t]
    if (age_val == 0) {
      age_lt1[i, t] <- 1
      age_ge1[i, t] <- 0
    } else {
      age_lt1[i, t] <- 0  
      age_ge1[i, t] <- 1
    }
  }
}

for (i in 1:nrow(obs_matrix)) {
  for (t in first[i]:K) {
    age_val <- age[i, t]
    age_ge2[i, t] <- as.integer(age_val > 1)
    age_1to2[i, t] <- as.integer(age_val == 1)
  }
}

my.constants <- list(
  N = nrow(obs_matrix),
  K = ncol(obs_matrix),
  first = first,
  age_std_mu2 = age_std_mu2,
  age_lt1  = age_lt1,
  age_ge1  = age_ge1,
  age_ge2  = age_ge2,
  age_1to2 = age_1to2
)

# data
my.data <- list(y = obs_matrix + 1)

# MCMC settings
parameters.to.save <- c("alpha_H_01",
                        "alpha_H_ge1",
                        "alpha_NH_01",
                        "alpha_NH_12",
                        "beta0_NH",
                        "beta1_NH",
                        "ll",
                        "lp_ageLE1",
                        "lp_ageGT1",
                        "p_ageLE1",
                        "p_ageGT1")

# MCMC settings
n.iter <- 10000
n.burnin <- 2500
n.chains <- 2

ilogit <- function(x) 1/(1+exp(-x))

viterbi_timeinhom <- function(obs, delta, gamma_list, omega_list, omega_init){
  Tlen <- length(obs); S <- length(delta); eps <- 1e-12
  logp <- function(x) log(pmax(x, eps))
  V <- matrix(-Inf, S, Tlen); ptr <- matrix(NA_integer_, S, Tlen)
  V[,1] <- logp(delta) + logp(omega_init[, obs[1]])
  ptr[,1] <- 0L
  for(t in 2:Tlen){
    gL <- logp(gamma_list[[t-1]]); oL <- logp(omega_list[[t-1]])
    for(s in 1:S){
      cand <- V[, t-1] + gL[, s]
      j <- which.max(cand)
      V[s, t] <- cand[j] + oL[s, obs[t]]
      ptr[s, t] <- j
    }
  }
  z <- integer(Tlen); z[Tlen] <- which.max(V[, Tlen])
  if(Tlen >= 2) for(t in (Tlen-1):1) z[t] <- ptr[z[t+1], t+1]
  z
}

# draw parameter starts
set.seed(1)
lp_ageLE1 <- rnorm(K, 0, 1)
lp_ageGT1 <- rnorm(K, 0, 1)
p_ageLE1  <- ilogit(lp_ageLE1)
p_ageGT1  <- ilogit(lp_ageGT1)

par0 <- list(
  alpha_H_01 = rnorm(1, 0, 1),
  alpha_H_ge1 = rnorm(1, 0, 1),
  alpha_NH_01 = rnorm(1, 0, 1),
  alpha_NH_12 = rnorm(1, 0, 1),
  beta0_NH = rnorm(1, 0, 1),
  beta1_NH = rnorm(1, 0, 1),
  ll = runif(1, 0, 1)
)

# fixed pieces
delta <- c(1,0,0,0)
omega_init <- rbind(
  c(0,1,0),
  c(0,0,1),
  c(1,0,0),
  c(1,0,0)
)

# Viterbi init for z
N <- nrow(obs_matrix); Y <- my.data$y
z_init <- matrix(NA_integer_, N, K)

for(i in 1:N){
  t0 <- first[i]
  if(t0 <= K){
    gamma_list <- vector("list", max(0, K - t0))
    omega_list <- vector("list", max(0, K - t0))
    if(t0 < K){
      for(t in t0:(K-1)){
        # transitions from A via multinomial logit
        etaH  <- par0$alpha_H_01 * age_lt1[i,t] + par0$alpha_H_ge1 * age_ge1[i,t]
        etaNH <- par0$alpha_NH_01 * age_lt1[i,t] + par0$alpha_NH_12 * age_1to2[i,t] +
          (par0$beta0_NH + par0$beta1_NH * age_std_mu2[i,t]) * age_ge2[i,t]
        den <- 1 + exp(etaH) + exp(etaNH)
        G <- matrix(0, 4, 4)
        G[1,1] <- 1/den; G[1,2] <- exp(etaH)/den; G[1,3] <- exp(etaNH)/den; G[1,4] <- 0
        G[2,4] <- 1; G[3,4] <- 1; G[4,4] <- 1
        gamma_list[[t - t0 + 1]] <- G
        
        # observations at time t+1 use slice t
        pA <- p_ageLE1[t] * age_lt1[i,t] + p_ageGT1[t] * age_ge1[i,t]
        O <- matrix(0, 4, 3)
        O[1,1] <- 1 - pA; O[1,2] <- pA; O[1,3] <- 0
        O[2,1] <- 1 - par0$ll; O[2,2] <- 0; O[2,3] <- par0$ll
        O[3,1] <- 1; O[4,1] <- 1
        omega_list[[t - t0 + 1]] <- O
      }
    }
    obs_i <- Y[i, t0:K]
    zi <- viterbi_timeinhom(obs_i, delta, gamma_list, omega_list, omega_init)
    z_init[i, t0:K] <- zi
  }
}

# pack inits for nimble
simple_inits <- list(
  alpha_H_01 = par0$alpha_H_01, alpha_H_ge1 = par0$alpha_H_ge1,
  alpha_NH_01 = par0$alpha_NH_01, alpha_NH_12 = par0$alpha_NH_12,
  beta0_NH = par0$beta0_NH, beta1_NH = par0$beta1_NH,
  ll = par0$ll,
  lp_ageLE1 = lp_ageLE1,
  lp_ageGT1 = lp_ageGT1,
  z = z_init
)

# run NIMBLE
out <- nimbleMCMC(code = hmm, 
                  constants = my.constants,
                  inits = simple_inits,
                  data = my.data,
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin, 
                  nchains = n.chains)

# get results
MCMCsummary(out, round = 2) 
# neff for ll, the alpha's and tne beta's are low
# would need to re-run w/ more iterations

# get MCMC draws
draws <- rbind(out[[1]], out[[2]])

# parameters
pars <- c("alpha_H_01","alpha_H_ge1","alpha_NH_01","alpha_NH_12","beta0_NH","beta1_NH")

# age grids
ages_plot <- 1:max(age)           
ages_plot <- 1:18

# standardization as in our script
age_std_grid     <- (ages_plot - age_mean) / age_sd
age_std_mu2_grid <- ifelse(ages_plot >= 3, age_std_grid, 0)

# extract draw vectors
aH01  <- draws[, "alpha_H_01"]
aHge1 <- draws[, "alpha_H_ge1"]
aNH01 <- draws[, "alpha_NH_01"]
aNH12 <- draws[, "alpha_NH_12"]
b0    <- draws[, "beta0_NH"]
b1    <- draws[, "beta1_NH"]

# prepare empty tibble to store results
df_results <- tibble(
  age = numeric(),
  metric = character(),
  median = numeric(),
  l95 = numeric(),
  u95 = numeric()
)

# loop over ages
for (k in seq_along(ages_plot)) {
  a <- ages_plot[k]
  
  # eta_H(age): with our age definition (≥1), this is alpha_H_ge1
  etaH <- aHge1
  
  # eta_NH(age): piecewise
  if (a < 1) {
    etaNH <- aNH01
  } else if (a == 2) {
    etaNH <- aNH12
  } else if (a >= 3) {
    etaNH <- b0 + b1 * age_std_mu2_grid[k]
  } else {
    # a == 1 → treat like 1–2 piece
    etaNH <- aNH12
  }
  
  # probabilities
  den  <- 1 + exp(etaH) + exp(etaNH)
  muA  <- 1 / den
  muH  <- exp(etaH)  / den
  muNH <- exp(etaNH) / den
  
  # summary function
  qfun <- function(x) c(median = median(x), l95 = quantile(x, 0.025), u95 = quantile(x, 0.975))
  sH   <- qfun(muH)
  sNH  <- qfun(muNH)
  
  # append rows for this age
  df_results <- bind_rows(
    df_results,
    tibble(
      age = a,
      metric = "Mortality: human",
      median = sH[1], l95 = sH[2], u95 = sH[3]
    ),
    tibble(
      age = a,
      metric = "Mortality: natural",
      median = sNH[1], l95 = sNH[2], u95 = sNH[3]
    )
  )
}

ggplot(df_results, aes(age, median)) +
  #  geom_ribbon(aes(ymin = l95, ymax = u95, fill = metric), alpha = 0.2) +
  geom_line(aes(linetype = metric)) +
  scale_y_continuous(limits = c(0,1), name = "Probability") +
  scale_x_continuous(breaks = pretty(ages_plot)) +
  labs(x = "Age (years)", linetype = NULL, fill = NULL,
       title = "Roe deer: age-specific and cause-specific mortalities") +
  theme_bw()


