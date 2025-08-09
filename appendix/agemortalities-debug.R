


# Read lines from file
lines <- readLines("dat/roedeer.txt")

# Remove trailing semicolon
lines <- gsub(";", "", lines)

# Split into components
split_lines <- strsplit(lines, " ")

# Number of rows (individuals)
n <- length(split_lines)

# Initialize ind vector
ind <- character(n)

# Get length of obs (assume all same length)
obs_length <- nchar(split_lines[[1]][1])

# Initialize matrix for obs
obs_matrix <- matrix(NA, nrow = n, ncol = obs_length)

# Initialize sex vector
# Fill matrix and sex vector
for (i in seq_along(split_lines)) {
  obs_str <- split_lines[[i]][1]

  # Fill observation matrix with digits
  obs_matrix[i, ] <- as.integer(strsplit(obs_str, "")[[1]])
  ind[i] <- as.integer(split_lines[[i]][6])

    # Determine sex
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

library(nimble)


hmm <- nimbleCode({

  # Priors
  mu[1:3] ~ ddirch(alpha[1:3])
  p ~ dunif(0, 1) # detection prob
  ll ~ dunif(0, 1) # recovery prob

  # vector of initial stats probs
  delta[1] <- 1
  delta[2] <- 0
  delta[3] <- 0
  delta[4] <- 0

  # probabilities of state z(t+1) given z(t)

  gamma[1,1] <- mu[3]
  gamma[1,2] <- mu[1]
  gamma[1,3] <- mu[2]
  gamma[1,4] <- 0
  gamma[2,1] <- 0
  gamma[2,2] <- 0
  gamma[2,3] <- 0
  gamma[2,4] <- 1
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 0
  gamma[3,4] <- 1
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- 1

    # probabilities of y(t) given z(t)

  omega[1,1] <- 1 - p
  omega[1,2] <- p
  omega[1,3] <- 0
  omega[2,1] <- 1 - ll
  omega[2,2] <- 0
  omega[2,3] <- ll
  omega[3,1] <- 1
  omega[3,2] <- 0
  omega[3,3] <- 0
  omega[4,1] <- 1
  omega[4,2] <- 0
  omega[4,3] <- 0

  omega.init[1,1] <- 0
  omega.init[1,2] <- 1
  omega.init[1,3] <- 0
  omega.init[2,1] <- 0
  omega.init[2,2] <- 0
  omega.init[2,3] <- 1
  omega.init[3,1] <- 1
  omega.init[3,2] <- 0
  omega.init[3,3] <- 0
  omega.init[4,1] <- 1
  omega.init[4,2] <- 0
  omega.init[4,3] <- 0

  # Likelihood
  for (i in 1:N){
    # Define latent state at first capture
    z[i,first[i]] ~ dcat(delta[1:4])
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]],1:3])
    for (t in (first[i]+1):K){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(omega[z[i,t],1:3])
    }
  }
})

# constants
get.first <- function(x) min(which(x != 0))
first <- apply(obs_matrix, 1, get.first)

my.constants <- list(first = first,
                     K = ncol(obs_matrix),
                     N = nrow(obs_matrix),
                     alpha = rep(1,3))

# data
my.data <- list(y = obs_matrix + 1)

# MCMC settings
parameters.to.save <- c("mu",
                        "p",
                        "ll")

n.iter <- 1500
n.burnin <- 500
n.chains <- 2

#-- initial values

# simple_inits <- list(
#   mu = c(0.30, 0.20, 0.50),  # mu[1]=Pr(A->H), mu[2]=Pr(A->NH), mu[3]=Pr(A->A)
#   p  = 0.60,                 # detection given A
#   ll = 0.80                  # recovery given H
# )
# 
# 
# z_init <- matrix(NA, nrow = my.constants$N, ncol = my.constants$K)
# Y <- my.data$y  # values in {1,2,3}; remember 1=not seen, 2=captured, 3=recovered dead (human)
# N <- my.constants$N
# K <- my.constants$K
# for (i in 1:nrow(obs_matrix)) {
#   t0 <- first[i]
#   z_init[i, t0] <- 1L  # start Alive at first capture (omega_init forces y=2 there)
#   dead <- FALSE
#   if (t0 < K) {
#     for (t in (t0+1):K) {
#       if (dead) {
#         z_init[i, t] <- 4L               # once dead -> remain Dead
#       } else if (Y[i, t] == 3L) {        # observed recovered dead
#         z_init[i, t] <- 2L               # H at the death occasion
#         if (t < K) z_init[i, (t+1):K] <- 4L
#         dead <- TRUE
#       } else if (Y[i, t] == 2L) {        # captured/recaptured
#         z_init[i, t] <- 1L               # must be Alive
#       } else {                           # not seen
#         z_init[i, t] <- 1L               # could be A/NH/D; pick A pre-death to avoid zeros
#       }
#     }
#   }
# }
# 
# simple_inits$z <- z_init

#-- initial values

start_params <- list(
  mu = c(0.30, 0.20, 0.50),  # mu[1]=Pr(A->H), mu[2]=Pr(A->NH), mu[3]=Pr(A->A)
  p  = 0.60,                 # detection given A
  ll = 0.80                  # recovery given H
)

build_mats_roedeer <- function(par){
  # --- transition gamma (4x4): states = A,H,NH,D = 1..4
  gamma <- matrix(0, 4, 4)
  gamma[1,1] <- par$mu[3]; gamma[1,2] <- par$mu[1]; gamma[1,3] <- par$mu[2]; gamma[1,4] <- 0
  gamma[2,4] <- 1
  gamma[3,4] <- 1
  gamma[4,4] <- 1

  # --- emission omega (4x3): obs = 1:not seen, 2:captured/recaptured, 3:found dead human-caused
  omega <- matrix(0, 4, 3)
  omega[1,1] <- 1 - par$p;  omega[1,2] <- par$p;    omega[1,3] <- 0      # A
  omega[2,1] <- 1 - par$ll; omega[2,2] <- 0;        omega[2,3] <- par$ll # H
  omega[3,1] <- 1;          omega[3,2] <- 0;        omega[3,3] <- 0      # NH
  omega[4,1] <- 1;          omega[4,2] <- 0;        omega[4,3] <- 0      # D

  # --- initial emission omega.init (4x3), if you want to keep it
  omega_init <- matrix(0, 4, 3)
  omega_init[1,2] <- 1  # first ever capture must be "captured" if A
  omega_init[2,3] <- 1  # if H at first time, must be "human-killed"
  omega_init[3,1] <- 1  # NH -> not seen at first time
  omega_init[4,1] <- 1  # D  -> not seen at first time

  # --- initial state probs delta (length 4): your model uses A with prob 1
  delta <- c(1,0,0,0)

  list(gamma = gamma, omega = omega, omega_init = omega_init, delta = delta)
}

viterbi_path <- function(obs_vec, delta, gamma, omega, omega_init){
  Tlen <- length(obs_vec); S <- nrow(gamma)
  eps <- 1e-12
  gL <- log(pmax(gamma, eps))
  oL <- log(pmax(omega, eps))
  dL <- log(pmax(delta, eps))
  o0L <- log(pmax(omega_init, eps))

  V   <- matrix(-Inf, nrow = S, ncol = Tlen)
  ptr <- matrix(NA_integer_, nrow = S, ncol = Tlen)

  # init emission
  if(is.null(omega_init)){
    V[,1] <- dL + oL[, obs_vec[1]]
  } else {
    V[,1] <- dL + o0L[, obs_vec[1]]
  }
  ptr[,1] <- 0L

  # forward
  for(t in 2:Tlen){
    for(s in 1:S){
      cand <- V[, t-1] + gL[, s]
      j <- which.max(cand)
      V[s, t] <- cand[j] + oL[s, obs_vec[t]]
      ptr[s, t] <- j
    }
  }
  # backtrack
  z <- integer(Tlen)
  z[Tlen] <- which.max(V[, Tlen])
  if(Tlen >= 2){
    for(t in (Tlen-1):1){
      z[t] <- ptr[z[t+1], t+1]
    }
  }
  z
}

generate_inits_viterbi_roedeer <- function(Y, first, par, use_init_emission = TRUE){
  mats <- build_mats_roedeer(par)
  gamma <- mats$gamma; omega <- mats$omega; omega_init <- mats$omega_init; delta <- mats$delta
  N <- nrow(Y); K <- ncol(Y)
  z <- matrix(NA_integer_, N, K)

  for(i in 1:N){
    t0 <- first[i]
    obs_i <- Y[i, t0:K]  # values in 1..3
    if(use_init_emission){
      zi <- viterbi_path(obs_i, delta, gamma, omega, omega_init)
    } else {
      zi <- viterbi_path(obs_i, delta, gamma, omega, omega_init = NULL)
    }
    z[i, t0:K] <- zi
  }

  list(
    mu = par$mu,
    p  = par$p,
    ll = par$ll,
    z  = z
  )
}

simple_inits <- generate_inits_viterbi_roedeer(my.data$y, my.constants$first, start_params,
                                        use_init_emission = TRUE)

model <- nimbleModel(code = hmm,
                     constants = my.constants,
                     inits = simple_inits,
                     data = my.data,
                     check = FALSE)

# quick sanity checks
rowSums(model$gamma[,])    # should be ~1
rowSums(model$omega[,])    # should be ~1
model$calculate()

out <- nimbleMCMC(code = hmm,
                  constants = my.constants,
                  inits = simple_inits,
                  data = my.data,
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin,
                  nchains = n.chains)

library(MCMCvis)

MCMCsummary(out, round = 2)


#----- sans valeurs initiales VITERBI

# > MCMCsummary(out, round = 2)
# mean   sd 2.5%  50% 97.5% Rhat n.eff
# ll    0.97 0.03 0.91 0.98  1.00 1.06   144
# mu[1] 0.01 0.00 0.00 0.01  0.01 1.01   893
# mu[2] 0.00 0.00 0.00 0.00  0.00 1.05   468
# mu[3] 0.99 0.00 0.99 0.99  1.00 1.02   870
# p     0.09 0.00 0.08 0.09  0.10 1.00   327
# 



#----- avec valeurs initiales VITERBI

# mean   sd 2.5%  50% 97.5% Rhat n.eff
# ll    0.57 0.23 0.23 0.54  0.96 1.35    13
# mu[1] 0.05 0.02 0.02 0.04  0.10 1.39    17
# mu[2] 0.25 0.02 0.20 0.25  0.29 1.39    17
# mu[3] 0.70 0.01 0.68 0.70  0.72 1.01   815
# p     0.52 0.01 0.49 0.52  0.55 1.01   365
# 
