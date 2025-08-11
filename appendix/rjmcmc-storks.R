
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
  for (j in 1:10){
    beta[j] ~ dnorm(0, sd = 1.5) # prior slopes
  }
  intercept ~ dnorm(0, sd = 1.5) # prior intercept
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
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
            control = list(mean = 0, scale = 2))

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

#----- post prob of presence / The Marginal Posterior Probability of Inclusion
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

#----- model-averaged effects, average beta * ksi across draws (you already compute betaksi in the model).

# combine chains
K <- 10

beta  <- out[, paste0("beta[", 1:K, "]"), drop = FALSE]
ksi   <- out[, paste0("ksi[", 1:K, "]"),  drop = FALSE]
bstar <- beta * ksi  # model-averaged betas on the logit scale

summ_bstar <- t(apply(bstar, 2, function(v) c(mean = mean(v),
                                              q2.5 = quantile(v, .025),
                                              q97.5 = quantile(v, .975),
                                              P_gt0 = mean(v > 0))))
round(summ_bstar, 3)

#----- model-averaged survival

invlogit <- function(x) 1 / (1 + exp(-x))

# design matrix used in the model
X <- rainfall[1:(16 - 1), 1:K, drop = FALSE]

# linear predictor for every time t and every draw
eta <- X %*% t(bstar)                     # (T-1) x ndraw
eta <- sweep(eta, 2, post[, "intercept"], "+")
phi_draws <- invlogit(eta)                # (T-1) x ndraw

# summarize by time
phi_summary <- t(apply(phi_draws, 1, function(v) c(mean = mean(v),
                                                   q2.5 = quantile(v, .025),
                                                   q97.5 = quantile(v, .975))))
round(phi_summary, 3)



#----- get model post prob
K <- 10
ksimat <- out[, paste0("ksi[", 1:K, "]")]

# encode each draw's model as a string like "0101001001"
model_id <- apply(ksimat, 1, paste0, collapse = "")

# posterior model probabilities (empirical frequencies)
post_model_prob <- sort(prop.table(table(model_id)), decreasing = TRUE)

# Top models
head(data.frame(post_model_prob), 10)

# model_id    Freq
# 1  0001000000 0.75470
# 2  0001010000 0.03965
# 3  0001100000 0.02910
# 4  0000000000 0.02310
# 5  0001000001 0.02055
# 6  0001000100 0.01710
# 7  0001000010 0.01660
# 8  1001000000 0.01575
# 9  0001001000 0.01320
# 10 0101000000 0.01045


save(out, file = "rjmcmc.RData")

# Bayes factor

a <- 1; b <- 1; p <- 10
prior_prob_for_k <- function(k, a=1, b=1, p=10) beta(k + a, p - k + b) / beta(a, b)

# vector of prior probs for each possible size k
prior_size <- sapply(0:p, prior_prob_for_k, a=a, b=b, p=p)

# prior prob for a specific model_id string
prior_model_prob <- function(model_id, a=1, b=1) {
  k <- sum(strsplit(model_id, "")[[1]] == "1")
  beta(k + a, length(model_id) - k + b) / beta(a, b)
}

BF_vs_null <- function(m, m0 = rep(0, K), a=1, b=1) {
  id  <- paste0(m, collapse = "")
  id0 <- paste0(m0, collapse = "")
  post <- as.numeric(post_model_prob[id])
  post0 <- as.numeric(post_model_prob[id0])
  prior <- prior_model_prob(id, a, b)
  prior0 <- prior_model_prob(id0, a, b)
  (post / post0) / (prior / prior0)
}

# example: BF of the top model vs. the null (all zeros)
top_id <- names(post_model_prob)[1]
BF_vs_null(strsplit(top_id, "")[[1]]) # best VS null model (intercept only)

