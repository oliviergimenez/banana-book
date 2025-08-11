

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


