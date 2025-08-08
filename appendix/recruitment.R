

# Read the file
lines <- readLines("dat/flamingo.txt")

# Keep only the first 8 characters and split them into individual bits with spaces
processed <- sapply(substr(lines, 1, 18), function(x) paste(strsplit(x, "")[[1]], collapse = " "))

# Save the result
writeLines(processed, "cleaned_data.txt")


#------------------------------------------------------------------------------------

# NO trap-dependence
# add transience w/ age effect on survival

# read in data
y <- as.matrix(read.table(here::here("dat", "calonectris_diomedea.txt")))

# OBS
# 1 not seen
# 2 seen

# STATES
# 1 alive
# 2 dead

library(nimble)

hmm <- nimbleCode({
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    for (t in first[i]:(K-1)){
      phi[i,t] <- beta[age[i,t]] # beta1 = phi1, beta2 = phi1+
      gamma[1,1,i,t] <- phi[i,t]      # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t]  # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0           # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1           # Pr(dead t -> dead t+1)
    }
  }
  beta[1] ~ dunif(0, 1) # phi1
  beta[2] ~ dunif(0, 1) # phi1+
  phi1 <- beta[1]
  phi2 <- beta[2]
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# age effect via an individual by time covariate and use nested indexing 
# to distinguish survival over the interval after first detection from survival afterwards: 
age <- matrix(NA, nrow = nrow(y), ncol = ncol(y) - 1)
for (i in 1:nrow(age)){
  for (j in 1:ncol(age)){
    if (j == first[i]) age[i,j] <- 1 # age = 1
    if (j > first[i]) age[i,j] <- 2  # age > 1
  }
}

first <- apply(y, 1, function(x) min(which(x != 0)))
my.constants <- list(N = nrow(y), K = ncol(y), first = first, age = age)
my.data <- list(y = y + 1)

zinits <- y
zinits[zinits == 0] <- 1
initial.values <- function() list(beta = runif(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

parameters.to.save <- c("phi1", "phi2", "p")

n.iter <- 2500
n.burnin <- 500
n.chains <- 2

out <- nimbleMCMC(code = hmm,
                  constants = my.constants,
                  data = my.data,
                  inits = initial.values,
                  monitors = parameters.to.save,
                  niter = n.iter,
                  nburnin = n.burnin,
                  nchains = n.chains)

library(MCMCvis)

MCMCsummary(out, round = 2)

MCMCtrace(out, pdf = FALSE)

save(out, file = "calonectris-wo-trap-dep.RData")



