


#------------------------------------------------------------------------------------

# stopover

# Read the file as a data frame
dat <- read.table("dat/newt-stopover.txt", header = TRUE, stringsAsFactors = FALSE)

# Get number of columns
ncol_total <- ncol(dat)

# Extract y matrix (all columns except the last)
y <- as.matrix(dat[, 1:(ncol_total - 1)])

# Convert y to integer (if read as character)
y <- apply(y, 2, as.integer)

# Extract sex vector (last column)
sex <- dat[[ncol_total]]

# Check result
y
sex

sex <- as.numeric(as.factor(sex))

dim(y)

# We have 12 occasions, including one first occasion with 
# only zeros to create a state in which all animals are outside the site. 

# To assess when individuals arrive at the pond, 
# the model is conditional to the first occasion, and not to the first capture. 
# 
# Three STATES 
# + not yet arrived in the pond (1)
# + arrived in the pond (2)
# + departed from the pond (3)]
# 
# Two OBS
# + not captured (1)
# + captured (2)] 
# 
# Individuals were considered initially out of the pond in the state 
# ‘not yet arrived’. Transition matrices included two probabilities: 
# the probability of arriving (r, and its complement: not arriving), 
# and the probability of staying (f, and its complement: leaving the pond). 
# The events matrix has only one probability as individuals are only 
# capturable when present in the pond, and are not capturable when they 
# have not yet arrived or have already departed.

library(nimble)

hmm <- nimbleCode({
  p ~ dunif(0, 1) # prior probability of being captured if present on the site
  
  for (j in 1:(K-1)){
    logit(r[1,j]) <- beta[1,1] + beta[1,2] * (j-4.5)/2.45 # prior probability of arriving in the pond / 1 = female
    logit(r[2,j]) <- beta[2,1] + beta[2,2] * (j-4.5)/2.45 # prior probability of arriving in the pond / 2 = male
    s[1,j] ~ dunif(0, 1) # prior retention probability (staying) / 1 = female
    s[2,j] ~ dunif(0, 1) # prior retention probability (staying) / 2 = male
  }

  beta[1,1] ~ dnorm(0, sd = 1.5)
  beta[1,2] ~ dnorm(0, sd = 1.5)
  beta[2,1] ~ dnorm(0, sd = 1.5)
  beta[2,2] ~ dnorm(0, sd = 1.5)
  
  for (i in 1:N){
    for (j in 1:(K-1)){
      gamma[1,1,i,j] <- 1 - r[sex[i],j]    # Pr(alive t -> non-detected t)
      gamma[1,2,i,j] <- r[sex[i],j]        # Pr(alive t -> detected t)
      gamma[1,3,i,j] <- 0        # Pr(alive t -> detected t)
      gamma[2,1,i,j] <- 0        # Pr(dead t -> non-detected t)
      gamma[2,2,i,j] <- s[sex[i],j]        # Pr(dead t -> detected t)
      gamma[2,3,i,j] <- 1 - s[sex[i],j]        # Pr(dead t -> detected t)
      gamma[3,1,i,j] <- 0        # Pr(dead t -> non-detected t)
      gamma[3,2,i,j] <- 0        # Pr(dead t -> detected t)
      gamma[3,3,i,j] <- 1        # Pr(dead t -> detected t)
    }
  }
  
  omega[1,1] <- 1      # Pr(alive t -> alive t+1)
  omega[1,2] <- 0  # Pr(alive t -> dead t+1)
  omega[2,1] <- 1 - p           # Pr(dead t -> alive t+1)
  omega[2,2] <- p           # Pr(dead t -> dead t+1)
  omega[3,1] <- 1           # Pr(dead t -> alive t+1)
  omega[3,2] <- 0           # Pr(dead t -> dead t+1)
  
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 1          # Pr(alive t = 1) = 1
  delta[3] <- 0          # Pr(dead t = 1) = 0
  
  # likelihood
  for (i in 1:N){
    z[i,1] ~ dcat(delta[1:3])
    for (j in 2:K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3, i, j - 1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

my.constants <- list(N = nrow(y), K = ncol(y), sex = sex)
my.data <- list(y = y + 1)

# Function to generate initial values for latent states z
generate_zinits_simple <- function(y) {
  N <- nrow(y)
  K <- ncol(y)
  z <- matrix(NA, nrow = N, ncol = K)
  
  for (i in 1:N) {
    # All start in state 1 (not yet arrived)
    z[i, 1] <- 1
    
    # Find capture occasions (y[i,j] == 1 in original 0/1 data)
    captures <- which(y[i, ] == 1)
    
    if (length(captures) > 0) {
      first_capture <- min(captures)
      last_capture <- max(captures)
      
      # Strategy: create a simple but valid trajectory
      # Before first capture: state 1 (not yet arrived)
      if (first_capture > 2) {
        z[i, 2:(first_capture-1)] <- 1
      }
      
      # From first to last capture: state 2 (present)
      # This ensures all captures occur when individual is present
      z[i, first_capture:last_capture] <- 2
      
      # After last capture: state 3 (departed) 
      if (last_capture < K) {
        z[i, (last_capture+1):K] <- 3
      }
    } else {
      # Never captured: stay in state 1 (not yet arrived)
      z[i, 2:K] <- 1
    }
  }
  
  return(z)
}
initial.values <- function() list(beta = matrix(rnorm(4,0,1), nrow = 2),
                                  s = matrix(runif(2*(ncol(y)-1),0,1), nrow = 2),
                                  p = runif(1,0,1),
                                  z = generate_zinits_simple(y))

parameters.to.save <- c("s", "r", "p")

n.iter <- 10000
n.burnin <- 5000
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

# compute stopover duration

save(out, file = "stopover.RData")

library(MCMCvis)
load(here::here("dat","stopover.RData"))

samps <- rbind(out[[1]], out[[2]])

library(tidyverse)

# samps: matrix from rbind(out[[1]], out[[2]])
df <- as.data.frame(samps) %>% mutate(iter = row_number())

# Helpers to tidy blocks like r[i,j] and s[i,j]
tidy_block <- function(df, prefix){
  df %>%
    select(iter, matches(paste0("^", prefix, "\\["))) %>%
    pivot_longer(-iter, names_to = "param", values_to = "value") %>%
    mutate(
      i = as.integer(str_match(param, paste0(prefix, "\\[(\\d+),"))[,2]),
      j = as.integer(str_match(param, ",\\s*(\\d+)\\]")[,2]),
      sex = ifelse(i == 1, "Female", "Male"),
      week = j
    ) %>%
    select(iter, sex, week, value)
}

r_long <- tidy_block(df, "r")
s_long <- tidy_block(df, "s")

# -------- r: arrival probability, both sexes on one figure --------
r_sum <- r_long %>%
  group_by(sex, week) %>%
  summarise(
    mean = mean(value),
    median = median(value),
    l95 = quantile(value, 0.025),
    u95 = quantile(value, 0.975),
    .groups = "drop"
  )

p_r <- ggplot(r_sum, aes(x = week, y = mean, color = sex, fill = sex, group = sex)) +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.20, linewidth = 0) +
  geom_line(size = 1) +
  labs(x = "Week", y = "Arrival probability (posterior mean & 95% CI)",
       title = "Arrival probability r by sex") +
  theme_bw() +
  theme(legend.position = "top")
p_r



# ---- s -> expected stopover duration from each arrival week ----
# duration[w] = 1 + sum_{k=1..K-w} prod_{u=w..w+k-1} s[u]
# (stops at K-1 so staying past last week doesn't add observable time)
compute_duration_curve <- function(s_vec){
  s_vec <- pmin(s_vec, 0.999999)        # guard against s=1 exactly
  K <- length(s_vec)
  dur <- numeric(K)
  for(w in 1:K){
    if(w == K){
      dur[w] <- 1
    } else {
      dur[w] <- 1 + sum(cumprod(s_vec[w:(K-1)]))
    }
  }
  dur
}

K_weeks <- max(s_long$week)

duration_draws <- s_long %>%
  arrange(iter, sex, week) %>%
  group_by(iter, sex) %>%
  summarise(
    duration = list(compute_duration_curve(value)),
    .groups = "drop"
  ) %>%
  mutate(week = list(seq_len(K_weeks))) %>%
  unnest(c(week, duration))

duration_sum <- duration_draws %>%
  group_by(sex, week) %>%
  summarise(
    mean = mean(duration),
    median = median(duration),
    l95 = quantile(duration, 0.025),
    u95 = quantile(duration, 0.975),
    .groups = "drop"
  )

p_dur <- ggplot(duration_sum, aes(x = week, y = mean, color = sex, fill = sex, group = sex)) +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.20, linewidth = 0) +
  geom_line(size = 1) +
  labs(x = "Week", y = "Stopover duration (weeks)",
       title = "Expected stopover duration from arrival week (posterior mean & 95% CI)") +
  theme_bw() +
  theme(legend.position = "top")
p_dur


# a linear date effect on arrival probability associated with a sex effect, 
# and a time and sex effect on retention probability
# constant detection


# 
# # Read lines from file
# lines <- readLines("dat/data_2013.txt")
# 
# # Remove trailing semicolon
# lines <- gsub(";", "", lines)
# 
# # Split into components
# split_lines <- strsplit(lines, " ")
# 
# # Number of rows (individuals)
# n <- length(split_lines)
# 
# # Get length of obs (assume all same length)
# obs_length <- nchar(split_lines[[1]][1])
# 
# # Initialize matrix for obs
# obs_matrix <- matrix(NA, nrow = n, ncol = obs_length)
# 
# # Initialize sex vector
# sex <- character(n)
# 
# # Fill matrix and sex vector
# for (i in seq_along(split_lines)) {
#   obs_str <- split_lines[[i]][1]
#   sex1 <- as.integer(split_lines[[i]][2])
#   sex2 <- as.integer(split_lines[[i]][3])
#   
#   # Fill observation matrix with digits
#   obs_matrix[i, ] <- as.integer(strsplit(obs_str, "")[[1]])
#   
#   # Determine sex
#   if (sex1 == 1 && sex2 == 0) {
#     sex[i] <- "female"
#   } else if (sex1 == 0 && sex2 == 1) {
#     sex[i] <- "male"
#   } else {
#     sex[i] <- "unknown"
#   }
# }
# 
# # Combine into final data.frame
# result <- data.frame(obs_matrix, sex = sex, stringsAsFactors = FALSE)
# 
# # View result
# print(result)
# 
# # Optionally write to file
# write.table(result, "newt-stopover.txt", quote = FALSE, row.names = FALSE)
# 
# sum(sex == "male")
# sum(sex == "female")
