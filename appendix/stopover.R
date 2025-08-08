


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
  delta[2] <- 0          # Pr(dead t = 1) = 0
  
  # likelihood
  for (i in 1:N){
    z[i,1] ~ dcat(delta[1:2])
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


samps <- rbind(out[[1]], out[[2]])

# Extract column names related to retention probabilities s[sex, time]
s_cols <- grep("^s\\[", colnames(samps), value = TRUE)
stopover_duration <- lapply(s_cols, function(col) 1 / (1 - samps[, col]))
names(stopover_duration) <- make.names(s_cols)  # becomes e.g., s.1..1.

# Combine into data.frame
stopover_duration_df <- as.data.frame(stopover_duration)

# Create a mapping of syntactic names → sex and time
raw_names <- names(stopover_duration_df)

# Extract sex and time from the original names
sex_time <- do.call(rbind, lapply(raw_names, function(n) {
  # Reverse the transformation done by make.names
  # from s.1..3. → 1, 3
  parts <- unlist(regmatches(n, regexec("s\\.(\\d+)\\.\\.(\\d+)\\.", n)))
  c(sex = as.integer(parts[2]), time = as.integer(parts[3]))
}))


# Add names back and convert to data.frame
colnames(sex_time) <- c("sex", "time")
sex_time <- as.data.frame(sex_time)
sex_time$param <- raw_names

# Reshape to long format
library(reshape2)
stopover_long <- melt(stopover_duration_df, variable.name = "param", value.name = "duration")

# Merge sex/time info
stopover_long <- merge(stopover_long, sex_time, by = "param")

# Compute summary statistics
stopover_summary <- aggregate(duration ~ sex + time, data = stopover_long, FUN = function(x) {
  c(mean = mean(x), median = median(x), lwr = quantile(x, 0.025), upr = quantile(x, 0.975))
})

stopover_summary <- do.call(data.frame, stopover_summary)
colnames(stopover_summary) <- c("sex", "time", "mean", "median", "lwr_95", "upr_95")

# Optional: turn sex into factor with labels
stopover_summary$sex <- factor(stopover_summary$sex, levels = 1:2, labels = c("female", "male"))

# View
print(stopover_summary)


# Now you can summarise or plot!
library(ggplot2)

# Plot
ggplot(stopover_summary, aes(x = time, y = mean, color = sex, fill = sex)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.2, color = NA) +
  geom_point(aes(y = median), shape = 21, size = 2, fill = "white") +
  labs(
    x = "Time",
    y = "Expected stopover duration (days)",
    color = "Sex",
    fill = "Sex",
    title = "Time-specific stopover duration by sex"
  ) +
  theme_minimal(base_size = 14)


ggplot(stopover_summary, aes(x = time, y = median)) +
  geom_line(color = "steelblue", size = 1.2) +
  #geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), fill = "steelblue", alpha = 0.2) +
  geom_point(aes(y = median), shape = 21, size = 2, fill = "white") +
  facet_wrap(~ sex) +
  labs(
    x = "Time",
    y = "Expected stopover duration (days)",
    title = "Time-specific stopover duration by sex"
  ) +
  theme_minimal(base_size = 14)


#------ cumulative survival

# Get all s columns
s_cols <- grep("^s\\[", colnames(samps), value = TRUE)

# Sort by sex and time for safety
s_cols_sorted <- s_cols[order(
  as.numeric(sub("s\\[(\\d+),.*", "\\1", s_cols)),
  as.numeric(sub("s\\[\\d+,\\s*(\\d+)\\]", "\\1", s_cols))
)]


# Separate columns by sex
s_cols_by_sex <- split(s_cols_sorted, 
                       sapply(s_cols_sorted, function(x) as.integer(sub("s\\[(\\d+),.*", "\\1", x))))
# Result: list(female = c(...), male = c(...))

cum_surv_list <- list()

for (sex in names(s_cols_by_sex)) {
  s_mat <- samps[, s_cols_by_sex[[sex]], drop = FALSE]
  cum_surv_sex <- t(apply(s_mat, 1, cumprod))  # one row per posterior sample, one col per time
  colnames(cum_surv_sex) <- paste0("cum_surv_", sex, "_t", seq_len(ncol(cum_surv_sex)))
  cum_surv_list[[sex]] <- cum_surv_sex
}

# Combine
cum_surv_all <- do.call(cbind, cum_surv_list)


# Create summary data.frame
summary_list <- list()

for (col in colnames(cum_surv_all)) {
  values <- cum_surv_all[, col]
  summary_list[[col]] <- data.frame(
    param = col,
    mean = mean(values),
    median = median(values),
    lwr_95 = quantile(values, 0.025),
    upr_95 = quantile(values, 0.975)
  )
}

cumulative_survival_summary <- do.call(rbind, summary_list)

# Extract sex and time from param name
cumulative_survival_summary$sex <- sub("cum_surv_(\\d).*", "\\1", cumulative_survival_summary$param)
cumulative_survival_summary$time <- as.integer(sub(".*_t(\\d+)", "\\1", cumulative_survival_summary$param))

# Optional: label sex
cumulative_survival_summary$sex <- factor(cumulative_survival_summary$sex, 
                                          levels = c("1", "2"),
                                          labels = c("female", "male"))


ggplot(cumulative_survival_summary, aes(x = time, y = mean)) +
  geom_line(aes(color = sex), size = 1.2) +
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95, fill = sex), alpha = 0.2) +
  labs(
    title = "Cumulative survival (stopover retention)",
    x = "Time",
    y = "Cumulative probability of still being present",
    color = "Sex", fill = "Sex"
  ) +
  theme_minimal(base_size = 14)


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
