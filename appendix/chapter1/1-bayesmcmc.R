library(tidyverse)
theme_set(theme_light(base_size = 14))

# Bayesian statistics & MCMC {#crashcourse}

## Introduction

## Bayes' theorem

## What is the Bayesian approach?	

## Approximating posteriors via numerical integration {#numerical-approx}

# Say we capture, mark and release $n = 57$ animals at the beginning of a winter, 
# out of which we recapture $y = 19$ animals alive 
y <- 19 # nb of success
n <- 57 # nb of attempts

# Assuming all animals are independent of each other and have the same survival 
# probability, then $y$ the number of alive animals at the end of the winter is 
# a binomial distribution with $n$ trials and $\theta$ the probability of success:
# This likelihood can be visualised in `R`: 
grid <- seq(0, 1, 0.01) # grid of values for survival
likelihood <- dbinom(y, n, grid) # compute binomial likelihood
df <- data.frame(survival = grid, likelihood = likelihood) 
df %>%
  ggplot() + 
  aes(x = survival, y = likelihood) + 
  geom_line(linewidth = 1.5)

a <- 1; b <- 1; grid <- seq(0,1,0.01); prior <- dbeta(grid,a,b)
dfprior <- data.frame(survival = grid, prior = prior) 

# Now we apply Bayes' theorem. We write a `R` function that computes the product 
# of the likelihood times the prior, or the numerator in Bayes' theorem
numerator <- function(theta) dbinom(y, n, theta) * dunif(theta, 0, 1)

# We write another function that calculates the denominator, the average likelihood
denominator <- integrate(numerator,0,1)$value

# Then we get a numerical approximation of the posterior of winter survival by 
# applying Bayes' theorem:
grid <- seq(0, 1, 0.01) # grid of values for theta
numerical_posterior <- data.frame(survival = grid, 
                                  posterior = numerator(grid)/denominator) # Bayes' theorem
numerical_posterior %>%
  ggplot() +
  aes(x = survival, y = posterior) + 
  geom_line(linewidth = 1.5)

# The distribution beta($a$,$b$) for different values of $a$ and $b$
x <- seq(0, 1, length=200)
par(mfrow = c(2,3))
# distribution a posteriori beta
plot(x,dbeta(x, 1, 1),type='l',xlab='',ylab='Density',main='beta(1,1)',lwd=3,col='black',ylim=c(0,1.5))
plot(x,dbeta(x, 2, 1),type='l',xlab='',ylab='',main='beta(2,1)',lwd=3,col='black',ylim=c(0,2))
plot(x,dbeta(x, 1, 2),type='l',xlab='',ylab='',main='beta(1,2)',lwd=3,col='black',ylim=c(0,2))
plot(x,dbeta(x, 2, 2),type='l',xlab='',ylab='Density',main='beta(2,2)',lwd=3,col='black',ylim=c(0,1.5))
plot(x,dbeta(x, 10, 10),type='l',xlab='',ylab='',main='beta(10,10)',lwd=3,col='black',ylim=c(0,3.5))
plot(x,dbeta(x, 0.8, 0.8),type='l',xlab='',ylab='',main='beta(0.8,0.8)',lwd=3,col='black',ylim=c(0.5,2.5))

#If the likelihood of the data $y$ is binomial with $n$ trials and probability 
#of success $\theta$, and the prior is a beta distribution with parameters 
#$a$ and $b$, then the posterior is a beta distribution with parameters $a + y$ 
#and $b + n - y$. In our example, we have $n = 57$ trials and $y = 19$ animals 
#that survived and a uniform prior between 0 and 1 or a beta distribution with 
#parameters $a = b = 1$, therefore survival has a beta posterior distribution 
#with parameters 20 and 39. Let's superimpose the exact posterior and the 
#numerical approximation:
explicit_posterior <- dbeta(grid, y + a, n - y + b)
dfexpposterior <- data.frame(survival = grid, explicit_posterior = explicit_posterior)
ggplot() + 
  geom_line(data = numerical_posterior, 
            aes(x = survival, y = posterior), 
            size = 1.5, 
            col = wesanderson::wes_palettes$Royal1[2],
            alpha = 0.5) + 
  geom_line(data = dfexpposterior, 
            aes(x = survival, y = explicit_posterior),
            size = 1.5, 
            col = wesanderson::wes_palettes$Royal1[3], 
            linetype = "dashed")

## Markov chain Monte Carlo (MCMC)

### Monte Carlo integration

# Posterior mean can be calculated with Monte Carlo integration:
sample_from_posterior <- rbeta(1000, 20, 39) # draw 1000 values from posterior survival beta(20,39)
mean(sample_from_posterior) # compute mean with Monte Carlo integration

# You may check that the mean we have just calculated matches closely the 
# expectation of a beta distribution:
20/(20+39) # expectation of beta(20,39)

# A 95$\%$ credible interval for winter survival can be obtained in `R` with:
quantile(sample_from_posterior, probs = c(2.5/100, 97.5/100))

### Markov chains {#markovmodelmcmc}

# In our weather example, let's run the Markov chain for 20 steps:
weather <- matrix(c(0.8, 0.2, 0.9, 0.1), nrow = 2, byrow = T) # transition matrix
steps <- 20
for (i in 1:steps){
  weather <- weather %*% weather # matrix multiplication
}
round(weather, 2) # matrix product after 20 steps

### Metropolis algorithm {#metropolis-algorithm}

# 19 animals recaptured alive out of 57 captured, marked and released
survived <- 19
released <- 57

# binomial log-likelihood function
loglikelihood <- function(x, p){
  dbinom(x = x, size = released, prob = p, log = TRUE)
}

# uniform prior density
logprior <- function(p){
  dunif(x = p, min = 0, max = 1, log = TRUE)
}

# posterior density function (log scale)
posterior <- function(x, p){
  loglikelihood(x, p) + logprior(p) # - log(Pr(data))
}

# Enough of the theory, let's implement the Metropolis algorithm in `R`. 
steps <- 100 # number of steps
theta.post <- rep(NA, steps) # vector to store samples
accept <- rep(NA, steps) # keep track of accept/reject
set.seed(1234) # for reproducibility

# First, we pick a starting value, and store it (step 1):
inits <- 0.5
theta.post[1] <- inits
accept[1] <- 1

# Then, we need a function to propose a candidate value: 
move <- function(x, away = 1){ # by default, standard deviation of the proposal distribution is 1
  logitx <- log(x / (1 - x)) # apply logit transform (-infinity,+infinity)
  logit_candidate <- logitx + rnorm(1, 0, away) # add a value taken from N(0,sd=away) to current value
  candidate <- plogis(logit_candidate) # back-transform (0,1)
  return(candidate)
}

# We add a value taken from a normal distribution with mean zero and standard 
# deviation we call *away*. We work on the logit scale to make sure the 
# candidate value for survival lies between 0 and 1.

# Now we're ready for steps 2, 3 and 4. We write a loop to take care of step 5. 
# We start at initial value 0.5 and run the algorithm for 100 steps or iterations:
for (t in 2:steps){ # repeat steps 2-4 (step 5)
  
  # propose candidate value for survival (step 2)
  theta_star <- move(theta.post[t-1])
  
  # calculate ratio R (step 3)
  pstar <- posterior(survived, p = theta_star)  
  pprev <- posterior(survived, p = theta.post[t-1])
  logR <- pstar - pprev # likelihood and prior are on the log scale
  R <- exp(logR)
  
  # accept candidate value or keep current value (step 4)
  X <- runif(1, 0, 1) # spin continuous spinner
  if (X < R){
    theta.post[t] <- theta_star # accept candidate value
    accept[t] <- 1 # accept
  }
  else{
    theta.post[t] <- theta.post[t-1] # keep current value
    accept[t] <- 0 # reject
  }
}

# We get the following values:
head(theta.post) # first values
tail(theta.post) # last values

# Visually, you may look at the chain:
df <- data.frame(x = 1:steps, y = theta.post)
df %>%
  ggplot() +
  geom_line(aes(x = x, y = y), size = 1.5, color = wesanderson::wes_palettes$Zissou1[1]) + 
  labs(x = "iterations", y = "samples") + 
  ylim(0.1, 0.6)

# The acceptance probability is the average number of times we accepted a candidate 
# value, which is `r mean(accept)` and almost satisfying. 

# To make our life easier and avoid repeating the same lines of code again and again, 
# let's make a function out of the code we have written so far:
metropolis <- function(steps = 100, inits = 0.5, away = 1){
  
  # pre-alloc memory
  theta.post <- rep(NA, steps)
  
  # start
  theta.post[1] <- inits
  
  for (t in 2:steps){
    
    # propose candidate value for prob of success
    theta_star <- move(theta.post[t-1], away = away)
    
    # calculate ratio R
    pstar <- posterior(survived, p = theta_star)  
    pprev <- posterior(survived, p = theta.post[t-1])
    logR <- pstar - pprev
    R <- exp(logR)
    
    # accept candidate value or keep current value (step 4)
    X <- runif(1, 0, 1) # spin continuous spinner
    if (X < R){
      theta.post[t] <- theta_star
    }
    else{
      theta.post[t] <- theta.post[t-1]
    }
  }
  theta.post
}

# Can we run another chain and start at initial value 0.2 this time? 
# Yes, just go through the same algorithm again, and visualise the results 
# with trace plot of survival for two chains starting at 0.2 (yellow) and 0.5 
# (blue) run for 100 steps:
theta.post2 <- metropolis(steps = 100, inits = 0.2)
df2 <- data.frame(x = 1:steps, y = theta.post2)
ggplot() +
  geom_line(data = df, aes(x = x, y = y), size = 1.5, color = wesanderson::wes_palettes$Zissou1[1]) + 
  geom_line(data = df2, aes(x = x, y = y), size = 1.5, color = wesanderson::wes_palettes$Zissou1[3]) + 
  labs(x = "iterations", y = "values from posterior distribution") + 
  ylim(0.1, 0.6)

# Let's increase the number of steps, start at 0.5 and run a chain with 5000 iterations:
steps <- 5000
set.seed(1234)
theta.post <- metropolis(steps = steps, inits = 0.5)
df <- data.frame(x = 1:steps, y = theta.post)
df %>%
  ggplot() +
  geom_line(aes(x = x, y = y), size = 1, color = wesanderson::wes_palettes$Zissou1[1]) + 
  labs(x = "iterations", y = "values from posterior distribution") + 
  ylim(0.1, 0.6) + 
  geom_hline(aes(yintercept = mean(theta.post), linetype = "posterior mean")) + 
  scale_linetype_manual(name = "", values = c(2,2)) 

# I find it informative to look at the animated version of the figure above, 
# it helps understanding the stochastic behavior of the algorithm, and also to 
# realise how the chains converge to their stationary distribution

library(gganimate)
library(magick)

# deer data, 19 "success" out of 57 "attempts"
survived <- 19
released <- 57

#---------- apply Metropolis

steps <- 1000
chain1 <- metropolis(steps = steps, inits = 0.2)
chain2 <- metropolis(steps = steps, inits = 0.5)
chain3 <- metropolis(steps = steps, inits = 0.7)

df <- data.frame(iter = rep(1:steps, 3), 
                 value = c(chain1, chain2, chain3),
                 chain = c(rep("chain1", steps), 
                           rep("chain2", steps), 
                           rep("chain3", steps)))

#---------- time series
static_tsplot <- df %>%
  mutate(posterior_mean = mean(value)) %>%
  ggplot(aes(x = iter, y = value, group = chain, color = chain)) +
  geom_line(size = 1, alpha = 0.5) + 
  geom_hline(aes(yintercept = posterior_mean, linetype = "posterior mean")) + 
  scale_linetype_manual(name = "", values = c(2,2)) + 
  labs(color = "", x = "iterations", y = "survival")
static_tsplot  

# animate
animated_tsplot <- static_tsplot +
  transition_reveal(along = iter, 
                    range = as.integer(c(1, max(df$iter) + 50))) # trick to pause
animated_tsplot  

# save
a_gif <- animate(animated_tsplot,
                 width = 6, 
                 height = 3,
                 res = 600,
                 units = "in")

# check out file69a2608670f3.gif

## Assessing convergence {#convergence-diag}

### Burn-in
  
# The simplest method to determine the length of the burn-in period is to look at 
# trace plots. Going back to our example, let's have a look to a trace plot of a 
# chain that starts at value 0.99. 

# set up the scene
steps <- 1000
theta.post <- metropolis(steps = steps, inits = 0.99)
df <- data.frame(x = 1:steps, y = theta.post)
df %>%
  ggplot() +
  geom_line(aes(x = x, y = y), size = 1.2, color = wesanderson::wes_palettes$Zissou1[1]) + 
  labs(x = "iterations", y = "survival") + 
  theme_light(base_size = 14) + 
  annotate("rect", 
           xmin = 0, 
           xmax = 100, 
           ymin = 0.1, 
           ymax = 1, 
           alpha = .3) +
  scale_y_continuous(expand = c(0,0))

# The BGR statistic asks whether there is a chain effect, and is very much alike 
# the $F$ test in an analysis of variance. Values below 1.1 indicate likely convergence.

simul.bgr <- function(steps, inits){
  
  nb.replicates <- length(inits)
  theta.post <- matrix(NA, nrow = nb.replicates, ncol = steps)
  for (i in 1:nb.replicates){
    theta.post[i,1:steps] <- metropolis(steps = steps, inits = inits[i])
  }
  
  df <- data.frame(x = rep(1:steps, nb.replicates), 
                   y = c(t(theta.post)), 
                   chain = paste0("chain ",gl(nb.replicates, steps))) %>%
    filter(x > round(steps/2)) # apply burnin (half number of iterations)

  # compute BGR (R-hat)
  num <- quantile(df$y, probs = c(20/100, 80/100))[2] - quantile(df$y, probs = c(20/100, 80/100))[1]
  den <- df %>%
    group_by(chain) %>%
    summarise(ci = quantile(y, probs = c(20/100, 80/100))) %>%
    mutate(diff = ci - lag(ci, default = ci[1])) %>%
    filter(diff != 0) %>%
    pull(diff) %>%
    mean()
  
  bgr <- round(num / den, 3)
  return(bgr)
}

set.seed(1234)
steps <- seq(100, 5000, 100)
bgr <- rep(NA, length(steps))
for (i in 1:length(steps)){
  bgr[i] <- simul.bgr(steps = steps[i], inits = c(0.2, 0.8))
}
df <- data.frame(iterations = steps, bgr = bgr)

# Back to our example, we run two Markov chains with starting values 0.2 and 0.8 
# using 100 up to 5000 iterations, and calculate the BGR statistic using half the 
# number of iterations as the length of the burn-in (code not shown): 

df %>%
  ggplot() + 
  geom_line(aes(x = iterations, y = bgr), size = 1.2) +
  labs(y = "BGR statistic")

### Chain length
  
# The figure below shows trace plots for different values of the standard deviation 
# (parameter *away*) of the normal proposal distribution we use to propose 
# a candidate value (Section \@ref(metropolis-algorithm)): 

# inspired from https://bookdown.org/content/3686/markov-chain-monte-carlo.html

n_steps <- 10000

d <-
  tibble(away = c(0.1, 1, 10)) %>% 
  mutate(accepted_traj = map(away, metropolis, steps = n_steps, inits = 0.1)) %>% 
  unnest(accepted_traj)

d <-
  d %>% 
  mutate(proposal_sd = str_c("Proposal SD = ", away),
         iter        = rep(1:n_steps, times = 3))

trace <- d %>% 
  ggplot(aes(y = accepted_traj, x = iter)) +
  geom_path(size = 1/4, color = "steelblue") +
  geom_point(size = 1/2, alpha = 1/2, color = "steelblue") +
  scale_y_continuous("survival", breaks = 0:5 * 0.1, limits = c(0.15, 0.5)) +
  scale_x_continuous("iterations", 
                     breaks = seq(n_steps-n_steps*10/100,n_steps,by = 600), 
                     limits = c(n_steps-n_steps*10/100, n_steps)) +
  facet_wrap(~proposal_sd, ncol = 3) +
  theme_light(base_size = 14)

trace

# We obtain the autocorrelation function plots for different values of the 
# standard deviation of the proposal distribution with the R `forecast::ggAcf()` function: 

library(forecast)
plot1 <- ggAcf(x = d$accepted_traj[d$proposal_sd=="Proposal SD = 0.1"]) + ggtitle("Proposal SD = 0.1")
plot2 <- ggAcf(x = d$accepted_traj[d$proposal_sd=="Proposal SD = 1"]) + ggtitle("Proposal SD = 1")
plot3 <- ggAcf(x = d$accepted_traj[d$proposal_sd=="Proposal SD = 10"]) + ggtitle("Proposal SD = 10")

library(patchwork)
(plot1 + plot2 + plot3)

# In the animal survival example, `n.eff` can be calculated with the R 
# `coda::effectiveSize()` function:

neff1 <- coda::effectiveSize(d$accepted_traj[d$proposal_sd=="Proposal SD = 0.1"])
neff2 <- coda::effectiveSize(d$accepted_traj[d$proposal_sd=="Proposal SD = 1"])
neff3 <- coda::effectiveSize(d$accepted_traj[d$proposal_sd=="Proposal SD = 10"])
df <- tibble("Proposal SD" = c(0.1, 1, 10),
             "n.eff" = round(c(neff1, neff2, neff3)))
df

### What if you have issues of convergence?
  
## Summary

## Suggested reading 
