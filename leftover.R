#------------------------------ AS model

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1, 2, label = 'non-detection'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1.5, label = 'detection in site A'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1, label = 'detection in site B'), nudge_x = -0.6, size = 7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(2, 2, label = 'alive in site A'), nudge_x = 0.5, size = 7) +
  geom_text(aes(2, 1.5, label = 'alive in site B'), nudge_x = 0.5, size = 7) +
  geom_text(aes(2, 1, label = 'dead'), nudge_x = 0.5, size = 7) +
  xlim(0, 3) +
  ylim(0.5, 3) +
  annotate('text', x = 1, y = 2.6, label = 'Observations', size = 10) +
  annotate('text', x = 2, y = 2.6, label = 'States', size = 10) +
  theme_void()
```

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1, 2, label = 'non-detection'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1.5, label = 'detection in site A'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1, label = 'detection in site B'), nudge_x = -0.6, size = 7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(2, 2, label = 'alive in site A'), nudge_x = 0.5, size = 7) +
  geom_text(aes(2, 1.5, label = 'alive in site B'), nudge_x = 0.5, size = 7) +
  geom_text(aes(2, 1, label = 'dead'), nudge_x = 0.5, size = 7) +
  xlim(0, 3) +
  ylim(0.5, 3) +
  annotate('text', x = 1, y = 2.6, label = 'Observations', size = 10) +
  annotate('text', x = 2, y = 2.6, label = 'States', size = 10) +
  geom_segment(aes(x = 1, y = 1.5, xend = 2, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  theme_void()
```

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1, 2, label = 'non-detection'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1.5, label = 'detection in site A'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1, label = 'detection in site B'), nudge_x = -0.6, size = 7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(2, 2, label = 'alive in site A'), nudge_x = 0.5, size = 7) +
  geom_text(aes(2, 1.5, label = 'alive in site B'), nudge_x = 0.5, size = 7) +
  geom_text(aes(2, 1, label = 'dead'), nudge_x = 0.5, size = 7) +
  xlim(0, 3) +
  ylim(0.5, 3) +
  annotate('text', x = 1, y = 2.6, label = 'Observations', size = 10) +
  annotate('text', x = 2, y = 2.6, label = 'States', size = 10) +
  geom_segment(aes(x = 1, y = 1.5, xend = 2, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 1, xend = 2, yend = 1.5), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  theme_void()
```



```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(2, 2, label = 'non-detection'), nudge_x = 0.5, size = 7) +
  geom_text(aes(2, 1.5, label = 'detection in site A'), nudge_x = 0.6, size = 7) +
  geom_text(aes(2, 1, label = 'detection in site B'), nudge_x = 0.6, size = 7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1, 2, label = 'alive in site A'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1.5, label = 'alive in site B'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1, label = 'dead'), nudge_x = -0.6, size = 7) +
  xlim(0, 3) +
  ylim(0.5, 3) +
  annotate('text', x = 1, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2, y = 2.6, label = 'Observations', size = 10) +
  theme_void()
```

### The model construction: How we should think.

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(2, 2, label = 'non-detection'), nudge_x = 0.5, size = 7) +
  geom_text(aes(2, 1.5, label = 'detection in site A'), nudge_x = 0.6, size = 7) +
  geom_text(aes(2, 1, label = 'detection in site B'), nudge_x = 0.6, size = 7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1, 2, label = 'alive in site A'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1.5, label = 'alive in site B'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1, label = 'dead'), nudge_x = -0.6, size = 7) +
  xlim(0, 3) +
  ylim(0.5, 3) +
  annotate('text', x = 1, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2, y = 2.6, label = 'Observations', size = 10) +
  geom_segment(aes(x = 1, y = 1, xend = 2, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  theme_void()
```

### The model construction: How we should think.

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(2, 2, label = 'non-detection'), nudge_x = 0.5, size = 7) +
  geom_text(aes(2, 1.5, label = 'detection in site A'), nudge_x = 0.6, size = 7) +
  geom_text(aes(2, 1, label = 'detection in site B'), nudge_x = 0.6, size = 7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1, 2, label = 'alive in site A'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1.5, label = 'alive in site B'), nudge_x = -0.6, size = 7) +
  geom_text(aes(1, 1, label = 'dead'), nudge_x = -0.6, size = 7) +
  xlim(0, 3) +
  ylim(0.5, 3) +
  annotate('text', x = 1, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2, y = 2.6, label = 'Observations', size = 10) +
  geom_segment(aes(x = 1, y = 2, xend = 2, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 2, xend = 2, yend = 1.5), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 1, xend = 2, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  theme_void()
```


```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1.5, 2, label = 'not encountered (0)'), nudge_x = 1, size = 7) +
  geom_text(aes(1.5, 1.5, label = 'found, ascertained as breeder (1)'), nudge_x = 1.5, size = 7) +
  geom_text(aes(1.5, 1, label = 'found, ascertained as non-breeder (2)'), nudge_x = 1.7, size = 7) +
  geom_text(aes(1.5, 0.5, label = 'found, status unknown (3)'), nudge_x = 1.2, size = 7) +
  geom_point(aes(1.5, 0.5), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 2), size = 2.5, alpha = .7) +
  geom_text(aes(.5, 2, label = 'breeding'), nudge_x = 0, size = 7) +
  geom_text(aes(.5, 1.5, label = 'non-breeding'), nudge_x = -0.2, size = 7) +
  geom_text(aes(.5, 1, label = 'dead'), nudge_x = 0.1, size = 7) +
  xlim(0, 4.5) +
  ylim(0.5, 3) +
  annotate('text', x = .5, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2.5, y = 2.6, label = 'Observations', size = 10) +
  theme_void()
```

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1.5, 2, label = 'not encountered (0)'), nudge_x = 1, size = 7) +
  geom_text(aes(1.5, 1.5, label = 'found, ascertained as breeder (1)'), nudge_x = 1.5, size = 7) +
  geom_text(aes(1.5, 1, label = 'found, ascertained as non-breeder (2)'), nudge_x = 1.7, size = 7) +
  geom_text(aes(1.5, 0.5, label = 'found, status unknown (3)'), nudge_x = 1.2, size = 7) +
  geom_point(aes(1.5, 0.5), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 2), size = 2.5, alpha = .7) +
  geom_text(aes(.5, 2, label = 'breeding'), nudge_x = 0, size = 7) +
  geom_text(aes(.5, 1.5, label = 'non-breeding'), nudge_x = -0.2, size = 7) +
  geom_text(aes(.5, 1, label = 'dead'), nudge_x = 0.1, size = 7) +
  xlim(0, 4.5) +
  ylim(0.5, 3) +
  annotate('text', x = .5, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2.5, y = 2.6, label = 'Observations', size = 10) +
  geom_segment(aes(x = 1, y = 1, xend = 1.5, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  theme_void()
```

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1.5, 2, label = 'not encountered (0)'), nudge_x = 1, size = 7) +
  geom_text(aes(1.5, 1.5, label = 'found, ascertained as breeder (1)'), nudge_x = 1.5, size = 7) +
  geom_text(aes(1.5, 1, label = 'found, ascertained as non-breeder (2)'), nudge_x = 1.7, size = 7) +
  geom_text(aes(1.5, 0.5, label = 'found, status unknown (3)'), nudge_x = 1.2, size = 7) +
  geom_point(aes(1.5, 0.5), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 2), size = 2.5, alpha = .7) +
  geom_text(aes(.5, 2, label = 'breeding'), nudge_x = 0, size = 7) +
  geom_text(aes(.5, 1.5, label = 'non-breeding'), nudge_x = -0.2, size = 7) +
  geom_text(aes(.5, 1, label = 'dead'), nudge_x = 0.1, size = 7) +
  xlim(0, 4.5) +
  ylim(0.5, 3) +
  annotate('text', x = .5, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2.5, y = 2.6, label = 'Observations', size = 10) +
  geom_segment(aes(x = 1, y = 2, xend = 1.5, yend = .5), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 2, xend = 1.5, yend = 1.5), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 2, xend = 1.5, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  
  theme_void()
```

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1.5, 2, label = 'not encountered (0)'), nudge_x = 1, size = 7) +
  geom_text(aes(1.5, 1.5, label = 'found, ascertained as breeder (1)'), nudge_x = 1.5, size = 7) +
  geom_text(aes(1.5, 1, label = 'found, ascertained as non-breeder (2)'), nudge_x = 1.7, size = 7) +
  geom_text(aes(1.5, 0.5, label = 'found, status unknown (3)'), nudge_x = 1.2, size = 7) +
  geom_point(aes(1.5, 0.5), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1.5, 2), size = 2.5, alpha = .7) +
  geom_text(aes(.5, 2, label = 'breeding'), nudge_x = 0, size = 7) +
  geom_text(aes(.5, 1.5, label = 'non-breeding'), nudge_x = -0.2, size = 7) +
  geom_text(aes(.5, 1, label = 'dead'), nudge_x = 0.1, size = 7) +
  xlim(0, 4.5) +
  ylim(0.5, 3) +
  annotate('text', x = .5, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2.5, y = 2.6, label = 'Observations', size = 10) +
  geom_segment(aes(x = 1, y = 1.5, xend = 1.5, yend = 1), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 1.5, xend = 1.5, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 1.5, xend = 1.5, yend = .5), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  theme_void()
```



#----------------------------------- random effects

#### Temporal

Include temporal covariates, say $x_t$ with $\text{logit}(\phi_t) = \beta_1 + \beta_2 x_t$. 
If temporal variation not fully explained by covariates, add random effects $\text{logit}(\phi_t) 
= \beta_1 + \beta_2 x_t + \varepsilon_t, \; \varepsilon_t \sim N(0,\sigma^2)$. 
We may wish to allow for extra variation in the survival vs. water flow relationship. 
To do so, we consider a yearly random effect. The prior on the standard deviation of 
the random effect is uniform between 0 and 10. **explain how to pick prior for SD**
  ```{r eval = FALSE}
hmm.phiflowREpt <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    logit(phi[t]) <- beta[1] + beta[2] * flow[t] + eps[t] # eps is random effect
    eps[t] ~ dnorm(0, sd = sdeps) 
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
    p[t] ~ dunif(0, 1)          # prior detection
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0           # Pr(dead t -> detected t)
  }
  beta[1] ~ dnorm(0, 1.5) # prior intercept
  beta[2] ~ dnorm(0, 1.5) # prior slope
  sdeps ~ dunif(0,10)
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
    }
  }
})
```

Initial values. 
```{r eval = FALSE}
initial.values <- function() list(beta = rnorm(2,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  sdeps = runif(1,0,3),
                                  z = zinits)
```

Parameters to be monitored. 
```{r eval = FALSE}
parameters.to.save <- c("beta", "p", "phi", "sdeps")
```

MCMC details. Note that we've increased the number of iterations and the length of the burn-in period.
```{r eval = FALSE}
n.iter <- 10000
n.burnin <- 5000
n.chains <- 2
```

Run NIMBLE:
```{r eval = FALSE}
mcmc.phiflowREpt <- nimbleMCMC(code = hmm.phiflowREpt, 
                               constants = my.constants,
                               data = my.data,              
                               inits = initial.values,
                               monitors = parameters.to.save,
                               niter = n.iter,
                               nburnin = n.burnin, 
                               nchains = n.chains)
```

Display outputs. Seems that the water flow effect is not so important anymore. 
```{r eval = FALSE}
MCMCsummary(object = mcmc.phiflowREpt, round = 2)
```
```{r echo = FALSE}
load(here::here("dat/phiflowREpt.RData"))
MCMCsummary(object = mcmc.phiflowREpt, round = 2)
```

Trace plots for the standard deviation of the random effect.
```{r}
MCMCtrace(object = mcmc.phiflowREpt, 
          params = "sdeps", 
          pdf = FALSE)
```

#### Individual


#---------------------------------- covariates

Sex effect

+ Let's use a covariate $\text{sex}$ that takes value 0 if male, and 1 if female

+ And write $\text{logit}(\phi_i) = \beta_1 + \beta_2 \; \text{sex}_i$ for bird $i$

+ Then male survival is

$$\text{logit}(\phi_i) = \beta_1$$

+ And female survival is

$$\text{logit}(\phi_i) = \beta_1 + \beta_2$$

Nimble implementation with sex as a covariate

```{r eval = FALSE}
hmm.phisexp <- nimbleCode({
...
  for (i in 1:N){
    logit(phi[i]) <- beta[1] + beta[2] * sex[i]
    gamma[1,1,i] <- phi[i]      # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i]  # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0        # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1        # Pr(dead t -> dead t+1)
  }
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  phi_male <- 1/(1+exp(-beta[1]))
  phi_female <- 1/(1+exp(-(beta[1]+beta[2])))
...
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})
```

```{r echo = FALSE}
load(here::here("dat/phisexp.RData"))
MCMCsummary(object = mcmc.phisexp, round = 2)
```


#----------------------------------- further examples multievent


## Animal epidemiology with uncertain disease states

Let's have a look to another example. Very similar to the previous example. We consider a system of an emerging pathogen *Mycoplasma gallisepticum* Edward and Kanarek and its host the house finch, *Carpodacus mexicanus* Müller.

```{r}
knitr::include_graphics("images/infectedhousefinch.jpg")
```

A house finch with a heavy infection (Jim Mondok).

We consider a system of an emerging pathogen *Mycoplasma gallisepticum* Edward and Kanarek and its host the house finch, *Carpodacus mexicanus* Müller. [Faustino et al. (2004)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.0021-8790.2004.00840.x) and [Conn & Cooch (2009)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2664.2008.01597.x) studied impact of pathogen on host demographic rates. Problem is true disease state for some encountered individuals is ambiguous because seen at distance. In this context, how to study the dynamics of the disease?

### States and observations

+ 3 states
    + healthy (H)
    + ill (I)
    + dead (D)

+ 4 observations
    + not seen (0)
    + captured healthy (1)
    + captured ill (2)
    + health status unknown, i.e. seen at distance (3)


### How states generate observations.

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1.5, 2, label = 'not seen (0)'), nudge_x = 1.2, size = 7) +
  geom_text(aes(1.5, 1.5, label = 'captured healthy (1)'), nudge_x = 1.5, size = 7) +
  geom_text(aes(1.5, 1, label = 'captured ill (2)'), nudge_x = 1.3, size = 7) +
  geom_text(aes(1.5, 0.5, label = 'status unknown (3)'), nudge_x = 1.5, size = 7) +
  geom_point(aes(2, 0.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(.5, 2, label = 'healthy'), nudge_x = 0, size = 7) +
  geom_text(aes(.5, 1.5, label = 'ill'), nudge_x = 0.2, size = 7) +
  geom_text(aes(.5, 1, label = 'dead'), nudge_x = 0.1, size = 7) +
  xlim(0, 4.5) +
  ylim(0.5, 3) +
  annotate('text', x = .5, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2.7, y = 2.6, label = 'Observations', size = 10) +

  theme_void()
```

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1.5, 2, label = 'not seen (0)'), nudge_x = 1.2, size = 7) +
  geom_text(aes(1.5, 1.5, label = 'captured healthy (1)'), nudge_x = 1.5, size = 7) +
  geom_text(aes(1.5, 1, label = 'captured ill (2)'), nudge_x = 1.3, size = 7) +
  geom_text(aes(1.5, 0.5, label = 'status unknown (3)'), nudge_x = 1.5, size = 7) +
  geom_point(aes(2, 0.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(.5, 2, label = 'healthy'), nudge_x = 0, size = 7) +
  geom_text(aes(.5, 1.5, label = 'ill'), nudge_x = 0.2, size = 7) +
  geom_text(aes(.5, 1, label = 'dead'), nudge_x = 0.1, size = 7) +
  xlim(0, 4.5) +
  ylim(0.5, 3) +
  annotate('text', x = .5, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2.7, y = 2.6, label = 'Observations', size = 10) +

  geom_segment(aes(x = 1, y = 1, xend = 2, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +

  theme_void()
```

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1.5, 2, label = 'not seen (0)'), nudge_x = 1.2, size = 7) +
  geom_text(aes(1.5, 1.5, label = 'captured healthy (1)'), nudge_x = 1.5, size = 7) +
  geom_text(aes(1.5, 1, label = 'captured ill (2)'), nudge_x = 1.3, size = 7) +
  geom_text(aes(1.5, 0.5, label = 'status unknown (3)'), nudge_x = 1.5, size = 7) +
  geom_point(aes(2, 0.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(.5, 2, label = 'healthy'), nudge_x = 0, size = 7) +
  geom_text(aes(.5, 1.5, label = 'ill'), nudge_x = 0.2, size = 7) +
  geom_text(aes(.5, 1, label = 'dead'), nudge_x = 0.1, size = 7) +
  xlim(0, 4.5) +
  ylim(0.5, 3) +
  annotate('text', x = .5, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2.7, y = 2.6, label = 'Observations', size = 10) +


  geom_segment(aes(x = 1, y = 1.5, xend = 2, yend = 1), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 1.5, xend = 2, yend = 2), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 1.5, xend = 2, yend = .5), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +

  theme_void()
```

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1.5, 2, label = 'not seen (0)'), nudge_x = 1.2, size = 7) +
  geom_text(aes(1.5, 1.5, label = 'captured healthy (1)'), nudge_x = 1.5, size = 7) +
  geom_text(aes(1.5, 1, label = 'captured ill (2)'), nudge_x = 1.3, size = 7) +
  geom_text(aes(1.5, 0.5, label = 'status unknown (3)'), nudge_x = 1.5, size = 7) +
  geom_point(aes(2, 0.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(.5, 2, label = 'healthy'), nudge_x = 0, size = 7) +
  geom_text(aes(.5, 1.5, label = 'ill'), nudge_x = 0.2, size = 7) +
  geom_text(aes(.5, 1, label = 'dead'), nudge_x = 0.1, size = 7) +
  xlim(0, 4.5) +
  ylim(0.5, 3) +
  annotate('text', x = .5, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2.7, y = 2.6, label = 'Observations', size = 10) +


  geom_segment(aes(x = 1, y = 2, xend = 2, yend = .5), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 2, xend = 2, yend = 1.5), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 2, xend = 2, yend = 2), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +

  theme_void()
```

```{r, echo = FALSE}
ggplot() +
  geom_point(aes(1, 1), size = 2.5, alpha = .7) +
  geom_point(aes(1, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(1, 2), size = 2.5, alpha = .7) +
  geom_text(aes(1.5, 2, label = 'not seen (0)'), nudge_x = 1.2, size = 7) +
  geom_text(aes(1.5, 1.5, label = 'captured healthy (1)'), nudge_x = 1.5, size = 7) +
  geom_text(aes(1.5, 1, label = 'captured ill (2)'), nudge_x = 1.3, size = 7) +
  geom_text(aes(1.5, 0.5, label = 'status unknown (3)'), nudge_x = 1.5, size = 7) +
  geom_point(aes(2, 0.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1), size = 2.5, alpha = .7) +
  geom_point(aes(2, 1.5), size = 2.5, alpha = .7) +
  geom_point(aes(2, 2), size = 2.5, alpha = .7) +
  geom_text(aes(.5, 2, label = 'healthy'), nudge_x = 0, size = 7) +
  geom_text(aes(.5, 1.5, label = 'ill'), nudge_x = 0.2, size = 7) +
  geom_text(aes(.5, 1, label = 'dead'), nudge_x = 0.1, size = 7) +
  xlim(0, 4.5) +
  ylim(0.5, 3) +
  annotate('text', x = .5, y = 2.6, label = 'States', size = 10) +
  annotate('text', x = 2.7, y = 2.6, label = 'Observations', size = 10) +

  geom_segment(aes(x = 1, y = 1, xend = 2, yend = 2), alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +

  geom_segment(aes(x = 1, y = 1.5, xend = 2, yend = 1), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 1.5, xend = 2, yend = 2), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 1.5, xend = 2, yend = .5), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +

  geom_segment(aes(x = 1, y = 2, xend = 2, yend = .5), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 2, xend = 2, yend = 1.5), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_segment(aes(x = 1, y = 2, xend = 2, yend = 2), lty = 2, alpha = 0.7, arrow = arrow(length = unit(0.02, "npc"))) +

  theme_void()
```

Vector of initial state probabilities

$$
\begin{matrix}
& \\
\mathbf{\delta} =
    \left ( \vphantom{ \begin{matrix} 12 \end{matrix} } \right .
\end{matrix}
\hspace{-1.2em}
\begin{matrix}
    z_t=H & z_t=I & z_t=D \\ \hdashline
\pi_H & 1 - \pi_{H} & 0\\
\end{matrix}
\hspace{-0.2em}
\begin{matrix}
& \\
\left . \vphantom{ \begin{matrix} 12 \end{matrix} } \right )
    \begin{matrix}
    \end{matrix}
\end{matrix}
$$

$\pi_H$ is the probability that a newly encountered individual is healthy. $\pi_{I} = 1 - \pi_H$ is the probability that a newly encountered individual is ill.

### HMM model for disease states with uncertainty

Transition matrix

$$
\begin{matrix}
& \\
\mathbf{\Gamma} =
    \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right .
\end{matrix}
\hspace{-1.2em}
\begin{matrix}
    z_t=H & z_t=I & z_t=D \\ \hdashline
\phi_H (1-\psi_{HI}) & \phi_H \psi_{HI} & 1 - \phi_H\\
\phi_{I} \psi_{IH} & \phi_{I} (1-\psi_{IH}) & 1 - \phi_{I}\\
0 & 0 & 1
\end{matrix}
\hspace{-0.2em}
\begin{matrix}
& \\
\left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
    \begin{matrix}
    z_{t-1}=H \\ z_{t-1}=I \\ z_{t-1}=D
    \end{matrix}
\end{matrix}
$$

$\phi_H$ is the survival probability of healthy individuals, $\phi_I$ that of ill individuals. $\psi_{HI}$ is the probability of getting sick, $\psi_{IH}$ that of recovering from the disease.

Transition matrix, incurable disease

$$
\begin{matrix}
& \\
\mathbf{\Gamma} =
    \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right .
\end{matrix}
\hspace{-1.2em}
\begin{matrix}
    z_t=H & z_t=I & z_t=D \\ \hdashline
\phi_H (1-\psi_{HI}) & \phi_H \psi_{HI} & 1 - \phi_H\\
0 & \phi_{I}  & 1 - \phi_{I}\\
0 & 0 & 1
\end{matrix}
\hspace{-0.2em}
\begin{matrix}
& \\
\left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
    \begin{matrix}
    z_{t-1}=H \\ z_{t-1}=I \\ z_{t-1}=D
    \end{matrix}
\end{matrix}
$$

No possibility of recovering from the disease, that is $\psi_{IH} = 0$. Once you get sick, you remain sick $\psi_{II} = 1 - \psi_{IH} = 1$. For analysing the house finch data, we allow recovering from the disease, and we will use transition matrix from previous slide.

Observation matrix

$$
\begin{matrix}
& \\
\mathbf{\Omega} =
    \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12\end{matrix} } \right .
\end{matrix}
\hspace{-1.2em}
\begin{matrix}
    y_t=0 & y_t=1 & y_t=2 & y_t=3\\ \hdashline
1-p_H & p_H \beta_H & 0 & p_H (1-\beta_H)\\
1-p_I & 0 & p_{I} \beta_{I} & p_{I} (1-\beta_{I})\\
1 & 0 & 0 & 0
\end{matrix}
\hspace{-0.2em}
\begin{matrix}
& \\
\left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12\end{matrix} } \right )
    \begin{matrix}
    z_{t}=H \\ z_{t}=I \\ z_{t}=D
    \end{matrix}
\end{matrix}
$$

$\beta_H$ is the probability to assign a healthy individual to state H. $\beta_{I}$ is the probability to assign a sick individual to state I. $p_H$ is the detection probability of healthy individuals, $p_I$ that of sick individuals.

### Results

```{r echo = FALSE}
library(MCMCvis)
load(here::here("dat","disease.RData"))
MCMCsummary(out, round = 2)
```

Healthy individuals are correctly assigned, while infected individuals are difficult to ascertain. Sounds like being infected has an effect on detection and survival. Run models without effects and compare with WAIC for formal testing. Infection rate is 22%, recovery rate is 46%.

```{r, echo = FALSE}
library(MCMCvis)
load(here::here("dat","disease.RData"))
MCMCplot(out)
```

## Individual heterogeneity with finite mixtures.

Our last example is about individual heterogeneity and how to account for it with HMMs. Gray wolf is a social species with hierarchy in packs which may reflect in species demography. As an example, we'll work with gray wolves.

```{r}
knitr::include_graphics("images/wolfdominance.jpg")
```

Gray wolf is a social species with hierarchy in packs which may reflect in demography. Shirley Pledger in a series of papers developed heterogeneity models in which individuals are assigned in two or more classes with class-specific survival/detection probabilities. [Cubaynes et al. (2010)](https://conbio.onlinelibrary.wiley.com/doi/abs/10.1111/j.1523-1739.2009.01431.x) used HMMs to account for heterogeneity in the detection process due to social status, see also [Pradel et al. (2009)](https://link.springer.com/chapter/10.1007%2F978-0-387-78151-8_36). Dominant individuals tend to use path more often than others, and these paths are where we look for scats.

### Individual heterogeneity

+ 3 states
+ alive in class 1 (A1)
+ alive in class 2 (A2)
+ dead (D)

+ 4 observations
+ not captured (0)
+ captured (1)

### HMM model for individual heterogeneity

Vector of initial state probabilities

$$
  \begin{matrix}
& \\
\mathbf{\delta} =
  \left ( \vphantom{ \begin{matrix} 12 \end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          z_t=A1 & z_t=A2 & z_t=D \\ \hdashline
          \pi & 1 - \pi & 0\\
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \end{matrix} } \right )
\begin{matrix}
\end{matrix}
\end{matrix}
$$
  
  $\pi$ is the probability of being alive in class 1. $1 - \pi$ is the probability of being in class 2.


### HMM model for individual heterogeneity

Transition matrix

$$
  \begin{matrix}
& \\
\mathbf{\Gamma} =
  \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          z_t=A1 & z_t=A2 & z_t=D \\ \hdashline
          \phi  & 0 & 1 - \phi\\
          0 & \phi & 1 - \phi\\
          0 & 0 & 1
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
\begin{matrix}
z_{t-1}=A1 \\ z_{t-1}=A2 \\ z_{t-1}=D
\end{matrix}
\end{matrix}
$$
  
  $\phi$ is the survival probability, which could be made heterogeneous.

Transition matrix, with change in heterogeneity class

$$
  \begin{matrix}
& \\
\mathbf{\Gamma} =
  \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          z_t=A1 & z_t=A2 & z_t=D \\ \hdashline
          \phi (1-\psi_{12}) & \phi \psi_{12} & 1 - \phi\\
          \phi \psi_{21} & \phi (1-\psi_{21}) & 1 - \phi\\
          0 & 0 & 1
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
\begin{matrix}
z_{t-1}=A1 \\ z_{t-1}=A2 \\ z_{t-1}=D
\end{matrix}
\end{matrix}
$$
  
  $\psi_{12}$ is the probability for an individual to change class of heterogeneity, from 1 to 2. $\psi_{21}$ is the probability for an individual to change class of heterogeneity, from 2 to 1.

Observation matrix

$$
  \begin{matrix}
& \\
\mathbf{\Omega} =
  \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12\end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          y_t=0 & y_t=1\\ \hdashline
          1 - p_1 & p_1\\
          1 - p_2 & p_2\\
          1 & 0
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12\end{matrix} } \right )
\begin{matrix}
z_{t}=A1 \\ z_{t}=A2 \\ z_{t}=D
\end{matrix}
\end{matrix}
$$
  
  $p_1$ is detection for individuals in class 1, and $p_2$ that of individuals in class 2.

### Results

```{r echo = FALSE}
library(MCMCvis)
load(here::here("dat","wolf_het.RData"))
MCMCsummary(mcmc.phipmix, round = 2)
```

We have lowly detectable individuals (class A1 with $p_1$) in proportion 62%. And highly (or so) detectable individuals (class A2 with $p_2$) in proportion 38%. Note that interpretation of classes is made a posteriori. Survival is 81%.

```{r, echo = FALSE}
library(MCMCvis)
load(here::here("dat","wolf_het.RData"))
MCMCplot(mcmc.phipmix)
```

You may consider more classes, and select among models, see [Cubaynes et al. (2012)](https://oliviergimenez.github.io/pubs/Cubaynesetal2011MEE.pdf). You may also go for a non-parametric approach and let the data tell you how many classes you need. This is relatively easy to do in Nimble, see [Turek et al. (2021)](https://arxiv.org/abs/2007.10163). More about individual heterogeneity in [Gimenez et al. (2018)](https://oliviergimenez.github.io/pubs/GimenezCamGaillard2017Oikos.pdf).

## HMMs to analyse capture-recapture data

With the same data, ask further questions, just consider different states.

### How to make our models remember?

So far, the dynamics of the states are first-order Makovian. The site where you will be depends only on the site where you are, and not on the sites you were previously. How to relax this assumption, and go second-order Markovian?
  
  Memory models were initially proposed by [Hestbeck et al. (1991)](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.2307/2937193) and [Brownie et al. (1993)](https://www.jstor.org/stable/2532259?origin=crossref&seq=1#metadata_info_tab_contents), then formulated as HMMs in [Rouan et al. (2009)](https://link.springer.com/article/10.1198/jabes.2009.06108). See also [Cole et al. (2014)](https://onlinelibrary.wiley.com/doi/10.1002/ece3.1037).
                                                                                                                                                                     
                                                                                                                                                                     ### Remember HMM model for dispersal between 2 sites
                                                                                                                                                                     
                                                                                                                                                                     Transition matrix
                                                                                                                                                                     
                                                                                                                                                                     $$
                                                                                                                                                                       \begin{matrix}
                                                                                                                                                                     & \\
                                                                                                                                                                     \mathbf{\Gamma} =
                                                                                                                                                                       \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right .
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-1.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               z_t=A & z_t=B & z_t=D \\ \hdashline
                                                                                                                                                                               \phi_A (1-\psi_{AB}) & \phi_A \psi_{AB} & 1 - \phi_A\\
                                                                                                                                                                               \phi_B \psi_{BA} & \phi_B (1-\psi_{BA}) & 1 - \phi_B\\
                                                                                                                                                                               0 & 0 & 1
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-0.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               & \\
                                                                                                                                                                               \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
                                                                                                                                                                     \begin{matrix}
                                                                                                                                                                     z_{t-1}=A \\ z_{t-1}=B \\ z_{t-1}=D
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     $$
                                                                                                                                                                       
                                                                                                                                                                       Observation matrix
                                                                                                                                                                     
                                                                                                                                                                     $$
                                                                                                                                                                       \begin{matrix}
                                                                                                                                                                     & \\
                                                                                                                                                                     \mathbf{\Omega} =
                                                                                                                                                                       \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right .
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-1.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               y_t=0 & y_t=1 & y_t=2 \\ \hdashline
                                                                                                                                                                               1 - p_A & p_A & 0\\
                                                                                                                                                                               1 - p_B & 0 & p_B\\
                                                                                                                                                                               1 & 0 & 0
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-0.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               & \\
                                                                                                                                                                               \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
                                                                                                                                                                     \begin{matrix}
                                                                                                                                                                     z_{t}=A \\ z_{t}=B \\ z_{t}=D
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     $$
                                                                                                                                                                       
                                                                                                                                                                       ### HMM formulation of the memory model
                                                                                                                                                                       
                                                                                                                                                                       To keep track of the sites previously visited, the trick is to consider states as being pairs of sites occupied
                                                                                                                                                                     
                                                                                                                                                                     + States
                                                                                                                                                                     + AA is for alive in site A at $t$ and alive in site A at $t-1$
                                                                                                                                                                       + AB is for alive in site A at $t$ and alive in site B at $t-1$
                                                                                                                                                                       + BA is for alive in site B at $t$ and alive in site A at $t-1$
                                                                                                                                                                       + BB is for alive in site B at $t$ and alive in site B at $t-1$
                                                                                                                                                                       + D is for dead
                                                                                                                                                                     
                                                                                                                                                                     + Observations
                                                                                                                                                                     + 0 not captured
                                                                                                                                                                     + 1 captured at site A
                                                                                                                                                                     + 2 captured at site B
                                                                                                                                                                     
                                                                                                                                                                     Vector of initial state probabilities
                                                                                                                                                                     
                                                                                                                                                                     $$
                                                                                                                                                                       \begin{matrix}
                                                                                                                                                                     & \\
                                                                                                                                                                     \mathbf{\delta} =
                                                                                                                                                                       \left ( \vphantom{ \begin{matrix} 12 \end{matrix} } \right .
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-1.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               z_t=AA & z_t=AB & z_t=BA & z_t=BB &z_t=D \\ \hdashline
                                                                                                                                                                               \pi_{AA} & \pi_{AB} & \pi_{BA} & \pi_{BB} & 0\\
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-0.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               & \\
                                                                                                                                                                               \left . \vphantom{ \begin{matrix} 12 \end{matrix} } \right )
                                                                                                                                                                     \begin{matrix}
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     $$
                                                                                                                                                                       
                                                                                                                                                                       where $\pi_{BB} = 1 - (\pi_{AA} + \pi_{AB} + \pi_{BA})$, and $\pi_{ij}$ at site $j$ when first captured at $t$ and site $i$ at $t - 1$.
                                                                                                                                                                     
                                                                                                                                                                     Transition matrix
                                                                                                                                                                     
                                                                                                                                                                     $$
                                                                                                                                                                       \begin{matrix}
                                                                                                                                                                     & \\
                                                                                                                                                                     \mathbf{\Gamma} =
                                                                                                                                                                       \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\ 12\end{matrix} } \right .
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-1.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               z_t=AA & z_t=AB & z_t=BA & z_t=BB & z_t=D \\ \hdashline
                                                                                                                                                                               \phi_{AAA} & \phi_{AAB} & 0 & 0 & 1 - \phi_{AAA} - \phi_{AAB}\\
                                                                                                                                                                               0 & 0 & \phi_{ABA} & \phi_{ABB} & 1 - \phi_{ABA} - \phi_{ABB}\\
                                                                                                                                                                               \phi_{BAA} & \phi_{BAB} & 0 & 0 & 1 - \phi_{BAA} - \phi_{BAB}\\
                                                                                                                                                                               0 & 0 & \phi_{BBA} & \phi_{BBB} & 1 - \phi_{BBA} - \phi_{BBB}\\
                                                                                                                                                                               0 & 0 & 0 & 0 & 1
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-0.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               & \\
                                                                                                                                                                               \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\ 12\end{matrix} } \right )
                                                                                                                                                                     \begin{matrix}
                                                                                                                                                                     z_t=AA \\ z_t=AB \\ z_t=BA \\ z_t=BB \\ z_t=D
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     $$
                                                                                                                                                                       
                                                                                                                                                                       $\phi_{ijk}$ is probability to be in site $k$ at time $t + 1$ for an individual
                                                                                                                                                                     present in site $j$ at $t$ and in site $i$ at $t - 1$
                                                                                                                                                                       
                                                                                                                                                                       Transition matrix, alternate parameterization
                                                                                                                                                                     
                                                                                                                                                                     $$
                                                                                                                                                                       \begin{matrix}
                                                                                                                                                                     & \\
                                                                                                                                                                     \mathbf{\Gamma} =
                                                                                                                                                                       \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\ 12\end{matrix} } \right .
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-1.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               z_t=AA & z_t=AB & z_t=BA & z_t=BB & z_t=D \\ \hdashline
                                                                                                                                                                               \phi \psi_{AAA} & \phi (1 - \psi_{AAA}) & 0 & 0 & 1 - \phi\\
                                                                                                                                                                               0 & 0 & \phi (1 - \psi_{ABB}) & \phi \psi_{ABB} & 1 - \phi\\
                                                                                                                                                                               \phi \psi_{BAA} & \phi (1 - \psi_{BAA}) & 0 & 0 & 1 - \phi\\
                                                                                                                                                                               0 & 0 & \phi (1-\psi_{BBB}) & \phi \psi_{BBB} & 1 - \phi\\
                                                                                                                                                                               0 & 0 & 0 & 0 & 1
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-0.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               & \\
                                                                                                                                                                               \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\ 12\end{matrix} } \right )
                                                                                                                                                                     \begin{matrix}
                                                                                                                                                                     z_t=AA \\ z_t=AB \\ z_t=BA \\ z_t=BB \\ z_t=D
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     $$
                                                                                                                                                                       
                                                                                                                                                                       $\phi$ is the probability of surviving from one occasion to the next. $\psi_{ijj}$ is the probability an animal stays at the same site $j$ given that it was at site $i$ on the previous occasion.
                                                                                                                                                                     
                                                                                                                                                                     Observation matrix
                                                                                                                                                                     
                                                                                                                                                                     $$
                                                                                                                                                                       \begin{matrix}
                                                                                                                                                                     & \\
                                                                                                                                                                     \mathbf{\Omega} =
                                                                                                                                                                       \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\ 12\end{matrix} } \right .
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-1.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               y_t=0 & y_t=1 & y_t=2 \\ \hdashline
                                                                                                                                                                               1 - p_A & p_A & 0\\
                                                                                                                                                                               1 - p_B & 0 & p_B\\
                                                                                                                                                                               1 - p_A & p_A & 0\\
                                                                                                                                                                               1 - p_B & 0 & p_B\\
                                                                                                                                                                               1 & 0 & 0
                                                                                                                                                                               \end{matrix}
                                                                                                                                                                               \hspace{-0.2em}
                                                                                                                                                                               \begin{matrix}
                                                                                                                                                                               & \\
                                                                                                                                                                               \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\ 12\end{matrix} } \right )
                                                                                                                                                                     \begin{matrix}
                                                                                                                                                                     z_t=AA \\ z_t=AB \\ z_t=BA \\ z_t=BB \\ z_t=D
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     \end{matrix}
                                                                                                                                                                     $$

#------------------------------------ flexibility of multistate

## Multistate models are very flexible

Access to reproduction

Temporary emigration

Combination of life and dead encounters

### Access to reproduction

Transition matrix:
  
  $$
  \begin{matrix}
& \\
\mathbf{\Gamma} =
  \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\ 12 \end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          z_t=J & z_t=1yNB & z_t=2yNB & z_t=B & z_t=D \\ \hdashline
          0 & \phi_1 (1-\alpha_1) & 0 & \phi_1 \alpha_1 & 1 - \phi_1\\
          0 & 0 & \phi_2 (1-\alpha_2) & \phi_2 \alpha_2 & 1 - \phi_2\\
          0 & 0 & 0 & \phi_3 & 1 - \phi_3\\
          0 & 0 & 0 & \phi_B & 1 - \phi_B\\
          0 & 0 & 0 & 0 & 1
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\ 12 \end{matrix} } \right )
\begin{matrix}
z_{t-1} = J \\ z_{t-1} = 1yNB \\ z_{t-1} = 2yNB \\ z_{t-1} = B \\ z_{t-1} = D
\end{matrix}
\end{matrix}
$$
  
  First-year and second-year individuals breed with probabilities $\alpha_1$ and $\alpha_2$.

Then, everybody breeds from age 3.

### Access to reproduction

Observation matrix:
  
  $$
  \begin{matrix}
& \\
\mathbf{\Omega} =
  \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\12 \end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          y_t = 0 & y_t = 1 & y_t = 2 & y_t = 3\\ \hdashline
          1 & 0 & 0 & 0\\
          1 - p_1 & p_1 & 0 & 0\\
          1 - p_2 & 0 & p_2 & 0\\
          1 - p_3 & 0 & 0 & p_3\\
          1 & 0 & 0 & 0
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \\ 12 \\12 \end{matrix} } \right )
\begin{matrix}
z_t = J \\ z_t = 1yNB \\ z_t = 2yNB \\ z_t = B \\ z_t = D
\end{matrix}
\end{matrix}
$$
  
  Juveniles are never detected.

## Temporary emigration

Transition matrix:
  
  $$
  \begin{matrix}
& \\
\mathbf{\Gamma} =
  \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          z_t=\text{in} & z_t=\text{out} & z_t=\text{D} \\ \hdashline
          \phi (1-\psi_{\text{in} \rightarrow \text{out}}) & \phi \psi_{\text{in} \rightarrow \text{out}} & 1 - \phi\\
          \phi \psi_{\text{out} \rightarrow \text{in}} & \phi (1-\psi_{\text{out} \rightarrow \text{in}}) & 1 - \phi\\
          0 & 0 & 1
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
\begin{matrix}
z_{t-1}=\text{in} \\ z_{t-1}=\text{out} \\ z_{t-1}=\text{D}
\end{matrix}
\end{matrix}
$$
  
  Observation matrix:
  
  $$
  \begin{matrix}
& \\
\mathbf{\Omega} =
  \left ( \vphantom{ \begin{matrix} 12 \\ 12 \end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          y_t=0 & y_t=1 \\ \hdashline
          1 - p & p\\
          1 & 0\\
          1 & 0
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
\begin{matrix}
z_{t}=\text{in} \\ z_{t}=\text{out} \\ z_{t}=\text{D}
\end{matrix}
\end{matrix}
$$
  
  ### Combination of life and dead encounters
  
  Transition matrix

$$
  \begin{matrix}
& \\
\mathbf{\Gamma} =
  \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          z_t=A & z_t=JD & z_t=D \\ \hdashline
          s & 1-s & 0\\
          0 & 0 & 1\\
          0 & 0 & 1
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
\begin{matrix}
z_{t-1}=\text{alive} \\ z_{t-1}=\text{just dead} \\ z_{t-1}=\text{dead for good}
\end{matrix}
\end{matrix}
$$
  
  Observation matrix

$$
  \begin{matrix}
& \\
\mathbf{\Omega} =
  \left ( \vphantom{ \begin{matrix} 12 \\ 12 \\ 12\end{matrix} } \right .
          \end{matrix}
          \hspace{-1.2em}
          \begin{matrix}
          y_t=0 & y_t=1 & y_t=2 \\ \hdashline
          1 - p & 0 & p\\
          1 - r & r & 0\\
          1 & 0 & 0
          \end{matrix}
          \hspace{-0.2em}
          \begin{matrix}
          & \\
          \left . \vphantom{ \begin{matrix} 12 \\ 12 \\ 12 \end{matrix} } \right )
\begin{matrix}
z_{t}=A \\ z_{t}=JD \\ z_{t}=D
\end{matrix}
\end{matrix}
$$
  


#------------------------------------- prior elicitation

### Prior elicitation

## Why Bayes? Incorporate prior information

So far, we have assumed a non-informative prior on survival $\text{Beta}(1,1) = \text{Uniform}(0,1)$. With this prior, mean posterior survival is $\phi = 0.56$ with credible interval $[0.52,0.62]$. Graphically we may represent the posterior distribution of survival obtained with two chains with different colors, and our prior in gray dashed line:
  ```{r, echo = FALSE}
load(here::here("dat","dipper.RData"))
PR <- runif(1500, 0, 1)
MCMCtrace(mcmc.phip,
          params = c('phi'),
          ISB = FALSE,
          exact = TRUE,
          priors = PR,
          pdf = FALSE,
          Rhat = FALSE,
          n.eff = FALSE,
          #          post_zm = TRUE,
          sz_txt = NULL,
          ind = TRUE,
          type = "density",
          lwd_den = 3,
          lwd_pr = 3,
          col_pr = "gray70",
          lty_pr = 2,
          main_den = "",
          xlab_den = "survival")
```

The thing is that we know a lot about passerines and it is a shame not to be able to use this information and act as if we have to start from scratch and know nothing. We illustrate how to incorporate prior information by acknowledging that species with similar body masses have similar survival. By gathering information on several other European passerines than the dipper, let's assume we have built a regression of survival vs. body mass -- allometric relationship. Knowing dippers weigh on average 59.8g, we're now in the position to build a prior for dipper survival probability by predicting its value using the regression. We obtain a predicted survival of 0.57 and a standard deviation of 0.075. Using an informative prior $\text{Normal}(0.57, sd = 0.073)$ in NIMBLE, we get a mean posterior of $0.56$ with credible interval $[0.52, 0.61]$. There's barely no difference with the non-informative prior, quite a disappointment. 

Now let's assume that we had only the three first years of data, what would have happened? We fit the model with constant parameters with both the non-informative and informative priors to the dataset from which we delete the final 4 years of data. Now the benefit of using the prior information becomes clear as the credible interval when prior information is ignored has a width of 0.53, which is more than twice as much as when prior information is used (0.24), illustrating the increased precision provided by the prior. We may assess visually this gain in precision by comparing the survival posterior distributions with and without informative prior:
  
  ```{r, echo = FALSE}
load(here::here("dat","phip3y.RData"))
phinoprior <- c(mcmc.phip$chain1[,"phi"], mcmc.phip$chain2[,"phi"])
load(here::here("dat","phipriorp3y.RData"))
phiprior <- c(mcmc.phip$chain1[,"phi"], mcmc.phip$chain2[,"phi"])
df <- data.frame(posterior = c(phinoprior, phiprior),
                 type = c(rep("w/ vague prior", length(phinoprior)),
                          rep("w/ informative prior", length(phiprior))))
df %>%
  ggplot() +
  aes(x = posterior, fill = type) +
  geom_density(aes(y = ..density..),
               bins = 40,
               color = "white",
               alpha = 0.6) +
  labs(x = "survival", fill = "") +
  scale_fill_manual(values = wesanderson::wes_palette("Royal1")[2:1])
```

In brief, if the aim is to get an estimate of survival, Gilbert did not have to conduct further data collection after 3 years, and he could have reached the same precision as with 7 years of data by using prior information derived from body mass. In brief, the prior information was worth 4 years of field data. Of course, this is assuming that the ecological question remains the same whether you have 3 or 7 years of data, which is unlikely to be the case, as with long-term data, there is so much we can ask, more than "just" what annual survival probability is. 

### Moment matching

The prior $\text{Normal}(0.57, sd = 0.073)$ is not entirely satisfying because it is not constrained to be positive or less than one, which is the minimum for a probability (of survival) to be well defined. In our specific example, the prior distribution is centered on positive values far from 0, and the sandard deviation is small enough so that the chances to get values smaller than 0 or higher than 1 are null (to convince yourself, `hist(rnorm(1000, mean = 0.57, sd = 0.073))`). Can we do better? The answer is yes. 

Remember the Beta distribution? Recall that the Beta distribution is a continuous distribution with values between 0 and 1. It is therefore convenient to specify priors for survival and detection probabilities. Plus we know everything about the Beta distribution, in particular its moments. If $X \sim Beta(\alpha,\beta)$, then the first (mean) and second moments (variance) of $X$ are $\mu = \text{E}(X) = \frac{\alpha}{\alpha + \beta}$ and $\sigma^2 = \text{Var}(X) = \frac{\alpha\beta}{(\alpha + \beta)^2 (\alpha + \beta + 1)}$. 

In the capture-recapture example, we know a priori that the mean of the probability we're interested in is $\mu = 0.57$ and its variance is $\sigma^2 = 0.073^2$. Parameters $\mu$ and $\sigma^2$ are seen as the moments of a $Beta(\alpha,\beta)$ distribution. Now we look for values of $\alpha$ and $\beta$ that match the observed moments of the Beta distribution $\mu$ and $\sigma^2$. We need another set of equations:
$$\alpha = \bigg(\frac{1-\mu}{\sigma^2}- \frac{1}{\mu} \bigg)\mu^2$$
$$\beta = \alpha \bigg(\frac{1}{\mu}-1\bigg)$$
For our model, that means:
```{r echo = TRUE}
(alpha <- ( (1 - 0.57)/(0.073*0.073) - (1/0.57) )*0.57^2)
(beta <- alpha * ( (1/0.57) - 1))
```

Now we simply have to use $\text{Beta}(\alpha = 25.6,\beta = 19.3)$ as a prior instead of our $\text{Normal}(0.57, sd = 0.073)$. 



#---------------------------------------------------- Prior predictive checks

### Linear regression

Unreasonable prior $\beta \sim N(0, 1000^2)$
  
  ```{r, echo = FALSE}
df <- rnorm(1000, 0, 1000)
df %>%
  as_tibble() %>%
  ggplot(aes(x = value)) +
  geom_density(size = 2) +
  labs(x = "Height (m)") +
  theme_light(base_size = 20)
```

Reasonable prior $\beta \sim N(2, 0.5^2)$
  
  ```{r, echo = FALSE}
df <- rnorm(1000, 2, 0.5)
df %>%
  as_tibble() %>%
  ggplot(aes(x = value)) +
  geom_density(size = 2) +
  labs(x = "Height (m)") +
  theme_light(base_size = 20)
```

### Logistic regression

Unreasonable prior $\text{logit}(\phi) = \beta \sim N(0, 10^2)$
  
  ```{r, echo = FALSE}
df <- plogis(rnorm(1000, 0, 10))
df %>%
  as_tibble() %>%
  ggplot(aes(x = value)) +
  geom_density(size = 2) +
  labs(x = "survival") +
  theme_light(base_size = 20)
```

Reasonable prior $\text{logit}(\phi) = \beta \sim N(0, 1.5^2)$
  
  ```{r, echo = FALSE}
df <- plogis(rnorm(1000, 0, 1.5))
df %>%
  as_tibble() %>%
  ggplot(aes(x = value)) +
  geom_density(size = 2) +
  labs(x = "survival") +
  theme_light(base_size = 20)
```




##------------------------------------------------- Parameter-redundancy issue

There are two potential issues, either intrinsic or extrinsic parameter redundacy. 
Intrinsic redundancy means that the model likelihood can be expressed by a smaller 
number of parameters; it is a feature of the model. Extrinsic redundancymeans that 
model structure is fine, but the lack of data makes a parameter non-estimable; 
this is a feature of the data.

```{r, echo = FALSE}
load(here::here("dat", "profiledeviance.RData"))
df <- data.frame(last_survival = grid_lastphi,
                 max_dev = devmax)
mytable <- df %>% slice(c(20, 27, 38, 44)) %>% round(2)
ggplot() +
  geom_line(data = df,
            aes(x = last_survival, y = max_dev),
            size = 1.5,
            color = "gray70") +
  geom_point(data = df %>% slice(c(20, 27, 38, 44)),
             aes(x = last_survival, y = max_dev),
             size = 3.5,
             pch = 16,
             color = "darkblue") +
  labs(x = "survival over last time interval", y = "-log-likelihood") +
  theme_light(base_size = 14) +
  annotation_custom(gridExtra::tableGrob(mytable, rows=NULL), xmin=0.4, xmax=0.9, ymin=350, ymax=380)
```

Last survival and recapture probabilities cannot be estimated separately.

Poor mixing of the chains.

Two issues

Intrinsic redundancy: Likelihood can be expressed by a smaller number of parameters; Feature of the model

Extrinsic redundancy: Model structure is fine, But lack of data makes a parameter non-estimable, Feature of the data.

### Prior-posterior overlap for $\phi_4$ and $\phi_6$

```{r, echo = FALSE}
load(here::here("dat","dipper.RData"))
PR <- runif(1500, 0, 1)
MCMCtrace(mcmc.phitpt,
          params = c('phi[4]'),
          ISB = FALSE,
          exact = TRUE,
          priors = PR,
          pdf = FALSE,
          Rhat = FALSE,
          n.eff = FALSE,
          ind = FALSE,
          type = "density",
          lwd_den = 3,
          lwd_pr = 3,
          col_pr = "gray70",
          lty_pr = 2,
          main_den = "",
          xlab_den = "survival prob. between years 1984 and 1985",
          sz_txt = 1.8,
          sz_ax = 1.8,
          sz_ax_txt = 1.8,
          sz_tick_txt = 1.8,
          sz_main_txt = 1.8)
```


```{r, echo = FALSE}
load(here::here("dat","dipper.RData"))
PR <- runif(1500, 0, 1)
MCMCtrace(mcmc.phitpt,
          params = c('phi[6]'),
          ISB = FALSE,
          exact = TRUE,
          priors = PR,
          pdf = FALSE,
          Rhat = FALSE,
          n.eff = FALSE,
          ind = FALSE,
          type = "density",
          lwd_den = 3,
          lwd_pr = 3,
          col_pr = "gray70",
          lty_pr = 2,
          main_den = "",
          xlab_den = "survival prob. between years 1986 and 1987",
          sz_txt = 1.8,
          sz_ax = 1.8,
          sz_ax_txt = 1.8,
          sz_tick_txt = 1.8,
          sz_main_txt = 1.8)

```

### Prior-posterior overlap for $p_3$ and $p_7$

```{r, echo = FALSE}
load(here::here("dat","dipper.RData"))
PR <- runif(1500, 0, 1)
MCMCtrace(mcmc.phitpt,
          params = c('p[2]'),
          ISB = FALSE,
          exact = TRUE,
          priors = PR,
          pdf = FALSE,
          Rhat = FALSE,
          n.eff = FALSE,
          ind = FALSE,
          type = "density",
          lwd_den = 3,
          lwd_pr = 3,
          col_pr = "gray70",
          lty_pr = 2,
          main_den = "",
          xlab_den = "recapture prob. at year 1983",
          sz_txt = 1.8,
          sz_ax = 1.8,
          sz_ax_txt = 1.8,
          sz_tick_txt = 1.8,
          sz_main_txt = 1.8)

```

```{r, echo = FALSE}
load(here::here("dat","dipper.RData"))
PR <- runif(1500, 0, 1)
MCMCtrace(mcmc.phitpt,
          params = c('p[6]'),
          ISB = FALSE,
          exact = TRUE,
          priors = PR,
          pdf = FALSE,
          Rhat = FALSE,
          n.eff = FALSE,
          ind = FALSE,
          type = "density",
          lwd_den = 3,
          lwd_pr = 3,
          col_pr = "gray70",
          lty_pr = 2,
          main_den = "",
          xlab_den = "recapture prob. at year 1987",
          sz_txt = 1.8,
          sz_ax = 1.8,
          sz_ax_txt = 1.8,
          sz_tick_txt = 1.8,
          sz_main_txt = 1.8)

```



#----------------------------------------------------- WAIC

Models with smaller DIC values represent more parsimonious descriptions 
of the data than models with larger DIC values. The relative parsimony 
of the models can be assessed by comparing the difference in the DIC values. 
The DIC is a Bayesian equivalent of Akaike's information criterion (AIC; 
Akaike 1973; Hilborn & Mangel 1997; Burnham & Anderson 2002) and the rules 
of thumb suggested by Burnham & Anderson (2002) for comparing models with 
AIC seem to apply for DIC (Spiegelhalter et al. 2002). Therefore, differences 
of less than 2 indicate that the two models are indistinguishable, differences 
of 4–7 that the poorer model has considerably less support, and differences of 
more than 10 that the poorer model 
has essentially no support (Burnham & Anderson 2002).

The effective number of parameterŝ pwaic can be used as measure of 
complexity of the model, but it should not be overinterpreted, as the 
original goal is to estimate the difference between lpd and elpd.

This selection is based on assessing the trade-off between the fit and 
complexity of the models, with the aim to find the most parsimonious model 
of set of models. This is what the Akaike information criterion (AIC) is 
trying to achieve, with $AIC = - 2 \log(L(\hat{\theta}_1,\ldots,\hat{\theta}_K)) 
+ 2 K$ where $L$ the likelihood and $K$ the number of parameters $\theta_i$. 
First term is a measure of goodness-of-fit of the model to the data: the more 
parameters you have, the smaller the deviance is (or the bigger the likelihood is). 
Second term is a penalty: twice the number of parameters $K$. AIC makes the 
balance between *quality of fit* and *complexity* of a model. Best model is 
the one with lowest AIC value.

In Bayesian statistics, the relative parsimony of the models can be compared 
using the deviance information criterion (DIC). DIC is a measure of the fit of 
the model to the data that is penalized for the model's complexity. The measure 
of fit is based on the likelihood of obtaining the observed data given the means 
of the posterior distribution of the parameters. Parameter values provide a better 
fit if they are more likely to have produced the observed data. The complexity 
of the model is measured by the effective number of estimated parameters. 
**difference between AIC and DIC, and problems with DIC.**
  
  Bayesian version exists with Watanabe-Akaike (Widely-Applicable) Information 
Criteria or WAIC (Widely Applicable Information Criterion) given by 
$\textrm{WAIC} = -2 \sum_{i = 1}^n \log E[\Pr(y_i \mid \theta)] + 2 p_\text{WAIC}$ 
  where $E[p(y_i \mid \theta)]$ is the posterior mean of the likelihood evaluated 
pointwise at each $i$th observation, and $p_\text{WAIC}$ is a penalty computed 
using the posterior variance of the likelihood.


Even if all of the models being considered have mismatches with the data 
(forward reference to gof section, say relative fit vs absolute fit), it can 
be informative to evaluate their predictive accuracy, compare them, and consider 
where to go next. The challenge then is to estimate predictive model accuracy, 
correcting for the bias inherent in evaluating a model’s predictions of the data 
that were used to fit it.

A more general summary of predictive fit is the log predictive density, log p(y|), 
which is proportional to the mean squared error if the model is normal with constant 
variance. The log predictive density is also sometimes called the log-likelihood. 
The log predictive density has an important role in statistical model comparison 
because of its connection to the Kullback-Leibler information measure (see Burnham 
and Anderson, 2002, and Robert, 1996).

In the limit of large sample sizes, the model with the lowest Kullback-Leibler 
information—and thus, the highest expected log predictive density—will have the 
highest posterior probability. Thus, it seems reasonable to use expected log 
predictive density as a measure of overall model fit.

AIC uses the maximum likelihood as a measure of goodness-of-fit, and the number 
of free parameters as a measure of flexibility, with more parameters resulting in 
harsher penalties. DIC uses the average log-likelihood over the posterior distribution 
as a measure of goodness-of-fit, and the difference between this average and the 
log-likelihood at some fixed, central point of the posterior as a measure of 
flexibility, with greater differences resulting in harsher penalties. Although 
the mean of the parameter values over the joint posterior is often used as the 
point estimate in this calculation (Spiegelhalter et al., 2002), I instead use 
the point of minimum deviance in the posterior for the point estimate (also 
recommended by Spiegelhalter et al., 2002), as the use of the mean results in 
the strong assumption that the joint posterior distribution is a multivariate 
normal, and can result in negative estimates of flexibility when this assumption 
is violated (Vehtari et al., 2017). WAIC uses a similar measure of goodness-of-fit 
as DIC, being the log of the average posterior likelihood for each data point, 
but uses the variance in log-likelihood over the posterior distribution as a 
measure of flexibility, with greater variances resulting in harsher penalties.

A natural way to estimate out-of-sample prediction error is cross-validation (see Geisser and
Eddy, 1979, and Vehtari and Lampinen, 2002, for a Bayesian perspective), but researchers have
always sought alternative measures, as cross-validation requires repeated model fits and can run
into trouble with sparse data. For practical reasons alone, there remains a place for simple bias
corrections such as AIC (Akaike, 1973), DIC (Spiegelhalter et al., 2002, van der Linde, 2005), and,
more recently, WAIC (Watanabe, 2010), and all these can be viewed as approximations to different
versions of cross-validation (Stone, 1977).

All the different measures discussed above are based on adjusting the log predictive density of the
observed data by subtracting an approximate bias correction. The measures differ both in their
starting points and in their adjustments.

These methods penalize models for overfitting, as over-fitting will 
result in poor prediction of future data.

In addition to varying in computational tractability, these methods vary on their 
theoretical basis for selecting between models. AIC, DIC, and WAIC all attempt 
to find the model with the best “predictive accuracy”, which is the model that 
is able to best predict future empirical data, given some fixed set of parameter 
values or distributions (Akaike, 1974; Spiegelhalter et al., 2002; Vehtari et al., 2017).  



#------------------------------------------ model validation

The fit of model to data can be assessed using posterior predictive checks (Rubin,
1984), prior predictive checks (when evaluating potential replications involving new parameter
values)


#------------------------------------------- covariates

Proportion of variance explained. Path analyses, structural equation
models. Splines (more about spatial stats? CAR model?). Imputation 
and multistate models to account for missing data. Explain basics 
of parametric statistical modeling (linear models, GLMs and random effects).

Talk about ANODEV?
  
