#--------------------------- chapitre Extensions dans Transitions part

## Hidden semi-Markov models

See @ChoquetEtAl2011, choquet2013stopover, @choquet2014stopover, @guerin2017stopover and @king2016semi. 
The latter illustrates the approach with the house finch data which we use in the book. 
Ruth and Roland say: "We extend the Arnason–Schwarz model by specifying a semi-Markov model for the state 
process, where the dwell-time distribution is specified more generally, using, for example, a shifted 
Poisson or negative binomial distribution. A state expansion technique is applied in order to represent 
the resulting semi-Markov Arnason–Schwarz model in terms of a simpler and computationally tractable 
hidden Markov model. Semi-Markov Arnason–Schwarz models come with only a very modest increase in the 
number of parameters, yet permit a significantly more flexible state process.". The expansion method 
is from @langrock2011expansion, and is very well explained in Chapter 12 of @ZucchiniEtAl2016. Totally 
doable I guess. I think to remember that in house finch app, only one state has dwell time non-geometric, 
and if so, see Section 12.3.3 in @ZucchiniEtAl2016. There is even an R function to convert the HSMM into HMM.

## Continuous-time HMM

Shouldn't be too difficult to code the matrix exponential, or to use `mexp()` with a `nimbleRcall()`. 
The likelihood is the same basically, with the distribution of the dwell-time in there.


#---------------------------- case studies à voir

## Dependence among individuals

@culina_multievent_2013 and @cubaynes_modeling_2021

## Others

Multispecies? Phylogeny? Social networks? Path analysis? Structural Equation Modelling? 

## Prevalence 

Prevalence estimation with hybrid (@SantostasiEtAl2019) or sex-ratio (@pradel2008sex) example. 
Insist on prop in newly marked, how to go to proportion in whole population. See formula derived for C. Duchamp.

## Covariate on multinomial logit link or Dirichlet

Example?
  

#---------------------------- Dirichelt

```{r dirichlet, echo = FALSE, fig.cap = "Dirichlet prior with parameter alpha"}
library(gtools) # to make the rdirichlet() function available
library(ggtern) # to visually represent multidim prob distribution 
set.seed(123)
n <- 500
alpha1 <- c(1, 1, 1)
p1 <- gtools::rdirichlet(n, alpha1)
alpha2 <- c(5, 5, 5)
p2 <- rdirichlet(n, alpha2)
alpha3 <- c(1, 2, 2)
p3 <- rdirichlet(n, alpha3)
alpha4 <- c(2, 4, 8)
p4 <- rdirichlet(n, alpha4)
df <- cbind(rbind(p1, p2, p3, p4), c(rep("alpha = c(1, 1, 1)", n),
                                     rep("alpha = c(5, 5, 5)", n),
                                     rep("alpha = c(1, 2, 2)", n),
                                     rep("alpha = c(2, 4, 8)", n))) %>%
  as_tibble() %>%
  mutate(x = as.numeric(V1),
         y = as.numeric(V2),
         z = as.numeric(V3),
         alpha = V4)

df %>%
  ggtern(aes(x = x, y = y, z = z)) +
  stat_density_tern(aes(fill=..level.., alpha=..level..),
                    geom = 'polygon',
                    bdl = 0.005) + # a 2D kernel density estimation of the distribution
  scale_fill_viridis_b() +
  #  geom_point(alpha = 0.3, pch = "+") +
  theme_light(base_size = 14) +
  #  theme_showarrows() +
  #  scale_T_continuous(breaks = seq(0, 1, by = 0.2),
  #                     labels = seq(0, 1, by = 0.2)) +
  #  scale_L_continuous(breaks = seq(0, 1, by = 0.2),
  #                     labels = seq(0, 1, by = 0.2)) +
  #  scale_R_continuous(breaks = seq(0, 1, by = 0.2),
  #                     labels = seq(0, 1, by = 0.2)) +
  #  labs(x = "",
  #       y = "",
  #       z = "",
  #       Tarrow = "psiAA",
#       Larrow = "psiAB",
#       Rarrow = "psiAC") +
guides(color = "none", fill = "none", alpha = "none") +
  facet_wrap(~alpha)
```


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



#------------------------------------ flexibility of multistate



#------------------------------------- prior elicitation


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
  
