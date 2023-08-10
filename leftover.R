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


#------------------------------------------ model validation

The fit of model to data can be assessed using posterior predictive checks (Rubin,
1984), prior predictive checks (when evaluating potential replications involving new parameter
values)
