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
