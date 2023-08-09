## Prior predictive checks

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
