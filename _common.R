# packages
library(tidyverse)
theme_set(theme_light())
library(nimble)
library(MCMCvis)
library(magick)
library(pdftools)
library(wesanderson)
library(RColorBrewer)
library(hrbrthemes)
#library(nimbleEcology)
#library(basicMCMCplots)

# R options
options(width = 60)

# chunk options
knitr::opts_chunk$set(
  comment = "##",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
  )
