# packages
library(tidyverse)
theme_set(theme_light(base_size = 14))
library(nimble)
library(MCMCvis)
library(magick)
library(pdftools)
library(wesanderson)
library(RColorBrewer)
library(patchwork)
library(emo)
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
