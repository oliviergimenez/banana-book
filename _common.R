# packages
library(tidyverse)
theme_set(theme_light())
library(nimble)
library(MCMCvis)
library(magick)
library(pdftools)
library(wesanderson)

# example R options set globally
options(width = 60)

# example chunk options set globally
knitr::opts_chunk$set(
  comment = "##",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
  )
