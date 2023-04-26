---
title: "Bayesian Analysis of Capture-Recapture Data with Hidden Markov Models"
subtitle: "Theory and Case Studies in R"
author: "Olivier Gimenez"
date: "2023-04-26"
documentclass: krantz
bibliography: [book.bib]
#biblio-style: apalike
link-citations: yes
colorlinks: yes
lot: yes
lof: yes
fontsize: 12pt
site: bookdown::bookdown_site
description: "This is a textbook on the analysis of capture-recapture data with hidden Markov models (HMM) implemented in the Bayesian framework with R."
url: 'https\://oliviergimenez.github.io/bayesian-cr-workshop/'
github-repo: oliviergimenez/banana-book
header-includes: 
  - \usepackage{tikz}
  - \usepackage{pgfplots}
  - \usepackage{blkarray}
---




# Welcome {-}

<!-- bookdown::render_book("index.Rmd", "bookdown::pdf_book") -->

Welcome to the online version of the book *Bayesian Analysis of Capture-Recapture Data with Hidden Markov Models – Theory and Case Studies in R*. <!-- The book is also available in [PDF format](https://github.com/oliviergimenez/banana-book/raw/master/docs/bayesHMMcapturerecapture.pdf). -->

The HMM framework has gained much attention in the ecological literature over the last decade, and has been suggested as a general modelling framework for the demography of plant and animal populations. In particular, HMMs are increasingly used to analyse capture-recapture data and estimate key population parameters (e.g., survival, dispersal, recruitment or abundance) with applications in all fields of ecology. 

In parallel, Bayesian statistics is well established and fast growing in ecology and related disciplines, because it resonates with scientific reasoning and allows accommodating uncertainty smoothly. The popularity of Bayesian statistics also comes from the availability of free pieces of software (WinBUGS, OpenBUGS, JAGS, Stan, NIMBLE) that allow practitioners to code their own analyses.

This book offers a Bayesian treatment of HMMs applied to capture-recapture data. You will learn to use the R package NIMBLE which is seen by many as the future of Bayesian statistical ecology to deal with complex models and/or big data. An important part of the book consists in case studies presented in a tutorial style to abide by the “learning by doing” philosophy.

I'm currently writing this book, and I welcome any feedback. You may raise an issue [here](https://github.com/oliviergimenez/banana-book/issues), amend directly the R Markdown file that generated the page you're reading by clicking on the 'Edit this page' icon in the right panel, or [email me](mailto:olivier.gimenez@cefe.cnrs.fr). Many thanks!

Olivier Gimenez, Montpellier, France  
Last updated: April 26, 2023

## License {-}

The online version of this book is licensed under the [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](http://creativecommons.org/licenses/by-nc-nd/4.0/). 

The code is public domain, licensed under [Creative Commons CC0 1.0 Universal (CC0 1.0)](https://creativecommons.org/publicdomain/zero/1.0/).


<!--chapter:end:index.Rmd-->


# Preface {-}

## Why this book? {-}

**To be completed.** Why and what of capture-recapture data and models, with fields of application.^[Watch out nice Johnny Ball's video https://www.youtube.com/watch?v=tyX79mPm2xY.] Brief history of capture-recapture, with switch to state-space/hidden Markov model (HMM) formulation. Flexibility of HMM to decompose complex problems in smaller pieces that are easier to understand, model and analyse. From satellite guidance to conservation of endangered species. Why Bayes? Also three of my fav research topics -- capture-recapture, HMM and Bayes statistics -- let's enjoy this great cocktail together. 

## Who should read this book? {-}

This book is aimed at beginners who're comfortable using R and write basic code (including loops), as well as connoisseurs of capture-recapture who'd like to tap into the power of the Bayesian side of statistics. For both audiences, thinking in the HMM framework will help you in confidently building models and make the most of your capture-recapture data. 

## What will you learn? {-}

The book is divided into five parts. The first part is aimed at getting you up-to-speed with Bayesian statistics, NIMBLE, and hidden Markov models. The second part will teach you all about capture-recapture models for open populations, with reproducible R code to ease the learning process. In the third part, we will focus on issues in inferring states (dealing with uncertainty in assignment, modelling waiting time distribution). The fourth part provides real-world case studies from the scientific literature that you can reproduce using material covered in previous chapters. These problems can either i) be used to cement and deepen your understanding of methods and models, ii) be adapted for your own purpose, or iii) serve as teaching projects. The fifth and last chapter closes the book with take-home messages and recommendations, a list of frequently asked questions and references cited in the book. **Likely to be amended after feedbacks.**

## What won't you learn? {-}

There is hardly any maths in this book. The equations I use are either simple enough to be understood without a background in maths, or can be skipped without prejudice. I do not cover Bayesian statistics or even hidden Markov models fully, I provide just what you need to work with capture-recapture data. If you are interested in knowing more about these topics, hopefully the section Suggested reading at the end of each chapter will put you in the right direction. There are also a number of important topics specific to capture-recapture that I do not cover, including closed-population capture-recapture models [@WilliamsEtAl2002], and spatial capture-recapture models [@RoyleEtAl2013book]. These models can be treated as HMMs, but for now the usual formulation is just fine.  **There will be spatial considerations in the Covariates chapter w/ splines and CAR. I'm not sure yet about SCR models (R. Glennie's Biometrics paper on HMMs and open pop SCR will not be easy to Bayes transform and implement in NIMBLE).**

## Prerequisites {-}

This book uses primarily the R package NIMBLE, so you need to install at least R and NIMBLE. A bunch of other R packages are used. You can install them all at once by running:




```r
install.packages(c(
  "magick", "MCMCvis", "nimble", "pdftools", 
  "tidyverse", "wesanderson" 
))
```

## Acknowledgements {-}

**To be completed.**

## How this book was written {-}

I am writing this book in [RStudio](http://www.rstudio.com/ide/) using [bookdown](http://bookdown.org/). The [book website](https://oliviergimenez.github.io/banana-book) is hosted with [GitHub Pages](https://pages.github.com/), and automatically updated after every push by [Github Actions](https://github.com/features/actions). The source is available from [GitHub](https://github.com/oliviergimenez/banana-book).

The version of the book you're reading was built with R version 4.2.3 (2023-03-15) and the following packages:


|package     |version |source                                                               |
|:-----------|:-------|:--------------------------------------------------------------------|
|magick      |2.7.3   |CRAN (R 4.2.0)                                                       |
|MCMCvis     |0.15.5  |CRAN (R 4.2.0)                                                       |
|nimble      |0.12.3  |Github (nimble-dev/nimble\@6992a0db8d4dca99b85fb6865094c6ac38ff3160) |
|pdftools    |3.3.2   |CRAN (R 4.2.0)                                                       |
|tidyverse   |1.3.2   |CRAN (R 4.2.0)                                                       |
|wesanderson |0.3.6   |CRAN (R 4.2.0)                                                       |





<!--chapter:end:preface.Rmd-->


# About the author {-}

My name is Olivier Gimenez (https://oliviergimenez.github.io/). I am a senior (euphemism for not so young anymore) scientist at the National Centre for Scientific Research (CNRS) in the beautiful city of Montpellier, France. 

I struggled studying maths, obtained a PhD in applied statistics a long time ago in a galaxy of wine and cheese. I was awarded my habilitation (https://en.wikipedia.org/wiki/Habilitation) in ecology and evolution so that I could stop pretending to understand what my colleagues were talking about. More recently I embarked in sociology studies because hey, why not. 

Lost somewhere at the interface of animal ecology, statistical modeling and social sciences, my so-called expertise lies in population dynamics and species distribution modeling to address questions in ecology and conservation biology about the impact of human activities and the management of large carnivores. I would be nothing without the students and colleagues who are kind enough to bear with me.

You may find me on Twitter (https://twitter.com/oaggimenez), GitHub (https://github.com/oliviergimenez), or get in touch [by email](mailto:olivier.gimenez@cefe.cnrs.fr).

<!--chapter:end:author.Rmd-->


\mainmatter

# (PART) I. Foundations {-}

# Introduction {-}


<!--chapter:end:introductionpartone.Rmd-->


# NIMBLE tutorial {#intronimble}

## Introduction

In this second chapter, you will get familiar with NIMBLE, an R package that implements up-to-date MCMC algorithms for fitting complex models. NIMBLE spares you from coding the MCMC algorithms by hand, and requires only the specification of a likelihood and priors for model parameters. We will illustrate NIMBLE main features with a simple example, but the ideas hold for other problems.

## What is NIMBLE?

NIMBLE stands for **N**umerical **I**nference for statistical **M**odels using **B**ayesian and **L**ikelihood **E**stimation. Briefly speaking, NIMBLE is an R package that implements for you MCMC algorithms to generate samples from the posterior distribution of model parameters. Freed from the burden of coding your own MCMC algorithms, you only have to specify a likelihood and priors to apply the Bayes theorem. To do so, NIMBLE uses a syntax very similar to the R syntax, which should make your life easier. This so-called BUGS language is also used by other programs like WinBUGS, OpenBUGS, and JAGS. 

So why use NIMBLE you may ask? The short answer is that NIMBLE is capable of so much more than just running MCMC algorithms! First, you will work from within R, but in the background NIMBLE will translate your code in C++ for (in general) faster computation. Second, NIMBLE extends the BUGS language for writing new functions and distributions of your own, or borrow those written by others. Third, NIMBLE gives you full control of the MCMC samplers, and you may pick other algorithms than the defaults. Fourth, NIMBLE comes with a library of numerical methods other than MCMC algorithms, including sequential Monte Carlo (for particle filtering) and Monte Carlo Expectation Maximization (for maximum likelihood). Last but not least, the development team is friendly and helpful, and based on users' feedbacks, NIMBLE folks work constantly at improving the package capabilities. 

<img src="images/nimble-icon.png" alt="Logo of the NIMBLE R package designed by Luke Larson. **Ask Perry for context and meaning.**" width="50%" style="display: block; margin: auto;" />

<!-- Why NIMBLE over Stan? i) The BUGS language is cool, ii) discrete latent states easier to deal with NIMBLE, no need to marginalise like with Stan (ref forward to relevant section of the book for marginalization in NIMBLE), iii) also HMC is on its way in NIMBLE, so NIMBLE includes STAN and has so much more to offer the users.  -->

## Getting started {#start-nimble}





































































































































































