# Generalized Linear Mixed Model Tutorial in R

This repository contains a (relatively) brief tutorial on generalized linear mixed models (GLMMs) using R to fit and compare models. The general content of the tutorial was inspired by Richard McElreath's excellent statistics course, Statistical Rethinking. The most current take on this material can be found in Richard's [textbook of the same name](http://xcelab.net/rm/statistical-rethinking/). In particular, I wrote this script borrowing ideas from a series of problems that appeared on the course's final exam. These exercises seemed particularly enlightening to me because they illustrate that inclusion of random effects (aka varying effects) can not only change relative model rankings but also emphasize that adding in random effects can drastically alter our estimates of fixed effects (i.e., the things we usually care about most in our models). This tutorial uses the R packages `lme4`, `AICcmodavg`, and `rethinking`. Akaike's Information Criterion (AIC) is used to compare fit models.

### Repository Contents

- The `glmm_tutorial_script.R` file contains my code and tutorial commentary
- The `glmm_tutorial_data.csv` file contains the example data I created for use in this tutorial
