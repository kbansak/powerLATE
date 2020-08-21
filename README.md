
<!-- README.md is generated from README.Rmd. Please edit that file -->

# powerLATE

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/powerLATE?color=green)](https://cran.r-project.org/package=powerLATE)
[![Travis build
status](https://travis-ci.com/kbansak/powerLATE.svg?branch=master)](https://travis-ci.com/kbansak/powerLATE)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/kbansak/powerLATE?branch=master&svg=true)](https://ci.appveyor.com/project/kbansak/powerLATE)
[![Codecov test
coverage](https://codecov.io/gh/kbansak/powerLATE/branch/master/graph/badge.svg)](https://codecov.io/gh/kbansak/powerLATE?branch=master)
[![](http://cranlogs.r-pkg.org/badges/grand-total/powerLATE?color=blue)](https://cran.r-project.org/package=powerLATE)
<!-- badges: end -->

powerLATE implements the generalized power analysis for the local
average treatment effect (LATE), proposed by [Bansak
(2020)](https://projecteuclid.org/download/pdfview_1/euclid.ss/1591171230).

Power analysis is in the context of estimating the LATE (also known as
the complier average causal effect, or CACE), with calculations based on
a test of the null hypothesis that the LATE equals 0 with a two-sided
alternative. The method uses standardized effect sizes to place a
conservative bound on the power under minimal assumptions. powerLATE
allows users to recover power, sample size requirements, or minimum
detectable effect sizes. It also allows users to work with absolute
effects rather than effect sizes, to specify an additional assumption to
narrow the bounds, and to incorporate covariate adjustment.

## Installation

You can install the released version of powerLATE from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("powerLATE")
```

Or the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("kbansak/powerLATE")
```

## Getting started

``` r
library(powerLATE)
```

powerLATE provides two main functions:

  - `powerLATE`, which computes the power of the Wald IV estimator, or
    determines parameters to obtain a target power.

  - `powerLATE.cov`, which is similar to `powerLATE` but additioanlly
    allows the inclusion of covariates.

## Navigating main functions

![](https://github.com/kbansak/powerLATE/blob/master/aux/powerLATE_tree.png)

For examples on how to use the package, see
[here](https://htmlpreview.github.io/?https://github.com/kbansak/powerLATE/blob/master/aux/powerLATE_Examples.html)

## Reference

For a detailed description of the method see:

  - [Bansak, K. (2020). A Generalized Approach to Power Analysis for
    Local Average Treatment
    Effects](https://projecteuclid.org/download/pdfview_1/euclid.ss/1591171230)

## Maintainer

  - [Eddie Yang](https://github.com/EddieYang211)
