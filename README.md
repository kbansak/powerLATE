
<!-- README.md is generated from README.Rmd. Please edit that file -->

# powerLATE

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/powerLATE?color=green)](https://cran.r-project.org/package=powerLATE)
[![Codecov test
coverage](https://codecov.io/gh/kbansak/powerLATE/branch/master/graph/badge.svg)](https://codecov.io/gh/kbansak/powerLATE?branch=master)
[![R-CMD-check](https://github.com/kbansak/powerLATE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kbansak/powerLATE/actions/workflows/R-CMD-check.yaml)
[![](http://cranlogs.r-pkg.org/badges/grand-total/powerLATE?color=blue)](https://cran.r-project.org/package=powerLATE)
<!-- badges: end -->

powerLATE implements the generalized power analysis for the local
average treatment effect (LATE), proposed by [Bansak
(2020)](https://projecteuclid.org/download/pdfview_1/euclid.ss/1591171230).
A comprehensive tutorial on using this package [can be found
here](https://github.com/kbansak/powerLATE_tutorial).

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
#> powerLATE: Generalized Power Analysis for LATE
#> Version: 0.1.2
#> Reference: Bansak, K. (2020). A Generalized Approach to Power Analysis for Local Average Treatment Effects. Statistical Science, 35(2), 254-271.
```

powerLATE provides two main functions:

- `powerLATE()`, which computes the power of the Wald IV estimator, or
  determines parameters (e.g. required sample size) to obtain a target
  power.

- `powerLATE.cov()`, which is similar to `powerLATE()` but additionally
  allows the inclusion of covariates.

## Navigating main functions

![](https://github.com/EddieYang211/powerLATE_aux/blob/master/powerLATE_tree.png?raw=true)

## Tutorial

For a comprehensive tutorial on conducting a LATE power analysis with
this package, see [here](https://github.com/kbansak/powerLATE_tutorial).

## More Examples

For more examples on how to use the package, see
[here](https://htmlpreview.github.io/?https://github.com/EddieYang211/powerLATE_aux/blob/master/powerLATE_Examples.html).

## Reference

For a detailed description of the method see:

- [Bansak, K. (2020). A Generalized Approach to Power Analysis for Local
  Average Treatment
  Effects](https://projecteuclid.org/download/pdfview_1/euclid.ss/1591171230)

## Maintainer

- [Eddie Yang](https://github.com/EddieYang211)
