
<!-- README.md is generated from README.Rmd. Please edit that file -->

# powerLATE

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/EddieYang211/powerLATE.svg?branch=master)](https://travis-ci.com/EddieYang211/powerLATE)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/EddieYang211/powerLATE?branch=master&svg=true)](https://ci.appveyor.com/project/EddieYang211/powerLATE)
[![Codecov test
coverage](https://codecov.io/gh/EddieYang211/powerLATE/branch/master/graph/badge.svg)](https://codecov.io/gh/EddieYang211/powerLATE?branch=master)
<!-- badges: end -->

The goal of powerLATE is to provide a generalized approach to power
analysis in the context of estimating a local average treatment effect
(LATE), where the study subjects exhibit noncompliance with treatment
assignment.

powerLATE uses standardized effect sizes to place bounds on the power
for the most commonly used estimator of the LATE, the Wald IV estimator,
whereby variance terms and distributional parameters need not be
specified nor assumed. Instead, in addition to the effect size, sample
size, and error tolerance parameters, the only other parameter that must
be specified by the researcher is the compliance rate. Additional
conditions can also be introduced to further narrow the bounds on the
power calculation.

## Installation

You can install the released version of powerLATE from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("powerLATE")
```

Or the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("EddieYang211/powerLATE")
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

## Examples

Suppose we are interested in the power of an experiment, in which
exactly half of the experimental subjects will be assigned to treatment
and the other half to control. Further suppose the sample size is 800,
the compliance rate is 0.35 and the treatment effect size is 0.8. We can
compute the power of this experiment by the following:

``` r
library(powerLATE)
res <- powerLATE(pZ = 0.5, N = 800, pi = 0.35, kappa =0.8)
#> Power analysis for two-sided test that LATE equals zero
#> 
#>  pZ = 0.5
#>  pi = 0.35
#>  N = 800
#>  kappa = 0.8
#>  sig.level  = 0.05
#> 
#> Given these parameter values, the conservative bound for Power:
#>     power 
#> 0.8213495 
#> 
#> NOTE: The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.
```

If absolute effect is desired instead of effect size, simply set
`effect.size=FALSE` and supply `omega`, the pooled standard deviation
across treatment groups, and `tau`, the absolute effect.

``` r
res <- powerLATE(pZ=0.5, N=800, pi=0.35, omega = 0.6, tau = 0.4, effect.size = FALSE)
#> Power analysis for two-sided test that LATE equals zero
#> 
#>  pZ = 0.5
#>  pi = 0.35
#>  N = 800
#>  tau = 0.4
#>  omega = 0.6
#>  sig.level  = 0.05
#> 
#> Given these parameter values, the conservative bound for Power:
#>     power 
#> 0.7104445 
#> 
#> NOTE: The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.
```

We can also input a vector of values for any one parameter in the set
{N, kappa, power, pi, tau, omega}.

``` r
res <- powerLATE(pZ = 0.5, N = 800, pi = 0.35, kappa = seq(0.5, 1.0, 0.1))
#> Power analysis for two-sided test that LATE equals zero
#> 
#>  pZ = 0.5
#>  pi = 0.35
#>  N = 800
#>  kappa = Multiple values inputted (see table below)
#>  sig.level  = 0.05
#> 
#> Given these parameter values, the conservative bound for Power:
#>   power     User-inputted kappa
#> 1 0.5181033 0.5                
#> 2 0.6399776 0.6                
#> 3 0.7419496 0.7                
#> 4 0.8213495 0.8                
#> 5 0.8797640 0.9                
#> 6 0.9208686 1.0                
#> 
#> NOTE: The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.
```

With covariates,

``` r
res <- powerLATE.cov(pZ = 0.5, N = 800, pi = 0.35, kappa = seq(0.5, 1.0, 0.1), r2dw = 0.2, r2yw = 0.3)
#> Power analysis for two-sided test that LATE equals zero
#> 
#>  pZ = 0.5
#>  pi = 0.35
#>  N = 800
#>  kappa = Multiple values inputted (see table below)
#>  r2dw = 0.2
#>  r2yw = 0.3
#>  sig.level  = 0.05
#> 
#> Given these parameter values, the conservative bound for Power:
#>   power     User-inputted kappa
#> 1 0.5522840 0.5                
#> 2 0.6757927 0.6                
#> 3 0.7755263 0.7                
#> 4 0.8502388 0.8                
#> 5 0.9030156 0.9                
#> 6 0.9386345 1.0                
#> 
#> NOTE: The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.
```

In addition to the message that is automatically printed out,
`powerLATE` and `powerLATE.cov` also return the output parameter values
as an invisible object, simply query

``` r
res <- powerLATE(pZ = 0.5, N = 800, pi = 0.35, kappa = seq(0.5, 1.0, 0.1))
#> Power analysis for two-sided test that LATE equals zero
#> 
#>  pZ = 0.5
#>  pi = 0.35
#>  N = 800
#>  kappa = Multiple values inputted (see table below)
#>  sig.level  = 0.05
#> 
#> Given these parameter values, the conservative bound for Power:
#>   power     User-inputted kappa
#> 1 0.5181033 0.5                
#> 2 0.6399776 0.6                
#> 3 0.7419496 0.7                
#> 4 0.8213495 0.8                
#> 5 0.8797640 0.9                
#> 6 0.9208686 1.0                
#> 
#> NOTE: The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.
res$output.parameter
#>            [,1]      [,2]      [,3]      [,4]     [,5]      [,6]
#> power 0.5181033 0.6399776 0.7419496 0.8213495 0.879764 0.9208686
```

## Reference

For a detailed description of the method see:

  - [Bansak, K. (2020). A Generalized Approach to Power Analysis for
    Local Average Treatment
    Effects](https://projecteuclid.org/download/pdfview_1/euclid.ss/1591171230)

## Maintainer

  - [Eddie Yang](https://github.com/EddieYang211)
