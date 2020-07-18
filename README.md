
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

The goal of powerLATE is to provide a generalized approach power
analysis in the context of estimating a local average treatment effect
(LATE), where the study subjects exhibit noncompliance with treatment
assignment. powerLATE uses standardized effect sizes to place bounds on
the power for the most commonly used estimator of the LATE, the Wald IV
estimator, whereby variance terms and distributional parameters need not
be specified nor assumed. Instead, in addition to the effect size,
sample size, and error tolerance parameters, the only other parameter
that must be specified by the researcher is the compliance rate.
Additional conditions can also be introduced to further narrow the
bounds on the power calculation.

## Installation

You can install the released version of powerLATE from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("powerLATE")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(powerLATE)
## basic example code
res <- powerLATE.cov(pZ=0.5, pi=0.35, kappa=seq(0.4, 1.0, 0.1), power=0.8,
                     r2dw=0.2, r2yw=0.15)
#> Power analysis for two-sided test that LATE equals zero
#> 
#>  pZ = 0.5
#>  pi = 0.35
#>  kappa = Multiple values inputted (see table below)
#>  Power = 0.8
#>  r2dw = 0.2
#>  r2yw = 0.15
#>  sig.level  = 0.05
#> 
#> Given these parameter values, the conservative bound for N (required sample size):
#>   N         User-inputted kappa
#> 1 2201.1650 0.4                
#> 2 1521.2035 0.5                
#> 3 1137.4857 0.6                
#> 4  897.4863 0.7                
#> 5  736.1270 0.8                
#> 6  621.6712 0.9                
#> 7  537.0657 1.0                
#> 
#> NOTE: The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.

# also returns an invisible object: output.parameter
res$output.parameter
#>    output.1  output.2  output.3  output.4  output.5  output.6  output.7
#> N 2201.1650 1521.2035 1137.4857  897.4863  736.1270  621.6712  537.0657

res <- powerLATE.cov(pZ=0.67, pi=0.35, kappa=seq(0.4, 1.0, 0.1), 
                power=0.8, assume.ord.means=TRUE, r2dw=0.2, r2yw=0.15)
#> Power analysis for two-sided test that LATE equals zero
#> 
#>  pZ = 0.67
#>  pi = 0.35
#>  kappa = Multiple values inputted (see table below)
#>  Power = 0.8
#>  r2dw = 0.2
#>  r2yw = 0.15
#>  sig.level  = 0.05
#> 
#> Given these parameter values, the conservative bound for N (required sample size):
#>   N         User-inputted kappa
#> 1 1839.9830 0.4                
#> 2 1202.6270 0.5                
#> 3  856.4089 0.6                
#> 4  647.6503 0.7                
#> 5  512.1579 0.8                
#> 6  419.2648 0.9                
#> 7  352.8189 1.0                
#> 
#> NOTE: The Ordered-Means assumption is being employed. User should confirm that the assumption is reasonable in the context of interest. The Homoskedasticity assumption is currently being made because pZ does not equal 0.5.
```
