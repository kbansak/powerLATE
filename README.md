
<!-- README.md is generated from README.Rmd. Please edit that file -->

# powerLATE

<!-- badges: start -->

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
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
