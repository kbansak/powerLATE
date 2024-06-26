% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unequal_ordered_cov.R
\name{unequal.ordered.cov}
\alias{unequal.ordered.cov}
\title{Subsidiary powerLATE Function}
\usage{
unequal.ordered.cov(
  power = NULL,
  sig.level = NULL,
  pi = NULL,
  kappa = NULL,
  N = NULL,
  pZ = NULL,
  r2dw = NULL,
  r2yw = NULL
)
}
\arguments{
\item{power}{power of test (1 minus Type II error probability)}

\item{sig.level}{significance level (Type I error probability).}

\item{pi}{compliance rate. Equivalently, average causal effect of Z on D.}

\item{kappa}{effect size}

\item{N}{total number of observations}

\item{pZ}{probability of being assigned to treatment.}

\item{r2dw}{proportion of variation in D left unexplained by Z that is explained by W.}

\item{r2yw}{proportion of variation in Y left unexplained by Z that is explained by W.}
}
\value{
A vector of values for one in \{kappa, N, power\} that is not supplied by the user.
}
\description{
Subsidiary function to perform power calculation with covariates under unequal assignment probability with ordered mean assumption.
}
\note{
This function is called internally and thus should not be used directly.
}
\references{
Bansak, K. (2020). A Generalized Approach to Power Analysis for Local Average Treatment Effects. Statistical Science, 35(2), 254-271.
}
\seealso{
\code{\link{equal.unordered.cov}}, \code{\link{equal.ordered.cov}}, \code{\link{unequal.unordered.cov}}.
}
\author{
Kirk Bansak and Eddie Yang
}
