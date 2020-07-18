#' @title Subsidiary powerLATE Function
#' @description Subsidiary function to perform power calculation with covariates under equal assignment probability with ordered mean assumption.
#' @param pi            compliance rate. Equivalently, average causal effect of Z on D.
#' @param N             total number of observations
#' @param kappa         effect size
#' @param sig.level     significance level (Type I error probability).
#' @param power         power of test (1 minus Type II error probability)
#' @param r2dw			proportion of variation in D left unexplained by Z that is explained by W.
#' @param r2yw 			proportion of variation in Y left unexplained by Z that is explained by W.
#' @return A vector of values for one in {kappa, N, power} that is not supplied by the user.
#' @note This function is called internally and thus should not be used directly.
#' @author Kirk Bansak
#' @seealso \code{\link{equal.unordered.cov}}, \code{\link{unequal.unordered.cov}}, \code{\link{unequal.ordered.cov}}.
#' @importFrom stats pnorm qnorm
#' @references Bansak, K. (2020). A Generalized Approach to Power Analysis for Local Average Treatment Effects. Statistical Science, 35(2), 254-271.

equal.ordered.cov <- function(power = NULL,
								sig.level = NULL,
								pi = NULL,
								kappa = NULL,
								N = NULL,
								r2dw = NULL,
								r2yw = NULL){

	G <- (0.5-(pi/2))*(0.5+(pi/2))
	S <- 1-r2yw^2
	T <- 1-r2dw^2

	# mdes
	if (!is.null(power) && !is.null(N)){
		beta <- 1 - power
		M <- qnorm(1-(sig.level/2)) + qnorm(1-beta)
		numer <- 2*M*sqrt(S)
		denom <- sqrt(N*pi^2 - 4*M^2*T*G)
		upper.effect <- numer/denom
		return(upper.effect)
	}

	# sample size
	if (!is.null(power) && !is.null(kappa)){
		beta <- 1 - power
		M <- qnorm(1-(sig.level/2)) + qnorm(1-beta)
		numer <- 4*M^2*(kappa^2*T*G + S)
		denom <- kappa^2*pi^2
		upper.N <- numer/denom
		return(upper.N)
	}

	# power
	if (!is.null(N) && !is.null(kappa)){
		c.val <- qnorm(1-(sig.level/2)) # 1-beta
		effect.bound <- (0.5*kappa*pi*sqrt(N))/sqrt(S + kappa^2*T*G)
		power <- pnorm(-c.val + effect.bound) + pnorm(-c.val - effect.bound)
		return(power)
	}
}