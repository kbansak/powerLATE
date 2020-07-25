#' @title Subsidiary powerLATE Function
#' @description Subsidiary function to perform power calculation with covariates under unequal assignment probability with ordered mean assumption.
#' @param pZ            probability of being assigned to treatment.
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
#' @seealso \code{\link{equal.unordered.cov}}, \code{\link{equal.ordered.cov}}, \code{\link{unequal.unordered.cov}}.
#' @importFrom stats pnorm qnorm
#' @references Bansak, K. (2020). A Generalized Approach to Power Analysis for Local Average Treatment Effects. Statistical Science, 35(2), 254-271.

unequal.ordered.cov <- function(
	power = NULL,
	sig.level = NULL,
	pi = NULL,
	kappa = NULL,
	N = NULL,
	pZ = NULL,
	r2dw = NULL,
	r2yw = NULL){

	S <- 1-r2yw^2
	T <- 1-r2dw^2
	J <- pZ*(1-pZ)

	# mdes
	if (!is.null(power) && !is.null(N)){
		beta <- 1 - power
		M <- qnorm(1-(sig.level/2)) + qnorm(1-beta)
		numer <- 2*M*sqrt(S)
		denom <- sqrt(4*pi^2*N*J - M^2*T)
		upper.effect <- numer/denom
		return(upper.effect)
	}

	# sample size
	if (!is.null(power) && !is.null(kappa)){
		beta <- 1 - power
		M <- qnorm(1-(sig.level/2)) + qnorm(1-beta)
		numer <- M^2*(4*S + kappa^2*T)
		denom <- 4*kappa^2*pi^2*J
		upper.N <- numer/denom
		return(upper.N)
	}

	# power
	if (!is.null(N) && !is.null(kappa)){
		c.val <- qnorm(1-(sig.level/2)) # 1-beta
		effect.bound <- (kappa*pi*sqrt(J*N))/sqrt(S + 0.25*kappa^2*T)
		power <- pnorm(-c.val + effect.bound) + pnorm(-c.val - effect.bound)
		return(power)
	}
}