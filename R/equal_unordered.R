#' @title Subsidiary PowerLATE Function
#' @description Subsidiary function to perform power calculation under equal assignment probability and without ordered mean assumption.
#' @param pi            compliance rate. Equivalently, average causal effect of Z on D.
#' @param N             total number of observations
#' @param kappa         effect size
#' @param sig.level     significance level (Type I error probability).
#' @param power         power of test (1 minus Type II error probability)
#' @return A vector of values for one in \{kappa, N, power\} that is not supplied by the user.
#' @note This function is called internally and thus should not be used directly.
#' @author Kirk Bansak and Eddie Yang
#' @seealso \code{\link{equal.ordered}}, \code{\link{unequal.ordered}}, \code{\link{unequal.unordered}}.
#' @importFrom stats pnorm qnorm
#' @references Bansak, K. (2020). A Generalized Approach to Power Analysis for Local Average Treatment Effects. Statistical Science, 35(2), 254-271.

equal.unordered <- function(
	power = NULL,
	sig.level = NULL,
	pi = NULL,
	kappa = NULL,
	N = NULL){

	# mdes
	if (!is.null(power) && !is.null(N)){
		beta <- 1 - power
		M <- qnorm(1-(sig.level/2)) + qnorm(1-beta)
		upper.effect <- (2*M)/(pi*sqrt(N)-2*M*sqrt((0.5-(pi/2))*(0.5+(pi/2))))
		return(upper.effect)
	}

	# sample size
	if (!is.null(power) && !is.null(kappa)){
		beta <- 1 - power
		M <- qnorm(1-(sig.level/2)) + qnorm(1-beta)
		upper.N <- (4*M^2*(1+kappa*sqrt((0.5-(pi/2))*(0.5+(pi/2))))^2)/(kappa^2*pi^2)
		return(upper.N)
	}

	# power
	if (!is.null(N) && !is.null(kappa)){
		c.val <- qnorm(1-(sig.level/2)) # 1-alpha/2
		effect.bound <- (0.5*kappa*pi*sqrt(N))/(1+kappa*sqrt((0.5-(pi/2))*(0.5+(pi/2))))
		power <- pnorm(-c.val + effect.bound) + pnorm(-c.val - effect.bound)
		return(power)
	}
}