#' @title A Generalized Approach to Power Analysis for Local Average Treatment Effects (with covariates)
#' @aliases powerLATE.cov
#' @description Function to perform generalized power analysis for the LATE (i.e. under noncompliance with treatment assignment), 
#' allowing for covariate adjustment. Function allows for user to work with either standardized effect sizes or absolute effects. 
#' The results provided presume a test of the null hypothesis that the LATE equals 0 with a two-sided alternative.
#' @usage powerLATE.cov(pZ = 0.5, pi, N, kappa,
#' 	sig.level = 0.05, power,
#' 	effect.size = TRUE, tau = NULL, omega = NULL,
#' 	assume.ord.means = FALSE, r2dw, r2yw)
#' @param pZ            	probability of being assigned to treatment. Default is 0.5, i.e. equal assignment probability.
#' @param pi            	compliance rate. Equivalently, average causal effect of treatment assignment on treatment uptake.
#' @param N             	total sample size.
#' @param kappa         	effect size.
#' @param sig.level     	significance level (Type I error probability). Default is 0.05.
#' @param power         	power of test (1 minus Type II error probability).
#' @param effect.size   	whether effect size (kappa) rather than absolute effect (tau) is used in power calculations. Default is \code{TRUE}.
#' @param tau           	absolute effect. Must only be supplied if \code{effect.size = FALSE}.
#' @param omega             within-group standard deviation of the outcome. Must be supplied if \code{effect.size = FALSE}. See Details.
#' @param assume.ord.means	whether ordered means assumption is made. Default is \code{FALSE}. See Details.
#' @param r2dw				proportion of variation in D left unexplained by Z that is explained by W.
#' @param r2yw 				proportion of variation in Y left unexplained by Z that is explained by W.
#' @details If \code{effect.size = TRUE} (the default setting), exactly two of the parameters \{\code{kappa, N, power}\} must be supplied, 
#' from which the third (target) parameter will be calculated. If \code{effect.size = FALSE}, \code{omega} must be supplied, and exactly two of 
#' the parameters \{\code{tau, N, power}\} must be supplied. \code{pi} must always be supplied, and the user can change \code{pZ} and \code{sig.level} 
#' from their default values.
#' 
#' Values between 0 and 1 must also be supplied for r2dw and r2yw. See "Power with Covariates) section in Bansak (2020) for more information and guidance.
#' 
#' The user may also supply one of \{\code{kappa, N, power, pi, tau, r2dw, r2yw}\} as a vector of values to perform multiple power calculations at a time, 
#' in which case the target parameter will be calculated for that entire vector.
#' 
#' If \code{effect.size = FALSE}, \code{omega} represents the reference within-assignment-group standard deviation of the outcome. The user may wish to 
#' use an estimate of the standard deviation of the outcome prior to the intervention (i.e. in the absence of the treatment). 
#' See "Discussion on Effect Sizes" section in Bansak (2020) for more information and guidance.
#' 
#' The \code{assume.ord.means} argument allows the user to choose whether or not to make the ordered means assumption, presented and described in Bansak (2020). 
#' Users should only make this assumption (i.e. set \code{assume.ord.means = TRUE}) if they are reasonably confident that it will be met in their context of interest. 
#' See "Narrowing the Bounds" section in Bansak (2020) for more information and guidance.
#' 
#' @return A list that includes the values of the input parameters supplied by the user (\code{input.parameter}) and the corresponding output value(s) 
#' of the target parameter (\code{output.parameter}).
#' 
#' Note also that the results along with additional information will be displayed in the console.
#' @author Kirk Bansak and Eddie Yang
#' @references Bansak, K. (2020). A Generalized Approach to Power Analysis for Local Average Treatment Effects. Statistical Science, 35(2), 254-271.
#' @export
#' @examples
#' 
#' #######################################################
#' ### Example 1: 
#' ### LATE w/ equal assigment w/o ordered mean assumption
#' #######################################################
#' 
#' res <- powerLATE.cov(pZ=0.5, pi=0.35, kappa=seq(0.4, 1.0, 0.1), power=0.8,
#'                      r2dw=0.2, r2yw=0.15)
#' 
#' # also returns an invisible object: output.parameter
#' res$output.parameter
#' 
#' #######################################################
#' ### Example 2: 
#' ### LATE w/o equal assigment w/ ordered mean assumption
#' #######################################################
#' 
#' res <- powerLATE.cov(pZ=0.67, pi=0.35, kappa=seq(0.4, 1.0, 0.1), 
#'                 power=0.8, assume.ord.means=TRUE, r2dw=0.2, r2yw=0.15)
#' 
#' # also returns an invisible object: output.parameter
#' res$output.parameter


powerLATE.cov <- function(
	pZ = 0.5,
	pi,
	N,
	kappa,
	sig.level = 0.05,
	power,
	effect.size = TRUE,
	tau = NULL,
	omega = NULL,
	assume.ord.means = FALSE,
	r2dw,
	r2yw){

	# checks
	if (missing(pi)) stop("pi (compliance rate) needs to be specified")
	if (pZ <= 0 | pZ >= 1) stop("pZ (assignemnt probability) needs to be (0, 1)")
	if (any(pi <= 0 | pi > 1)) stop("pi (compliance rate) needs to be (0, 1]")
	if (sig.level <= 0 | sig.level >= 1) stop("sig.level needs to be between (0, 1)")
	if (!missing(power) && any(power <= 0 | power > 1)) stop("power needs to be between (0, 1]")
	if (!missing(power) && any(power < 0.5)) stop("Supplied power is lower than 0.5. Results will be unreliable")
	if (!missing(kappa) && any(kappa<=0)) stop("kappa needs to be positive")
	if (!is.null(tau) && any(tau<=0)) stop("tau needs to be positive")
	if (!is.null(omega) && any(omega<=0)) stop("omega needs to be positive")
	if (!missing(N) && any(N<=0)) stop("N needs to be positive")
	if (!missing(kappa) && !is.null(tau)) stop("kappa and tau cannot be both supplied")
	if (effect.size && !is.null(tau)) stop("tau is supplied instead of kappa while effect.size = FALSE. User should set effect.size = TRUE in order to work with tau in place of kappa")
	if (effect.size && !is.null(omega)) stop("omega is supplied while effect.size = TRUE.")
	if (!effect.size && is.null(omega)) stop("omega must be supplied when effect.size=FALSE")
	if (effect.size && sum(missing(kappa), missing(N), missing(power))!=1){
		stop("two of args {kappa, N, power} need to be specified")
	}
	if (!effect.size && sum(missing(tau), missing(N), missing(power))!=1){
		stop("two of args {N, power, tau} need to be specified")
	}
	#if (effect.size && !is.null(tau)) warning("tau is supplied while effect.size==TRUE. tau will be ignored")
	#if (effect.size && !is.null(omega)) warning("omega is supplied while effect.size==TRUE. omega will be ignored")
	if (missing(r2dw)) stop("r2dw must be supplied")
	if (missing(r2yw)) stop("r2yw must be supplied")
	if (any(r2dw < 0) | any(r2dw >= 1)) stop("r2dw must be in the range [0, 1)")
	if (any(r2yw < 0) | any(r2yw >= 1)) stop("r2yw must be in the range [0, 1)")

	# set up elements
	if (!effect.size && !is.null(tau)){
		kappa <- tau/omega
	}

	if(missing(kappa)) kappa <- NULL
	if(missing(N)) N <- NULL
	if(missing(power)) power <- NULL

	if (sort(lengths(list(N, kappa, power, pi, tau, r2dw, r2yw)), decreasing = TRUE)[2]>1){
		stop("Only one in {N, kappa, power, pi, tau, r2dw, r2yw} can be a vector")
	}
	if(effect.size){
		input.para <- list(pZ=pZ, pi=pi, N=N, kappa=kappa, sig.level=sig.level, power=power, effect.size=effect.size, assume.ord.means=assume.ord.means,
			r2dw=r2dw, r2yw=r2yw)
	}else{
		input.para <- list(pZ=pZ, pi=pi, N=N, tau=tau, omega=omega, sig.level=sig.level, power=power, effect.size=effect.size, assume.ord.means=assume.ord.means,
			r2dw=r2dw, r2yw=r2yw)
	}

	# main
	if (pZ==0.5 && !assume.ord.means){
		out <- equal.unordered.cov(
			pi = pi,
			power = power,
			N = N,
			kappa = kappa,
			sig.level = sig.level,
			r2dw=r2dw,
			r2yw=r2yw)
	}

	if (pZ==0.5 && assume.ord.means){
		out <- equal.ordered.cov(
			pi = pi,
			power = power,
			N = N,
			kappa = kappa,
			sig.level = sig.level,
			r2dw=r2dw,
			r2yw=r2yw)
	}

	if (pZ!=0.5 && !assume.ord.means){
		out <- unequal.unordered.cov(
			pZ = pZ,
			pi = pi,
			power = power,
			N = N,
			kappa = kappa,
			sig.level = sig.level,
			r2dw=r2dw,
			r2yw=r2yw)
	}

	if (pZ!=0.5 && assume.ord.means){
		out <- unequal.ordered.cov(
			pZ = pZ,
			pi = pi,
			power = power,
			N = N,
			kappa = kappa,
			sig.level = sig.level,
			r2dw=r2dw,
			r2yw=r2yw)
	}

	#names <- c("Compliance Rate", "Effect Size", "N", "Power")

	# output
	if (is.null(kappa)){
		target <- 1
	}else if (is.null(N)){
		target <- 2
	}else if (is.null(power)){
		target <- 3
	}

	input <- list(main = paste0("Power analysis for two-sided test that LATE equals zero", "\n\n"),
		pZ = paste0("pZ = ", as.character(checkVec(pZ)), "\n"),
		pi = paste0("pi = ", as.character(checkVec(pi)), "\n"),
		N = paste0("N = ", as.character(checkVec(N)), "\n"),
		kappa = paste0("kappa = ", as.character(checkVec(kappa)), "\n"),
		tau = paste0("tau = ", as.character(checkVec(tau)), "\n"),
		omega = paste0("omega = ", as.character(checkVec(omega)), "\n"),
		power = paste0("Power = ", as.character(checkVec(power)), "\n"),
		r2dw = paste0("r2dw = ", as.character(checkVec(r2dw)), "\n"),
		r2yw = paste0("r2yw = ", as.character(checkVec(r2yw)), "\n"),
		sig.level = paste0("sig.level  = ", as.character(checkVec(sig.level)), "\n\n"))

	if (!effect.size) input$kappa <- NULL
	input <- input[!grepl(pattern = "= \n", x=unlist(input))]
	cat(unlist(input))

	multiple.input <- names(input[grepl(pattern = "= Multiple", x=unlist(input))])
	message.target <- c("kappa (minimum detectable effect size)", "N (required sample size)", "Power")
	output.name <- c(c("kappa", "N", "power")[target], paste0("User-inputted ", multiple.input))

	if (length(out)==1 && any(out<0, !is.finite(out))){
		stop("The returned value is invalid. Results are not feasible with given parameters. You must increase the values of pi and/or N")
	}
	if (length(out)>1 && any(out<0, !is.finite(out))){
		if (sum(out<0, !is.finite(out))==length(out)){
			stop("The returned values are invalid. Results are not feasible with given parameters. You must increase the values of pi and/or N")
		}else{
			out[c(which(out<0), which(!is.finite(out)))] <- NA
			warning("Some of the returned values are invalid, indicating that results are not feasible for those parameters, coercing into NAs")
		}		
	}

	if (!effect.size && target==1){
		out <- out*omega
		output.name[1]  <- "tau"
		message.target[1] <- "tau (minimum detectable effect)"
	}

	cat(paste0("Given these parameter values, the conservative bound for ", message.target[target], ":\n"))

	if(length(multiple.input)!=0){
		res <- structure(list(
        	target = out,
        	multiple.input = eval(parse(text=multiple.input))),
        	.Names = output.name,
        	row.names = c(NA, length(out)),
        	class = "data.frame")
		output <- t(res[,output.name[1]])
		rownames(output) <- output.name[1]
	} else{
		res <- out
		names(res) <- output.name[1]
		output <- res[output.name[1]]
	}

	#if (!all(apply(res, 2, is.finite))){
	#	warning("Some returned values are not finite. Experiment may not be feasible with given parameters")
	#}
	print(res, right=F)

	if (pZ == 0.5 && !assume.ord.means){
		cat("\nNOTE: The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE.", fill=TRUE)
	}
	if (pZ == 0.5 && assume.ord.means){
		cat("\nNOTE: The Ordered-Means assumption is being employed. User should confirm that the assumption is reasonable in the context of interest.", fill=TRUE)
	}
	if (pZ != 0.5 && !assume.ord.means){
		cat("\nNOTE: The Ordered-Means assumption is not being employed. If the user would like to make this assumption to narrow the bounds, set the argument assume.ord.means to TRUE. The Homoskedasticity assumption is currently being made because pZ does not equal 0.5.", fill=TRUE)
	}
	if (pZ != 0.5 && assume.ord.means){
		cat("\nNOTE: The Ordered-Means assumption is being employed. User should confirm that the assumption is reasonable in the context of interest. The Homoskedasticity assumption is currently being made because pZ does not equal 0.5.", fill=TRUE)
	}

	#output.para.name <- colnames(res)[1]
	#output.para <- structure(
	#	list(output = output),
	#	row.names = output.para.name,
	#	class = "data.frame")
	out <- list(input.parameter=input.para, output.parameter=output)
	return(invisible(out))
}


