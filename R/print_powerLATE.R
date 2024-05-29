#' @title Print Function for powerLATE
#' @method print powerLATE
#' @description Print output for powerLATE and powerLATE.cov.
#' @param x List of message.input, message.output, res, note to be printed
#' @param ... Further arguments to be passed to \code{print.powerLATE()}.
#' @return strings and a dataframe for output.
#' @note This function is called internally and should not be used directly.
#' @author Kirk Bansak and Eddie Yang

print.powerLATE <- function(x, ...){
	cat(x$message.input, x$message.output)
	print(x$res, right=F)
	cat("\nNOTE: ", x$note, "\n", fill=TRUE)
}