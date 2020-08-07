#' @title Print Function for powerLATE
#' @description Print output for powerLATE and powerLATE.cov.
#' @param message.input   vector of strings specifying input parameters.
#' @param message.output  string prompt specifying the target parameter.
#' @param res             dataframe of output parameters.
#' @param note            note specifying assumptions made and additional details.
#' @return strings and a dataframe for output.
#' @note This function is called internally and should not be used directly.
#' @author Kirk Bansak and Eddie Yang

print.powerLATE <- function(message.input, message.output, res, note){
	cat(message.input, message.output)
	print(res, right=F)
	cat("\nNOTE: ", note, "\n", fill=TRUE)
}