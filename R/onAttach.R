#' @title Subsidiary PowerLATE Function
#' @description Subsidiary function to perform power calculation under equal assignment probability and ordered mean assumption.
#' @param lib   libname
#' @param pkg   package name
#' @importFrom utils packageDescription
#' @references Bansak, K. (2020). A Generalized Approach to Power Analysis for Local Average Treatment Effects. Statistical Science, 35(2), 254-271.


.onAttach <- function(lib, pkg) {
  mylib <- dirname(system.file(package = pkg))
  title <- packageDescription(pkg, lib.loc = mylib)$Title
  ver <- packageDescription(pkg, lib.loc = mylib)$Version
  #author <- packageDescription(pkg, lib.loc = mylib)$Author
  packageStartupMessage(pkg, ": ", title, "\nVersion: ", ver, "\nReference: ", 
  	"Bansak, K. (2020). A Generalized Approach to Power Analysis for Local Average Treatment Effects. Statistical Science, 35(2), 254-271.", "\n")
}