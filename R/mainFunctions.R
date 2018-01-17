#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib RMAPI
#' @importFrom Rcpp evalCpp
NULL

#------------------------------------------------
#' Dummy function
#'
#' This is a dummy function
#'
#' @param x Some parameter
#'
#' @export
#' @examples
#' dummy1()

dummy1 <- function() {
	
	# convert R objects into a list of arguments
	args <- list(foo=99, bar=c(1,3,6))

	# run efficient c++ code
	output_raw <- dummy1_cpp(args)
	
	# process raw output
	ret <- output_raw$bar
	
	print(output_raw)
}