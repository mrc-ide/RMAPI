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

# Node and other parameters as R objects

# convert R objects into a list of arguments
args <- list(xnode=c(1.1,2.1,3.9,2.8,3.6,4.5,0.7,1.8,2.1,3.3),ynode=c(1.2,3.2,4.3,1.1,2.2,3.3,4.5,4.6,3.2,2.6),vnode=c(4.1,2.1,1.4,0.9,1.1,2.2,3.6,4.1,2.1,1.1),a_multiplier=-0.45)

# run efficient c++ code
output_raw <- dummy1_cpp(args)

# process raw output
ret <- output_raw

vmap=matrix(ret[["matrix_values"]],ret[["dim_matrix"]])

#print(ret)
#print(vmap)
filled.contour(ret[["xpoints"]], ret[["ypoints"]], vmap, color = terrain.colors, plot.axes = { points(args[["xnode"]],args[["ynode"]]);axis(1, at = ret[["xtick"]], label = ret[["xtick"]]); axis(2, at = ret[["ytick"]], label = ret[["ytick"]]) })

}
