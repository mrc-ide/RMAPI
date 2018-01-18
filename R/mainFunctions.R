#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib RMAPI
#' @importFrom Rcpp evalCpp
NULL

#--------------------------------------------------------------------------------
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

# Node and other parameters as R objects -----------------------------------------

xnode_i=c(1.1,2.1,3.9,2.8,3.6,4.5,0.7,1.8,2.1,3.3)	#Positions of data nodes on x-axis
ynode_i=c(1.2,3.2,4.3,1.1,2.2,3.3,4.5,4.6,3.2,2.6)	#Positions of data nodes on y-axis
vnode_i=c(4.1,2.1,1.4,0.9,1.1,2.2,3.6,4.1,2.1,1.1)	#Values at data nodes (of whatever type - calculation of ellipse 
									#values may have to change depending on type of values used)
									#NOTE: Above three vectors should have same size. Size of xnode used as number of nodes.
a_multiplier_i=-0.45						#Controls relationship between a (ellipse long radius) and c (ellipse short radius equal to distance between foci): a = c*(1 + a_multiplier)

# Set up arguments for input into C++ ---------------------------------------------

args <- list(xnode=xnode_i,ynode=ynode_i,vnode=vnode_i,a_multiplier=a_multiplier_i) 

# Carry out calculations in C++ to generate map data
output_raw <- dummy1_cpp(args)

# process raw output
ret <- output_raw

vmap=matrix(ret[["matrix_values"]],ret[["dim_matrix"]])

#print(ret)
#print(vmap)
filled.contour(ret[["xpoints"]], ret[["ypoints"]], vmap, color = terrain.colors, plot.axes = { points(args[["xnode"]],args[["ynode"]]);axis(1, at = ret[["xtick"]], label = ret[["xtick"]]); axis(2, at = ret[["ytick"]], label = ret[["ytick"]]) })

}







