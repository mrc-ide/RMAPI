
#------------------------------------------------
#' Simple filled contour plot of RMAPI output
#'
#' TODO - some help text here.
#'
#' @param proj the current RMAPI project
#'
#' @export

RMAPI_plot1 <- function(proj) {
	
	# get output
	vmap <- matrix(proj$output[["matrix_values"]], proj$output[["dim_matrix"]])
	
	xlabels <- round(proj$output[["xtick"]],2)
	ylabels <- round(proj$output[["ytick"]],2)
	
	# filled contour
	filled.contour(proj$output[["xpoints"]], proj$output[["ypoints"]], vmap, color = terrain.colors, plot.axes = FALSE)
	
	# add node points
	#points(proj$output[["xnode"]], proj$output[["ynode"]])
	
	# add aces
	axis(1, at = xlabels, label = xlabels)
	axis(2, at = ylabels, label = ylabels)
	
}







