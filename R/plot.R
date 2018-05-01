
#------------------------------------------------
#' Simple filled contour plot of RMAPI output
#'
#' TODO - some help text here.
#'
#' @param proj the current RMAPI project
#' @param variable which element of the project output to use as map colours
#'
#' @export

plot_map <- function(proj, variable = "map_values") {
	
  # check inputs
  assert_that( is.rmapi_project(proj) )
  assert_that( variable %in% names(proj$output))
  
	# get output
	plot_vals <- proj$output[[variable]]
	
	plot(proj$map$hex, col = smooth_cols(plot_vals), border = NA, main = variable)
	
}

#------------------------------------------------
#' Add lines to a plot representing barriers
#'
#' TODO - some help text here.
#'
#' @param barrier_x TODO
#' @param barrier_y TODO
#' @param barrier_angle TODO
#' @param ... other arguments that get passed directly to \code{abline()}
#'
#' @export

barrier_lines <- function(barrier_x = 0, barrier_y = 0, barrier_angle = 0, ...) {
  
  # recycle inputs if different lengths
  n <- max(mapply(length, list(barrier_x, barrier_y, barrier_angle)))
  barrier_x <- rep(barrier_x, n)[1:n]
  barrier_y <- rep(barrier_y, n)[1:n]
  barrier_angle <- rep(barrier_angle, n)[1:n]
  
  # add barrier lines
  for (i in 1:n) {
    if (barrier_angle[i]==0 | barrier_angle[i]==180) {
      abline(v=barrier_x[i], ...)
    } else if (barrier_angle[i]==90 | barrier_angle[i]==270) {
      abline(h=barrier_y[i], ...)
    } else {
      theta <- barrier_angle[i]/360*2*pi
      abline(a=barrier_y[i] + barrier_x[i]/tan(theta), b=1/tan(theta), ...)
    }
  }
}

#------------------------------------------------
# produce smooth colours directly from numeric values
# (not exported)

smooth_cols <- function (x, xmin = min(x, na.rm = T), xmax = max(x, na.rm = TRUE), 
          n = 1000, raw_cols = NULL) {
  
  # scale to between 0 and 1
  x <- (x - xmin)/(xmax - xmin)
  
  # define default colours
  if (is.null(raw_cols)) {
    rawCols <- c("#00008F", "#00009F", "#0000AF", "#0000BF", 
                 "#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0010FF", 
                 "#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF", 
                 "#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF", 
                 "#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", 
                 "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", 
                 "#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", 
                 "#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", 
                 "#FFFF00", "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", 
                 "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", 
                 "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", 
                 "#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", 
                 "#BF0000", "#AF0000", "#9F0000", "#8F0000", "#800000")
  }
  
  # get colours from colour palette
  myPal <- colorRampPalette(rawCols)
  ret <- myPal(n + 1)[floor(x*n) + 1]
  
  return(ret)
}
