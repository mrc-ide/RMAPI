
#------------------------------------------------
#' @title Plot pairwise spatial distance against loaded statistic
#'
#' @description Plot pairwise spatial distance against loaded statistic.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param col the colour of points.
#' 
#' @export

plot_dist <- function(proj, col = "#00000050") {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_length(col, 1)
  
  # create basic plot
  df_plot <- data.frame(spatial = as.vector(proj$data$spatial_dist), stat = as.vector(proj$data$stat_dist))
  plot1 <- ggplot(df_plot) + theme_bw()
  
  # add points
  plot1 <- plot1 + geom_point(aes(x = spatial, y = stat), col = col)
  
  # titles etc
  plot1 <- plot1 + xlab("spatial distance") + ylab("statistical distance")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Simple filled contour plot of RMAPI output
#'
#' @description Simple filled contour plot of RMAPI output.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param variable which element of the project output to use as map colours.
#' 
#' @export

plot_map <- function(proj, variable = NULL) {
	
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  if (is.null(variable)) {
    if ("hex_values" %in% names(proj$output)) {
      col_vec <- proj$output$hex_values
    } else {
      col_vec <- rep(0, length(proj$map$hex))
    }
  } else {
    assert_in(variable, names(proj$output))
    col_vec <- proj$output[[variable]]
  }
  
  # get hex data into dataframe
  hex_df <- ggplot2::fortify(proj$map$hex)
  hex_df$col <- col_vec[as.numeric(hex_df$group)]
  
  # produce plot
  plot1 <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank())
  
  # add hexs
  plot1 <- plot1 + geom_polygon(aes(x = long, y = lat, group = group, fill = col), data = hex_df)
  
  # add points
  plot1 <- plot1 + geom_point(aes(x = long, y = lat), size = 0.5, data = proj$data$coords)
  
  # titles and legends
  plot1 <- plot1 + scale_fill_gradientn(colours = c("#4575B4", "#91BFDB", "#E0F3F8", "#FEE090", "#FC8D59", "#D73027"), name = "hex_value")
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  
  
	# return plot object
  return(plot1)
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

