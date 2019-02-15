
#------------------------------------------------
#' @title Plot pairwise spatial distance against loaded statistic
#'
#' @description Plot pairwise spatial distance against loaded statistic.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param col the colour of points.
#' @param overlay_model whether to overlay the model fit, if \code{fit_model()}
#'   has been called.
#' 
#' @import ggplot2
#' @export

plot_dist <- function(proj, col = "#00000050", overlay_model = TRUE) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_length(col, 1)
  assert_single_logical(overlay_model)
  
  # create basic plot
  df_plot <- data.frame(spatial = as.vector(proj$data$spatial_dist), stat = as.vector(proj$data$stat_dist))
  plot1 <- ggplot(df_plot) + theme_bw()
  
  # add points
  plot1 <- plot1 + geom_point(aes(x = spatial, y = stat), col = col)
  
  # add model fit
  if (overlay_model) {
    df_predict <- data.frame(spatial = df_plot$spatial, stat = predict(proj$model$model_fit))
    plot1 <- plot1 + geom_line(aes(x = spatial, y = stat), col = "red", data = df_predict)
  }
  
  # titles etc
  plot1 <- plot1 + xlab("spatial distance") + ylab("statistical distance")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot hex map of RMAPI output
#'
#' @description Plot hex map of RMAPI output.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param variable which element of the project output to use as map colours.
#' 
#' @import ggplot2
#' @export

plot_map <- function(proj, variable = NULL) {
	
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  if (!is.null(variable)) {
    assert_in(variable, names(proj$output))
  }
  
  # produce default colours
  add_legend <- TRUE
  if (is.null(variable)) {
    if ("hex_values" %in% names(proj$output)) {
      col_vec <- proj$output$hex_values
    } else {
      add_legend <- FALSE
      col_vec <- rep(0, length(proj$map$hex))
    }
  } else {
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
  if (add_legend) {
    plot1 <- plot1 + scale_fill_gradientn(colours = c("#4575B4", "#91BFDB", "#E0F3F8", "#FEE090", "#FC8D59", "#D73027"), name = variable)
  } else {
    plot1 <- plot1 + guides(fill = FALSE)
  }
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  
	# return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Interactive hex map of RMAPI output
#'
#' @description Interactive hex map of RMAPI output.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param variable which element of the project output to use as map colours.
#' @param fill_opacity opacity of fill; 0 = fully transparent, 1 = fully opaque.
#' @param legend_opacity opacity of legend; 0 = fully transparent, 1 = fully
#'   opaque.
#' 
#' @import leaflet
#' @export

plot_leaflet <- function(proj, variable = NULL, fill_opacity = 1, legend_opacity = 1) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  if (!is.null(variable)) {
    assert_in(variable, names(proj$output))
  }
  
  # produce default colours
  add_legend <- TRUE
  if (is.null(variable)) {
    if ("hex_values" %in% names(proj$output)) {
      x <- proj$output$hex_values
    } else {
      add_legend <- FALSE
      x <- rep(0, length(proj$map$hex))
    }
  } else {
    x <- proj$output[[variable]]
  }
  
  # combine chosen variable with hex map to produce SpatialPolygonsDataFrame object 
  poly_df <- data.frame(col = x)
  rownames(poly_df) <- sapply(slot(proj$map$hex, "polygons"), function(x) slot(x, "ID"))
  hex_spdf <- SpatialPolygonsDataFrame(p$map$hex, poly_df)
  
  # define colour ramp and palette
  col_ramp <- colorRamp(c("#4575B4", "#91BFDB", "#E0F3F8", "#FEE090", "#FC8D59", "#D73027"))
  if (all(x == 0)) {
    col_ramp <- colorRamp(grey(0.5))
  }
  pal <- colorNumeric(col_ramp, domain = x)
  
  # produce basic leaflet plot
  plot1 <- leaflet(hex_spdf)
  plot1 <-  addProviderTiles(plot1, leaflet::providers[[97]])
  
  # add hex polygons
  plot1 <- addPolygons(plot1, color = NA, fillColor = ~pal(col), fillOpacity = fill_opacity)
  
  # add legend
  if (add_legend) {
    plot1 <- addLegend(plot1, position = "bottomright", pal = pal, values = ~col,
                       title = variable, opacity = legend_opacity)
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Add points to dynamic map
#'
#' @description Add points to dynamic map
#'
#' @param myplot dynamic map produced by \code{plot_leaflet()} function
#' @param lon longitude of points
#' @param lat latitude of points
#' @param col colour of points
#' @param size size of points
#' @param opacity opacity of points
#'
#' @import leaflet
#' @export

overlay_points <- function(myplot, lon, lat, col = "black", size = 2, opacity = 1.0) {
  
  # check inputs
  assert_custom_class(myplot, "leaflet")
  assert_numeric(lon)
  assert_vector(lon)
  assert_numeric(lat)
  assert_vector(lat)
  assert_same_length(lon, lat)
  assert_single_string(col)
  assert_single_pos(size, zero_allowed = FALSE)
  assert_single_pos(opacity, zero_allowed = TRUE)
  assert_bounded(opacity)
  
  # add circle markers
  myplot <- addCircleMarkers(myplot, lng = lon, lat = lat, radius = size,
                             fillColor = col, stroke = FALSE, fillOpacity = opacity)
  
  # return plot object
  return(myplot)
}
