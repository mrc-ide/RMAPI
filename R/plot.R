
#------------------------------------------------
#' @title Red to blue colours
#'
#' @description Simple sequence of red-to-blue colours.
#'
#' @param n the number of colours.
#'
#' @importFrom grDevices colorRampPalette
#' @export

col_hotcold <- function(n = 6) {
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_pal <- colorRampPalette(raw_cols)
  return(my_pal(n))
}

#------------------------------------------------
# a series of internally-used colours
#' @noRd
daily_cols <- function() {
  c("firebrick1", "chartreuse3", "dodgerblue", "dodgerblue4", "purple", "darkorange", "firebrick4")
}

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
  df_plot <- data.frame(spatial <- as.vector(proj$data$spatial_dist), stat <- as.vector(proj$data$stat_dist))
  plot1 <- ggplot(df_plot) + theme_bw()
  
  # add points
  plot1 <- plot1 + geom_point(aes(x = spatial, y = stat), data=df_plot, col = col)
  
  # add model fit
  if (overlay_model & !is.null(proj$model)) {
    if(is.null(proj$model$model_fit_pred)) { df_predict <- data.frame(spatial <- df_plot$spatial, stat <- stats::predict(proj$model$model_fit)) }
    else { df_predict <- data.frame(spatial <- df_plot$spatial, stat <- proj$model$model_fit_pred) }
    plot1 <- plot1 + geom_line(aes(x = spatial, y = stat), col = "red", data = df_predict)
  }
  
  # titles etc
  plot1 <- plot1 + xlab("spatial distance") + ylab("statistical distance")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot hex coverage
#'
#' @description Given an RMAPI project with a hap already assigned, plot the
#'   coverage (number of intersecting ellipses) of each hex. Useful as a
#'   diagnostic plot, as hexes with low coverage will conform less well to the
#'   assumptions of the permutation test (see \emph{detials}).
#'
#' @details Good coverage is needed to ensure the validity of the statistical
#'   procedure. For each hex, the observed value is the mean of the (normalised)
#'   values of the edges that intersect it. Under the null model that these
#'   normalised values have no systematic bias, the distribution of this test
#'   statistic is approximately normal as a consequence of the central limit
#'   theorem. Hence, we can use a permutation test to characterise the mean and
#'   standard deviation of this null distribution, then we can quantify how
#'   extreme the observed value is in terms of its z-score. If coverage is too
#'   low then the null distribution will not be normally distributed, and hence
#'   the z-score will not be an accurate description of how extreme the observed
#'   data really is.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param return_type format of return object: 1 = a \code{gridExtra} object, 2
#'   = a list of \code{ggplot2} objects.
#' 
#' @import sf
#' @import ggplot2
#' @import gridExtra
#' @export

plot_coverage <- function(proj, return_type = 1) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  rmapi_proj.check_coords_loaded(proj)
  rmapi_proj.check_map_assigned(proj)
  assert_in(return_type, 1:2)
  
  # get number of intersections per hex
  intersect_vec <- mapply(length, p$map$hex_edges)
  
  # bin intersection counts
  cut_breaks <- c(0,50,100,200,300,400,600,800,1000,2000,Inf)
  intersect_bin <- cut(intersect_vec, breaks = cut_breaks)
  
  # -------- Map --------
  
  # basic plot
  plot1 <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank())
  
  # add hexs
  plot1 <- plot1 + geom_sf(aes(fill = intersect_bin), data = proj$map$hex)
  
  # add points
  coords <- data.frame(long <- proj$coords$long, lat <- proj$coords$lat)
  plot1 <- plot1 + geom_point(aes(x = long, y = lat), data = coords, size = 0.5)
  
  # titles and legends
  plot1 <- plot1 + scale_fill_manual(values = col_hotcold(10), name = "no.intersections")
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  plot1 <- plot1 + guides(fill = guide_legend(reverse = T))
  
  # -------- Histogram --------
  
  # basic plot
  plot2 <- ggplot() + theme_bw()
  
  # add histogram
  data_intersect <- data.frame(x <- intersect_vec)
  plot2 <- plot2 + geom_histogram(aes(x = x), color = "black", fill = grey(0.2),
                                  bins = 50, data = data_intersect)
  
  # titles, legends, limits
  plot2 <- plot2 + xlab("no.intersections")
  plot2 <- plot2 + scale_y_continuous(expand = c(0, 0))
  
  # ----------------
  
  # arrange in grid
  plot_grid <- grid.arrange(plot1, plot2, ncol = 2)
  
  # return
  if (return_type == 1) {
    invisible(plot_grid)
  } else {
    return(list(plot1, plot2))
  }
}

#------------------------------------------------
#' @title Plot hex map of RMAPI output
#'
#' @description Plot hex map of RMAPI output.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param col_scale the colour scale to use.
#' @param plot_sampling_points whether to overlay sampling locations.
#' @param plot_hex_grid whether to plot the hex map. If the project contains
#'   output then hexes will be coloured based on z-scores, otherwise a flat
#'   colour scheme will be used.
#' @param plot_significance whether to outline areas that were identified as
#'   significant outliers.
#' @param barrier_list optional list of polygon coordinates that are added to
#'   plot.
#' 
#' @import ggplot2
#' @importFrom viridisLite magma
#' @export

plot_map <- function(proj,
                     col_scale = viridisLite::magma(100),
                     plot_sampling_points = TRUE,
                     plot_hex_grid = TRUE,
                     plot_significance = TRUE,
                     empirical_tail = "both",
                     alpha_raw = 0.05,
                     barrier_list = list()) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_single_logical(plot_sampling_points)
  assert_single_logical(plot_hex_grid)
  assert_single_logical(plot_significance)
  assert_single_string(empirical_tail)
  assert_in(empirical_tail, c("left", "right", "both"))
  assert_single_bounded(alpha_raw)
  assert_list(barrier_list)
  nb <- length(barrier_list)
  if (nb > 0) {
    for (i in 1:nb) {
      assert_dataframe(barrier_list[[i]])
      assert_in("long", names(barrier_list[[i]]))
      assert_in("lat", names(barrier_list[[i]]))
      assert_eq(barrier_list[[i]][1,], barrier_list[[i]][nrow(barrier_list[[i]]),], 
                message = "barrier polygons must be closed, i.e. the last node coordinate equals the first")
    }
  }
  
  # determine which aspects can/should be plotted
  if (plot_sampling_points) {
    plot_sampling_points <- !is.null(proj$data$coords)
  }
  if (plot_hex_grid) {
    plot_hex_grid <- !is.null(proj$map$hex)
  }
  plot_hex_values <- !is.null(proj$output$hex_values)
  
  # determine hex colours
  if (plot_hex_values) {
    add_legend <- TRUE
    col_vec <- proj$output$hex_values
  } else {
    add_legend <- FALSE
    col_vec <- rep(0, length(proj$map$hex))
    plot_significance <- FALSE
  }
  
  # produce basic plot
  plot1 <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank())
  
  # add hexs
  if (plot_hex_grid) {
    plot1 <- plot1 + geom_sf(aes(fill = col_vec),
                             color = NA,
                             data = proj$map$hex)
  }
  
  # outline significance
  if (plot_significance) {
    
    # get significant hexes
    hex_signif <- get_significant_hexes(proj,
                                        empirical_tail = empirical_tail,
                                        alpha_raw = alpha_raw)
    
    # outline low values
    w <- hex_signif$which_lower
    if (length(w) != 0) {
      merged_poly_lower <- get_merged_poly(proj$map$hex[w], d = proj$map$hex_size/10)
      plot1 <- plot1 + geom_sf(color = "white", fill = NA, data = merged_poly_lower)
    }
    
    w <- hex_signif$which_upper
    if (length(w) != 0) {
      merged_poly_upper <- get_merged_poly(proj$map$hex[w], d = proj$map$hex_size/10)
      plot1 <- plot1 + geom_sf(color = "black", fill = NA, data = merged_poly_upper)
    }
    
  }
  
  # add points
  if (plot_sampling_points) {
    plot1 <- plot1 + geom_point(aes(x = long, y = lat),
                                shape = 21, color = "white", fill = "black", size = 1,
                                data = proj$data$coords)
  }
  
  # titles and legends
  if (add_legend) {
    plot1 <- plot1 + scale_fill_gradientn(colours = col_scale, name = "z-score")
  } else {
    plot1 <- plot1 + guides(fill = FALSE)
  }
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  
  # add barrier polygons
  if (nb > 0) {
    for (i in 1:nb) {
      plot1 <- plot1 + geom_polygon(aes(x = long, y = lat),
                                    col = "white", fill = NA, linetype = "dotted",
                                    data = as.data.frame(barrier_list[[i]]))
    }
  }
  
  # return plot object
  return(plot1)
}


#------------------------------------------------
#' title Plot hex map of RMAPI output with hexes considered statistically significant highlighted
#'
#' description Plot hex map of RMAPI output.
#'
#' param proj object of class \code{rmapi_project}.
#' param variable which element of the project output to use as map colours.
#' param col_scale the colour scale to use.
#' param barrier_list optional list of polygon coordinates that are added to
#'   plot.
#' param tails boundaries of statistical significance
#' 
#' import ggplot2
#' importFrom viridisLite magma
#' export
#' @noRd

plot_map2 <- function(proj, variable = NULL, col_scale = magma(100), barrier_list = list(), tails=c(0,1)) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  type=typeof(proj$map$hex)
  assertthat::assert_that(type=="list" || type=="S4")
  if(type=="list") {map_type=1 } else { map_type=2 }
  assert_list(barrier_list)
  assertthat::assert_that(length(tails)==2)
  nb <- length(barrier_list)
  if (nb > 0) {
    for (i in 1:nb) {
      assert_dataframe(barrier_list[[i]])
      assert_in("long", names(barrier_list[[i]]))
      assert_in("lat", names(barrier_list[[i]]))
      assert_eq(barrier_list[[i]][1,], barrier_list[[i]][nrow(barrier_list[[i]]),])
    }
  }
  
  # produce default colours
  add_legend <- TRUE
  if (is.null(variable)) {
    if ("hex_values" %in% names(proj$output)) {
      variable <- "hex_values"
      col_vec <- proj$output$hex_values
    } else {
      add_legend <- FALSE
      col_vec <- rep(0, length(proj$map$hex))
    }
  } else {
    col_vec <- proj$output[[variable]]
  }
  
  # Create plot
  plot1 <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank())
  
  # Hexes without borders with fill colour corresponding to value
  nhexes=length(col_vec)
  hex_df <- list()
  for(i in 1:nhexes){
    if(map_type==1){hex_coords=proj$map$hex[[i]][1][[1]]} else { hex_coords=proj$map$hex[i][1]@polygons[[1]]@Polygons[[1]]@coords }
    hex_df[[i]]=data.frame(long=hex_coords[,1],lat=hex_coords[,2],value=rep(col_vec[i],7))
    plot1 <- plot1 + geom_polygon(aes(x = long, y = lat, fill=value), data = hex_df[[i]], colour = NA, size = 0)
  }
  
  # Unfilled hexes with white borders showing which hexes are statistically significant
  for(i in 1:nhexes){
    if(is.na(proj$output$hex_ranks[i])==FALSE && proj$output$hex_ranks[i]>=0.1 && proj$output$hex_ranks[i]<=1){ 
      plot1 <- plot1 + geom_polygon(aes(x = long, y = lat), data = hex_df[[i]], colour = "white", fill = NA, size = 0.5)
    }
  }
  
  # Add points
  coords <- data.frame(long <- proj$data$coords$long, lat <- proj$data$coords$lat)
  plot1 <- plot1 + geom_point(aes(x = long, y = lat), shape = 21, color = "white", fill = "black", size = 1, data = coords)
  plot1
  
  # Titles and legends
  if (add_legend) {
    plot1 <- plot1 + scale_fill_gradientn(colours = col_scale, name = "hex_values")
  } else {
    plot1 <- plot1 + guides(fill = FALSE)
  }
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  
  # add barrier polygons
  if (nb > 0) {
    for (i in 1:nb) {
      #plot1 <- plot1 + geom_polygon(aes(x = long, y = lat), col = "black", fill = "#00000050", data = as.data.frame(barrier_list[[i]]))
      plot1 <- plot1 + geom_polygon(aes(x = long, y = lat), col = "white", fill = NA, data = as.data.frame(barrier_list[[i]]))
    }
  }
  
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
#' @param map_type an index from 1 to 137 indicating the type of base map. The
#'   map types are taken from \code{leaflet::providers}, see
#'   \href{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}{here}
#'   for an interactive gallary.
#' @param col_scale the colour scale to use.
#' 
#' @import leaflet
#' @importFrom viridisLite magma
#' @importFrom grDevices colorRamp
#' @importFrom methods slot
#' @export

plot_leaflet <- function(proj, variable = NULL, fill_opacity = 0.8, legend_opacity = 1,
                         map_type = 97, col_scale = viridisLite::magma(100)) {
  
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
  
  # add colour variable to hex map
  hex_sf <- st_sf(geom = proj$map$hex, col = x)
  
  # define colour ramp and palette
  col_ramp <- colorRamp(col_scale)
  if (all(x == 0)) {
    col_ramp <- colorRamp(grey(0.5))
  }
  pal <- colorNumeric(col_ramp, domain = x)
  
  # produce basic leaflet plot
  plot1 <- leaflet(hex_sf)
  plot1 <- addProviderTiles(plot1, leaflet::providers[[map_type]])
  
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

#------------------------------------------------
#' @title Plot daily counts of each host state from simulation
#'
#' @description For a set of simulation output produced by the function
#'   \code{sim_falciparum()}, plots daily simulated values in a given deme.
#'
#' @param x object of class \code{rmapi_sim}.
#' @param deme which deme to plot.
#' @param states which states to plot. Can be any subset of \code{c("S", "E",
#'   "I", "Sv", "Ev", "Iv")}.
#'
#' @importFrom grDevices grey
#' @import tidyr
#' @export

plot_daily_states <- function(x, deme = 1, states = c("S", "E", "I")) {
  
  # check inputs
  assert_custom_class(x, "rmapi_sim")
  assert_leq(deme, length(x$daily_values), message = "deme not found within simulation output")
  assert_in(states, c("S", "E", "I", "Sv", "Ev", "Iv"))
  
  # subset to desired rows and columns
  df_wide <- x$daily_values[[deme]][, c("time", states), drop = FALSE]
  
  # get to long format
  df_long <- tidyr::gather(df_wide, states, factor_key = TRUE)
  
  # choose plotting colours
  raw_cols <- daily_cols()
  plot_cols <- c(S = grey(0.5), E = raw_cols[2], I = raw_cols[1],
                 Sv = grey(0.8), Ev = raw_cols[6], Iv = raw_cols[7],
                 EIR = grey(0.0))
  
  # produce plot
  ggplot2::ggplot(df_long) + ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes_(x = ~time, y = ~count, color = ~state)) +
    ggplot2::scale_color_manual(values = plot_cols) +
    ggplot2::xlab("time (days)") + ggplot2::ggtitle(sprintf("Deme %s", deme))
}