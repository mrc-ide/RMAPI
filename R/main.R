
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib RMAPI
#' @import assertthat
#' @import graphics
#' @import stats
#' @import utils
#' @importFrom grDevices chull colorRampPalette
#' @importFrom Rcpp evalCpp
NULL

#------------------------------------------------
#' @title Check that RMAPI package has loaded successfully
#'
#' @description Simple function to check that RMAPI package has loaded
#'   successfully. Prints "RMAPI loaded successfully!" if so.
#'
#' @export

check_RMAPI_loaded <- function() {
  message("RMAPI loaded successfully!")
}

#------------------------------------------------
#' @title Get great circle distance between spatial points
#'
#' @description Get great circle distance between spatial points.
#' 
#' @param long vector of longitudes.
#' @param lat vector of latitudes.
#'
#' @export

get_spatial_distance <- function(long, lat) {
  
  # check inputs
  assert_vector(long)
  assert_numeric(long)
  assert_vector(lat)
  assert_numeric(lat)
  assert_same_length(long, lat)
  
  # calculate distance matrix
  ret <- apply(cbind(long, lat), 1, function(y) {lonlat_to_bearing(long, lat, y[1], y[2])$gc_dist})
  ret <- as.dist(ret, upper = TRUE)
  
  return(ret)
}

#------------------------------------------------
#' @title Calculate great circle distance and bearing between coordinates
#'
#' @description Calculate great circle distance and bearing between spatial
#'   coordinates.
#'
#' @param origin_lon the origin longitude
#' @param origin_lat the origin latitude
#' @param dest_lon the destination longitude
#' @param dest_lat the destination latitude
#'
#' @export
#' @examples
#' # one degree longitude should equal approximately 111km at the equator
#' lonlat_to_bearing(0, 0, 1, 0)

lonlat_to_bearing <- function(origin_lon, origin_lat, dest_lon, dest_lat) {
  
  # check inputs
  assert_vector(origin_lon)
  assert_numeric(origin_lon)
  assert_vector(origin_lat)
  assert_numeric(origin_lat)
  assert_vector(dest_lon)
  assert_numeric(dest_lon)
  assert_vector(dest_lat)
  assert_numeric(dest_lat)
  
  # convert input arguments to radians
  origin_lon <- origin_lon*2*pi/360
  origin_lat <- origin_lat*2*pi/360
  dest_lon <- dest_lon*2*pi/360
  dest_lat <- dest_lat*2*pi/360
  
  # get change in lon
  delta_lon <- dest_lon - origin_lon
  
  # calculate bearing
  bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  
  # calculate great circle angle. Use temporary variable to avoid acos(>1) or 
  # acos(<0), which can happen due to underflow issues
  tmp <- sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon)
  tmp <- ifelse(tmp > 1, 1, tmp)
  tmp <- ifelse(tmp < 0, 0, tmp)
  gc_angle <- acos(tmp)
  
  # convert bearing from radians to degrees measured clockwise from due north,
  # and convert gc_angle to great circle distance via radius of earth (km)
  bearing <- bearing*360/(2*pi)
  bearing <- (bearing+360)%%360
  earth_rad <- 6371
  gc_dist <- earth_rad*gc_angle
  
  # return list
  ret <-list(bearing = bearing,
             gc_dist = gc_dist)
  return(ret)
}

#------------------------------------------------
#' @title Load data into RMAPI project
#'
#' @description Load data into RMAPI project.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param long vector of node longitudes.
#' @param lat vector of node latitudes.
#' @param stat_dist matrix of pairwise statistics between nodes.
#' @param check_delete_output if \code{TRUE} (the default) then check before
#'   overwriting any existing data loaded into a project.
#'
#' @export

bind_data <- function(proj, long, lat, stat_dist, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_vector(long)
  assert_numeric(long)
  assert_vector(lat)
  assert_numeric(lat)
  assert_same_length(long, lat)
  assert_matrix(stat_dist)
  assert_numeric(stat_dist)
  assert_nrow(stat_dist, length(long))
  assert_ncol(stat_dist, length(long))
  assert_single_logical(check_delete_output)
  
  # check whether there is data loaded into project already
  if (!is.null(proj$data$coords)) {
    
    # return existing project if user not happy to continue
    if (check_delete_output) {
      if (!user_yes_no("All existing output for this project will be lost. Continue? (Y/N): ")) {
        message("returning original project\n")
        return(proj)
      }
    }
    
    # delete old output
    message("overwriting data\n")
  }
  
  # update project with new data
  proj$data$coords <- data.frame(long = long, lat = lat)
  proj$data$stat_dist <- as.dist(stat_dist, upper = TRUE)
  proj$data$spatial_dist <- get_spatial_distance(long, lat)
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Fit a simple model to data
#'
#' @description Fit a simple linear or exponential model to the distribution of
#'   the observed statistic as a function of distance. Can be used in the
#'   simulation step to generate a null model against which to compare.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param type which model to fit to the data: 1 = linear model, 2 =
#'   exponential model.
#'
#' @export

fit_model <- function(proj, type = 1) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_single_pos_int(type)
  assert_in(type, 1:2)
  
  # fit model
  if (type == 1) {
    model_fit <- lm(proj$data$stat_dist ~ proj$data$spatial_dist)
  } else {
    df <- data.frame(x = as.vector(proj$data$spatial_dist), y = as.vector(proj$data$stat_dist))
    model_fit <- nls(y ~ SSasymp(x, alpha, beta, log_lambda), data = df)
  }
  
  # save model
  proj$model <- list(type = type,
                     model_fit = model_fit)
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Create mapping space
#'
#' @description Create mapping space.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param hex_size size of hexagons.
#' @param buffer size of buffer zone around the data.
#'
#' @import rgeos
#' @import sp
#' @export

create_map <- function(proj, hex_size = 1, buffer = 2*hex_size) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_single_pos(hex_size)
  assert_single_pos(buffer)
  
  message("Creating hex map")
  
  # unfortunately we have to go through a long process to get hexs that cover
  # all the nodes, the reason being that a raw call to sp::spsample() only
  # creates hexs whose centroid is fully within the bounding poly, which can
  # leave some nodes outside. The solution implemented here is to 1) create a
  # bounding poly from the convex hull of the data, 2) apply a large buffer to
  # the bounding poly, 3) generate hexs from the buffered poly, 4) subset to
  # hexs that intersect the original poly, 5) create a new bounding poly from
  # the convex hull of the centroids of the remaining hexs, 5) this new bounding
  # poly is used to create the hex map, with optional buffer applied by the
  # user.
  
  # get convex hull of data
  ch_data <- chull(proj$data$coords[,c("long", "lat")])
  ch_data_coords <- proj$data$coords[c(ch_data, ch_data[1]), c("long", "lat")]
  
  # get convex hull in SpatialPolygons format and expand by fixed buffer of two hexs
  bounding_poly_original_raw <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(ch_data_coords)), ID = 1)))
  bounding_poly_original <- rgeos::gBuffer(bounding_poly_original_raw, width = 2*hex_size)
  
  # get hex centroids and polygons
  hex_pts_original <- sp::spsample(bounding_poly_original, type = "hexagonal", cellsize = hex_size, offset = c(0,0))
  hex_pts_original_df <- as.data.frame(hex_pts_original)
  names(hex_pts_original_df) <- c("long", "lat")
  hex_polys_original <- sp::HexPoints2SpatialPolygons(hex_pts_original)
  
  # convert original bounding poly and hexs to sf format
  bounding_poly_original_raw_sfc <- sf::st_as_sfc(bounding_poly_original_raw)
  hex_polys_original_sfc <- sf::st_as_sfc(hex_polys_original)
  
  # subset hex centroids and polys to those that intersect original bounding poly
  intersect_vec <- as.matrix(sf::st_intersects(hex_polys_original_sfc, bounding_poly_original_raw_sfc))[,1]
  hex_pts_original_df <- hex_pts_original_df[which(intersect_vec),]
  hex_polys_original <- hex_polys_original[which(intersect_vec)]
  
  # get convex hull of hex centroids
  ch_hex <- chull(hex_pts_original_df)
  ch_hex_coords <- hex_pts_original_df[c(ch_hex, ch_hex[1]),]
  
  # get convex hull in SpatialPolygons format and expand by user-defined buffer
  bounding_poly_raw <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(ch_hex_coords)), ID = 1)))
  bounding_poly <- rgeos::gBuffer(bounding_poly_raw, width = buffer)
  
  # get hex centroids and polygons
  hex_pts <- sp::spsample(bounding_poly, type = "hexagonal", cellsize = hex_size, offset = c(0,0))
  hex_pts_df <- as.data.frame(hex_pts)
  names(hex_pts_df) <- c("long", "lat")
  hex_polys <- sp::HexPoints2SpatialPolygons(hex_pts)
  nhex <- length(hex_polys)
  
  message(sprintf("%s hexagons created", nhex))
  
  # add to project
  proj$map$hex <- hex_polys
  proj$map$hex_centroid <- hex_pts_df
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Perform RMAPI simulation
#'
#' @description Perform RMAPI simulation.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param eccentricity eccentricity of ellipses, defined as half the distance
#'   between foci divided by the semi-major axis. We can say \eqn{e = sqrt{1 -
#'   b^2/a^2}}, where \eqn{e} is the eccentricity, \eqn{a} is the length of the
#'   semi-major axis, and \eqn{b} is the length of the semi-minor axis.
#'   Eccentricity ranges between 0 (perfect circle) and 1 (straight line between
#'   foci).
#' @param null_method what method to use for the null model:
#'   \enumerate{
#'     \item permutation test on raw statistic. If \code{n_breaks = 1} then this
#'     is identical to the original MAPI method, if \code{n_breaks > 1} then a
#'     spatial permutation test is implemented in which edges are binned based
#'     on spatial distance and permutation only occurs within a bin.
#'     \item permutation test on residuals. The fitted model from a previously
#'     call of \code{fit_model()} is used to compute residuals, and the
#'     permutation test method is carried out on residual values. As with method
#'     1, this permutation test can be spatially binned.
#'   }
#' @param n_perms number of permutations to run when checking statistical
#'   significance. Set to 0 to skip this step.
#' @param n_breaks alternative to defining sequence of \code{dist_breaks}. If
#'   \code{dist_breaks == NULL} then the full spatial range of the data is split
#'   into \code{n_breaks} equal groups.
#' @param dist_breaks vector of breaks that divide spatial distances into
#'   permutation groups. Defaults to the full spatial range of the data, i.e.
#'   all edges are included in permutation testing.
#' @param empirical_tail whether to do calculate empirical p-values using a
#'   one-sided test (\code{empirical_tail = "left"} or \code{empirical_tail =
#'   "right"}) or a two-sided test (\code{empirical_tail = "both"}).
#' @param min_intersections minimum number of ellipses that must intersect a hex
#'   for it to be included in the final map, otherwise \code{NA}.
#' @param report_progress if \code{TRUE} then a progress bar is printed to the
#'   console during the permutation testing procedure.
#'
#' @export

run_sims <- function(proj, eccentricity = 0.5, null_method = 1,
                     n_perms = 1e2, n_breaks = 1, dist_breaks = NULL,
                     empirical_tail = "both", min_intersections = 5, report_progress = TRUE) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_single_pos(eccentricity, zero_allowed = TRUE)
  assert_bounded(eccentricity, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = FALSE)
  assert_single_pos_int(null_method)
  assert_in(null_method, 1:2)
  assert_single_pos_int(n_perms, zero_allowed = TRUE)
  assert_single_pos_int(min_intersections, zero_allowed = FALSE)
  assert_single_pos_int(n_breaks, zero_allowed = FALSE)
  if (is.null(dist_breaks)) {
    dist_breaks <- seq(min(proj$data$spatial_dist, na.rm = TRUE), max(proj$data$spatial_dist, na.rm = TRUE), l = n_breaks + 1)
  }
  assert_vector(dist_breaks)
  assert_numeric(dist_breaks)
  assert_greq(length(dist_breaks), 2)
  assert_eq(dist_breaks, sort(dist_breaks), message = "dist_breaks must be in increasing order")
  assert_leq(min(dist_breaks), min(proj$data$spatial_dist))
  assert_greq(max(dist_breaks), max(proj$data$spatial_dist))
  assert_in(empirical_tail, c("left", "right", "both"))
  assert_single_logical(report_progress)
  
  # ---------------------------------------------
  # Process data differently depending on null method
  
  # get x-values, y-values predicted values
  x <- proj$data$spatial_dist
  y <- proj$data$stat_dist
  y_pred <- predict(proj$model$model_fit)
  
  # subtract model fit under null_method = 2
  if (null_method == 2) {
    y <- y - y_pred
  }
  
  # break spatial distances into groups
  edge_group <- as.numeric(cut(x, breaks = dist_breaks, include.lowest = TRUE))
  edge_group_list <- mapply(function(x) which(edge_group == x) - 1, 1:n_breaks, SIMPLIFY = FALSE)
  
  # ---------------------------------------------
  # Set up arguments for input into C++
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  pb <- txtProgressBar(0, n_perms, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # create argument list
  args <- list(node_long = proj$data$coords$long,
               node_lat = proj$data$coords$lat,
               edge_value = y,
               edge_value_pred = y_pred,
               edge_group_list = edge_group_list,
               hex_long = proj$map$hex_centroid$long,
               hex_lat = proj$map$hex_centroid$lat,
               eccentricity = eccentricity,
               n_perms = n_perms,
               min_intersections = min_intersections,
               report_progress = report_progress)
  
  # ---------------------------------------------
  # Carry out simulations in C++ to generate map data
  
  output_raw <- run_sims_cpp(args, args_functions, args_progress)
  
  # ---------------------------------------------
  # Process raw output
  
  # get hex values
  hex_values <- output_raw$hex_values
  Nintersections <- output_raw$Nintersections
  hex_values[Nintersections < min_intersections] <- NA
  
  # get hex weights
  hex_weights <- output_raw$hex_weights
  hex_weights[Nintersections < min_intersections] <- NA
  
  # get hex value rank proportions
  hex_ranks <- NULL
  if (n_perms > 0) {
    hex_ranks <- (output_raw$hex_ranks+1)/(n_perms+1)
    hex_ranks[Nintersections < min_intersections] <- NA
  }
  
  # save output as list
  proj$output <- list(hex_values = hex_values,
                      hex_weights = hex_weights,
                      hex_ranks = hex_ranks,
                      n_intersections = output_raw$Nintersections)
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Calculate ellipse polygon coordinates from foci and eccentricity
#'
#' @description Calculate ellipse polygon coordinates from foci and eccentricity.
#'
#' @param f1 x- and y-coordinates of the first focus.
#' @param f2 x- and y-coordinates of the first focus.
#' @param ecc eccentricity of the ellipse, defined as half the distance between
#'   foci divided by the semi-major axis. We can say \eqn{e = sqrt{1 -
#'   b^2/a^2}}, where \eqn{e} is the eccentricity, \eqn{a} is the length of the
#'   semi-major axis, and \eqn{b} is the length of the semi-minor axis.
#'   Eccentricity ranges between 0 (perfect circle) and 1 (straight line between
#'   foci).
#' @param n number of points in polygon.
#'
#' @export

get_ellipse <- function(f1 = c(-3,-2), f2 = c(3,2), ecc = 0.8, n = 100) {
  
  # check inputs
  assert_vector(f1)
  assert_length(f1, 2)
  assert_numeric(f1)
  assert_vector(f2)
  assert_length(f2, 2)
  assert_numeric(f2)
  assert_single_pos(ecc)
  assert_bounded(ecc, inclusive_left = FALSE)
  assert_single_pos_int(n)
  
  # define half-distance between foci (c), semi-major axis (a) and semi-minor
  # axis(b)
  c <- 0.5*sqrt(sum((f2-f1)^2))
  a <- c/ecc
  b <- sqrt(a^2-c^2)
  
  # define slope of ellipse (alpha) and angle of points from centre (theta)
  alpha <- atan2(f2[2]-f1[2], f2[1]-f1[1])
  theta <- seq(0, 2*pi, l = n+1)
  
  # define x and y coordinates
  x <- (f1[1]+f2[1])/2 + a*cos(theta)*cos(alpha) - b*sin(theta)*sin(alpha)
  y <- (f1[2]+f2[2])/2 + a*cos(theta)*sin(alpha) + b*sin(theta)*cos(alpha)
  
  # ensure ellipse closes perfectly
  x[n+1] <- x[1]
  y[n+1] <- y[1]
  
  # return as dataframe
  return(data.frame(x = x, y = y))
}

