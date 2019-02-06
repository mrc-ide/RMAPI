
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
#' @examples
#' # TODO

create_map <- function(proj, hex_size = 1, buffer = 2*hex_size) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_single_pos(hex_size)
  assert_single_pos(buffer)
  
  message("Creating hex map")
  
  # get convex hull of data
  ch <- chull(proj$data$coords[,c("long", "lat")])
  ch_coords <- proj$data$coords[c(ch, ch[1]), c("long", "lat")]
  
  # get convex hull in SpatialPolygons format and expand by buffer
  sp_poly_raw <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(ch_coords)), ID = 1)))
  sp_poly <- rgeos::gBuffer(sp_poly_raw, width = buffer)
  
  # get hex centre points and polygons
  hex_pts <- sp::spsample(sp_poly, type = "hexagonal", cellsize = hex_size, offset = c(0,0))
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
#'   between foci divided by the semi-major axis. Ranges between 0 (perfect
#'   circle) and 1 (straight line between foci).
#' @param n_perms number of permutations to run when checking statistical
#'   significance. Set to 0 to skip this step.
#' @param min_intersections minimum number of ellipses that must intersect a hex
#'   for it to be included in the final map.
#' @param dist_breaks vector of breaks that divide spatial distances into
#'   permutation groups. Defaults to the full spatial range of the data, i.e.
#'   all edges are included in permutation testing.
#' @param n_breaks alternative to defining sequence of \code{dist_breaks}. If
#'   \code{dist_breaks == NULL} then the full spatial range of the data is split
#'   into \code{n_breaks} equal groups.
#'
#' @export
#' @examples
#' # TODO

run_sims <- function(proj, eccentricity = 0.5, n_perms = 1e2,
                     min_intersections = 5, dist_breaks = NULL, n_breaks = 1,
                     divide_weights = TRUE) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_single_pos(eccentricity, zero_allowed = TRUE)
  assert_bounded(eccentricity, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = FALSE)
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
  
  # ---------------------------------------------
  # Set up arguments for input into C++
  
  # break spatial distances into groups
  edge_group <- as.numeric(cut(proj$data$spatial_dist, breaks = dist_breaks, include.lowest = TRUE))
  edge_group_list <- mapply(function(x) which(edge_group == x) - 1, 1:n_breaks, SIMPLIFY = FALSE)
  
  # create argument list
  args <- list(node_long = proj$data$coords$long,
               node_lat = proj$data$coords$lat,
               edge_value = proj$data$stat_dist,
               #edge_group = edge_group,
               edge_group_list = edge_group_list,
               hex_long = proj$map$hex_centroid$long,
               hex_lat = proj$map$hex_centroid$lat,
               n_perms = n_perms,
               min_intersections = min_intersections,
               eccentricity = eccentricity,
               divide_weights = divide_weights)
  
  # ---------------------------------------------
  # Carry out simulations in C++ to generate map data
  
  output_raw <- run_sims_cpp(args)
  
  # ---------------------------------------------
  # Process raw output
  
  # get hex values
  hex_values <- output_raw$hex_values
  Nintersections <- output_raw$Nintersections
  hex_values[Nintersections < min_intersections] <- NA
  
  # get hex weights
  hex_weights <- output_raw$hex_weights
  hex_weights[Nintersections < min_intersections] <- NA
  
  # get empirical p-values
  empirical_p <- NULL
  if (n_perms > 0) {
    empirical_p <- output_raw$empirical_p #/ n_perms
    empirical_p[Nintersections < min_intersections] <- NA
  }
  
  # save output as list
  proj$output <- list(hex_values = hex_values,
                      hex_weights = hex_weights,
                      empirical_p = empirical_p,
                      n_intersections = output_raw$Nintersections)
  
  # return invisibly
  invisible(proj)
}

