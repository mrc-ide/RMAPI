
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
#' @title Load data into RMAPI project
#'
#' @description Load data into RMAPI project.
#'
#' @param proj object of class \code{rmapi_project}.
#' @param long vector of node longitudes.
#' @param lat vector of node latitudes.
#' @param pairwise_dist matrix of pairwise distances between nodes.
#' @param check_delete_output if \code{TRUE} (the default) then check before
#'   overwriting any existing data loaded into a project.
#'
#' @export

bind_data <- function(proj, long, lat, pairwise_dist, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_vector(long)
  assert_numeric(long)
  assert_vector(lat)
  assert_numeric(lat)
  assert_same_length(long, lat)
  assert_matrix(pairwise_dist)
  assert_numeric(pairwise_dist)
  assert_nrow(pairwise_dist, length(long))
  assert_ncol(pairwise_dist, length(long))
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
  proj$data$pairwise_dist <- pairwise_dist
  
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
#' @param n_perms number of permutations to run when checking statistical
#'   significance. Set to 0 to skip this step.
#' @param eccentricity eccentricity of ellipses, defined as half the distance
#'   between foci divided by the semi-major axis. Ranges between 0 (perfect
#'   circle) and 1 (straight line between foci).
#'
#' @export
#' @examples
#' # TODO

run_sims <- function(proj, n_perms = 1e2, eccentricity = 0.5, flag_nullmap = 0, dist_model = 0) {
  
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_single_pos_int(n_perms, zero_allowed = TRUE)
  assert_single_pos(eccentricity, zero_allowed = TRUE)
  assert_bounded(eccentricity, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = FALSE)
  
  # ---------------------------------------------
  # Set up arguments for input into C++
  
  args <- list(long_node = proj$data$coords$long,
               lat_node = proj$data$coords$lat,
               vnode = mat_to_rcpp(proj$data$pairwise_dist),
               long_hex = proj$map$hex_centroid$long,
               lat_hex = proj$map$hex_centroid$lat,
               n_perms = n_perms,
               eccentricity = eccentricity,
	             flag_nullmap = flag_nullmap,
	             dist_model = dist_model) 
  
  # ---------------------------------------------
  # Carry out simulations in C++ to generate map data
  
  output_raw <- run_sims_cpp(args)
  return(output_raw)
  # ---------------------------------------------
  # Process raw output
  
  # get final map
  map_values1 <- output_raw$map_values1
  map_values2 <- output_raw$map_values2
  map_values3 <- output_raw$map_values3
  Nintersections <- output_raw$Nintersections
  map_values1[Nintersections==0] <- NA
  
  # get map weights
  map_weights <- output_raw$map_weights
  map_weights[Nintersections==0] <- NA
  
  # get empirical p-values
  empirical_p <- NULL
  if (Nperms>0) {
    empirical_p <- output_raw$empirical_p / Nperms
    empirical_p[Nintersections==0] <- NA
  }
  
  # save output as list
  proj[["output"]] <- list(map_values1 = map_values1, map_values2 = map_values2, map_values3 = map_values3, Nintersections = Nintersections, map_weights = map_weights, empirical_p = empirical_p)
  
  # return invisibly
  invisible(proj)
}

