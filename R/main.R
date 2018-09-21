
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib RMAPI
#' @import assertthat
#' @import sp
#' @import rgeos
#' @import graphics
#' @import stats
#' @import utils
#' @importFrom grDevices chull colorRampPalette
#' @importFrom Rcpp evalCpp
NULL

#------------------------------------------------
#' Load data into RMAPI project
#'
#' TODO - some help text here.
#'
#' @param proj the current RMAPI project
#' @param data a data frame, formatted into the correct RMAPI format (see Details)
#' @param check_delete_output whether to perform a check to see if project already contains data, in which case all old data and output will be lost
#'
#' @export

bind_data <- function(proj, data, check_delete_output = TRUE) {
  
  # check inputs
  assert_that( is.rmapi_project(proj) )
  assert_that( is.data.frame(data) )
  assert_that( ncol(data)>=4 )
  assert_that( nrow(data)==ncol(data)-3 )
  assert_that( is.logical(check_delete_output) )
  
  # check whether there is data loaded into project already
  if (!is.null(proj$data)) {
    
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
  proj[["data"]] <- data
  proj[["map"]] <- list()
  proj[["output"]] <- list()
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' Create mapping space
#'
#' TODO - some help text here.
#' TODO - move to hex grid
#'
#' @param proj the current RMAPI project
#' @param hex_size size of hexagons
#' @param buffer size of buffer zone around the data
#'
#' @export
#' @examples
#' # TODO

create_map <- function(proj, hex_size = 1, buffer = 2*hex_size) {
  
  # check inputs
  assert_that( is.rmapi_project(proj) )
  assert_that( is.pos_scalar(hex_size) )
  assert_that( is.pos_scalar(buffer) )
  
  # TODO - check project has correct elements
  
  message("Creating hex map")
  
  # get convex hull of data
  ch <- chull(proj$data[,2:3])
  ch_coords <- proj$data[c(ch, ch[1]), 2:3]
  
  # get convex hull in SpatialPolygons format and expand by buffer
  sp_poly_raw <- SpatialPolygons(list(Polygons(list(Polygon(ch_coords)), ID=1)))
  sp_poly <- gBuffer(sp_poly_raw, width = buffer)
  
  # get hex centre points and polygons
  hex_pts <- spsample(sp_poly, type="hexagonal", cellsize = hex_size, offset=c(0,0))
  hex_pts_df <- as.data.frame(hex_pts)
  names(hex_pts_df) <- c("long", "lat")
  hex_polys <- HexPoints2SpatialPolygons(hex_pts)
  nhex <- length(hex_polys)
  
  message(paste0(nhex, " hexagons"))
  
  # add to project
  proj[["map"]]$hex <- hex_polys
  proj[["map"]]$hex_centroid <- hex_pts_df
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' Perform RMAPI simulation
#'
#' TODO - some help text here.
#'
#' @param proj the current RMAPI project
#' @param Nperms number of permutations to run when checking statistical significance. Set to 0 to skip this step
#' @param eccentricity eccentricity of ellipses, defined as half the distance between foci divided by the semi-major axis. Ranges between 0 (perfect circle) and 1 (straight line between foci).
#'
#' @export
#' @examples
#' # TODO

run_sims <- function(proj, Nperms = 1e2, eccentricity = 0.5) {
  
  # check inputs
  assert_that( is.rmapi_project(proj) )
  assert_that( is.numeric(Nperms) )
  assert_that( length(Nperms)==1 )
  assert_that( is.int(Nperms) )
  assert_that( is.pos_scalar(eccentricity) )
  assert_that( eccentricity<1 )
  
  # TODO - check project has correct elements
  
  # Split inputs into components
  long_node <- proj$data$long
  lat_node <- proj$data$lat
  pairwise_stats <- as.matrix(proj$data[,4:ncol(proj$data)])
  long_hex <- proj$map$hex_centroid$long
  lat_hex <- proj$map$hex_centroid$lat
  
  # ---------------------------------------------
  # Set up arguments for input into C++
  
  args <- list(long_node = long_node,
               lat_node = lat_node,
               vnode = mat_to_rcpp(pairwise_stats),
               long_hex = long_hex,
               lat_hex = lat_hex,
               Nperms = Nperms,
               eccentricity = eccentricity) 
  
  # ---------------------------------------------
  # Carry out simulations in C++ to generate map data
  
  output_raw <- run_sims_cpp(args)
  
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

