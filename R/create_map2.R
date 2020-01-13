#------------------------------------------------
#' title Create mapping space
#'
#' description Create mapping space - alternate version where hex map is built around input region outline
#'
#' param proj object of class \code{rmapi_project}.
#' param hex_size size of hexagons.
#' param buffer size of buffer zone around the data. It is recommended to not
#'   use a buffer to avoid edge-effects.
#'
#' import rgeos
#' import sp
#' export
#' @noRd

create_map2 <- function(proj, hex_size = 1, buffer = 0,outline) {
   
  # check inputs
  assert_custom_class(proj, "rmapi_project")
  assert_single_pos(hex_size)
  assert_single_pos(buffer)
  assert_list(outline)
  
  # check hex size
  min_range <- min(apply(proj$data$coords, 2, function(x) diff(range(x))))
  if (hex_size > min_range/4) {
    stop(sprintf("hex_size too large for spatial range of data. Suggested size: %s", round(min_range/10, digits = 3)))
  }
  
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
  ch_data <- chull(x=outline$Long,y=outline$Lat)
  ch_data_coords <- outline[c(ch_data, ch_data[1]), c("Long", "Lat")]

  # get convex hull in SpatialPolygons format and expand by fixed buffer of two hexs
  bounding_poly_original_raw <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(outline)), ID = 1)))
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