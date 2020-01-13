
#------------------------------------------------
#' @title Simulate data from simple model
#'
#' @description Simulate data from simple model.
#'
#' @param node_long longitudes of nodes.
#' @param node_lat latitudes of nodes.
#' @param barrier_list list of polygons representing barriers. Each element of
#'   the list must be a dataframe with columns \code{long} and \code{lat}
#'   specifying the coordinates of points that make up the polygon. Polygons
#'   must be complete rings, meaning the final row of the dataframe must equal
#'   the first row.
#' @param barrier_penalty penalty values of each barrier.
#' @param barrier_method the method by which penalties are applied:
#'   \itemize{
#'     \item{bullet 1 compare line, apply penalty on intersection}
#'     \item{bullet 2 compare line, apply penalty per unit intersection}
#'     \item{bullet 3 compare ellipse, apply penalty per unit area intersection}
#'   }
#' @param eccentricity eccentricity of ellipses (only used under
#'   \code{barrier_method = 3}).
#' @param n_ell number of points that make up an ellipse (only used under
#'   \code{barrier_method = 3}).
#' @param dist_transform the method by which distances are transformed to
#'   produce the final statistic:
#'   \itemize{
#'     \item{bullet 1 linear}
#'     \item{bullet 2 exponential decay}
#'   }
#' @param lambda the rate of decay of the exponential function (only used under
#'   \code{dist_transform = 2}).
#' @param eps the standard deviation of white noise applied to final statistics.
#'
#' @import sf
#' @importFrom stats dist
#' @export

sim_simple <- function(node_long,
                       node_lat,
                       barrier_list = list(),
                       barrier_penalty = numeric(),
                       barrier_method = 1,
                       eccentricity = 0.9,
                       n_ell = 20,
                       dist_transform = 1,
                       lambda = 0.1,
                       eps = 0.1) {
  
  # check inputs
  assert_vector(node_long)
  assert_numeric(node_long)
  assert_vector(node_lat)
  assert_numeric(node_lat)
  assert_same_length(node_long, node_lat)
  assert_list(barrier_list)
  nb <- length(barrier_list)
  if (nb > 0) {
    for (i in 1:nb) {
      assert_dataframe(barrier_list[[i]])
      assert_in("long", names(barrier_list[[i]]))
      assert_in("lat", names(barrier_list[[i]]))
      assert_eq(barrier_list[[i]][1,], barrier_list[[i]][nrow(barrier_list[[i]]),])
    }
  }
  assert_vector(barrier_penalty)
  assert_numeric(barrier_penalty)
  assert_same_length(barrier_list, barrier_penalty)
  assert_single_pos_int(barrier_method)
  assert_in(barrier_method, 1:3)
  assert_single_numeric(eccentricity)
  assert_bounded(eccentricity, inclusive_left = FALSE)
  assert_single_pos_int(n_ell, zero_allowed = FALSE)
  assert_single_pos_int(dist_transform)
  assert_in(dist_transform, 1:2)
  assert_single_numeric(lambda)
  assert_single_pos(eps, zero_allowed = TRUE)
  
  # apply barrier penalties
  intersect_penalty <- 0
  if (nb > 0) {
    
    # convert barrier list to st_polygon
    poly_list <- list()
    for (i in 1:length(barrier_list)) {
      poly_list[[i]] <- sf::st_polygon(list(as.matrix(barrier_list[[i]])))
    }
    
    # get node coordinates in matrix
    node_mat <- cbind(node_long, node_lat)
    
    # if comparing lines
    if (barrier_method %in% c(1,2)) {
      
      # create all pairwise sf_linestring between nodes
      line_list <- list()
      n_node <- length(node_long)
      i2 <- 1
      for (i in 1:(n_node-1)) {
        for (j in (i+1):n_node) {
          line_list[[i2]] <- sf::st_linestring(node_mat[c(i,j),])
          i2 <- i2+1
        }
      }
      
      # convert lines and polys to st_sfc
      line_sfc <- sf::st_sfc(line_list)
      poly_sfc <- sf::st_sfc(poly_list)
      
      # get boolean intersection matrix
      intersect_mat <- as.matrix(sf::st_intersects(line_sfc, poly_sfc))
      
      # convert to length of intersection if using method 2
      if (barrier_method == 2) {
        intersect_mat[intersect_mat == TRUE] <- mapply(function(x) sf::st_length(x), sf::st_intersection(line_sfc, poly_sfc))
      }
    }
    
    # if comparing ellipse
    if (barrier_method == 3) {
      
      # create all pairwise ellipses between nodes
      ell_list <- list()
      n_node <- length(node_long)
      i2 <- 1
      for (i in 1:(n_node-1)) {
        for (j in (i+1):n_node) {
          ell_df <- get_ellipse(f1 = node_mat[i,], f2 = node_mat[j,], ecc = eccentricity, n = n_ell)
          ell_list[[i2]] <- sf::st_polygon(list(as.matrix(ell_df)))
          i2 <- i2+1
        }
      }
      
      # convert ellipses and polys to st_sfc
      ell_sfc <- sf::st_sfc(ell_list)
      poly_sfc <- sf::st_sfc(poly_list)
      
      # get boolean intersection matrix
      intersect_mat <- as.matrix(sf::st_intersects(ell_sfc, poly_sfc))
      
      # convert to area of intersection
      intersect_mat[intersect_mat == TRUE] <- mapply(function(x) sf::st_area(x), sf::st_intersection(ell_sfc, poly_sfc))
      
    }
    
    # apply penalty
    intersect_penalty <- rowSums(sweep(intersect_mat, 2, barrier_penalty, '*'))
    
  }  # end apply barrier penalties
  
  # get pairwise distance plus penalty
  d <- dist(cbind(node_long, node_lat)) + intersect_penalty
  
  # apply transformation
  if (dist_transform == 2) {
    d <- exp(-lambda*d)
  }
  
  # add noise
  d <- d + stats::rnorm(length(d), sd = eps)
  
  # return matrix
  return(as.matrix(d))
}
