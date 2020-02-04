
#------------------------------------------------
# if NULL then replace with chosen value, otherwise keep original value
#' @noRd
define_default <- function(x, default) {
  if (is.null(x)) {
    return(default)
  } else {
    return(x)
  }
}

#------------------------------------------------
# if a single value is provided then expand to a vector of length n
#' @noRd
force_vector <- function(x, n) {
  if (length(x) == 1) {
    return(rep(x,n))
  } else {
    return(x)
  }
}

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE.
#' @noRd
user_yes_no <- function(x = "continue? (Y/N): ") {
  user_choice <- NA
  while (!user_choice %in% c("Y", "y" ,"N", "n")) {
    user_choice <- readline(x)
  }
  return(user_choice %in% c("Y", "y"))
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
matrix_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_matrix <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}

# -----------------------------------
# takes list format returned from Rcpp and converts to three-dimensional array.
# Array indexing is in the same order as the underlying list, for example
# x[i,j,k] is equivalent to l[[i]][[j]][[k]]
#' @noRd
rcpp_to_array <- function(x) {
  ret <- array(unlist(x), dim = c(length(x[[1]][[1]]), length(x[[1]]), length(x)))
  ret <- aperm(ret, perm = c(3,2,1))
  return(ret)
}

#------------------------------------------------
# return 95% quantile
#' @importFrom stats quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs = c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#------------------------------------------------
# sum logged values without underflow, i.e. do log(sum(exp(x)))
#' @noRd
log_sum <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  x_max <- max(x, na.rm = TRUE)
  ret <- x_max + log(sum(exp(x - x_max)))
  return(ret)
}

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
#' @title Print unclassed object
#'
#' @description Print object after unclassing, thereby removing any custom print
#'   method.
#'
#' @param x object to print in full.
#'
#' @export

print_full <- function(x) {
  print(unclass(x))
}

#------------------------------------------------
#' @title Get great circle distance between spatial points
#'
#' @description Get great circle distance between spatial points, defined by a
#'   vector of longitudes and latitudes. Distances are returned in a pairwise
#'   distance matrix.
#' 
#' @param long,lat vector of longitudes and latitudes.
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
#'   coordinates, defined by longitude and latitude of both origin and
#'   destination points.
#'
#' @param origin_lon,origin_lat the origin longitude and latitude.
#' @param dest_lon,dest_lat the destination longitude and latitude
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
  bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat) - sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  
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
# pass in a set of hex centroids and polygons, and hex size. Return a list of
# coordinates representing all contigous blocks (concave hulls) that can be made
# from these hexes.
#' @noRd
get_hex_hulls <- function(hex_centroids, hex_polys, hex_size) {
  
  #w <- w_signif_high
  #hex_centroids <- proj$map$hex_centroid[w,]
  #hex_polys <- proj$map$hex[w]
  #hex_size <- proj$map$hex_size
  
  # get list of the neighbours of every hex. All adjacent hexes should be a
  # distance hex_size away from the centre (give or take floating point errors),
  # and the minimum distance to the nearest non-adjacent hex is hex_size*sqrt(3).
  # Hence, but using a distance 1.5*hex_size we ensure that we capture all
  # adjacent hexes.
  hex_adj_list <- get_neighbours(hex_centroids, 1.5*hex_size)
  
  # get list of contiguous blocks
  block_list <- get_contiguous_blocks(hex_adj_list)
  
  # apply method to each block
  conc_hull_list <- list()
  for (b in 1:length(block_list)) {
    
    # if a single poly then simply append these coordinates and skip ahead
    if (length(block_list[[b]]) == 1) {
      conc_hull_list[[b]] <- hex_polys[[ block_list[[b]][1] ]][[1]]
      next()
    }
    
    # get list of polygon coordinates, removing final coordinate which is
    # identical to first
    poly_list <- mapply(function(i) {
      ret <- hex_polys[[i]][[1]]
      ret[-1,]
    }, block_list[[b]], SIMPLIFY = FALSE)
    
    # convert list to dataframe over all polys
    poly_df <- as.data.frame(do.call(rbind, poly_list))
    names(poly_df) <- c("long", "lat")
    
    # get all possible unique long and lat coordinates of points
    unique_long_lat <- mapply(function(i) sort(unique(poly_df[,i])), 1:2, SIMPLIFY = FALSE)
    
    # determine the smallest possible distance between any two different points.
    # Anything smaller than this distance is actually on the same point (but may
    # not appear so due to floating point issues)
    smallest_diff <- hex_size/(2*sqrt(3))
    
    # drop from unique_long_lat any points that are within smallest_diff of
    # another point
    unique_long_lat <- mapply(function(x) {
      bad_diff <- which(diff(x) < 0.5*smallest_diff)
      if (length(bad_diff) == 0) {
        return(x)
      } else {
        return(x[-(bad_diff+1)])
      }
    }, unique_long_lat, SIMPLIFY = FALSE)
    
    # create an exhaustive list of all possible long/lat combinations
    egrid <- expand.grid(long = unique_long_lat[[1]], lat = unique_long_lat[[2]])
    
    # each point in poly_df must index a point in egrid. Find this index
    poly_df$index <- mapply(function(i) {
      which.min((poly_df$long[i] - egrid$long)^2 + (poly_df$lat[i] - egrid$lat)^2)
    }, 1:nrow(poly_df))
    
    # points on the concave hull are represented a maximum of twice in poly_df.
    # Use this fact to subset to unique points on the concave hull
    tab <- tabulate(poly_df$index)
    w <- which(tab == 1 | tab == 2)
    sub_df <- poly_df[poly_df$index %in% w,]
    sub_df <- sub_df[!duplicated(sub_df$index),]
    
    #plot(sub_df$long, sub_df$lat)
    
    # get all adjacent neighbours of each point
    point_adj_list <- get_neighbours(sub_df[,1:2], 1.1*hex_size/sqrt(3))
    
    # calculate the concave hull around the final points
    vertex_order <- concave_hull_from_neighbours(sub_df$long, sub_df$lat, point_adj_list)
    conc_hull <- sub_df[vertex_order, 1:2]
    
    # store in list
    conc_hull_list[[b]] <- conc_hull
  }
  
  # return list
  return(conc_hull_list)
}

#------------------------------------------------
# given a set of coordinates and a distance r, return a list of all neighbours
# of each coordinate, defined as points within a Euclidean distance r
#' @noRd
get_neighbours <- function(coords, r) {
  
  # get Euclidean distance between all coordinates
  d <- as.matrix(dist(coords))
  diag(d) <- Inf
  
  # find all other coordinates that are within a distance r
  adj_list <- mapply(function(i) which(d[i,] < r), 1:nrow(d), SIMPLIFY = FALSE)
  
  return(adj_list)
}

#------------------------------------------------
# given a list over containing the neighbours of each node, return a list of
# contiguous blocks such that all nodes have naighbours within a block and no
# nodes have neighbours between blocks
#' @noRd
get_contiguous_blocks <- function(neighbour_list) {
  
  # get basic properties
  n <- length(neighbour_list)
  
  # output will be a list of potentially multiple contiguous blocks
  block_list <- list()
  
  # repeatedly identify blocks, until a cutout limit is reached
  for (j in 1:1e3) {
    if (j == 1e3) {
      stop("error in get_contiguous_blocks(): reached maximum allowed j")
    }
    
    # start the current block with the first remaining value that has not been
    # allocated to any existing block...
    block <- setdiff(1:n, unlist(block_list))[1]
    block_size <- 1
    
    # ...if there is no such block then return
    if (is.na(block)) {
      return(block_list)
    }
    
    # repeatedly add to this block, until a cutout limit is reached
    for (i in 1:1e3) {
      if (i == 1e3) {
        stop("error in get_contiguous_blocks(): reached maximum allowed i")
      }
      
      # concatenate this block with all of its neighbours
      block <- unique(c(block, unlist(neighbour_list[block])))
      
      # if there was no change in block size then exit
      if (length(block) == block_size) {
        break()
      }
      
      # update block size
      block_size <- length(block)
    }
    
    # store current block
    block_list[[j]] <- block
  }
  
  # should return before reaching this point
  stop("error in get_contiguous_blocks(): reached end of function")
}

#------------------------------------------------
# given a series of coordinates x and y, and a list of the same length that
# contains the allowed neighbours of each point, calculate a concave hull around
# the points. Return the order in which the original x and y should be arranged
# to produce a polygon.
#' @noRd
concave_hull_from_neighbours <- function(x, y, neighbours) {
  
  #x <- sub_df$long
  #y <- sub_df$lat
  #neighbours <- point_adj_list
  
  # get basic properties
  n <- length(x)
  
  # start algorithm at node with greatest y-value
  j <- which.max(y)
  vertex_order <- j
  theta <- 0
  
  # loop through remaining nodes
  #for (i in 1:(n-1)) {
  for (i in 1:103) {
    print(i)
    
    # find all neighbours of the current point. Discard those that have already
    # been explored
    j_prop <- neighbours[[j]]
    j_prop <- setdiff(j_prop, vertex_order)
    
    # there should be at least one proposed point. If not we have a bug
    if (length(j_prop) == 0) {
      stop("error in concave_hull_from_neighbours(): no proposed points")
    }
    
    # calculate angle from current point j to all proposed points. Angle is
    # calculated relative to theta
    clockwise_angle <- atan2(x[j_prop] - x[j], y[j_prop] - y[j])
    angle_from_prev <- clockwise_angle - theta
    
    # ensure angle is between 0 and 2*pi radians
    angle_from_prev[angle_from_prev < 0] <- angle_from_prev[angle_from_prev < 0] + 2*pi
    angle_from_prev[angle_from_prev > 2*pi] <- angle_from_prev[angle_from_prev > 2*pi] - 2*pi
    
    # find the smallest angle
    w <- which.min(angle_from_prev)
    
    # update running values
    j <- j_prop[w]
    theta <- clockwise_angle[w] - pi
    
    # store vertices
    vertex_order <- c(vertex_order, j)
  }
  
  plot(x,y)
  lines(x[vertex_order], y[vertex_order])
  
  # complete ring by making final vertex equal to the first
  vertex_order <- c(vertex_order, vertex_order[1])
  
  return(vertex_order)
}