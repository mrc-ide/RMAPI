
#------------------------------------------------
#' example data set to test RMAPI functionality
#'
#' Made up data for testing - this data set should be deleted prior to release.
#' \cr
#' \cr
#' Spatial points corresponding to real DRC demographic and health survey cluster locations. These points can be used to generate simulated data under different spatial barriers using the \code{sim_custom()} function.
#'
#' @docType data
#'
#' @examples
#' # TODO
#'
#' @usage data(coords_DRC)
"coords_DRC"

#------------------------------------------------
#' Load "simulation" data from original MAPI paper

  load_old_data1 <- function()
{
  coord_data=read.table("C:/Users/Kjfras16/Desktop/NewGitHub/RMAPI/data/old_coord_data.txt")
  x=coord_data["V1"]$V1
  y=coord_data["V2"]$V2
  coords <- cbind(x, y)
  colnames(coords) <- c("long", "lat")
  v_data=read.table("C:/Users/Kjfras16/Desktop/NewGitHub/RMAPI/data/old_v_data.txt")
  v=v_data["V1"]$V1
  sim_data=matrix(v,nrow=200,byrow=TRUE)

  # produce final return object
  cluster_names <- paste0("cluster", 1:nrow(coords))
  ret <- cbind(data.frame(name = cluster_names, stringsAsFactors = FALSE), coords, sim_data)
  
  # return
  return(ret)
}
#------------------------------------------------
#' Simulate rectangular lattice of points with one or more barriers
#'
#' TODO - some help text here.
#'
#' @param x positions of points in the x-dimension
#' @param y positions of points in the y-dimension
#' @param barrier_x TODO
#' @param barrier_y TODO
#' @param barrier_angle TODO
#' @param barrier_penalty TODO
#'
#' @export

sim_custom <- function(x, y, barrier_x = 0, barrier_y = 0, barrier_angle = 0, barrier_penalty = 10) {
  
  # check inputs
  assert_that( all(is.numeric(x)) )
  assert_that( all(is.numeric(y)) )
  assert_that( length(x)==length(y) )
  assert_that( all(is.numeric(barrier_x)) )
  assert_that( all(is.numeric(barrier_y)) )
  assert_that( all(is.numeric(barrier_angle)) )
  assert_that( all(is.numeric(barrier_penalty)) )
  
  # recycle inputs if different lengths
  Nbarrier <- max(mapply(length, list(barrier_x, barrier_y, barrier_angle, barrier_penalty)))
  barrier_x <- rep(barrier_x, Nbarrier)[1:Nbarrier]
  barrier_y <- rep(barrier_y, Nbarrier)[1:Nbarrier]
  barrier_angle <- rep(barrier_angle, Nbarrier)[1:Nbarrier]
  barrier_penalty <- rep(barrier_penalty, Nbarrier)[1:Nbarrier]
  
  # calculate Euclidian distance between nodes
  coords <- cbind(x, y)
  colnames(coords) <- c("long", "lat")
  sim_stat <- as.matrix(dist(coords))
  
  # apply barriers
  for (i in 1:Nbarrier) {
    
    # get shortest distance to barrier
    theta_rad <- (barrier_angle[i] %% 360)/360*2*pi
    x_intercept <- barrier_x[i] - barrier_y[i]*tan(theta_rad)
    dist_barrier <- y*sin(theta_rad) - (x-x_intercept)*cos(theta_rad)
    
    # find which points lie each side of barrier
    w1 <- which(dist_barrier>=0)
    w2 <- which(dist_barrier<=0)
    
    # apply penalty across barrier
    sim_stat[w1,w2] <- sim_stat[w1,w2] + barrier_penalty[i]
    sim_stat[w2,w1] <- sim_stat[w2,w1] + barrier_penalty[i]
  }
  
  # keep upper diagonal only
  sim_stat[row(sim_stat)>=col(sim_stat)] <- NA
  
  # produce final return object
  cluster_names <- paste0("cluster", 1:nrow(coords))
  ret <- cbind(data.frame(name = cluster_names, stringsAsFactors = FALSE), coords, sim_stat)
  
  # return
  return(ret)
}

#------------------------------------------------
#' Simulate rectangular lattice of points with one or more barriers
#'
#' TODO - some help text here.
#'
#' @param Nx number of points in the x-dimension
#' @param Ny number of points in the y-dimension
#' @param sep TODO
#' @param barrier_x TODO
#' @param barrier_y TODO
#' @param barrier_angle TODO
#' @param barrier_penalty TODO
#'
#' @export

sim_rect <- function(Nx = 10, Ny = 10, sep = 1, barrier_x = 0, barrier_y = 0, barrier_angle = 0, barrier_penalty = 10) {
  
  # check inputs
  assert_that( is.pos_scalar(Nx) )
  assert_that( is.int(Nx) )
  assert_that( is.pos_scalar(Ny) )
  assert_that( is.int(Ny) )
  assert_that( is.pos_scalar(sep) )
  
  # generate coordinates
  x_vec <- 1:Nx*sep
  x_vec <- x_vec - mean(x_vec)
  y_vec <- 1:Ny*sep
  y_vec <- y_vec - mean(y_vec)
  coords <- as.matrix(expand.grid(x_vec, y_vec))
  x <- coords[,1]
  y <- coords[,2]
  
  # generate data
  ret <- sim_custom(x, y, barrier_x = barrier_x, barrier_y = barrier_y, barrier_angle = barrier_angle, barrier_penalty = barrier_penalty)
  
  # return
  return(ret)
}