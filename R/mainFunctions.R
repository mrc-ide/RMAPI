
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib RMAPI
#' @import assertthat
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

load_data <- function(proj, data, check_delete_output = TRUE) {
  
  # check inputs
  assert_that( is.rmapi_project(proj) )
  assert_that( is.data.frame(data) )
  assert_that( ncol(data)>=4 )
  assert_that( nrow(data)==ncol(data)-3 )
  
  # check whether there is data loaded already
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
    proj[["output"]] <- NULL
  }
  
  # update project with new data
  proj[["data"]] <- data
  
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
#' @param a_multiplier controls relationship between a (ellipse long radius) and c (ellipse short radius equal to distance between foci): a = c*(1 + a_multiplier)
#'
#' @export
#' @examples
#' run_sims()

run_sims <- function(proj, Nperms = 0, a_multiplier = 0.01, dim_matrix = 101) {
  
  # TODO - set default x and y limits
  
  # ---------------------------------------------
  # Set up arguments for input into C++
  
  pairwise_stats <- as.matrix(proj$data[,4:ncol(proj$data)])
  args <- list(xnode = proj$data$long,
               ynode = proj$data$lat,
               vnode = mat_to_rcpp(pairwise_stats),
               a_multiplier = a_multiplier,
               dim_matrix = dim_matrix,
               Nperms = Nperms) 
  
  # ---------------------------------------------
  # Carry out simulations in C++ to generate map data
  output_raw <- run_sims_cpp(args)
  
  # ---------------------------------------------
  # process raw output
  proj[["output"]] <- output_raw
  
  # return invisibly
  invisible(proj)
}

