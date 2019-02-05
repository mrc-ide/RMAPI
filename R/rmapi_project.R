
#' #------------------------------------------------
#' @title Define empty RMAPI project
#'
#' @description Define empty RMAPI project.
#'
#' @export

rmapi_project <- function() {
  
  # initialise project with default values
  ret <- list(data = list(coords = NULL,
                          pairwise_dist = NULL),
              map = NULL,
              output = NULL)
  
  # create class and return
  class(ret) <- "rmapi_project"
  return(ret)
}

#------------------------------------------------
# overload print() function
#' @noRd
print.rmapi_project <- function(proj, ...) {
  
  # print summary
  summary(proj)
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
# overload summary() function.
#' @noRd
summary.rmapi_project <- function(proj, ...) {
  
  # print data details
  message("# DATA:")
  if (is.null(proj$data$coords)) {
    message("no data loaded")
  } else {
    message(sprintf("n = %s samples", nrow(proj$data$coords)))
  }
  message("")
  
  # print map details
  message("# HEX MAP:")
  if (is.null(proj$map)) {
    message("no map generated")
  } else {
    message(sprintf("h = %s hexagons", length(proj$map$hex)))
  }
  message("")
  
  # print output details
  message("# OUTPUT:")
  if (length(proj$output) == 0) {
    message("no output")
  } else {
    message("TODO - some details of output")
  }
}

#------------------------------------------------
# make new print function for standard print
#' @noRd
print_full <- function(proj) {
  
  # check that viewing rmapi_project
  assert_custom_class(proj, "rmapi_project")
  
  # print raw list
  print(unclass(proj))
}

#------------------------------------------------
# determine if object is of class rmapi_project
#' @noRd
is.rmapi_project <- function(x) {
  inherits(x, "rmapi_project")
}

