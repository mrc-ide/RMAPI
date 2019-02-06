
#' #------------------------------------------------
#' @title Define empty RMAPI project
#'
#' @description Define empty RMAPI project.
#'
#' @export

rmapi_project <- function() {
  
  # initialise project with default values
  ret <- list(data = list(coords = NULL,
                          stat_dist = NULL,
                          spatial_dist = NULL),
              map = NULL,
              output = NULL)
  
  # create class and return
  class(ret) <- "rmapi_project"
  return(ret)
}

#------------------------------------------------
#' @title Overload print function for RMAPI project
#'
#' @description Overload print function for RMAPI project.
#'
#' @export

print.rmapi_project <- function(proj, ...) {
  
  # print summary
  summary(proj)
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Overload summary function for RMAPI project
#'
#' @description Overload summary function for RMAPI project.
#'
#' @export

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
#' @title Print unclassed RMAPI project
#'
#' @description Print unclassed RMAPI project.
#'
#' @export

print_full <- function(proj) {
  
  # check that viewing rmapi_project
  assert_custom_class(proj, "rmapi_project")
  
  # print raw list
  print(unclass(proj))
}

#------------------------------------------------
#' @title Overload plot function for RMAPI project
#'
#' @description Overload plot function for RMAPI project.
#'
#' @export

plot.rmapi_project <- function(x, y, ...) {
  plot_map(x)
}

#------------------------------------------------
#' @title Determine if object is of class rmapi_project
#'
#' @description Determine if object is of class rmapi_project.
#'
#' @export

is.rmapi_project <- function(x) {
  inherits(x, "rmapi_project")
}

