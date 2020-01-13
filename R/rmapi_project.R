
#### Member functions for objects of class rmapi_project

#------------------------------------------------
# Overload print()
#' @method print rmapi_project
#' @export
print.rmapi_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# Overload summary()
#' @method summary rmapi_project
#' @export
summary.rmapi_project <- function(object, ...) {
  
  # print data details
  message("# DATA:")
  if (is.null(object$data$coords)) {
    message("no data loaded")
  } else {
    message(sprintf("n = %s samples", nrow(object$data$coords)))
  }
  message("")
  
  # print map details
  message("# HEX MAP:")
  if (is.null(object$map)) {
    message("no map generated")
  } else {
    message(sprintf("h = %s hexagons", length(object$map$hex)))
  }
  message("")
  
  # print output details
  message("# OUTPUT:")
  if (length(object$output) == 0) {
    message("no output")
  } else {
    message("TODO - some details of output")
  }
}

#------------------------------------------------
# Overload plot()
#' @method plot rmapi_project
#' @export
plot.rmapi_project <- function(x, y, ...) {
  plot_map(x)
}

#------------------------------------------------
#' @title Determine if object is of class rmapi_project
#'
#' @description Determine if object is of class \code{rmapi_project}.
#'
#' @param x object to query.
#'
#' @export

is.rmapi_project <- function(x) {
  inherits(x, "rmapi_project")
}

