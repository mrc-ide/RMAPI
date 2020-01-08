
#### Member functions for objects of class rmapi_sim

#------------------------------------------------
# Overload print()
#' @method print rmapi_sim
#' @export
print.rmapi_sim <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# Overload summary()
#' @method summary rmapi_sim
#' @export
summary.rmapi_sim <- function(object, ...) {
  
  # print data details
  message("Simulation details:")
  message(sprintf("  demes = %s", length(object$daily_values)))
  message(sprintf("  max time = %s", nrow(object$daily_values[[1]])))
  
}

#------------------------------------------------
# Overload plot()
#' @method plot rmapi_project
#' @export
plot.rmapi_sim <- function(x, y, ...) {
  plot_daily_states(x)
}

#------------------------------------------------
#' @title Determine if object is of class rmapi_sim
#'
#' @description Determine if object is of class \code{rmapi_sim}.
#'
#' @param x object to query.
#'
#' @export

is.rmapi_sim <- function(x) {
  inherits(x, "rmapi_sim")
}

