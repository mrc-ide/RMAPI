
#------------------------------------------------
# define empty rmapiProject object
# TODO - documentation
#' @export

rmapi_project <- function() {
  
  # initialise project with default values
  data <- NULL
  parameters = list(alpha = 1,
                    beta = 2,
                    xlimits = NULL)
  output <- list()
  
  # create class and return
  ret <- list(data = data, parameters = parameters, output = output)
  class(ret) <- "rmapi_project"
  return(ret)
}

#------------------------------------------------
# overload print() function to print summary only
# (not exported)

print.rmapi_project <- function(proj, ...) {
  
  # print summary only
  summary(proj)
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
# overload summary() function.
# (not exported)

summary.rmapi_project <- function(proj, ...) {
  
  # print data details
  cat("### DATA:\n")
  if (is.null(proj$data)) {
    cat("no data loaded\n")
  } else {
    cat(paste0("n = ", nrow(proj$data), " samples\n"))
    cat(paste0("X% missing data\n"))	# TODO - count missing data
  }
  cat("\n")
  
  # print parameters
  cat("### PARAMETERS:\n")
  cat(paste0("alpha:\t", proj$parameters$alpha, "\n"))
  cat(paste0("beta:\t", proj$parameters$beta, "\n"))
  cat("\n")
  
  # print output details
  cat("### OUTPUT:\n")
  if (length(proj$output)==0) {
    cat("no output\n")
  } else {
    cat("TODO - some details of output\n")
  }
}

#------------------------------------------------
# make new print function for doing standard print
# (not exported)

print_full <- function(proj) {
  
  # check that viewing rmapi_project
  assert_that(is.rmapi_project(proj))
  
  # print raw list
  print(unclass(proj))
}

#------------------------------------------------
# overload list assignment symbol. Assignment occurs one level down in the list, in the parameters sublist
# (not exported)

`$<-.rmapi_project` <- function(proj, i, value) {
  
  # check that i is a valid parameter name
  if ( !(i %in% c("alpha", "beta", "xlimits")) ) {
    stop(paste(i,"is not a valid parameter name"))
  }
  
  # assign value within the parameters sublist
  proj[["parameters"]][[i]] <- value
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
# create function for determining if object is of class rmapiProject
# (not exported)

is.rmapi_project <- function(x) {
  inherits(x, "rmapi_project")
}

