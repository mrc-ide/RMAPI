
#------------------------------------------------
#' Define empty RMAPI project
#'
#' TODO - some text here.
#'
#' @export

rmapi_project <- function() {
  
  # initialise project with default values
  data <- NULL
  map <- list()
  output <- list()
  
  # create class and return
  ret <- list(data = data, map = map, output = output)
  class(ret) <- "rmapi_project"
  return(ret)
}

#------------------------------------------------
# overload print() function to print only certain elements
# (not exported)

print.rmapi_project <- function(proj, ...) {
  
  # print certain elements
  print_full(proj)
  
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
  
  # print map details
  cat("### HEX MAP:\n")
  cat("(TODO)\n")
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
# create function for determining if object is of class rmapiProject
# (not exported)

is.rmapi_project <- function(x) {
  inherits(x, "rmapi_project")
}

