
#------------------------------------------------
# test if integer
#' @noRd
is.int <- function(x) {
  as.integer(x)==x
}

#------------------------------------------------
# test if positive scalar value
#' @noRd
is.pos_scalar <- function(x) {
  is.numeric(x) & length(x)==1 & all(x>0)
}

#------------------------------------------------
# replace NULL value with default
#' @noRd
define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
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
mat_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}