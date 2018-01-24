
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib RMAPI
#' @importFrom Rcpp evalCpp
NULL

# -----------------------------------
# mat_to_Rcpp
# takes matrix as input, converts to list format for use within Rcpp code
# (not exported)

mat_to_Rcpp <- function(x) {
    return(split(x,f=1:nrow(x)))
}

# -----------------------------------
# Rcpp_to_mat
# Takes list format returned from Rcpp and converts to matrix.
# (not exported)

Rcpp_to_mat <- function(x) {
    ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
    return(ret)
}

# -----------------------------------
# user_yesNo
# Ask user a yes/no question. Return TRUE/FALSE.
# (not exported)

user_yesNo <- function(x="continue? (Y/N): ") {
    
    userChoice <- NA
    while (!userChoice %in% c("Y", "y" ,"N", "n")) {
        userChoice <- readline(x)
    }
    return(userChoice %in% c("Y", "y"))
}
