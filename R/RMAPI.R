#' @title RMAPI
#'
#' @description An Rcpp implementation of the MAPI software, which can be used
#'   to infer barriers/corridors of gene flow (as one example application).
#' 
#' @name RMAPI
#' @docType package
NULL

#------------------------------------------------
#' @useDynLib RMAPI, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("RMAPI", libpath)  # nocov
}