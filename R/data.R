
#------------------------------------------------
#' example data set to test RMAPI functionality
#'
#' Made up data for testing - this data set should be deleted prior to release.
#' \cr
#' \cr
#' Spatial points correspond to real DRC DHS cluster locations, giving some spatial realism. Values between points are completely made up, and represent a barrier to gene flow transecting the country from North to South. Hopefully RMAPI should be able to detect this.
#'
#' @docType data
#'
#' @examples
#' # create project and load in data
#' data(fakeData_DRC)
#' proj <- rmapi_project()
#' proj <- loadData(proj, fakeData_DRC)
#'
#' # define parameters
#' proj$alpha <- 2
#'
#' # run simulations
#' proj <- runSims(proj)
#'
#' # plot output
#' RMAPI_plot1(proj)
#'
#' @usage data(fakeData_DRC)
"fakeData_DRC"
