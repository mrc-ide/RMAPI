
#------------------------------------------------
#' Load data into RMAPI project
#'
#' TODO - some help text here.
#'
#' @param proj the current RMAPI project
#' @param data a data frame, formatted into the correct RMAPI format (see Details)
#' @param checkDeleteOutput whether to perform a check to see if project already contains data, in which case all old data and output will be lost
#'
#' @export

loadData <- function(proj, data, checkDeleteOutput=TRUE) {
    
    # check whether there is data loaded already
    if (!is.null(proj$data)) {
        
        # return existing project if user not happy to continue
        if (checkDeleteOutput) {
            if (!user_yesNo("All existing output for this project will be lost. Continue? (Y/N): ")) {
                cat("returning original project\n")
                return(proj)
            }
        }
        
        # delete old output
        cat("overwriting data\n")
        proj[["output"]] <- NULL
    }
    
    # perform checks on data
    stopifnot( is.data.frame(data) )
    stopifnot( ncol(data)>=4 )
    stopifnot( nrow(data)==ncol(data)-3 )
    
    # update project with new data
    proj[["data"]] <- data
    
    # return invisibly
    invisible(proj)
}

#------------------------------------------------
#' Perform RMAPI simulation
#'
#' TODO - some help text here.
#'
#' @param proj the current RMAPI project
#'
#' @export
#' @examples
#' runSims()

runSims <- function(proj, Nperms=1e2) {
	
	# Node and other parameters as R objects -----------------------------------------
	
	a_multiplier_i=-0.45						#Controls relationship between a (ellipse long radius) and c (ellipse short radius equal to distance between foci): a = c*(1 + a_multiplier)
	
	# set default x and y limits
	if (is.null(proj$parameters$xlimits)) {
		# TODO - set default limits
	}
	
	# Set up arguments for input into C++ ---------------------------------------------
	
	pairwiseStats <- as.matrix(proj$data[,4:ncol(proj$data)])
	args <- list(xnode=proj$data$long, ynode=proj$data$lat, vnode=mat_to_Rcpp(pairwiseStats), a_multiplier=a_multiplier_i, Nperms=Nperms) 
	
	# Carry out simulations in C++ to generate map data-------------------------------
	
	output_raw <- dummy1_cpp(args)	# TODO - rename C++ function to "runSims_cpp"
	
	# process raw output---------------------------------------------------------------	
	
	proj[["output"]] <- output_raw
	
    # return invisibly
    invisible(proj)
}







