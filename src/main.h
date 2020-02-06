
#pragma once

#include "misc_v4.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif


//------------------------------------------------
// assign edges to hexes based on intersection
// [[Rcpp::export]]
Rcpp::List assign_map_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);

//------------------------------------------------
// compute map and run permutation test
// [[Rcpp::export]]
Rcpp::List rmapi_analysis_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);

//------------------------------------------------
// simulate from simple individual-based model
// [[Rcpp::export]]
Rcpp::List sim_falciparum_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);