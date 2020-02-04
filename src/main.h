
#pragma once

#include "misc_v4.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// check if value is within ellipse
bool ellipse_check(double x, double y, double xf1, double yf1, double xf2, double yf2, double a);

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