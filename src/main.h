
#pragma once

//------------------------------------------------
// check if value is within ellipse
bool ellipse_check(const double x, const double y, const double xf1, const double yf1, const double xf2, const double yf2, const double a);

//------------------------------------------------
// compute map and run permutation test
Rcpp::List run_sims_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);
