
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// take in vector of integers (by reference). Resample these without replacement
void reshuffle(std::vector<int> &x);
