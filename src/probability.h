
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// take in vector of integers (by reference). Resample these without replacement
void reshuffle(std::vector<int> &x);

//------------------------------------------------
// resample a vector of integers without replacement while always sampling from within specified groups
void reshuffle_group(std::vector<int> &x, std::vector<std::vector<int>> &group_list);