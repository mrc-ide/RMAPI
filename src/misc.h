
#pragma once

#include <Rcpp.h>
#include <chrono>

//------------------------------------------------
// define very large/small numbers for catching overflow/underflow problems
#define OVERFLO   1e100
#define UNDERFLO   1e-100

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types).
template<class TYPE>
TYPE sum(const std::vector<TYPE> &x) {
  TYPE output = 0;
  for (int i=0; i<int(x.size()); i++) {
    output += x[i];
  }
  return output;
}

//------------------------------------------------
// mean of vector (templated for different data types)
template<class TYPE>
double mean(const std::vector<TYPE> &x) {
  return sum(x)/double(x.size());
}

//------------------------------------------------
// min of vector (templated for different data types)
template<class TYPE>
TYPE min(const std::vector<TYPE> x) {
  return *min_element(x.begin(), x.end());
}

//------------------------------------------------
// max of vector (templated for different data types)
template<class TYPE>
TYPE max(const std::vector<TYPE> x) {
  return *max_element(x.begin(), x.end());
}

//------------------------------------------------
// push back multiple values to vector
template<class TYPE>
void push_back_multiple(std::vector<TYPE> &lhs, const std::vector<TYPE> &rhs) {
  lhs.insert(lhs.end(), rhs.begin(), rhs.end());
}

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double log_sum(const double logA, const double logB);

//------------------------------------------------
// helper function for printing a single value or series of values (templated for different data types)
template<class TYPE>
void print(const TYPE x) {
  Rcpp::Rcout << x << "\n";
  R_FlushConsole();
}
template<class TYPE1, class TYPE2>
void print(const TYPE1 x1, const TYPE2 x2) {
  Rcpp::Rcout << x1 << " " << x2 << "\n";
  R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3>
void print(const TYPE1 x1, const TYPE2 x2, const TYPE3 x3) {
  Rcpp::Rcout << x1 << " " << x2 << " " << x3 << "\n";
  R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4>
void print(const TYPE1 x1, const TYPE2 x2, const TYPE3 x3, const TYPE4 x4) {
  Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << "\n";
  R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5>
void print(const TYPE1 x1, const TYPE2 x2, const TYPE3 x3, const TYPE4 x4, const TYPE5 x5) {
  Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << "\n";
  R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5, class TYPE6>
void print(const TYPE1 x1, const TYPE2 x2, const TYPE3 x3, const TYPE4 x4, const TYPE5 x5, const TYPE6 x6) {
  Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
template<class TYPE>
void print_vector(const std::vector<TYPE> &x) {
  for (int i=0; i<int(x.size()); i++) {
    Rcpp::Rcout << x[i] << " ";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
template<class TYPE>
void print_matrix(const std::vector< std::vector<TYPE> > &x) {
  for (int i=0; i<int(x.size()); i++) {
    for (int j=0; j<int(x[i].size()); j++) {
      Rcpp::Rcout << x[i][j] << " ";
    }
    Rcpp::Rcout << "\n";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
template<class TYPE>
void print_array(const std::vector< std::vector< std::vector<TYPE> > > &x) {
  for (int i=0; i<int(x.size()); i++) {
    Rcpp::Rcout << "--- slice " << i+1 << " ---\n";
    for (int j=0; j<int(x[i].size()); j++) {
      for (int k=0; k<int(x[i][j].size()); k++) {
        Rcpp::Rcout << x[i][j][k] << " ";
      }
      Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void print_stars(const std::string title = "", const int n = 10);

//------------------------------------------------
// print "foo", with option number e.g. "foo2"
void foo(const int n = 0);

//------------------------------------------------
// print "bar", with option number e.g. "bar2"
void bar(const int n = 0);

//------------------------------------------------
// print "foobar", with option number e.g. "foobar2"
void foobar(const int n = 0);

//------------------------------------------------
// analogue of R function seq() for integers
std::vector<int> seq_int(int from, const int to, const int by = 1);

//------------------------------------------------
// converts input from Rcpp::List format to bool format.
int rcpp_to_bool(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to int format.
int rcpp_to_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to double format.
double rcpp_to_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to string format.
std::string rcpp_to_string(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<int> format.
std::vector<int> rcpp_to_vector_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<double> format.
std::vector<double> rcpp_to_vector_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<string> format.
std::vector<std::string> rcpp_to_vector_string(SEXP x);

// converts input from Rcpp::List format to vector<vector<bool>> format.
std::vector< std::vector<bool> > rcpp_to_mat_bool(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<int>> format.
std::vector< std::vector<int> > rcpp_to_mat_int(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
std::vector< std::vector<double> > rcpp_to_mat_double(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
std::vector< std::vector< std::vector<double> > > rcpp_to_array_double(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
std::vector< std::vector< std::vector<int> > > rcpp_to_array_int(Rcpp::List x);

//------------------------------------------------
// return timer
void chrono_timer(std::chrono::high_resolution_clock::time_point &t0);

//------------------------------------------------
// check if value is within ellipse
bool ellipse_check(const double x, const double y, const double xf1, const double yf1, const double xf2, const double yf2, const double a);

//------------------------------------------------
// square function (templated for different data types)
template<class TYPE>
TYPE sq(const TYPE x) {
  return x*x;
}

//------------------------------------------------
// Euclidian distance between points in 2 dimensions (templated for different data types)
template<class TYPE>
double dist_euclid_2d(const TYPE x1, const TYPE y1, const TYPE x2, const TYPE y2) {
  double ret = sqrt(sq(x1 - x2) + sq(y1 - y2));
  return ret;
}