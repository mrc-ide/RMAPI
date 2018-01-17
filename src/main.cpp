
#include <Rcpp.h>

using namespace std;

// [[Rcpp::export]]
Rcpp::List dummy1_cpp(Rcpp::List args) {

	// convert Rcpp args to native c++ args
	int foo = Rcpp::as<int>(args["foo"]);
	vector<int> bar = Rcpp::as< vector<int> >(args["bar"]);

	//------------------------------------------------

	// do some really hard maths
	for (int i = 0; i < bar.size(); i++) {
		bar[i] *= foo;
	}

	//------------------------------------------------

	// return output as an Rcpp list
	Rcpp::List ret;
	ret.push_back(Rcpp::wrap(bar));

	Rcpp::StringVector ret_names;
	ret_names.push_back("bar");

	ret.names() = ret_names;
    return ret ;
}
