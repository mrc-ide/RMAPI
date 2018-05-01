
#include "probability.h"

using namespace std;

// take in vector of integers (by reference). Resample these without replacement
void reshuffle(vector<int> &x) {
  int rnd1, tmp1; // dummy variables
  int n = x.size();
  for (int i=0; i<n; i++) {
    rnd1 = R::runif(i, n);	// draw random index from i to end of vector. Note that although runif returns a double, by forcing to int we essentially round this value down to nearest int.
    tmp1 = x[rnd1];		    // temporarily store current value of vector at this position
    x[rnd1] = x[i];			// swap for value at position i
    x[i] = tmp1;			// complete the swap
  }
}
