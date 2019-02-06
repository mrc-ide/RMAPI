
#include "probability.h"

using namespace std;

//------------------------------------------------
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

//------------------------------------------------
// resample a vector of integers without replacement while always sampling from within specified groups
void reshuffle_group(vector<int> &x, vector<vector<int>> &group_list) {
  for (int j=0; j<int(group_list.size()); ++j) {
    int rnd1, tmp1, old_pos, new_pos;
    int n = group_list[j].size();
    for (int i=0; i<n; i++) {
      rnd1 = R::runif(i, n);
      old_pos = group_list[j][i];
      new_pos = group_list[j][rnd1];
      tmp1 = x[new_pos];
      x[new_pos] = x[old_pos];
      x[old_pos] = tmp1;
    }
  }
}