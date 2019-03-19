
#include "sim.Mosquito.h"
#include "probability.h"
#include "misc_v4.h"

using namespace std;

//------------------------------------------------
// constructors
Mosquito::Mosquito() {
  n_haplotypes = 0;
  n_products = 0;
}

//------------------------------------------------
// initialise with single haplotype of random SNPs
void Mosquito::denovo() {
  haplotype1 = vector<int>(L);
  for (int i=0; i<L; ++i) {
    haplotype1[i] = rbernoulli1(0.5) ? 1 : 0;
  }
  n_haplotypes = 1;
}

//------------------------------------------------
// death
void Mosquito::die() {
  haplotype1.clear();
  haplotype2.clear();
  n_haplotypes = 0;
  products.clear();
  n_products = 0;
}

//------------------------------------------------
// get product of recombination
vector<int> Mosquito::get_product() {
  if (n_haplotypes == 0) {
    Rcpp::stop("attempt to get product from uninfected mosquito");
  }
  
  // if single haplotype then return clone
  //if (n_haplotypes == 1) {
    return haplotype1;
  //}
  /*
  // determine whether this product is novel
  bool new_product = rbernoulli1(1 - 0.25*n_products);
  
  // if new product then 
  if (new_product) {
    
  } else {
    
  }
  */
}
