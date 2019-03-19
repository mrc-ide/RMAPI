
#pragma once

#include "sim.Parameters.h"

#include <vector>

//------------------------------------------------
// class defining mosquito
class Mosquito : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // genotypes
  std::vector<int> haplotype1;
  std::vector<int> haplotype2;
  int n_haplotypes;
  std::vector<std::vector<int>> products;
  int n_products;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito();
  
  // methods
  void denovo();
  void die();
  std::vector<int> get_product();
  
};
