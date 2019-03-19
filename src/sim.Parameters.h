
#pragma once

#include <Rcpp.h>

#include "sim.Sampler.h"

//------------------------------------------------
// parameters of individual-based simulation model
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // genetic parameters
  static int L;
  static double prob_cotransmission;
  
  // epidemiological parameters
  static double a;
  static double mu;
  static int u;
  static int v;
  static int g;
  static std::vector<double> prob_infection;
  static int n_prob_infection;
  static std::vector<double> duration_infection;
  static int n_duration_infection;
  static double infectivity;
  static int max_innoculations;
  
  // deme parameters
  static int H;
  static std::vector<int> seed_infections;
  static std::vector<int> M_vec;
  static int n_demes;
  
  // demography
  static std::vector<double> life_table;
  static std::vector<double> age_death;
  static std::vector<double> age_stable;
  static int n_age;
  
  // objects for sampling from probability distributions
  static Sampler sampler_age_stable;
  static Sampler sampler_age_death;
  static Sampler sampler_duration_infection;
  
  // migration
  // TODO
  
  // run parameters
  static std::vector<int> time_out;
  static int n_time_out;
  static int max_time;
  
  // misc parameters
  static double prob_v_death;  // daily probability of mosquito death
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters() {};
  Parameters(const Rcpp::List &args);
  
  // methods
  void print_summary();
};