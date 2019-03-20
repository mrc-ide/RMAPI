
#pragma once

#include <set>

#include "sim.Parameters.h"
#include "sim.Host.h"
#include "sim.Mosquito.h"
#include "misc_v4.h"
#include "array.h"

//------------------------------------------------
// class defining individual-based simulation model
class Dispatcher : public Parameters  {
  
public:
  
  // PUBLIC OBJECTS
  
  // scheduler objects
  std::vector<std::set<int>> schedule_death;
  std::vector<std::vector<std::pair<int, int>>> schedule_Eh_to_Ih;
  std::vector<std::vector<std::pair<int, int>>> schedule_Ih_to_Sh;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective_recovery;
  
  // counts of host types
  int H_total;
  std::vector<int> Sh;
  std::vector<int> Eh;
  std::vector<int> Ih;
  
  // population of human hosts
  std::vector<Host> host_pop;
  int next_host_ID;
  
  // store the integer index of hosts in each deme
  array_2d_int host_index;
  std::vector<std::vector<int>> host_infective_index;
  
  // counts of mosquito types
  int M_total;
  std::vector<int> Sv;
  std::vector<int> Ev;
  std::vector<int> Iv;
  
  // population of mosquitoes
  std::vector<Mosquito> mosq_pop;
  
  // store integer index of mosquitoes at various stages
  array_2d_int Sv_index;
  array_3d_int Ev_death;
  array_3d_int Ev_to_Iv;
  array_2d_int Iv_index;
  
  // objects for storing results
  array_3d_double daily_values;
  array_4d_int age_innoculations;
  array_5d_int genotypes;
  array_4d_int genotype_metadata;
  
  // misc
  std::vector<double> EIR;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher();
  
  // methods
  void simulate();
  
};