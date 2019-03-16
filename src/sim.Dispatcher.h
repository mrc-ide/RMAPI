
#pragma once

#include <set>

#include "sim.Parameters.h"
#include "sim.Host.h"
#include "sim.Mosquito.h"

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
  std::vector<std::vector<int>> host_index;
  std::vector<std::vector<int>> host_infective_index;
  
  // counts of mosquito types
  int M_total;
  std::vector<int> Sv;
  std::vector<int> Ev;
  std::vector<int> Iv;
  
  // population of mosquitoes
  std::vector<Mosquito> mosq_pop;
  
  // store integer index of mosquitoes at various stages
  std::vector<std::vector<int>> Sv_index;
  std::vector<std::vector<std::vector<int>>> Ev_death;
  std::vector<std::vector<std::vector<int>>> Ev_to_Iv;
  std::vector<std::vector<int>> Iv_index;
  
  // objects for storing results
  std::vector<std::vector<std::vector<int>>> daily_counts;
  
  // misc
  std::vector<double> EIR;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher();
  
  // methods
  void new_infection(int this_host, Mosquito &mosq, int t);
  void denovo_infection(int this_host);
  void simulate();
  
};