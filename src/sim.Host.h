
#pragma once

#include <vector>

#include "sim.Parameters.h"

//------------------------------------------------
// enumerate possible asexual and sexual innoculation status
enum Status_asexual {Inactive_asexual, Liverstage_asexual, Bloodstage_asexual};
enum Status_sexual {Inactive_sexual, Active_sexual};

//------------------------------------------------
// class defining host
class Host : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // unique ID and record of current deme
  int ID;
  int deme;
  
  // indices relating to global distributions
  int prob_infection_index;
  
  // dates of birth and death
  int birth_day;
  int death_day;
  
  // innoculation objects
  std::vector<bool> innoc_active;
  std::vector<Status_asexual> innoc_status_asexual;
  std::vector<Status_sexual> innoc_status_sexual;
  
  // innoculation counts
  int n_latent;
  int n_infected;
  int n_infective;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  
  // getters and setters
  int get_n_innoculations();
  int get_n_asexual();
  double get_prob_infection();
  
  // other methods
  void init(int &ID, int deme);
  void reset(int &ID, int birth_day);
  void update_prob_infection();
  
};
