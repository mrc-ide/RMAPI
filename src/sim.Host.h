
#pragma once

#include "sim.Parameters.h"
#include "sim.Mosquito.h"
#include "array.h"

#include <vector>

//------------------------------------------------
// enumerate possible asexual and sexual innoculation status
enum Status_asexual {Inactive_asexual, Liverstage_asexual, Bloodstage_asexual};
enum Status_sexual {Inactive_sexual, Active_sexual};

//------------------------------------------------
// class defining host
class Host : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // identifiers
  int index; // where in the population (vector of hosts) this host resides
  int ID;    // unique ID, incremented upon death
  int deme;  // deme in which this host resides
  
  // pointers to external objects
  std::vector<int>* Sh_ptr;
  std::vector<int>* Eh_ptr;
  std::vector<int>* Ih_ptr;
  std::vector<std::set<int>>* schedule_death_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_Eh_to_Ih_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_Ih_to_Sh_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_infective_ptr;
  std::vector<std::vector<std::pair<int, int>>>* schedule_infective_recovery_ptr;
  std::vector<std::vector<int>>* host_infective_index_ptr;
  
  // indices relating to global distributions
  int prob_infection_index;
  
  // dates of birth and death
  int birth_day;
  int death_day;
  
  // innoculation objects
  std::vector<bool> innoc_active;
  std::vector<Status_asexual> innoc_status_asexual;
  std::vector<Status_sexual> innoc_status_sexual;
  
  // haplotypes
  std::vector<std::vector<std::vector<int>>> haplotypes;
  std::vector<int> n_haplotypes;
  int n_haplotypes_total;
  
  // innoculation counts
  int n_latent;
  int n_infected;
  int n_infective;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  
  // other methods
  void init(int index, int &ID, int deme,
            std::vector<int> &Sh, std::vector<int> &Eh, std::vector<int> &Ih,
            std::vector<std::vector<int>> &host_infective_index,
            std::vector<std::set<int>> &schedule_death,
            std::vector<std::vector<std::pair<int, int>>> &schedule_Eh_to_Ih,
            std::vector<std::vector<std::pair<int, int>>> &schedule_Ih_to_Sh,
            std::vector<std::vector<std::pair<int, int>>> &schedule_infective,
            std::vector<std::vector<std::pair<int, int>>> &schedule_infective_recovery);
  void death(int &ID, int birth_day);
  void new_infection(Mosquito &mosq, int t);
  void denovo_infection();
  void Eh_to_Ih(int this_slot);
  void Ih_to_Sh(int this_slot);
  void begin_infective(int this_slot);
  void end_infective(int this_slot);
  void update_prob_infection();
  
  // getters and setters
  int get_n_innoculations();
  int get_n_asexual();
  double get_prob_infection();
  
};
