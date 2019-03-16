
#include "sim.Host.h"
#include "probability.h"
#include "misc_v4.h"

using namespace std;

//------------------------------------------------
// get total number of innoculations
int Host::get_n_innoculations() {
  return sum_bool(innoc_active);
}

//------------------------------------------------
// get total number of asexual stage innoculations
int Host::get_n_asexual() {
  return n_latent + n_infected;
}

//------------------------------------------------
// get current probability of infection
double Host::get_prob_infection() {
  return prob_infection[prob_infection_index];
}

//------------------------------------------------
// initialise host
void Host::init(int &ID, int deme) {
  
  // unique ID and record of current deme
  this->ID = ID++;
  this->deme = deme;
  
  // indices relating to global distributions
  prob_infection_index = 0;
  
  // draw age from demography distribution
  int age_years = sampler_age_stable.draw() - 1;
  int extra_days = sample2(1, 365);
  int age_days = age_years*365 + extra_days;
  
  // draw duration of life from demography distribution looking forward from
  // current age. This is tricky, as we must account for the fact that if we
  // are already part way into an age group we have a reduced probability of
  // dying within that age group.
  int life_days = 0;
  double prop_year_remaining = 1 - extra_days/365.0;
  double prob_die_this_year = life_table[age_years]*prop_year_remaining;
  if (rbernoulli1(prob_die_this_year) || age_years == (n_age-1)) {
    life_days = age_years*365 + sample2(extra_days, 365);
  } else {
    for (int i=(age_years+1); i<n_age; ++i) {
      if (rbernoulli1(life_table[i])) {
        life_days = i*365 + sample2(1, 365);
        break;
      }
    }
  }
  
  // convert to final birth and death days
  birth_day = -age_days;
  death_day = life_days - age_days;
  if (death_day == 0) {  // in the unlikely even that due to die on day zero, delay death by one day
    death_day++;
  }
  
  // initialise innoculation objects
  innoc_active = vector<bool>(max_innoculations, false);
  innoc_status_asexual = vector<Status_asexual>(max_innoculations, Inactive_asexual);
  innoc_status_sexual = vector<Status_sexual>(max_innoculations, Inactive_sexual);
  
  // initialise innoculation counts
  n_latent = 0;
  n_infected = 0;
  n_infective = 0;
  
}

//------------------------------------------------
// reset host, e.g. upon death
void Host::reset(int &ID, int birth_day) {
  
  // new unique ID
  this->ID = ID++;
  
  // reset indices relating to global distributions
  prob_infection_index = 0;
  
  // date of birth
  this->birth_day = birth_day;
  
  // draw life duration from demography distribution
  int life_years = sampler_age_death.draw() - 1;
  int life_days = life_years*365 + sample2(1, 365);
  death_day = birth_day + life_days;
  
  // reset innoculation objects
  fill(innoc_active.begin(), innoc_active.end(), false);
  fill(innoc_status_asexual.begin(), innoc_status_asexual.end(), Inactive_asexual);
  fill(innoc_status_sexual.begin(), innoc_status_sexual.end(), Inactive_sexual);
  
  // reset innoculation counts
  n_latent = 0;
  n_infected = 0;
  n_infective = 0;
  
}

//------------------------------------------------
// update probabilty of infection
void Host::update_prob_infection() {
  if (prob_infection_index < (n_prob_infection-1)) {
    prob_infection_index++;
  }
}
