
#include "sim.Dispatcher.h"
#include "misc_v4.h"
#include "probability.h"

using namespace std;


//------------------------------------------------
// constructor
Dispatcher::Dispatcher() {
  
  // events are enacted using scheduler objects. New events (e.g. infection) are
  // generated in the current time step, and future events (e.g. transition to
  // blood-stage) are scheduled for future time steps using these objects. This
  // avoids the need to loop through every host in every time step, as we only
  // need to modify the hosts for which we have scheduled events.
  schedule_death = vector<set<int>>(max_time+1);
  schedule_Eh_to_Ih = vector<vector<pair<int, int>>>(max_time+1);
  schedule_Ih_to_Sh = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective_recovery = vector<vector<pair<int, int>>>(max_time+1);
  
  // counts of host types
  H_total = n_demes*H;
  Sh = vector<int>(n_demes, H);
  Eh = vector<int>(n_demes);
  Ih = vector<int>(n_demes);
  
  // initialise single population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single static population
  // we can simply change the "deme" attribute of a host to represent migration
  host_pop = vector<Host>(H_total);
  next_host_ID = 0;
  
  // for each deme, store the integer index of all hosts in that deme, and the
  // integer index of infective hosts only
  host_index = vector<vector<int>>(n_demes);
  host_infective_index = vector<vector<int>>(n_demes);
  for (int k=0; k<n_demes; ++k) {
    host_index[k] = seq_int(k*H, (k+1)*H-1);
    reshuffle(host_index[k]);
  }
  
  // initialise hosts
  for (int k=0; k<n_demes; ++k) {
    for (int i=0; i<H; ++i) {
      int this_host = host_index[k][i];
      
      // initialise host
      host_pop[this_host].init(next_host_ID, k);
      
      // add death_day to scheduler
      int death_day = host_pop[this_host].death_day;
      if (death_day <= max_time) {
        schedule_death[death_day].insert(this_host);
      }
    }
  }
  
  // seed initial infections
  for (int k=0; k<n_demes; ++k) {
    for (int i=0; i<seed_infections[k]; i++) {
      denovo_infection(host_index[k][i]);
    }
  }
  
  // counts of mosquito types
  M_total = sum(M_vec);
  Sv = M_vec;
  Ev = vector<int>(n_demes);
  Iv = vector<int>(n_demes);
  
  // population of mosquitoes
  mosq_pop = vector<Mosquito>(M_total);
  
  // store integer index of mosquitoes at various stages
  Sv_index = vector<vector<int>>(n_demes);
  int M_cum = 0;
  for (int k=0; k<n_demes; ++k) {
    Sv_index[k] = seq_int(M_cum, M_cum+M_vec[k]-1);
    M_cum += M_vec[k];
  }
  Ev_death = vector<vector<vector<int>>>(n_demes, vector<vector<int>>(v));
  Ev_to_Iv = vector<vector<vector<int>>>(n_demes, vector<vector<int>>(v));
  Iv_index = vector<vector<int>>(n_demes);
  
  // objects for storing results
  daily_counts = vector<vector<vector<int>>>(n_demes, vector<vector<int>>(max_time));
  
  // misc
  EIR = vector<double>(n_demes);
  
}

//------------------------------------------------
// new human infection by mosquito
void Dispatcher::new_infection(int this_host, Mosquito &mosq, int t) {
  
  // return if already at max_innoculations
  if (host_pop[this_host].get_n_innoculations() == max_innoculations) {
    return;
  }
  
  // update deme counts
  int this_deme = host_pop[this_host].deme;
  if (host_pop[this_host].get_n_asexual() == 0) {
    Sh[this_deme]--;
    Eh[this_deme]++;
  }
  
  // update host counts
  host_pop[this_host].n_latent++;
  
  // get next free innoculation slot
  int this_slot = 0;
  for (int i=0; i<max_innoculations; ++i) {
    if (!host_pop[this_host].innoc_active[i]) {
      break;
    }
    this_slot++;
  }
  if (this_slot == max_innoculations) {
    print(host_pop[this_host].get_n_innoculations(), max_innoculations);
    Rcpp::stop("could not find free innoculation slot");
  }
  
  // add new innoculation
  host_pop[this_host].innoc_active[this_slot] = true;
  host_pop[this_host].innoc_status_asexual[this_slot] = Liverstage_asexual;
  
  // draw duration of infection
  int duration_infection = sampler_duration_infection.draw() + 1;
  
  // get times of future events
  int t1 = t + u;                           // begin bloodstage
  int t2 = t + u + duration_infection;      // end bloodstage
  int t3 = t + u + g;                       // begin infective
  int t4 = t + u + g + duration_infection;  // end infective
  int this_death_day = host_pop[this_host].death_day;
  
  // schedule move to Ih
  if (t1 < this_death_day && t1 <= max_time) {
    schedule_Eh_to_Ih[t1].emplace_back(this_host, this_slot);
  }
  
  // schedule bloodstage recovery
  if (t2 < this_death_day && t2 <= max_time) {
    schedule_Ih_to_Sh[t2].emplace_back(this_host, this_slot);
  }
  
  // schedule begin infective
  if (t3 < this_death_day && t3 <= max_time) {
    schedule_infective[t3].emplace_back(this_host, this_slot);
  }
  
  // schedule end infective
  if (t4 < this_death_day && t4 <= max_time) {
    schedule_infective_recovery[t4].emplace_back(this_host, this_slot);
  }
  
}

//------------------------------------------------
// de-novo human infection, without mosquito
void Dispatcher::denovo_infection(int this_host) {
  
  // generating starting genotype in a dummy mosquito
  Mosquito dummy_mosquito;
  
  // infect human host
  new_infection(this_host, dummy_mosquito, 0);
  
}

//------------------------------------------------
// run simulation
void Dispatcher::simulate() {
  
  // define ring buffer indices
  int ringtime = 0;
  
  // loop through daily time steps
  for (int t=1; t<=max_time; t++) {
    
    // loop through demes
    for (int k=0; k<n_demes; ++k) {
      
      
      //-------- MOSQUITO EVENTS --------
      
      // update ring buffer index
      ringtime = (ringtime == v-1) ? 0 : ringtime+1;
      
      // carry out Ev death
      int Ev_death_size = Ev_death[k][ringtime].size();
      if (Ev_death_size > 0) {
        push_back_multiple(Sv_index[k], Ev_death[k][ringtime]);
        Ev_death[k][ringtime].clear();
        Ev[k] -= Ev_death_size;
        Sv[k] += Ev_death_size;
      }
      
      // move Ev into Iv
      int Ev_to_Iv_size = Ev_to_Iv[k][ringtime].size();
      if (Ev_to_Iv_size) {
        push_back_multiple(Iv_index[k], Ev_to_Iv[k][ringtime]);
        Ev_to_Iv[k][ringtime].clear();
        Ev[k] -= Ev_to_Iv_size;
        Iv[k] += Ev_to_Iv_size;
      }
      
      // draw number of new bites on infective hosts. Note, the rate is already
      // scaled by infectivity, meaning bites that do not lead to mosquito
      // infection are ignored automatically
      double rate_v_infected = a*infectivity*host_infective_index[k].size()/double(H); // rate of mosquito biting infective host and becoming infected
      double prob_v_infected_or_death = 1 - exp(-(rate_v_infected + mu));              // probability of mosquito becoming infected or dying (competing hazards)
      double prob_v_infected = rate_v_infected/(rate_v_infected + mu);                 // relative probability of mosquito becoming infected vs. dying
      int v_infected_or_death = rbinom1(Sv[k], prob_v_infected_or_death);              // number of susceptible mosquitoes becoming infected or dying
      int v_infected = rbinom1(v_infected_or_death, prob_v_infected);                  // number of susceptible mosquitoes becoming infected
      
      // loop through infective bites
      for (int i=0; i<v_infected; ++i) {
        
        // choose mosquito at random from susceptibles
        int rnd1 = sample2(0, Sv[k]-1);
        int this_mosq = Sv_index[k][rnd1];
        
        // drop from Sv
        quick_erase(Sv_index[k], rnd1);
        
        // the majority of new mosquito infections will die in lag phase. If so
        // then no need to store genotype, instead simply add to Ev_death object
        // to schedule move back to susceptible state (equivalent to death) at
        // future time point.
        int v_time_death = rgeom1(prob_v_death) + 1;
        if (v_time_death <= v) {
          
          // schedule death
          Ev_death[k][(ringtime+v_time_death) % v].push_back(this_mosq);
          
        } else {
          
          // choose host at random from infectives
          //int rnd2 = sample2(0, host_infective_index[k].size()-1);
          //int this_host = host_infective_index[k][rnd2];
          
          // TODO - copy genotypes
          
          // schedule move to Iv
          Ev_to_Iv[k][(ringtime+v) % v].push_back(this_mosq);
        }
        
      } // end loop through infective bites
      
      // update deme counts
      Sv[k] -= v_infected;
      Ev[k] += v_infected;
      
      
      //-------- NEW HUMAN EVENTS --------
      
      // get number of new infectious bites on humans
      EIR[k] = a*Iv[k]/double(H);
      double prob_h_infectious_bite = 1 - exp(-EIR[k]);            // probability of new infectious bite per host
      int h_infectious_bite = rbinom1(H, prob_h_infectious_bite);  // total number of new infectious bites
      
      // apply new infectious bites
      for (int i=0; i<h_infectious_bite; i++) {
        
        // choose host at random
        int rnd1 = sample2(0, H-1);
        int this_host = host_index[k][rnd1];
        
        // determine whether infectious bite is successful
        if (rbernoulli1(host_pop[this_host].get_prob_infection())) {
          
          // choose mosquito at random and carry out infection
          int rnd2 = sample2(0, Iv[k]-1);
          int this_mosq = Iv_index[k][rnd2];
          new_infection(this_host, mosq_pop[this_mosq], t);
        }
        
        // update prob_infection_index irrespective of whether infection takes hold
        host_pop[this_host].update_prob_infection();
        
      }
      
    } // end loop through demes
    
    
    //-------- SCHEDULED HUMAN EVENTS --------
    
    // scheduled deaths
    for (auto it = schedule_death[t].begin(); it != schedule_death[t].end(); ++it) {
      int this_host = *it;
      int this_deme = host_pop[this_host].deme;
      
      // update deme counts
      if (host_pop[this_host].n_infected > 0) {
        Ih[this_deme]--;
        Sh[this_deme]++;
      } else if (host_pop[this_host].n_latent > 0) {
        Eh[this_deme]--;
        Sh[this_deme]++;
      }
      
      // drop from infective list if necessary
      if (host_pop[this_host].n_infective > 0) {
        erase_remove(host_infective_index[this_deme], this_host);
      }
      
      // reset host
      host_pop[this_host].reset(next_host_ID, t);
      
      // add new death_day to scheduler
      int death_day = host_pop[this_host].death_day;
      if (death_day <= max_time) {
        schedule_death[death_day].insert(this_host);
      }
    }
    
    // scheduled Eh to Ih
    for (auto it = schedule_Eh_to_Ih[t].begin(); it != schedule_Eh_to_Ih[t].end(); ++it) {
      int this_host = it->first;
      int this_slot = it->second;
      int this_deme = host_pop[this_host].deme;
      
      // update deme counts
      if (host_pop[this_host].n_infected == 0) {
        Eh[this_deme]--;
        Ih[this_deme]++;
      }
      
      // update host status
      host_pop[this_host].innoc_status_asexual[this_slot] = Bloodstage_asexual;
      
      // update host counts
      host_pop[this_host].n_latent--;
      host_pop[this_host].n_infected++;
    }
    
    // scheduled Ih to Sh
    for (auto it = schedule_Ih_to_Sh[t].begin(); it != schedule_Ih_to_Sh[t].end(); ++it) {
      int this_host = it->first;
      int this_slot = it->second;
      int this_deme = host_pop[this_host].deme;
      
      // update deme counts
      if (host_pop[this_host].n_infected == 1) {
        Ih[this_deme]--;
        if (host_pop[this_host].n_latent == 0) {
          Sh[this_deme]++;
        } else {
          Eh[this_deme]++;
        }
      }
      
      // update host status
      host_pop[this_host].innoc_status_asexual[this_slot] = Inactive_asexual;
      
      // update host counts
      host_pop[this_host].n_infected--;
    }
    
    // scheduled become infective
    for (auto it = schedule_infective[t].begin(); it != schedule_infective[t].end(); ++it) {
      int this_host = it->first;
      int this_slot = it->second;
      int this_deme = host_pop[this_host].deme;
      
      // update host status
      host_pop[this_host].innoc_status_sexual[this_slot] = Active_sexual;
      
      // update host counts
      host_pop[this_host].n_infective++;
      
      // if newly infective then add to infectives list
      if (host_pop[this_host].n_infective == 1) {
        host_infective_index[this_deme].push_back(this_host);
      }
    }
    
    // scheduled infective recovery
    for (auto it = schedule_infective_recovery[t].begin(); it != schedule_infective_recovery[t].end(); ++it) {
      int this_host = it->first;
      int this_deme = host_pop[this_host].deme;
      int this_slot = it->second;
      
      // update host status
      host_pop[this_host].innoc_status_sexual[this_slot] = Inactive_sexual;
      host_pop[this_host].innoc_active[this_slot] = false;
      
      // update host counts
      host_pop[this_host].n_infective--;
      
      // if no longer infective then drop from infectives list
      if (host_pop[this_host].n_infective == 0) {
        erase_remove(host_infective_index[this_deme], this_host);
      }
    }
    
    
    //-------- STORE RESULTS --------
    
    // store daily counts
    for (int k=0; k<n_demes; ++k) {
      vector<int> tmp1 = {Sh[k], Eh[k], Ih[k]};
      push_back_multiple(daily_counts[k][t-1], tmp1);
    }
    
  } // end time loop
  
}
