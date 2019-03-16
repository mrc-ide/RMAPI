
#include "sim.Parameters.h"

#include "misc_v4.h"

using namespace std;

//------------------------------------------------
// declare static member variables

// epidemiological parameters
double Parameters::a;
double Parameters::mu;
int Parameters::u;
int Parameters::v;
int Parameters::g;
vector<double> Parameters::prob_infection;
int Parameters::n_prob_infection;
vector<double> Parameters::duration_infection;
int Parameters::n_duration_infection;
double Parameters::infectivity;
int Parameters::max_innoculations;

// deme parameters
int Parameters::H;
vector<int> Parameters::seed_infections;
vector<int> Parameters::M_vec;
int Parameters::n_demes;

// demography
vector<double> Parameters::life_table;
vector<double> Parameters::age_death;
vector<double> Parameters::age_stable;
int Parameters::n_age;

// objects for sampling from probability distributions
Sampler Parameters::sampler_age_stable;
Sampler Parameters::sampler_age_death;
Sampler Parameters::sampler_duration_infection;

// migration
// TODO

// run parameters
vector<int> Parameters::time_out;
int Parameters::n_time_out;
int Parameters::max_time;

// misc parameters
double Parameters::prob_v_death;

//------------------------------------------------
// constructor
Parameters::Parameters(const Rcpp::List &args) {
  
  // epidemiological parameters
  a = rcpp_to_double(args["a"]);
  mu = rcpp_to_double(args["mu"]);
  u = rcpp_to_int(args["u"]);
  v = rcpp_to_int(args["v"]);
  g = rcpp_to_int(args["g"]);
  prob_infection = rcpp_to_vector_double(args["prob_infection"]);
  n_prob_infection = int(prob_infection.size());
  duration_infection = rcpp_to_vector_double(args["duration_infection"]);
  n_duration_infection = int(duration_infection.size());
  infectivity = rcpp_to_double(args["infectivity"]);
  max_innoculations = rcpp_to_int(args["max_innoculations"]);
  
  // deme parameters
  H = rcpp_to_int(args["H"]);
  seed_infections = rcpp_to_vector_int(args["seed_infections"]);
  M_vec = rcpp_to_vector_int(args["M"]);
  n_demes = int(M_vec.size());
  
  // demography
  life_table = rcpp_to_vector_double(args["life_table"]);
  age_death = rcpp_to_vector_double(args["age_death"]);
  age_stable = rcpp_to_vector_double(args["age_stable"]);
  n_age = int(age_stable.size());
  
  // objects for sampling from probability distributions
  sampler_age_stable = Sampler(age_stable, 1000);
  sampler_age_death = Sampler(age_death, 1000);
  sampler_duration_infection = Sampler(duration_infection, 1000);
  
  // migration
  // (TODO)
  
  // run parameters
  time_out = rcpp_to_vector_int(args["time_out"]);
  n_time_out = int(time_out.size());
  max_time = max(time_out);
  
  // misc parameters
  prob_v_death = 1 - exp(-mu);  // daily probability of mosquito death
  
}

//------------------------------------------------
// print summary
void Parameters::print_summary() {
  
  // epidemiological parameters
  print("-- epidemiological parameters --");
  print("a: ", a);
  print("mu: ", mu);
  print("u: ", u);
  print("v: ", v);
  print("g: ", g);
  print("prob_infection: ");
  print_vector(prob_infection);
  print("duration_infection: ");
  print_vector(duration_infection);
  print("infectivity: ", infectivity);
  print("max_innoculations: ", max_innoculations);
  print("");
  
  // deme parameters
  print("-- deme parameters --");
  print("H: ", H);
  print("seed_infections: ");
  print_vector(seed_infections);
  print("M_vec: ");
  print_vector(M_vec);
  print("n_demes: ", n_demes);
  print("");
  
  // demography
  print("-- demography --");
  print("life_table: ");
  print_vector(life_table);
  print("age_death: ");
  print_vector(age_death);
  print("age_stable: ");
  print_vector(age_stable);
  print("n_age: ", n_age);
  print("");
  
  // migration
  print("-- migration --");
  print("migration: (TODO)");
  print("");
  
  // run parameters
  print("-- run parameters --");
  print("t_out: ");
  print_vector(time_out);
  print("n_time_out: ", n_time_out);
  print("max_time: ", max_time);
  print("");
  
  // misc parameters
  print("-- misc parameters --");
  print("prob_v_death: ", prob_v_death);
  print("");
}