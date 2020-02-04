
#include "main.h"
#include "probability.h"
#include "sim.Parameters.h"
#include "sim.Dispatcher.h"

#include <chrono>
#include <vector>

using namespace std;

//------------------------------------------------
// check if value is within ellipse
// x and y are the coordinates of the query point. f1 and f2 are the two foci of
// the ellipse. a is the semi-major axis of the ellipse.
bool ellipse_check(double x, double y,
                   double xf1, double yf1, double xf2, double yf2,
                   double a) {
  bool ret = (dist_euclid_2d(x, y, xf1, yf1) + dist_euclid_2d(x, y, xf2, yf2) <= 2*a);
  return ret;
}

//------------------------------------------------
// assign edges to hexes based on intersection
Rcpp::List assign_map_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  print("Assigning edges to hexes");
  
  // load data and parameters
  vector<double> node_long = rcpp_to_vector_double(args["node_long"]);                  //Longitude of data nodes
  vector<double> node_lat = rcpp_to_vector_double(args["node_lat"]);                    //Latitude of data nodes
  vector<double> hex_long = rcpp_to_vector_double(args["hex_long"]);                    //Longitude of hex cells
  vector<double> hex_lat = rcpp_to_vector_double(args["hex_lat"]);                      //Latitude of hex cells
  double eccentricity = rcpp_to_double(args["eccentricity"]);                           //Eccentricity of ellipses (see help for details)
  bool report_progress = rcpp_to_bool(args["report_progress"]);                         //Whether to update progress bar
  Rcpp::Function update_progress = args_functions["update_progress"];                   //R function for updating progress bar
  
  // get basic properties
  int n_node = node_long.size();
  int n_hex = hex_long.size();
  double inv_eccentricity = 1.0 / eccentricity;
  
  // store list of which edges intersect each hex
  vector<vector<int>> hex_edges(n_hex);
  
  // loop through hexes
  for (int hex = 0; hex < n_hex; ++hex) {
    
    // report progress
    if (report_progress) {
      update_progress(args_progress, "pb", hex+1, n_hex);
    }
    
    // loop through pairwise nodes
    int i = 0;
    for (int node1 = 0; node1 < (n_node-1); ++node1) {
      for (int node2 = (node1+1); node2 < n_node; ++node2) {
        i++;
        
        // determine whether ellipse intersects centroid of this hex
        double dist = dist_euclid_2d(node_long[node1], node_lat[node1], node_long[node2], node_lat[node2]);
        double linear_eccentricity = 0.5 * dist;
        double semi_major = linear_eccentricity*inv_eccentricity;
        bool intersects = ellipse_check(hex_long[hex], hex_lat[hex],
                                        node_long[node1], node_lat[node1],
                                        node_long[node2], node_lat[node2],
                                        semi_major);
        
        // push back edge index if intersects
        if (intersects) {
          hex_edges[hex].push_back(i);
        }
      }
    }
  }
  
  // return list
  return Rcpp::List::create(Rcpp::Named("hex_edges") = hex_edges);
}

//------------------------------------------------
// compute map and run permutation test
Rcpp::List rmapi_analysis_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
	
	// start timer
	//chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
	// ------------------------------------------------------------------------------------------------
	// Convert Rcpp arguments to native c++ arguments
  
  vector<int> perm_group = rcpp_to_vector_int(args["perm_group"]);              // The permutation group of each observed edge
	vector<vector<double>> perm_list = rcpp_to_matrix_double(args["perm_list"]);  // The set of values in each permutation group
	vector<vector<int>> hex_edges = rcpp_to_matrix_int(args["hex_edges"]);        // The edges that intersect each hex
	int n_perms = rcpp_to_int(args["n_perms"]);                                   // Number of permutations to run
	bool report_progress = rcpp_to_bool(args["report_progress"]);                 // Whether to update progress bar
	Rcpp::Function update_progress = args_functions["update_progress"];           // R function for updating progress bar
  
  
  // ------------------------------------------------------------------------------------------------
  // Derived values
  
  // numbers of edges and hexes
  int n_edge = int(perm_group.size());
  int n_hex = int(hex_edges.size());
  int n_breaks = int(perm_list.size());
  
  // size of list elements
  vector<int> perm_list_size(n_breaks);
  for (unsigned int i = 0; i < n_breaks; ++i) {
    perm_list_size[i] = perm_list[i].size();
  }
  vector<int> hex_edges_size(n_hex);
  for (unsigned int i = 0; i < n_hex; ++i) {
    hex_edges_size[i] = hex_edges[i].size();
  }
  
  // objects for storing results
  vector<double> edge_values(n_edge);
  vector<double> hex_values(n_hex);
  vector<double> ret_sum(n_hex);
  vector<vector<double>> ret_sum_sq(n_hex, vector<double>(n_hex));
  
  
  // ------------------------------------------------------------------------------------------------
  // Carry out permutation test
  
  print("Carrying out permutation test");
  
  // loop through permutations
  for (unsigned int perm = 0; perm < n_perms; ++perm) {
    
    // report progress
    if (report_progress) {
      update_progress(args_progress, "pb", perm+1, n_perms);
    }
    
    // resample edge values
    for (unsigned int i = 0; i < n_edge; ++i) {
      int pg = perm_group[i] - 1;
      int rnd_index = sample2(1, perm_list_size[pg]) - 1;
      edge_values[i] = perm_list[pg][rnd_index];
    }
    
    // recalculate hex values
    fill(hex_values.begin(), hex_values.end(), 0.0);
    for (unsigned int h = 0; h < n_hex; ++h) {
      
      // recalculate hex value
      for (unsigned int i = 0; i < hex_edges_size[h]; ++i) {
        hex_values[h] += edge_values[hex_edges[h][i] - 1];
      }
      hex_values[h] /= double(hex_edges_size[h]);
      
      // update running sums
      ret_sum[h] += hex_values[h];
      for (unsigned int j = 0; j < (h+1); ++j) {
        ret_sum_sq[h][j] += hex_values[j]*hex_values[h];
      }
      
    }
    
  }  // end loop over permutations
  
  // return list
  return Rcpp::List::create(Rcpp::Named("ret_sum") = ret_sum,
                            Rcpp::Named("ret_sum_sq") = ret_sum_sq);
}

//------------------------------------------------
// simulate from simple individual-based model
Rcpp::List sim_falciparum_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract model parameters into separate class
  Parameters parameters(args);
  //parameters.print_summary();
  
  // R functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // create simulation dispatcher object
  Dispatcher dispatcher(parameters, update_progress, args_progress);
  
  // carry out simulation
  dispatcher.simulate();
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  print("simulation completed in", time_span.count(), "seconds\n");
  
  // return list
  return Rcpp::List::create(Rcpp::Named("daily_values") = dispatcher.daily_values.arr,
                            Rcpp::Named("genotypes") = dispatcher.genotypes.arr,
                            Rcpp::Named("indlevel_data") = dispatcher.indlevel_data.arr);
}