
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
// the ellipse. a is the linear eccentricity of the ellipse.
bool ellipse_check(const double x, const double y, const double xf1, const double yf1, const double xf2, const double yf2, const double a) {
  return dist_euclid_2d(x, y, xf1, yf1) + dist_euclid_2d(x, y, xf2, yf2) <= 2*a;
}

//------------------------------------------------
// compute map and run permutation test
// [[Rcpp::export]]
Rcpp::List rmapi_analysis_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress)
{
	int hex, node1, node2, ell, perm, this_ellipse;
	double dist, linear_eccentricity, semi_minor;
	
	// start timer
	//chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
	//Convert Rcpp arguments to native c++ arguments ------------------------------------------------------------------------------------------------
  
	// load data and parameters
	print("Loading data and parameters");
	vector<double> node_long = rcpp_to_vector_double(args["node_long"]);                  //Longitude of data nodes
	vector<double> node_lat = rcpp_to_vector_double(args["node_lat"]);                    //Latitude of data nodes
	vector<double> edge_value = rcpp_to_vector_double(args["edge_value"]);                //Values of edges
	vector<double> edge_value_pred = rcpp_to_vector_double(args["edge_value_pred"]);      //Values of edges predicted from model fit
	vector<vector<int>> edge_group_list = rcpp_to_matrix_int(args["edge_group_list"]);    //List of which edges belong to each group
	vector<double> hex_long = rcpp_to_vector_double(args["hex_long"]);                    //Longitude of hex cells
	vector<double> hex_lat = rcpp_to_vector_double(args["hex_lat"]);                      //Latitude of hex cells
	int Nperms = rcpp_to_int(args["n_perms"]);                                            //Number of permutations to run (if 0, no permutation)
	int min_intersections = rcpp_to_int(args["min_intersections"]);                       //Minimum number of ellipses required to insect a hex
	double eccentricity = rcpp_to_double(args["eccentricity"]);                           //Eccentricity of ellipses (see help for details)
	double inv_eccentricity = 1.0 / eccentricity;
	bool report_progress = rcpp_to_bool(args["report_progress"]);                         //Whether to update progress bar
	Rcpp::Function update_progress = args_functions["update_progress"];                   //R function for updating progress bar
  
  
	//Derived values/fixed constants ------------------------------------------------------------------------------------------------
  
	int Nnodes = node_long.size();                  //Number of nodes
	int Nells = edge_value.size();                  //Number of ellipses (edges)
	int Nhex = hex_long.size();                     //Number of hex cells
  
	//Create ellipses ------------------------------------------------------------------------------------------------
  
	print("Setting up ellipses");
	
	// ellipse properties
	vector<double> semi_major(Nells, 0.0);            //Semi-major axis of ellipse
	vector<double> xfocus1(Nells, 0.0);               //x-coordinate of focus 1 of each ellipse
	vector<double> yfocus1(Nells, 0.0);               //y-coordinate of focus 1 of each ellipse
	vector<double> xfocus2(Nells, 0.0);               //x-coordinate of focus 2 of each ellipse
	vector<double> yfocus2(Nells, 0.0);               //y-coordinate of focus 2 of each ellipse
	vector<double> area_inv(Nells, 0.0);              //Inverse of area of each ellipse
	vector<double> edge_weighted(Nells, 0.0);         //Weighted metric value of each ellipse
  
    // loop through ellipses
	ell = 0;
	for (node1 = 0; node1 < Nnodes; node1++)
	{
		for (node2 = node1 + 1; node2 < Nnodes; node2++)
		{
			// store foci and centre of ellipse
			xfocus1[ell] = node_long[node1];
			yfocus1[ell] = node_lat[node1];
			xfocus2[ell] = node_long[node2];
			yfocus2[ell] = node_lat[node2];
      
			// store long radius and inverse area of ellipse
			dist = dist_euclid_2d(node_long[node1], node_lat[node1], node_long[node2], node_lat[node2]);
			linear_eccentricity = 0.5 * dist;
			semi_major[ell] = linear_eccentricity*inv_eccentricity;
			semi_minor = sqrt(sq(semi_major[ell]) - sq(linear_eccentricity));
			area_inv[ell] = 1.0 / (M_PI * semi_major[ell] * semi_minor);
      
			// store weighted value of ellipse
			edge_weighted[ell] = edge_value[ell] * area_inv[ell];
			
			ell++;
		}
	}
  
	//chrono_timer(t0);
  
	//Create map by summation of ellipses intersecting each hex ------------------------------------------------------------------------------------------------

	print("Computing map");
  
    // hex properties
	vector<double> hex_values(Nhex, 0.0);		//Final value of hex
	vector<double> hex_weights(Nhex, 0.0);		//Sum of weights of hex
	vector<double> inv_hex_weights(Nhex, 0.0);	//Inverse sum of weights of hex
	vector<int> Nintersections(Nhex, 0);		//Number of ellipses intersecting hex
	vector<vector<int>> intersections(Nhex);	//List of ellipses intersecting each hex
	
	// loop through hexs
	for (hex = 0; hex < Nhex; hex++)
	{
		// test every ellipse for intersection with this hex
		for (ell = 0; ell < Nells; ell++)
		{
			if (ellipse_check(hex_long[hex], hex_lat[hex], xfocus1[ell], yfocus1[ell], xfocus2[ell], yfocus2[ell], semi_major[ell]))
			{
				
				// add to hex value and weights
				hex_values[hex] += edge_value[ell]*area_inv[ell];
				hex_weights[hex] += area_inv[ell];
				
				// store this intersection
				intersections[hex].push_back(ell);
				Nintersections[hex]++;
			}
		}
		// divide hex value by weight
		if (Nintersections[hex] >= min_intersections) 
		{
			inv_hex_weights[hex] = 1.0 / hex_weights[hex];
			hex_values[hex] *= inv_hex_weights[hex]; 
		}
		else
		{
			inv_hex_weights[hex] = 1.0;
		}
	}
  
	//chrono_timer(t0);
  
	//Permute data to check statistical significance-------------------------------------------------------------------------------------------------------------------------------------------
  
	// hex rankings calculated by permutation test
	vector<int> hex_ranks(Nhex);
	
	// skip over if Nperms == 0
	if (Nperms > 0)
	{
		print("Running permutation test");
		
		// permuted hex values
		vector<double> hex_values_perm(Nhex, 0.0);
    
		// create vector for indexing random permutations
		vector<int> perm_vec = seq_int(0, Nells - 1);
		
		// loop through permutations
		for (perm = 0; perm < Nperms; perm++)
		{
			// report progress
			if (report_progress) { update_progress(args_progress, "pb", perm, Nperms-1); }
			
			// new permutation
			reshuffle_group(perm_vec, edge_group_list);
		  
			for (hex = 0; hex < Nhex; hex++)
			{
				hex_values_perm[hex] = 0;
				if (Nintersections[hex] >= min_intersections)
				{
					// compute hex value
					for (ell = 0; ell < Nintersections[hex]; ell++)
					{
						this_ellipse = intersections[hex][ell];
						hex_values_perm[hex] += edge_value[perm_vec[this_ellipse]] * area_inv[this_ellipse];
					}
					hex_values_perm[hex] *= inv_hex_weights[hex];
				  
					// compare with unpermuted value
					if (hex_values_perm[hex] < hex_values[hex]) { hex_ranks[hex]++; }
				}
			}
			
		} // end loop over Nperms
    
		//chrono_timer(t0);
    
	} // end Nperms > 0 condition
	
	// return list
	return Rcpp::List::create(Rcpp::Named("hex_values") = hex_values,
                           Rcpp::Named("hex_weights") = hex_weights,
                           Rcpp::Named("hex_ranks") = hex_ranks,
                           Rcpp::Named("Nintersections") = Nintersections);
}

//------------------------------------------------
// simulate from simple individual-based model
// [[Rcpp::export]]
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