
#include "main.h"
#include "probability.h"
#include "sim.Parameters.h"
#include "sim.Dispatcher.h"

#include <chrono>
#include <vector>

using namespace std;

//------------------------------------------------
// compute map and run permutation test
// [[Rcpp::export]]
Rcpp::List calc_hex_values_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress)
{
	int hex, ell, perm, this_ellipse, Nint;
	
	// start timer
	//chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
	//Convert Rcpp arguments to native c++ arguments ------------------------------------------------------------------------------------------------
  
	// load data and parameters
	print("Loading data and parameters");
	vector<double> edge_value = rcpp_to_vector_double(args["edge_value"]);                //Values of edges
	vector<vector<int>> edge_group_list = rcpp_to_matrix_int(args["edge_group_list"]);    //List of which edges belong to each group
	vector<int> Nintersections = rcpp_to_vector_int(args["Nintersections"]);			  //Number of ellipses intersecting hex
	vector<vector<int>> intersections = rcpp_to_matrix_int(args["intersections"]);	  //List of ellipses intersecting each hex
	vector<double> area_inv = rcpp_to_vector_double(args["area_inv"]);					  //Inverse of area of each ellipse
	vector<double> inv_hex_weights = rcpp_to_vector_double(args["inv_hex_weights"]);		  //Inverse sum of weights of hex
	int Nperms = rcpp_to_int(args["n_perms"]);                                            //Number of permutations to run (if 0, no permutation)
	int min_intersections = rcpp_to_int(args["min_intersections"]);                       //Minimum number of ellipses required to insect a hex
	bool report_progress = rcpp_to_bool(args["report_progress"]);                         //Whether to update progress bar
	Rcpp::Function update_progress = args_functions["update_progress"];                   //R function for updating progress bar  
  
	//Derived values/fixed constants ------------------------------------------------------------------------------------------------
  
	int Nells = edge_value.size();                  //Number of ellipses (edges)
	int Nhex = inv_hex_weights.size();                     //Number of hex cells
    
	//chrono_timer(t0);
  
	//Create map by summation of ellipses intersecting each hex ------------------------------------------------------------------------------------------------

	print("Computing map");
  
    // hex properties
	vector<double> hex_values(Nhex, 0.0);		//Final value of hex
	
	// loop through hexes
	for (hex = 0; hex < Nhex; hex++)
	{
		Nint = Nintersections[hex];
		// test every ellipse for intersection with this hex
		for (ell = 0; ell < Nint; ell++)
		{		
			// add to hex value
			this_ellipse = intersections[hex][ell];
			hex_values[hex] += edge_value[this_ellipse]*area_inv[this_ellipse];
		}
		// divide hex value by weight
		if (Nint >= min_intersections)
		{
			hex_values[hex] *= inv_hex_weights[hex]; 
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
				Nint = Nintersections[hex];
				if (Nint >= min_intersections)
				{
					// compute hex value
					for (ell = 0; ell < Nint; ell++)
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
                              Rcpp::Named("hex_ranks") = hex_ranks);
}