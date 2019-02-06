
#include "misc_v3.h"
#include "probability.h"

#include <Rcpp.h>
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
// [[Rcpp::export]]
Rcpp::List run_sims_cpp(Rcpp::List args)
{
	
	// start timer
	chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
	//Convert Rcpp arguments to native c++ arguments ------------------------------------------------------------------------------------------------
  
	// load data and parameters
	print("Loading data and parameters");
	vector<double> node_long = rcpp_to_vector_double(args["node_long"]);                  //Longitude of data nodes
	vector<double> node_lat = rcpp_to_vector_double(args["node_lat"]);                    //Latitude of data nodes
	vector<double> edge_value = rcpp_to_vector_double(args["edge_value"]);                //Values of edges
	//vector<int> edge_group = rcpp_to_vector_int(args["edge_group"]);                      //Spatial grouping of edges
	vector<vector<int>> edge_group_list = rcpp_to_matrix_int(args["edge_group_list"]);    //List of which edges belong to each group
	vector<double> hex_long = rcpp_to_vector_double(args["hex_long"]);                    //Longitude of hex cells
	vector<double> hex_lat = rcpp_to_vector_double(args["hex_lat"]);                      //Latitude of hex cells
	int Nperms = rcpp_to_int(args["n_perms"]);                                            //Number of permutations to run (if 0, no permutation)
	int min_intersections = rcpp_to_int(args["min_intersections"]);                       //Minimum number of ellipses required to insect a hex
	double eccentricity = rcpp_to_double(args["eccentricity"]);                           //Eccentricity of ellipses (see help for details)
	bool divide_weights = rcpp_to_bool(args["divide_weights"]);
  
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
	int ell = 0;
	for (int node1 = 0; node1 < Nnodes; node1++)
	{
		for (int node2 = node1 + 1; node2 < Nnodes; node2++)
		{
			// store foci and centre of ellipse
			xfocus1[ell] = node_long[node1];
			yfocus1[ell] = node_lat[node1];
			xfocus2[ell] = node_long[node2];
			yfocus2[ell] = node_lat[node2];
      
			// store long radius and inverse area of ellipse
			double dist = dist_euclid_2d(node_long[node1], node_lat[node1], node_long[node2], node_lat[node2]);
			double linear_eccentricity = 0.5 * dist;
			semi_major[ell] = linear_eccentricity / eccentricity;
			double semi_minor = sqrt(sq(semi_major[ell]) - sq(linear_eccentricity));
			area_inv[ell] = 1.0 / (M_PI * semi_major[ell] * semi_minor);
      
			// store weighted value of ellipse
			edge_weighted[ell] = edge_value[ell] * area_inv[ell];
			
			ell++;
		}
	}
  
	chrono_timer(t0);
  
	//Create map by summation of ellipses intersecting each hex ------------------------------------------------------------------------------------------------

	print("Computing map");
  
  // hex properties
	vector<double> hex_values(Nhex, 0.0);		  //Final value of hex
	vector<double> hex_weights(Nhex, 0.0);		//Sum of weights of hex
	vector<int> Nintersections(Nhex, 0);		  //Number of ellipses intersecting hex
	vector<vector<int>> intersections(Nhex);	//List of ellipses intersecting each hex
	
	// loop through hexs
	for (int hex = 0; hex < Nhex; hex++)
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
		  if (divide_weights) {
		    hex_values[hex] /= hex_weights[hex];
		  } else {
		    hex_values[hex] /= double(Nintersections[hex]);
		  }
		}
	}
  
	chrono_timer(t0);
  
	//Permute data to check statistical significance-------------------------------------------------------------------------------------------------------------------------------------------
  
  // empirical p-value calculated by permutation test
	vector<int> empirical_p(Nhex);
	
	// skip over if Nperms == 0
	if (Nperms > 0)
	{
    
		print("Running permutation test");
    
    // permuted hex values
    vector<double> hex_values_perm(Nhex, 0.0);
    
		// create vector for indexing random permutations
		vector<int> perm_vec = seq_int(0, Nells - 1);
		
		// loop through permutations
		for (int perm = 0; perm < Nperms; perm++)
		{
		  // new permutation
			reshuffle_group(perm_vec, edge_group_list);
		  
			for (int hex = 0; hex < Nhex; hex++)
			{
			  hex_values_perm[hex] = 0;
				if (Nintersections[hex] >= min_intersections)
				{
					for (int j = 0; j < Nintersections[hex]; j++)
					{
						int this_ellipse = intersections[hex][j];
					  hex_values_perm[hex] += edge_value[perm_vec[this_ellipse]] * area_inv[this_ellipse];
					}
					if (divide_weights) {
					  hex_values_perm[hex] /= hex_weights[hex];
					} else {
					  hex_values_perm[hex] /= double(Nintersections[hex]);
					}
				  if (hex_values_perm[hex] < hex_values[hex])
				  {
				    empirical_p[hex]++;
				  }
				}
			}
			//hex_values = hex_values_perm;
			/*
			double vmax1 = max(hex_values_perm[hex]);
			double vmin1 = min(hex_values_perm[hex]);
			double vinv1 = 1.0 / (vmax1 - vmin1);
			for (int hex = 0; hex < Nhex; hex++)
			{
				if (Nintersections[hex] > 0)
				{
					map_values1p[hex] = (map_values1p[hex] - vmin1)* vinv1;
					if (map_values1p[hex] - map_values2[hex] < map_values3[hex]) { empirical_p[hex] += dp; }
				}
			}
			*/
		} // end loop over Nperms

		chrono_timer(t0);
    /*
		for (int hex = 0; hex < Nhex; hex++)
		{
			if (Nintersections[hex] > 0)
			{
				if (empirical_p[hex] < 0.975 && empirical_p[hex] > 0.025) { empirical_p[hex] = 0.0; }
				else { empirical_p[hex] = 1.0; }
			}
			else
			{
				empirical_p[hex] = 0.0;
			}
		}
    */
	} // end Nperms > 0 condition

	//Normalize map values----------------------------------------------------------------------------------
  /*
	vmax3 = -1.0e99;
	vmin3 = 1.0e99;
	for (int hex = 0; hex < Nhex; hex++)
	{
		if (map_values3[hex] > vmax3) { vmax3 = map_values3[hex]; }
		if (map_values3[hex] < vmin3) { vmin3 = map_values3[hex]; }
	}
	vinv3 = vmax3 == vmin3 ? 1.0 : (1.0 / (vmax3 - vmin3));
	for (int hex = 0; hex < Nhex; hex++)
	{
		map_values3[hex] = (map_values3[hex] - vmin3) * vinv3;
		//Rcpp::Rcout << "\nMap1 value at hex " << hex << ": " << map_values1[hex] << " Map2 value: " << map_values2[hex] << " Map3 value: " << map_values3[hex] << " EP value: " << empirical_p[hex] << " Ints: " << Nintersections[hex] << " Wt: " << map_weights[hex];
	}

	//Return output as an Rcpp list ------------------------------------------------------------------------------------------------
	Rcpp::List ret;
	ret.push_back(Rcpp::wrap(Nintersections));
	ret.push_back(Rcpp::wrap(map_weights));
	ret.push_back(Rcpp::wrap(map_values1));
	ret.push_back(Rcpp::wrap(map_values2));
	ret.push_back(Rcpp::wrap(map_values3));
	ret.push_back(Rcpp::wrap(empirical_p));

	Rcpp::StringVector ret_names;
	ret_names.push_back("Nintersections");
	ret_names.push_back("map_weights");
	ret_names.push_back("map_values1");
	ret_names.push_back("map_values2");
	ret_names.push_back("map_values3");
	ret_names.push_back("empirical_p");

	ret.names() = ret_names;
	return ret;
  */
	
	return Rcpp::List::create(Rcpp::Named("hex_values") = hex_values,
                           Rcpp::Named("hex_weights") = hex_weights,
                           Rcpp::Named("empirical_p") = empirical_p,
                           Rcpp::Named("Nintersections") = Nintersections);
}