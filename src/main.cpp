
#include <Rcpp.h>
#include <chrono>
#include <vector>
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
//Global variables
double pi = 3.1415926536;

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List run_sims_cpp(Rcpp::List args)
{
	// start timer
	chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();

	//Convert Rcpp arguments to native c++ arguments ------------------------------------------------------------------------------------------------

	// load data and parameters
	print("Loading data and parameters");
	vector<double> long_node = rcpp_to_vector_double(args["long_node"]);    //Longitude of data nodes
	vector<double> lat_node = rcpp_to_vector_double(args["lat_node"]);      //Latitude of data nodes
	vector< vector<double> > vnode = rcpp_to_mat_double(args["vnode"]);     //Pairwise metric data
	vector<double> long_hex = rcpp_to_vector_double(args["long_hex"]);      //Longitude of hex cells
	vector<double> lat_hex = rcpp_to_vector_double(args["lat_hex"]);        //Latitude of hex cells
	int Nperms = rcpp_to_int(args["Nperms"]);                               //Number of permutations to run (if 0, no permutation)
	double eccentricity = rcpp_to_double(args["eccentricity"]);             //Eccentricity of ellipses (see help for details)

	int Nnodes = long_node.size();              //Number of nodes
	int Nhex = long_hex.size();                 //Number of hex cells
	
	double dist_min = 0.0;						//Minimum distance to nearest node for hex to be considered (set to 0 to disable)
	double dist_minsq = sq(dist_min);

	//Create ellipses ------------------------------------------------------------------------------------------------

	print("Setting up ellipses");
	print("Nperms\teccentricity\tdist_min");
	print(Nperms, eccentricity, dist_min);

	//Set up ellipses
	int Nells = ((Nnodes - 1) * Nnodes) / 2;          //Number of ellipses
	vector<double> linear_eccentricity(Nells, 0.0);   //Linear eccentricity of ellipse, i.e. half the distance between foci
	vector<double> semi_major(Nells, 0.0);            //Semi-major axis of ellipse
	vector<double> xfocus1(Nells, 0.0);               //x-coordinate of focus 1 of each ellipse
	vector<double> yfocus1(Nells, 0.0);               //y-coordinate of focus 1 of each ellipse
	vector<double> xfocus2(Nells, 0.0);               //x-coordinate of focus 2 of each ellipse
	vector<double> yfocus2(Nells, 0.0);               //y-coordinate of focus 2 of each ellipse
	vector<double> xcentre(Nells, 0.0);               //x-coordinate of centre of each ellipse
	vector<double> ycentre(Nells, 0.0);               //y-coordinate of centre of each ellipse
	vector<double> area_inv(Nells, 0.0);              //Inverse of area of each ellipse
	vector<double> v_ell(Nells, 0.0);                 //Pairwise metric value of each ellipse
	vector<double> vw_ell(Nells, 0.0);                //Weighted metric value of each ellipse

	int ell = 0;
	for (int node1 = 0; node1 < Nnodes; node1++)
	{
		for (int node2 = node1 + 1; node2 < Nnodes; node2++)
		{
			// store foci and centre of ellipse
			xfocus1[ell] = long_node[node1];
			yfocus1[ell] = lat_node[node1];
			xfocus2[ell] = long_node[node2];
			yfocus2[ell] = lat_node[node2];
			xcentre[ell] = 0.5 * (long_node[node1] + long_node[node2]);
			ycentre[ell] = 0.5 * (lat_node[node1] + lat_node[node2]);

			// store long radius and area of ellipse
			linear_eccentricity[ell] = 0.5 * dist_euclid_2d(long_node[node1], lat_node[node1], long_node[node2], lat_node[node2]);
			semi_major[ell] = linear_eccentricity[ell] / eccentricity;
			double semi_minor = sqrt(sq(semi_major[ell]) - sq(linear_eccentricity[ell]));
			area_inv[ell] = 1.0 / (pi * semi_major[ell] * semi_minor);

			// store original value attributed to ellipse
			v_ell[ell] = vnode[node1][node2];
			vw_ell[ell] = v_ell[ell] * area_inv[ell];
			ell++;
		}
	}

	//Create map by summation of ellipses intersecting each hex ------------------------------------------------------------------------------------------------

	print("Computing basic map");

	vector<double> map_values(Nhex, 0.0);       //Map final values
	vector<double> map_weights(Nhex, 0.0);      //Sum of weights
	vector<int> Nintersections(Nhex, 0);        //Number of ellipses intersecting each hex
	vector<vector<int>> intersections(Nhex);    //List of ellipses intersecting each hex

	for (int hex = 0; hex < Nhex; hex++)
	{
		if (dist_min > 0.0)
		{
			for (int node1 = 0; node1 < Nnodes; node1++)
			{
				if (sq(long_hex[hex] - long_node[node1]) + sq(lat_hex[hex] - lat_node[node1]) < dist_minsq) { goto start; }
			}
			goto skip;
		}
		
		start:
		for (int ell = 0; ell < Nells; ell++)
		{
			if (ellipse_check(long_hex[hex], lat_hex[hex], xfocus1[ell], yfocus1[ell], xfocus2[ell], yfocus2[ell], semi_major[ell]))
			{
				intersections[hex].push_back(ell);
				Nintersections[hex]++;
				// TODO - option to downweight as get further from centre of ellipse
				//double d1 = dist_euclid_2d(long_hex[hex], lat_hex[hex], xfocus1[ell], yfocus1[ell]);
				//double d2 = dist_euclid_2d(long_hex[hex], lat_hex[hex], xfocus2[ell], yfocus2[ell]);
				//map_values[hex] += (vw_ell[ell] * min(d1, d2)) / (d1 + d2);
				map_values[hex] += vw_ell[ell];
				map_weights[hex] += area_inv[ell];
			}
		}
		skip:
		map_values[hex] = Nintersections[hex] > 0 ? map_values[hex] / map_weights[hex] : 0.0;
	}
	chrono_timer(t0);

	//Permute data to check statistical significance-------------------------------------------------------------------------------------------------------------------------------------------

	vector<double> empirical_p; //Empirical p-value, calculated by permutation test
	if (Nperms > 0) 
	{

		print("Running permutation");

		// create vector for storing empirical p-values
		empirical_p = vector<double>(Nhex);

		// create vector for indexing random permutations
		vector<int> perm_vec = seq_int(0, Nells - 1);
		double dp = 1.0 / Nperms;

		// loop through permutations
		for (int perm = 0; perm < Nperms; perm++)
		{
			reshuffle(perm_vec);  // new permutation
			for (int hex = 0; hex < Nhex; hex++)
			{
				if (Nintersections[hex] == 0) {continue;}
				double vw_sum = 0.0;
				for (int j = 0; j < Nintersections[hex]; j++)
				{
					int this_ellipse = intersections[hex][j];
					vw_sum += v_ell[perm_vec[this_ellipse]] * area_inv[this_ellipse];
				}
				if (vw_sum / map_weights[hex] < map_values[hex]) {
					empirical_p[hex] +=dp;
				}
			}
		} // end loop over Nperms

		chrono_timer(t0);

	} // end if Nperms > 0

	//Return output as an Rcpp list ------------------------------------------------------------------------------------------------

	double vmax = 0.0;
	double vmin = 1.0e99;
	for (int hex = 0; hex < Nhex; hex++)
	{
		if (map_values[hex] > vmax) { vmax = map_values[hex]; }
		if (map_values[hex] < vmin) { vmin = map_values[hex]; }
	}
	for (int hex = 0; hex < Nhex; hex++)
	{
		map_values[hex] = (map_values[hex] - vmin) / (vmax - vmin);
		//if (empirical_p[hex] > 0.95) { empirical_p[hex] = 1.0; }
		//else { empirical_p[hex] = 0.0; }
	}

	Rcpp::List ret;
	ret.push_back(Rcpp::wrap(map_values));
	ret.push_back(Rcpp::wrap(Nintersections));
	ret.push_back(Rcpp::wrap(map_weights));
	ret.push_back(Rcpp::wrap(empirical_p));

	Rcpp::StringVector ret_names;
	ret_names.push_back("map_values");
	ret_names.push_back("Nintersections");
	ret_names.push_back("map_weights");
	ret_names.push_back("empirical_p");

	ret.names() = ret_names;
	return ret;
}
