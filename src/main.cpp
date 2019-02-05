
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
	int ell, Nells, hex, node1, node2;
	double dist, vmin1, vmax1, vinv1, vmin2, vmax2, vinv2, vmin3, vmax3, vinv3, semi_minor;
	
	// start timer
	chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();

	//Convert Rcpp arguments to native c++ arguments ------------------------------------------------------------------------------------------------

	// load data and parameters
	print("Loading data and parameters");
	vector<double> long_node = rcpp_to_vector_double(args["long_node"]);    //Longitude of data nodes
	vector<double> lat_node = rcpp_to_vector_double(args["lat_node"]);      //Latitude of data nodes
	vector< vector<double> > vnode = rcpp_to_matrix_double(args["vnode"]);     //Pairwise metric data
	vector<double> long_hex = rcpp_to_vector_double(args["long_hex"]);      //Longitude of hex cells
	vector<double> lat_hex = rcpp_to_vector_double(args["lat_hex"]);        //Latitude of hex cells
	int Nperms = rcpp_to_int(args["n_perms"]);                               //Number of permutations to run (if 0, no permutation)
	double eccentricity = rcpp_to_double(args["eccentricity"]);             //Eccentricity of ellipses (see help for details)
	int flag_nullmap = rcpp_to_int(args["flag_nullmap"]);					//Flag indicating whether to create null map to subtract from data map (0-No, 1-Yes)
	int dist_model = rcpp_to_int(args["dist_model"]);						//Mathematical model governing relationship between distance and pairwise data (0-linear, 1-?)

	//Derived values/fixed constants ------------------------------------------------------------------------------------------------

	int Nnodes = long_node.size();                                          //Number of nodes
	int Nhex = long_hex.size();                                             //Number of hex cells
																			//vector< vector<double> > distances(Nnodes, vector<double>(Nnodes,0.0);  //Euclidean distances between points //TODO - make an import from R?

	//double dist_min = 0.0;                                                  //Minimum distance to nearest node for hex to be considered (set to 0 to disable)
	//double dist_minsq = sq(dist_min);                                       //TODO: Remove this feature or make an import from R?

	//Create ellipses ------------------------------------------------------------------------------------------------

	Rcpp::Rcout << "\nNperms: " << Nperms << "\nEccentricity: " << eccentricity << "\nFlag_nullmap: " << flag_nullmap << "\ndist_model: " << dist_model;
	print("Setting up ellipses");

	//Set up ellipses
	Nells = ((Nnodes - 1) * Nnodes) / 2;          //Number of ellipses
	vector<double> linear_eccentricity(Nells, 0.0);   //Linear eccentricity of ellipse, i.e. half the distance between foci
	vector<double> semi_major(Nells, 0.0);            //Semi-major axis of ellipse
	vector<double> xfocus1(Nells, 0.0);               //x-coordinate of focus 1 of each ellipse
	vector<double> yfocus1(Nells, 0.0);               //y-coordinate of focus 1 of each ellipse
	vector<double> xfocus2(Nells, 0.0);               //x-coordinate of focus 2 of each ellipse
	vector<double> yfocus2(Nells, 0.0);               //y-coordinate of focus 2 of each ellipse
	vector<double> xcentre(Nells, 0.0);               //x-coordinate of centre of each ellipse
	vector<double> ycentre(Nells, 0.0);               //y-coordinate of centre of each ellipse
	vector<double> area_inv(Nells, 0.0);              //Inverse of area of each ellipse
	vector<double> v_ell(Nells, 0.0);                 //Pairwise metric value of each ellipse (supplied data)
	vector<double> v_ell_null(Nells, 0.0);            //Pairwise metric value of each ellipse (null - distance only)
	vector<double> vw_ell(Nells, 0.0);                //Weighted metric value of each ellipse (pairwise data)
	vector<double> vw_ell_null(Nells, 0.0);           //Weighted metric value of each ellipse (null - distance only)

	ell = 0;
	for (node1 = 0; node1 < Nnodes; node1++)
	{
		for (node2 = node1 + 1; node2 < Nnodes; node2++)
		{
			// store foci and centre of ellipse
			xfocus1[ell] = long_node[node1];
			yfocus1[ell] = lat_node[node1];
			xfocus2[ell] = long_node[node2];
			yfocus2[ell] = lat_node[node2];
			xcentre[ell] = 0.5 * (long_node[node1] + long_node[node2]);
			ycentre[ell] = 0.5 * (lat_node[node1] + lat_node[node2]);

			// store long radius and area of ellipse
			dist = dist_euclid_2d(long_node[node1], lat_node[node1], long_node[node2], lat_node[node2]); //TODO - Remove this calculation if distances are imported
																										 //distances[node1][node2] = dist;
			linear_eccentricity[ell] = 0.5 * dist;
			semi_major[ell] = linear_eccentricity[ell] / eccentricity;
			semi_minor = sqrt(sq(semi_major[ell]) - sq(linear_eccentricity[ell]));
			area_inv[ell] = 1.0 / (M_PI * semi_major[ell] * semi_minor);
      
			// store original value attributed to ellipse
			v_ell[ell] = vnode[node1][node2];
			vw_ell[ell] = v_ell[ell] * area_inv[ell];
			if (flag_nullmap == 1)
			{
				switch (dist_model)
				{
				case 0: {v_ell_null[ell] = dist; }
						break;
				default: {}
				}
				vw_ell_null[ell] = v_ell_null[ell] * area_inv[ell];
			}
			ell++;
		}
	}

	chrono_timer(t0);

	//Create map by summation of ellipses intersecting each hex ------------------------------------------------------------------------------------------------

	print("Computing map");

	vector<double> map_values1(Nhex, 0.0);		//Map final values (from pairwise data)
	vector<double> map_values2(Nhex, 0.0);		//Map final values (null map based on distance)
	vector<double> map_values3(Nhex, 0.0);		//Map final values (map1-map2)
	vector<double> map_weights(Nhex, 0.0);		//Sum of weights (inverse area values)
	vector<int> Nintersections(Nhex, 0);		//Number of ellipses intersecting each hex
	vector<vector<int>> intersections(Nhex);	//List of ellipses intersecting each hex

	for (hex = 0; hex < Nhex; hex++)
	{
		//Rcpp::Rcout << "\nhex " << hex << ": lat=" << lat_hex[hex] << " long=" << long_hex[hex];
		/*if (dist_min > 0.0)
		{
		for (node1 = 0; node1 < Nnodes; node1++)
		{
		if (sq(long_hex[hex] - long_node[node1]) + sq(lat_hex[hex] - lat_node[node1]) > dist_minsq) { goto start; }
		}
		goto skip;
		}*/

		//start:
		for (ell = 0; ell < Nells; ell++)
		{
			if (ellipse_check(long_hex[hex], lat_hex[hex], xfocus1[ell], yfocus1[ell], xfocus2[ell], yfocus2[ell], semi_major[ell]))
			{
				intersections[hex].push_back(ell);
				Nintersections[hex]++;
				// TODO - option to downweight as get further from centre of ellipse
				//double d1 = dist_euclid_2d(long_hex[hex], lat_hex[hex], xfocus1[ell], yfocus1[ell]);
				//double d2 = dist_euclid_2d(long_hex[hex], lat_hex[hex], xfocus2[ell], yfocus2[ell]);
				//map_values1[hex] += (vw_ell[ell] * min(d1, d2)) / (d1 + d2);
				map_values1[hex] += vw_ell[ell];
				if (flag_nullmap == 1) { map_values2[hex] += vw_ell_null[ell]; }
				map_weights[hex] += area_inv[ell];
			}
		}
		//skip:
		if (Nintersections[hex] > 0)
		{
			map_values1[hex] = map_values1[hex] / map_weights[hex];
			map_values2[hex] = map_values2[hex] / map_weights[hex];
		}
		else
		{
			map_values1[hex] = 0.0;
			map_values2[hex] = 0.0;
		}
	}


	vmax1 = -1.0e99;
	vmin1 = 1.0e99;
	for (int hex = 0; hex < Nhex; hex++)
	{
		if (map_values1[hex] > vmax1) { vmax1 = map_values1[hex]; }
		if (map_values1[hex] < vmin1) { vmin1 = map_values1[hex]; }
	}
	vinv1 = 1.0 / (vmax1 - vmin1);

	if (flag_nullmap == 1) //Construct map_values3 by comparing normalized map_values1 and map_values2
	{
		vmax2 = -1.0e99;
		vmin2 = 1.0e99;
		for (int hex = 0; hex < Nhex; hex++)
		{
			if (map_values2[hex] > vmax2) { vmax2 = map_values2[hex]; }
			if (map_values2[hex] < vmin2) { vmin2 = map_values2[hex]; }
		}
		vinv2 = 1.0 / (vmax2 - vmin2);
		for (int hex = 0; hex < Nhex; hex++)
		{
			map_values1[hex] = (map_values1[hex] - vmin1)* vinv1;
			map_values2[hex] = (map_values2[hex] - vmin2)* vinv2;
			map_values3[hex] = map_values1[hex] - map_values2[hex];
		}
	}
	else //Construct map_values3 by simply normalizing map_values1
	{
		for (int hex = 0; hex < Nhex; hex++)
		{
			map_values3[hex] = (map_values1[hex] - vmin1)* vinv1;
		}
	}

	chrono_timer(t0);

	//Permute data to check statistical significance-------------------------------------------------------------------------------------------------------------------------------------------

	vector<double> empirical_p; //Empirical p-value, calculated by permutation test
	vector<double> map_values1p(Nhex, 0.0);		//Randomized version of map1
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
				if (Nintersections[hex] > 0)
				{
					for (int j = 0; j < Nintersections[hex]; j++)
					{
						int this_ellipse = intersections[hex][j];
						map_values1p[hex] += v_ell[perm_vec[this_ellipse]] * area_inv[this_ellipse];
					}
				}
			}
			vmax1 = -1.0e99;
			vmin1 = 1.0e99;
			for (int hex = 0; hex < Nhex; hex++)
			{
				if (Nintersections[hex] > 0)
				{
					map_values1p[hex] = map_values1p[hex] / map_weights[hex];
					if (map_values1p[hex] > vmax1) { vmax1 = map_values1p[hex]; }
					if (map_values1p[hex] < vmin1) { vmin1 = map_values1p[hex]; }
				}
			}
			vinv1 = 1.0 / (vmax1 - vmin1);
			for (int hex = 0; hex < Nhex; hex++)
			{
				if (Nintersections[hex] > 0)
				{
					map_values1p[hex] = (map_values1p[hex] - vmin1)* vinv1;
					if (map_values1p[hex] - map_values2[hex] < map_values3[hex]) { empirical_p[hex] += dp; }
				}
			}
		} // end loop over Nperms

		chrono_timer(t0);

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

	} // end if Nperms > 0

	  //Normalize map values----------------------------------------------------------------------------------

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
}