
#include <Rcpp.h>
#include <chrono>
#include <vector>
#include "misc.h"
#include "probability.h"

using namespace std;

//Global variables-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double pi = 3.1415926536;
Rcpp::List dummy2_cpp(Rcpp::List args);

// [[Rcpp::export]]
Rcpp::List run_sims_cpp(Rcpp::List args)
{
	// start timer
	chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();

	// initialise variables
	int Nnodes, ell, Nells, nx, ny, ni, coord, dim_matrix, area_matrix, Nperms, int_total;
	double a_multiplier, x, y, psep, long_min, long_max, lat_min, lat_max, matrix_value, n, ell_sum, vmax, d1, d2, dist_min, dist_min2, linear_eccentricity, semi_major_axis, semi_minor_axis;
	vector<double> long_node;
	vector<double> lat_node;
	vector< vector<double> > vnode;

	// convert Rcpp arguments to native c++ arguments ------------------------------------------------------------------------------------------------------------------------------------------

	// load data and parameters
	print("Loading data.\n");
	long_node = rcpp_to_vector_double(args["xnode"]);		//Positions of data nodes on longitude (x-) axis
	lat_node = rcpp_to_vector_double(args["ynode"]);		//Positions of data nodes on latitude (y-) axis
	vnode = rcpp_to_mat_double(args["vnode"]);				//Pairwise metric data
	a_multiplier = rcpp_to_double(args["a_multiplier"]);	//Controls relationship between a (ellipse semi-major axis) and c (linear eccentricity): a = c*(1 + a_multiplier)

	Nnodes = long_node.size();								//Number of nodes
	Nperms = rcpp_to_int(args["Nperms"]);					//Number of permutations to run (If 0, no permutation)
	dim_matrix = rcpp_to_int(args["dim_matrix"]);			//Dimensions of matrix of map points
	dist_min = 0.5;											//Minimum distance to nearest node for point to be considered	
	dist_min2 = sq(dist_min);
	
	print("a_multiplier =", a_multiplier, "\tNperms =", Nperms);

	//Perform computations on imported data ---------------------------------------------------------------------------------------------------------------------------------------------------

	chrono_timer(t0);
	print("Establishing map boundaries.\n");

	long_min = min(long_node);
	long_max = max(long_node);
	lat_min = min(lat_node);
	lat_max = max(lat_node);

	if ((long_max - long_min) > (lat_max - lat_min))
	{
		psep = (long_max - long_min) / (dim_matrix - 3);
		long_max += 1.0 * psep;
		lat_max = lat_min + ((dim_matrix - 1) * psep);
	}
	else
	{
		psep = (lat_max - lat_min) / (dim_matrix - 3);
		lat_max += 1.0 * psep;
		long_max = long_min + ((dim_matrix - 1) * psep);
	}
	long_min -= 1.0 * psep;
	lat_min -= 1.0 * psep;

	chrono_timer(t0);
	print("Setting up ellipses.\n");

	//Set up ellipses
	Nells = ((Nnodes - 1) * Nnodes) / 2;      //Number of ellipses
	vector<double> a_ell(Nells, 0.0);         //Long radius of each ellipse
	vector<double> asq_ell(Nells, 0.0);        //Square of long radius of each ellipse
	vector<double> xf1_ell(Nells, 0.0);       //x-coordinate of focus 1 of each ellipse
	vector<double> yf1_ell(Nells, 0.0);       //y-coordinate of focus 1 of each ellipse
	vector<double> xf2_ell(Nells, 0.0);       //x-coordinate of focus 2 of each ellipse
	vector<double> yf2_ell(Nells, 0.0);       //y-coordinate of focus 2 of each ellipse
	vector<double> xc_ell(Nells, 0.0);        //x-coordinate of centre of each ellipse
	vector<double> yc_ell(Nells, 0.0);        //y-coordinate of centre of each ellipse
	vector<double> area_inv_ell(Nells, 0.0);  //Inverse of area of each ellipse
	vector<double> v_ell(Nells, 0.0);         //Pairwise metric value of each ellipse
	vector<double> vw_ell(Nells, 0.0);        //Weighted metric value of each ellipse

	ell = 0;
	for (int node1 = 0; node1 < Nnodes; node1++)
	{
		for (int node2 = node1 + 1; node2 < Nnodes; node2++)
		{
		  // store foci and centre of ellipse
			xf1_ell[ell] = long_node[node1];
			yf1_ell[ell] = lat_node[node1];
			xf2_ell[ell] = long_node[node2];
			yf2_ell[ell] = lat_node[node2];
			xc_ell[ell] = 0.5 * (long_node[node1] + long_node[node2]);
			yc_ell[ell] = 0.5 * (lat_node[node1] + lat_node[node2]);
			
			// store long radius and area of ellipse
			linear_eccentricity = 0.5 * dist_euclid_2d(long_node[node1], lat_node[node1], long_node[node2], lat_node[node2]);
			semi_major_axis = linear_eccentricity * (1.0 + a_multiplier);
			semi_minor_axis = sqrt(sq(semi_major_axis) - sq(linear_eccentricity));
			a_ell[ell] = semi_major_axis;
			asq_ell[ell] = sq(a_ell[ell]);
			area_inv_ell[ell] = 1.0 / (pi * semi_major_axis * semi_minor_axis);
			
			// store original value attributed to ellipse
			v_ell[ell] = vnode[node1][node2];
			vw_ell[ell] = v_ell[ell] * area_inv_ell[ell];
			ell++;
		}
	}

	//Create map by summation of ellipses intersecting each point

	chrono_timer(t0);
	print("Setting up map.\n");

	area_matrix = sq(dim_matrix);
	vector<double> matrix_values(area_matrix, 0.0);                         //Map data points
	vector<int> nintersections(area_matrix, 0);                             //Number of ellipses intersecting each point
	vector<vector<int>> intersections(area_matrix, vector<int>(Nells, 0));  //List of ellipses intersecting each point
	
	coord = 0;
	int_total = 0;
	y = lat_min;
	
	for (ny = 0; ny < dim_matrix; ny++)
	{
		x = long_min;
		for (nx = 0; nx < dim_matrix; nx++)
		{
			ell_sum = 0.0;
			matrix_value = 0.0;
			for (int node1 = 0; node1 < Nnodes; node1++)
			{
				if (sq(x - long_node[node1]) + sq(y - lat_node[node1]) < dist_min2) { goto start; }
			}
			goto skip;
			start:
			for (ell = 0; ell < Nells; ell++)
			{
				if (circle_check(x, y, xc_ell[ell], yc_ell[ell], asq_ell[ell]))
				{
					if (ellipse_check(x, y, a_ell[ell], xf1_ell[ell], yf1_ell[ell], xf2_ell[ell], yf2_ell[ell]))
					{
						intersections[coord][nintersections[coord]] = ell;
						nintersections[coord]++;
						int_total++;
						//d1 = dist_euclid_2d(x, y, xf1_ell[ell], yf1_ell[ell]);
						//d2 = dist_euclid_2d(x, y, xf2_ell[ell], yf2_ell[ell]);
						matrix_value += vw_ell[ell];
						//matrix_value += (vw_ell[ell] * min(d1, d2)) / (d1 + d2);
						ell_sum += area_inv_ell[ell];
					}
				}
			}
			skip:
			matrix_values[coord] = ell_sum > 0.0 ? matrix_value / ell_sum : -1.0 / 0.0;
			coord++;
			x += psep;
		}
		y += psep;
	}

	//Permute data to check statistical significance-------------------------------------------------------------------------------------------------------------------------------------------

	chrono_timer(t0);
	vector<double> matrix_values2(area_matrix, 0.0);								//Data points after any alterations
	vector<int> intersections2(int_total, 0);										//Condensed list of which ellipses intersect with each map point
	Rcpp::List to_perm;
	Rcpp::List from_perm;
	Rcpp::StringVector to_perm_names;
	int ip_total = 0;
	if (Nperms == 0)
	{
		matrix_values2 = matrix_values;
		goto finish;
	}

	int_total = 0;
	for (coord = 0; coord < area_matrix; coord++)
	{
		if (nintersections[coord] > 0)
		{
			ip_total++;
			for (ni = 0; ni < nintersections[coord]; ni++)
			{
				intersections2[int_total] = intersections[coord][ni];
				int_total++;
			}
		}
	}
	Rprintf("\nNumber of points with ellipse intersections = %i / %i\n", ip_total, area_matrix);

	to_perm.push_back(Rcpp::wrap(Nperms));
	to_perm.push_back(Rcpp::wrap(Nells));
	to_perm.push_back(Rcpp::wrap(dim_matrix));
	to_perm.push_back(Rcpp::wrap(v_ell));
	to_perm.push_back(Rcpp::wrap(area_inv_ell));
	to_perm.push_back(Rcpp::wrap(nintersections));
	to_perm.push_back(Rcpp::wrap(intersections2));
	to_perm.push_back(Rcpp::wrap(matrix_values));

	to_perm_names.push_back("Nperms");
	to_perm_names.push_back("Nells");
	to_perm_names.push_back("dim_matrix");
	to_perm_names.push_back("v_ell");
	to_perm_names.push_back("area_inv_ell");
	to_perm_names.push_back("nintersections");
	to_perm_names.push_back("intersections2");
	to_perm_names.push_back("matrix_values");
	to_perm.names() = to_perm_names;

	from_perm = dummy2_cpp(to_perm);
	matrix_values2 = rcpp_to_vector_double(from_perm["matrix_values2"]);

	//Set up final data to be output to R-----------------------------------------------------------------------------------------------------------------------------------------------------
finish:

	print("Setting up output data.\n");

	vector<double> xpoints(dim_matrix, 0.0);
	vector<double> ypoints(dim_matrix, 0.0);
	for (nx = 0; nx < dim_matrix; nx++)
	{
		xpoints[nx] = long_min + (nx*psep);
		ypoints[nx] = lat_min + (nx*psep);
	}
	n = (dim_matrix - 2) / 5.0;
	vector<double> xtick(6, 0.0);
	vector<double> ytick(6, 0.0);
	xtick[0] = long_min + (0.5 * psep);
	ytick[0] = lat_min + (0.5 * psep);
	for (nx = 1; nx < 6; nx++)
	{
		xtick[nx] = xtick[nx - 1] + (n * psep);
		ytick[nx] = ytick[nx - 1] + (n * psep);
	}
	vmax = -1.0e99;
	for (coord = 0; coord < area_matrix; coord++)
	{
		if (matrix_values2[coord] > vmax) { vmax = matrix_values2[coord]; }
	}
	for (coord = 0; coord < area_matrix; coord++)
	{
		matrix_values2[coord] = matrix_values2[coord] / vmax;
	}

	// return output as an Rcpp list-----------------------------------------------------------------------------------------------------------------------------------------------------------

	Rcpp::List ret;
	ret.push_back(Rcpp::wrap(dim_matrix));
	ret.push_back(Rcpp::wrap(long_node));
	ret.push_back(Rcpp::wrap(lat_node));
	ret.push_back(Rcpp::wrap(xpoints));
	ret.push_back(Rcpp::wrap(ypoints));
	ret.push_back(Rcpp::wrap(xtick));
	ret.push_back(Rcpp::wrap(ytick));
	ret.push_back(Rcpp::wrap(matrix_values2));

	Rcpp::StringVector ret_names;
	ret_names.push_back("dim_matrix");
	ret_names.push_back("xnode");
	ret_names.push_back("ynode");
	ret_names.push_back("xpoints");
	ret_names.push_back("ypoints");
	ret_names.push_back("xtick");
	ret_names.push_back("ytick");
	ret_names.push_back("matrix_values");
	ret.names() = ret_names;

	print("Output data set up.\n");

	return ret;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Rcpp::List dummy2_cpp(Rcpp::List args)
{
	int Nells, coord, ni, dim_matrix, area_matrix, perm, Nperms, this_ell, int_total;
	double matrix_value, ell_sum, Cprob;

	// start timer
	chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
	print("Loading data for statistical significance check.\n");

	Nperms = Rcpp::as<int>(args["Nperms"]);															//Number of permutations
	Nells = Rcpp::as<int>(args["Nells"]);															//Number of ellipses
	dim_matrix = Rcpp::as<int>(args["dim_matrix"]);													//Dimensions of map
	area_matrix = dim_matrix * dim_matrix;
	vector<double> v_ell = rcpp_to_vector_double(args["v_ell"]);									//Pairwise metric values
	vector<double> area_inv_ell = rcpp_to_vector_double(args["area_inv_ell"]);						//Inverse ellipse areas
	vector<int> nintersections = rcpp_to_vector_int(args["nintersections"]);						//Number of ellipses intersecting each point
	vector<int> intersections2 = rcpp_to_vector_int(args["intersections2"]);						//List of ellipses intersecting each point
	vector<double> matrix_values = rcpp_to_vector_double(args["matrix_values"]);					//Original map data

	vector<int> perm_index = seq_int(1, Nells);														//Index of which ellipse value to take
	vector<vector<double>> matrix_values_p(area_matrix, vector<double>(Nperms, 0.0));				//Map data for permutations
	vector<double> permvalues(Nperms, 0.0);
	vector<double> matrix_values2(area_matrix, 0.0);												//Output map data
	matrix_values2 = matrix_values;

	chrono_timer(t0);
	print("Starting permutation.\n");

	for (perm = 0; perm < Nperms; perm++)
	{
		// report progress periodically
		Rprintf("p=%i\n", perm);

		// reshuffle perm_index
		reshuffle(perm_index);

		//Calculate new map data
		int_total = 0;
		for (coord = 0; coord < area_matrix; coord++)
		{
			if (nintersections[coord] > 0)
			{
				ell_sum = 0.0;
				matrix_value = 0.0;
				for (ni = 0; ni < nintersections[coord]; ni++)
				{
					this_ell = intersections2[int_total];
					int_total++;
					matrix_value += v_ell[perm_index[this_ell]] * area_inv_ell[this_ell];
					ell_sum += area_inv_ell[this_ell];
				}
				matrix_values_p[coord][perm] = matrix_value / ell_sum;
			}
		}
	}

	chrono_timer(t0);
	print("Permutation values generated. Calculating cumulative probability distribution at each point\n");

	//Calculate cumulative probability distribution for each map point and determine whether observed value is statistically significant

	for (coord = 0; coord < area_matrix; coord++)
	{
		for (perm = 0; perm < Nperms; perm++)
		{
			permvalues[perm] = matrix_values_p[coord][perm];
		}
		sort(permvalues.begin(), permvalues.begin() + Nperms);
		Cprob = 1.0;
		for (perm = 0; perm < Nperms; perm++)
		{
			if (permvalues[perm] >= matrix_values[coord])
			{
				Cprob = (perm + 0.5) / Nperms;
				goto done;
			}
		}
	done:
		if (Cprob < 0.95 && Cprob > 0.05) { matrix_values2[coord] = -1.0 / 0.0; }
	}

	chrono_timer(t0);
	print("Statistical significance check complete. Outputting data back to main function.\n");

	Rcpp::List from_perm;
	from_perm.push_back(Rcpp::wrap(matrix_values2));

	Rcpp::StringVector from_perm_names;
	from_perm_names.push_back("matrix_values2");
	from_perm.names() = from_perm_names;

	chrono_timer(t0);
	print("Statistical significance check data output complete\n");

	return from_perm;
}