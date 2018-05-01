
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
<<<<<<< HEAD
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
=======
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
  // TODO - keep dist_min?
  //double dist_min = 100;                            //Minimum distance to nearest node for point to be considered	
  //double dist_minsq = sq(dist_min);
  
  //Create ellipses ------------------------------------------------------------------------------------------------
  
	print("Setting up ellipses");
  
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
>>>>>>> b767aaa46ea5a96739687c3e26e63ca382597b04
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
<<<<<<< HEAD
			linear_eccentricity = 0.5 * dist_euclid_2d(long_node[node1], lat_node[node1], long_node[node2], lat_node[node2]);
			semi_major_axis = linear_eccentricity * (1.0 + a_multiplier);
			semi_minor_axis = sqrt(sq(semi_major_axis) - sq(linear_eccentricity));
			a_ell[ell] = semi_major_axis;
			asq_ell[ell] = sq(a_ell[ell]);
			area_inv_ell[ell] = 1.0 / (pi * semi_major_axis * semi_minor_axis);
=======
			linear_eccentricity[ell] = 0.5 * dist_euclid_2d(long_node[node1], lat_node[node1], long_node[node2], lat_node[node2]);
			semi_major[ell] = linear_eccentricity[ell] / eccentricity;
			double semi_minor = sqrt(sq(semi_major[ell]) - sq(linear_eccentricity[ell]));
			area_inv[ell] = 1.0 / (pi * semi_major[ell] * semi_minor);
>>>>>>> b767aaa46ea5a96739687c3e26e63ca382597b04
			
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
		/*
		for (int node1 = 0; node1 < Nnodes; node1++)
		{
			if (sq(x - long_node[node1]) + sq(y - lat_node[node1]) < dist_min2) {
			  goto start;
			}
<<<<<<< HEAD
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
=======
>>>>>>> b767aaa46ea5a96739687c3e26e63ca382597b04
		}
		goto skip;
		start:
		*/
		for (int ell = 0; ell < Nells; ell++)
		{
			if (ellipse_check(long_hex[hex], lat_hex[hex], xfocus1[ell], yfocus1[ell], xfocus2[ell], yfocus2[ell], semi_major[ell]))
			{
			  intersections[hex].push_back(ell);
				Nintersections[hex]++;
				// TODO - option to downweight as get further from centre of ellipse
				//d1 = dist_euclid_2d(x, y, xf1_ell[ell], yf1_ell[ell]);
				//d2 = dist_euclid_2d(x, y, xf2_ell[ell], yf2_ell[ell]);
				map_values[hex] += vw_ell[ell];
				//matrix_value += (vw_ell[ell] * min(d1, d2)) / (d1 + d2);
				map_weights[hex] += area_inv[ell];
			}
		}
		map_values[hex] = Nintersections[hex] > 0 ? map_values[hex] / map_weights[hex] : 0;
	}
<<<<<<< HEAD
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

=======
	chrono_timer(t0);
	
	//Permute data to check statistical significance-------------------------------------------------------------------------------------------------------------------------------------------
	
	vector<double> empirical_p; //Empirical p-value, calculated by permutation test
	if (Nperms>0) {
	  
	  print("Permutation testing");
	  
	  // create vector for storing empirical p-values
	  empirical_p = vector<double>(Nhex);
	  
	  // create vector for indexing random permutations
	  vector<int> perm_vec = seq_int(0, Nells-1);
	  
	  // loop through permutations
	  for (int perm = 0; perm < Nperms; perm++)
	  {
	    reshuffle(perm_vec);  // new permutation
	    for (int hex = 0; hex < Nhex; hex++)
	    {
	      if (Nintersections[hex]==0) {
	        continue;
	      }
	      double vw_sum = 0.0;
	      for (int j = 0; j < Nintersections[hex]; j++)
	      {
	        int this_ellipse = intersections[hex][j];
	        vw_sum += v_ell[perm_vec[this_ellipse]] * area_inv[this_ellipse];
	      }
	      if (vw_sum / map_weights[hex] < map_values[hex]) {
	        empirical_p[hex] ++;
	      }
	    }
	  } // end loop over Nperms
	  
	  chrono_timer(t0);
	  
	} // end if Nperms > 0
	
  //Return output as an Rcpp list ------------------------------------------------------------------------------------------------
  
>>>>>>> b767aaa46ea5a96739687c3e26e63ca382597b04
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
