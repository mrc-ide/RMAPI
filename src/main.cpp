
#include <Rcpp.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include "misc.h"

using namespace std;

//Global variables-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double pi = 3.1415926536;

//Additional functions-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double positive(double value);																												//Find absolute value of a number
int ellipse_check(double x, double y, double a, double xf1, double yf1, double xf2, double yf2);	//Check whether point lies within ellipse

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List dummy1_cpp(Rcpp::List args)
{
	int node, node1, node2, Nnodes, ellipse, Nellipses, nx, ny, nxy, ni, dim_matrix, area_matrix, perm, Nperms, rnd1, rnd2;
	double a_multiplier, x, y, dx, xmin, xmax, ymin, ymax, c, n, ellipse_sum, tmp, Cprob, vmax;
	vector<double> xnode;
	vector<double> ynode;
	vector< vector<double> > vnode;
    
	// convert Rcpp arguments to native c++ arguments------------------------------------------------------------------------------------------------------------------------------------------

	Rprintf("\nLoading data.\n");
	xnode = Rcpp_to_vector_double(args["xnode"]);			//Positions of data nodes on x-axis
	ynode = Rcpp_to_vector_double(args["ynode"]);			//Positions of data nodes on y-axis
	vnode = Rcpp_to_mat_double(args["vnode"]);				//Pairwise metric data
	a_multiplier = Rcpp::as<double>(args["a_multiplier"]);	//Controls relationship between a (ellipse long radius) and c (ellipse short radius equal to distance between foci): a = c*(1 + a_multiplier)
	Nnodes = xnode.size();									//Number of nodes
	
    //printMatrix(vnode);
    //return Rcpp::List::create(Rcpp::Named("foo")=9);
    
	//Perform computations on imported data---------------------------------------------------------------------------------------------------------------------------------------------------

	dim_matrix = 101;	//Dimensions of matrix of map points

	//Calculate separation between points and x, y limits
	
	Rprintf("\nEstablishing map boundaries.\n");

	xmin = 1.0e99;		
	xmax = -1.0e99;
	ymin = 1.0e99;
	ymax = -1.0e99;
	for (node = 0; node < Nnodes; node++)
	{
		if (xnode[node] < xmin) { xmin = xnode[node]; }
		else { if (xnode[node] > xmax) { xmax = xnode[node]; } }
		if (ynode[node] < ymin) { ymin = ynode[node]; }
		else { if (ynode[node] > ymax) { ymax = ynode[node]; } }
	}
	if (xmax - xmin > ymax - ymin)
	{
		dx = (xmax - xmin) / (dim_matrix - 10);
		xmax += 5.0*dx;
		ymax = ymin + ((dim_matrix - 1)*dx);
	}
	else
	{
		dx = (ymax - ymin) / (dim_matrix - 10);
		ymax += 5.0*dx;
		xmax = xmin + ((dim_matrix - 1)*dx);
	}
	xmin -= 5.0*dx;
	ymin -= 5.0*dx;

	Rprintf("\nSetting up ellipses.\n");

	//Set up ellipses
	Nellipses = ((Nnodes - 1)*Nnodes) / 2;				//Number of ellipses
	vector<double> a_ellipse(Nellipses, 0.0);			//Long radius of each ellipse
	vector<double> xf1_ellipse(Nellipses, 0.0);			//x-coordinate of focus 1 of each ellipse
	vector<double> yf1_ellipse(Nellipses, 0.0);			//y-coordinate of focus 1 of each ellipse
	vector<double> xf2_ellipse(Nellipses, 0.0);			//x-coordinate of focus 2 of each ellipse
	vector<double> yf2_ellipse(Nellipses, 0.0);			//y-coordinate of focus 2 of each ellipse
	vector<double> v_ellipse(Nellipses, 0.0);			//Metric value of each ellipse
	vector<double> area_ellipse(Nellipses, 0.0);		//Area of each ellipse
	vector<double> vw_ellipse(Nellipses, 0.0);			//Weighted metric value of each ellipse

	ellipse = 0;
	ellipse_sum = 0.0;
	for (node1 = 0; node1 < Nnodes; node1++)
	{
		for (node2 = node1 + 1; node2 < Nnodes; node2++)
		{
			v_ellipse[ellipse] = vnode[node1][node2];
			xf1_ellipse[ellipse] = xnode[node1];
			yf1_ellipse[ellipse] = ynode[node1];
			xf2_ellipse[ellipse] = xnode[node2];
			yf2_ellipse[ellipse] = ynode[node2];
			c = pow(((xnode[node1] - xnode[node2])*(xnode[node1] - xnode[node2])) + ((ynode[node1] - ynode[node2])*(ynode[node1] - ynode[node2])), 0.5);
			a_ellipse[ellipse] = c*(1 + a_multiplier);
			area_ellipse[ellipse] = pi*a_ellipse[ellipse]*c;
			ellipse_sum += area_ellipse[ellipse];
			ellipse++;
		}
	}

	for (ellipse = 0; ellipse < Nellipses; ellipse++)
	{
		vw_ellipse[ellipse] = (v_ellipse[ellipse] * ellipse_sum) / area_ellipse[ellipse];
	}

	//Create map by summation of ellipses intersecting each point

	Rprintf("\nSetting up map.\n");

	area_matrix = dim_matrix*dim_matrix;
	vector<double> matrix_values(area_matrix, 0.0);								//Map data points
	vector<int> nintersections(area_matrix, 0);									//Number of ellipses intersecting each point
	vector<vector<int>> intersections(area_matrix, vector<int>(Nellipses,0));	//List of ellipses intersecting each point
	nxy = 0;
	y=ymin;
	for (ny = 0; ny < dim_matrix; ny++)
	{
		x = xmin;
		for (nx = 0; nx < dim_matrix; nx++)
		{
			for (ellipse = 0; ellipse < Nellipses; ellipse++)
			{
				if (ellipse_check(x, y, a_ellipse[ellipse], xf1_ellipse[ellipse], yf1_ellipse[ellipse], xf2_ellipse[ellipse], yf2_ellipse[ellipse]) == 1) 
				{
					intersections[nxy][nintersections[nxy]] = ellipse;
					nintersections[nxy]++;
					matrix_values[nxy] += vw_ellipse[ellipse]; 
				}
			}
			nxy++;
			x += dx;
		}
		y += dx;
	}

	//Permute data to check statistical significance-------------------------------------------------------------------------------------------------------------------------------------------
	    
	Nperms = 1000;
	vector< vector<double> > vnode2;															//Re-ordered node values
	vector<double> vw_ellipse2(Nellipses, 0.0);													//Weighted metric value of each ellipse after permutations
	//vector<double> matrix_values_p(area_matrix*Nperms, 0.0);									//Map data for permutations
	vector<vector<double>> matrix_values_p(area_matrix, vector<double>(Nperms, 0.0));				//Map data for permutations

	Rprintf("\nStarting permutation.\n");

	vnode2 = vnode;
	for (perm = 0; perm < Nperms; perm++)
	{
		//Rprintf("p=%i\n", perm);
		//Rprintf("\nnode1\tnode2\tV1\trnd1\trnd2\tV2");
		for (node1 = 0; node1 < Nnodes; node1++)
		{
			for (node2 = node1 + 1; node2 < Nnodes; node2++)
			{				
				rnd1 = R::runif(0, Nnodes); // draw random index from node to end of vector. Note that although runif returns a double, by forcing to int we essentially round this value down to nearest int.
				rnd2 = R::runif(rnd1 + 1, Nnodes);
				tmp = vnode2[rnd1][rnd2];   // temporarily store current value of vector at this position 
				vnode2[rnd1][rnd2] = vnode2[node1][node2]; // swap for value at position i
				vnode2[node1][node2] = tmp;  // complete the swap
				//Rprintf("\n%i\t%i\t%.2e\t%i\t%i\t%.2e", node1, node2, vnode[node1][node2], rnd1, rnd2, vnode2[node1][node2]);
			}
		}

		//Calculate new ellipse values
		ellipse = 0;
		for (node1 = 0; node1 < Nnodes; node1++)
		{
			for (node2 = node1 + 1; node2 < Nnodes; node2++)
			{
				vw_ellipse2[ellipse] = (vnode2[node1][node2] * ellipse_sum) / area_ellipse[ellipse];
				ellipse++;
			}
		}

		//Calculate new map data
		nxy = 0;
		for (ny = 0; ny < dim_matrix; ny++)
		{
			for (nx = 0; nx < dim_matrix; nx++)
			{
				matrix_values_p[nxy][perm] = 0.0;
				if (nintersections[nxy] > 0)
				{
					for (ni = 0; ni < nintersections[nxy]; ni++)
					{
						matrix_values_p[nxy][perm] += vw_ellipse2[intersections[nxy][ni]];
					}
				}
				nxy++;
			}
		}
	}

	Rprintf("\nNew values generated. Calculating cumulative probability distribution at each point.\n");

	//Calculate cumulative probability distribution for each map point and determine whether observed value is statistically significant
	vector<double> permvalues(Nperms, 0.0);	
	nxy = 0;
	for (ny = 0; ny < dim_matrix; ny++)
	{
		for (nx = 0; nx < dim_matrix; nx++)
		{
			for (perm = 0; perm < Nperms; perm++)
			{
				permvalues[perm] = matrix_values_p[nxy][perm];
			}
			/*if (nxy == 200)
			{
				Rprintf("\nObserved:\t%.3f\nPerm\tCprob\tValue1", matrix_values[nxy]);
				for (perm = 0; perm < Nperms; perm++)
				{
					Rprintf("\n%i\t%.3f\t%.1f", perm, (perm + 1.0) / Nperms, permvalues[perm]);
				}
			}*/
			sort(permvalues.begin(), permvalues.begin() + Nperms);
			Cprob = 1.0;
			for (perm = 0; perm < Nperms; perm++)
			{
				if(permvalues[perm] >= matrix_values[nxy])
				{ 
					Cprob = (perm + 1.0) / Nperms; 
					goto done;
				}
			}
			done:	
			if (Cprob <= 0.95 && Cprob >= 0.05) { matrix_values[nxy] = 0.0; }
			nxy++;
		}
	}
	
	Rprintf("\nPermutation complete.\n");

	//Set up final data to be output to R-----------------------------------------------------------------------------------------------------------------------------------------------------
	vector<double> xpoints(dim_matrix, 0.0);
	vector<double> ypoints(dim_matrix, 0.0);
	for (nx = 0; nx < dim_matrix; nx++)
	{ 
		xpoints[nx] = xmin + (nx*dx);
		ypoints[nx] = ymin + (nx*dx);
	}
	n = (dim_matrix - 2) / 5.0;
	vector<double> xtick(6, 0.0);
	vector<double> ytick(6, 0.0);
	xtick[0] = xmin + (0.5*dx);
	ytick[0] = ymin + (0.5*dx);
	for (nx = 1; nx < 6; nx++)
	{
		xtick[nx] = xtick[nx - 1] + (n*dx);
		ytick[nx] = ytick[nx - 1] + (n*dx);
	}
	vmax = -1.0e99;
	for (nxy = 0; nxy < area_matrix; nxy++)
	{
		matrix_values[nxy] = log(matrix_values[nxy]);
		//if (matrix_values[nxy] > vmax) { vmax = matrix_values[nxy]; }
	}
	/*for (nxy = 0; nxy < area_matrix; nxy++)
	{
		matrix_values[nxy] *= 1.0 / vmax;
	}*/

	// return output as an Rcpp list-----------------------------------------------------------------------------------------------------------------------------------------------------------
	    
	Rcpp::List ret;
	ret.push_back(Rcpp::wrap(dim_matrix));
	ret.push_back(Rcpp::wrap(xnode));
	ret.push_back(Rcpp::wrap(ynode));
	ret.push_back(Rcpp::wrap(xpoints));
	ret.push_back(Rcpp::wrap(ypoints));
	ret.push_back(Rcpp::wrap(xtick));
	ret.push_back(Rcpp::wrap(ytick));
	ret.push_back(Rcpp::wrap(matrix_values));

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
    return ret ;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double positive(double value)
{
	double output = max(value, -value);
	return output;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int ellipse_check(double x, double y, double a, double xf1, double yf1, double xf2, double yf2)
{
	int output = 0;
	if (0.5*(pow(((x - xf1)*(x - xf1)) + ((y - yf1)*(y - yf1)), 0.5) + pow(((x - xf2)*(x - xf2)) + ((y - yf2)*(y - yf2)), 0.5)) <= a) { output = 1; }
	return output;
}
