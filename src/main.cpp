
#include <Rcpp.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

//Global variables-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double pi = 3.1415926536;

//Additional functions-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double positive(double value);																												//Find absolute value of a number
int ellipse_check(double x, double y, double a_ellipse, double xf1_ellipse, double yf1_ellipse, double xf2_ellipse, double yf2_ellipse);	//Check whether point lies within ellipse

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List dummy1_cpp(Rcpp::List args)
{

	int node, node1, node2, Nnodes, ellipse, Nellipses, nx, ny, nxy, dim_matrix, area_matrix, perm, Nperms, rnd;
	double a_multiplier, x, y, dx, xmin, xmax, ymin, ymax, c, n, ellipse_sum, tmp, Cprob;
	vector<double> xnode;
	vector<double> ynode;
	vector<double> vnode;

	// convert Rcpp arguments to native c++ arguments------------------------------------------------------------------------------------------------------------------------------------------

	xnode = Rcpp::as<vector<double>>(args["xnode"]);		//Positions of data nodes on x-axis
	ynode = Rcpp::as<vector<double>>(args["ynode"]);		//Positions of data nodes on y-axis
	vnode = Rcpp::as<vector<double>>(args["vnode"]);		//Values at data nodes (of whatever type - calculation of ellipse values may have to change depending on type of values used)
	a_multiplier = Rcpp::as<double>(args["a_multiplier"]);	//Controls relationship between a (ellipse long radius) and c (ellipse short radius equal to distance between foci): a = c*(1 + a_multiplier)
	//xnode = { 1.0,2.0,3.0,4.0,5.0,6.0,1.0,2.0,3.0,4.0,5.0,6.0,1.0,2.0,3.0,4.0,5.0,6.0,1.0,2.0,3.0,4.0,5.0,6.0,1.0,2.0,3.0,4.0,5.0,6.0,1.0,2.0,3.0,4.0,5.0,6.0 };
	//ynode = { 1.0,1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,2.0,2.0,3.0,3.0,3.0,3.0,3.0,3.0,4.0,4.0,4.0,4.0,4.0,4.0,5.0,5.0,5.0,5.0,5.0,5.0,6.0,6.0,6.0,6.0,6.0,6.0 };
	//vnode = { 6.0,5.0,4.0,3.0,2.0,1.0,6.0,5.0,4.0,3.0,2.0,1.0,6.0,5.0,4.0,3.0,2.0,1.0,6.0,5.0,4.0,3.0,2.0,1.0,6.0,5.0,4.0,3.0,2.0,1.0,6.0,5.0,4.0,3.0,2.0,1.0 };
	//a_multiplier = -0.45;
	Nnodes = xnode.size();									//Number of nodes
	
	//Perform computations on imported data---------------------------------------------------------------------------------------------------------------------------------------------------

	dim_matrix = 101;	//Dimensions of matrix of map points

	//Calculate separation between points and x, y limits
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
			v_ellipse[ellipse] = positive(vnode[node1] - vnode[node2]);
			xf1_ellipse[ellipse] = xnode[node1];
			yf1_ellipse[ellipse] = ynode[node1];
			xf2_ellipse[ellipse] = xnode[node2];
			yf2_ellipse[ellipse] = ynode[node2];
			c = pow(pow(xnode[node1] - xnode[node2], 2.0) + pow(ynode[node1] - ynode[node2], 2.0), 0.5);
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
	area_matrix = dim_matrix*dim_matrix;
	vector<double> matrix_values(area_matrix, 0.0);						//Map data points
	vector<int> intersections(area_matrix*Nellipses, 0);				//Y/N data for which ellipses intersect each point
	nxy = 0;
	for (ny = 0; ny < dim_matrix; ny++)
	{
		y = ymin + (ny*dx);
		for (nx = 0; nx < dim_matrix; nx++)
		{
			x = xmin + (nx*dx);
			for (ellipse = 0; ellipse < Nellipses; ellipse++)
			{
				if (ellipse_check(x, y, a_ellipse[ellipse], xf1_ellipse[ellipse], yf1_ellipse[ellipse], xf2_ellipse[ellipse], yf2_ellipse[ellipse]) == 1) 
				{ 
					intersections[(ellipse*area_matrix) + nxy] = 1;
					matrix_values[nxy] += vw_ellipse[ellipse]; 
				}
			}
			nxy++;
		}
	}

	//Permute data to check statistical significance-------------------------------------------------------------------------------------------------------------------------------------------

	Nperms = 100;
	vector<double> vnode2(Nnodes, 0.0);									//Re-ordered node values
	vector<double> vw_ellipse2(Nellipses, 0.0);							//Weighted metric value of each ellipse after permutations
	vector<double> matrix_values_p(area_matrix*Nperms, 0.0);			//Map data for permutations

	vnode2 = vnode;
	for (perm = 0; perm < Nperms; perm++)
	{
		for (node = 0; node < Nnodes; node++)
		{
			rnd = R::runif(node, Nnodes); // draw random index from node to end of vector. Note that although runif returns a double, by forcing to int we essentially round this value down to nearest int.
			tmp = vnode2[rnd];   // temporarily store current value of vector at this position
			vnode2[rnd] = vnode2[node]; // swap for value at position i
			vnode2[node] = tmp;  // complete the swap
		}

		//Calculate new ellipse values
		ellipse = 0;
		for (node1 = 0; node1 < Nnodes; node1++)
		{
			for (node2 = node1 + 1; node2 < Nnodes; node2++)
			{
				vw_ellipse2[ellipse] = (positive(vnode2[node1] - vnode2[node2]) * ellipse_sum) / area_ellipse[ellipse];
				ellipse++;
			}
		}

		//Calculate new map data
		nxy = 0;
		for (ny = 0; ny < dim_matrix; ny++)
		{
			for (nx = 0; nx < dim_matrix; nx++)
			{
				matrix_values_p[(perm*area_matrix) + nxy] = 0.0;
				for (ellipse = 0; ellipse < Nellipses; ellipse++)
				{
					if (intersections[(ellipse*area_matrix) + nxy] == 1)
					{
						matrix_values_p[(perm*area_matrix) + nxy] += vw_ellipse2[ellipse];
					}
				}
				nxy++;
			}
		}
	}

	//Calculate cumulative probability distribution for each map point and determine whether observed value is statistically significant
	vector<double> permvalues(Nperms, 0.0);	
	nxy = 0;
	for (ny = 0; ny < dim_matrix; ny++)
	{
		for (nx = 0; nx < dim_matrix; nx++)
		{
			for (perm = 0; perm < Nperms; perm++)
			{
				permvalues[perm] = matrix_values_p[(perm*area_matrix) + nxy];
			}
			/*if (nxy == 5000)
			{
				Rprintf("\nPerm\tCprob\tValue1");
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
			if (Cprob < 0.95) { matrix_values[nxy] = 0.0; }
			nxy++;
		}
	}
	
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

	// return output as an Rcpp list-----------------------------------------------------------------------------------------------------------------------------------------------------------

	Rcpp::List ret;
	ret.push_back(Rcpp::wrap(dim_matrix));
	ret.push_back(Rcpp::wrap(xpoints));
	ret.push_back(Rcpp::wrap(ypoints));
	ret.push_back(Rcpp::wrap(xtick));
	ret.push_back(Rcpp::wrap(ytick));
	ret.push_back(Rcpp::wrap(matrix_values));

	Rcpp::StringVector ret_names;
	ret_names.push_back("dim_matrix");
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

int ellipse_check(double x, double y, double a_ellipse, double xf1_ellipse, double yf1_ellipse, double xf2_ellipse, double yf2_ellipse)
{
	int output = 0;
	if (0.5*(pow(pow(x - xf1_ellipse, 2.0) + pow(y - yf1_ellipse, 2.0), 0.5) + pow(pow(x - xf2_ellipse, 2.0) + pow(y - yf2_ellipse, 2.0), 0.5)) <= a_ellipse) { output = 1; }
	return output;
}