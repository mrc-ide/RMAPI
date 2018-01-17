
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

double positive(double value);
double ellipse_check(double x, double y, double a_ellipse, double xf1_ellipse, double yf1_ellipse, double xf2_ellipse, double yf2_ellipse, double v_ellipse);

// [[Rcpp::export]]
Rcpp::List dummy1_cpp(Rcpp::List args) {

	int node, node1, node2, Nnodes, ellipse, Nellipses, nx, ny, nxy, dim_matrix;
	double a_multiplier, x, y, dx, xmin, xmax, ymin, ymax, da, c, n;
	vector<double> xnode;
	vector<double> ynode;
	vector<double> vnode;

	// convert Rcpp args to native c++ args
	xnode = Rcpp::as<vector<double>>(args["xnode"]);
	ynode = Rcpp::as<vector<double>>(args["ynode"]);
	vnode = Rcpp::as<vector<double>>(args["vnode"]);
	a_multiplier = Rcpp::as<double>(args["a_multiplier"]);
	Nnodes = xnode.size();

	//------------------------------------------------

	// do some really hard maths

	xmin = 1.0e99;
	xmax = -1.0e99;
	ymin = 1.0e99;
	ymax = -1.0e99;
	dim_matrix = 31;

	for (node = 0; node < Nnodes; node++)
	{
		if (xnode[node] < xmin) { xmin = xnode[node]; }
		if (ynode[node] < ymin) { ymin = ynode[node]; }
		if (xnode[node] > xmax) { xmax = xnode[node]; }
		if (ynode[node] > ymax) { ymax = ynode[node]; }
		//Rprintf("xnode[%i]=%.2f\tynode[%i]=%.2f\tvnode[%i]=%.2f\n", node,xnode[node], node,ynode[node], node,vnode[node]);
	}


	if (xmax - xmin > ymax - ymin)
	{
		dx = (xmax - xmin) / (dim_matrix - 10);
		xmin -= 5.0*dx;
		xmax += 5.0*dx;
		ymin -= 5.0*dx;
		ymax = ymin + ((dim_matrix - 1)*dx);
	}
	else
	{
		dx = (ymax - ymin) / (dim_matrix - 10);
		ymin -= 5.0*dx;
		ymax += 5.0*dx;
		xmin -= 5.0*dx;
		xmax = xmin + ((dim_matrix - 1)*dx);
	}
	
	Nellipses = ((Nnodes - 1)*Nnodes) / 2;
	vector<double> a_ellipse(Nellipses, 0.0);
	vector<double> xf1_ellipse(Nellipses, 0.0);
	vector<double> yf1_ellipse(Nellipses, 0.0);
	vector<double> xf2_ellipse(Nellipses, 0.0);
	vector<double> yf2_ellipse(Nellipses, 0.0);
	vector<double> v_ellipse(Nellipses, 0.0);

	ellipse = 0;
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
			da = c*a_multiplier;
			a_ellipse[ellipse] = c + da;
			ellipse++;
		}
	}

	vector<double> matrix_values(dim_matrix*dim_matrix, 0.0);
	nxy = 0;
	for (ny = 0; ny < dim_matrix; ny++)
	{
		y = ymin + (ny*dx);
		for (nx = 0; nx < dim_matrix; nx++)
		{
			x = xmin + (nx*dx);
			for (ellipse = 0; ellipse < Nellipses; ellipse++)
			{
				matrix_values[nxy] += ellipse_check(x, y, a_ellipse[ellipse], xf1_ellipse[ellipse], yf1_ellipse[ellipse], xf2_ellipse[ellipse], yf2_ellipse[ellipse], v_ellipse[ellipse]);
			}
			nxy++;
		}
	}

	vector<double> xpoints(dim_matrix, 0.0);
	vector<double> ypoints(dim_matrix, 0.0);
	for (nx = 0; nx < dim_matrix; nx++)
	{ 
		xpoints[nx] = xmin + (nx*dx);
		ypoints[nx] = ymin + (nx*dx);
	}
	n = (dim_matrix - 1) / 5.0;
	vector<double> xtick(6, 0.0);
	vector<double> ytick(6, 0.0);
	for (nx = 0; nx < 6; nx++)
	{
		xtick[nx] = xmin + (nx*n*dx);
		ytick[nx] = ymin + (nx*n*dx);
	}

	//------------------------------------------------

	// return output as an Rcpp list
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

//-----------------------------------------------------

double positive(double value)
{
	double output = max(value, -value);
	return output;
}

//-----------------------------------------------------

double ellipse_check(double x, double y, double a_ellipse, double xf1_ellipse, double yf1_ellipse, double xf2_ellipse, double yf2_ellipse, double v_ellipse)
{
	double output = 0.0;
	double dist = 0.5*(pow(pow(x - xf1_ellipse, 2.0) + pow(y - yf1_ellipse, 2.0), 0.5) + pow(pow(x - xf2_ellipse, 2.0) + pow(y - yf2_ellipse, 2.0), 0.5));
	if (dist <= a_ellipse) { output = v_ellipse; }
	return output;
}