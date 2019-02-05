
#include "misc_v3.h"

#include <Rcpp.h>
#include <chrono>
#include <vector>

using namespace std;

//Sub-function(s)
double v_intersect(double xl1, double yl1, double xl2, double yl2, double xp, double yp, double rbarrier);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List hexbarrier01(Rcpp::List args_h)
{
	int barr, node1, node2;
	double x1, y1, x2, y2, vb;

	// start timer
	chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();

	//Convert Rcpp arguments to native c++ arguments ------------------------------------------------------------------------------------------------

	vector<double> long_node = rcpp_to_vector_double(args_h["long_node"]);			//Longitude of data nodes
	vector<double> lat_node = rcpp_to_vector_double(args_h["lat_node"]);			//Latitude of data nodes
	vector<double> long_barrier = rcpp_to_vector_double(args_h["long_barrier"]);    //Longitude of barrier cells
	vector<double> lat_barrier = rcpp_to_vector_double(args_h["lat_barrier"]);      //Latitude of barrier cells
	int Nnodes = long_node.size();													//Number of nodes
	int Nbarriers = long_barrier.size();                                            //Number of barrier cells
	double rbarrier = rcpp_to_double(args_h["rbarrier"]);							//Barrier cell radius
	vector<double> vbarrier = rcpp_to_vector_double(args_h["vbarrier"]);   		//Barrier cell distance multiplier

																					//Calculate pairwise data ------------------------------------------------------------------------------------------------

	print("Calculating distance modifiers due to barriers");

	vector<double> vbsum(Nnodes*Nnodes, 0.0);

	for (node1 = 0; node1 < Nnodes; node1++)
	{
		x1 = long_node[node1];
		y1 = lat_node[node1];
		for (node2 = node1 + 1; node2 < Nnodes; node2++)
		{
			x2 = long_node[node2];
			y2 = lat_node[node2];
			vb = 0.0;
			for (barr = 0; barr < Nbarriers; barr++)
			{
				vb += (vbarrier[barr] - 1.0)*v_intersect(x1, y1, x2, y2, long_barrier[barr], lat_barrier[barr], rbarrier);
			}

			vbsum[(node1*Nnodes) + node2] = vb;
		}
	}

	chrono_timer(t0);

	//Return output as an Rcpp list ------------------------------------------------------------------------------------------------

	Rcpp::List ret;
	ret.push_back(Rcpp::wrap(vbsum));
	Rcpp::StringVector ret_names;
	ret_names.push_back("vbsum");
	ret.names() = ret_names;

	return ret;
}

//------------------------------------------------------------------------------------------------------------------------------------

double v_intersect(double xl1, double yl1, double xl2, double yl2, double xp, double yp, double rbarrier)
{
	double value = 0.0;
	double dist, numerator, denominator;

	if (xl1 > xp && xl2 > xp) { goto skip; }
	if (xl1 < xp && xl2 < xp) { goto skip; }
	if (yl1 > yp && yl2 > yp) { goto skip; }
	if (yl1 < yp && yl2 < yp) { goto skip; }

	numerator = ((yl2 - yl1)*xp) - ((xl2 - xl1)*yp) + (xl2*yl1) - (yl2*xl1);
	denominator = dist_euclid_2d(xl1, yl1, xl2, yl2);
	dist = max(numerator, -numerator) / denominator;
	if (dist < rbarrier)
	{
		value = sqrt(1.0 - sq(dist / rbarrier));
	}

skip:
	return value;
}
