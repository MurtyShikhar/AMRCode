#include "math.h"
#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/hankel.hpp>

#define pi 3.14159
dcomp isq = -1;
dcomp jota = sqrt(isq);
double k_b;
const int M_mls = 100;
int M;

using namespace boost::math;


bool sort_complex_abs(dcomp a, dcomp b){
  return abs(a) < abs(b);
}

void multiply(comp_matrix& A, comp_matrix& B, comp_matrix& C)
{
	int sz1 = A.size();
	int sz2 = B[0].size();
	int sz3 = B.size();
	for (int i = 0; i < sz1; ++i) {
		for (int k = 0; k < sz3; ++k) {
			dcomp r = A[i][k];
			for (int j = 0; j < sz2; ++j){
				C[i][j] += r*B[k][j];
			}
		}
	}
		
	
}

//Populate the coordinates of source, g
void populate_source_coords(vector<pair<double, double> >& source_coords, double diameter, int M)
{
	double theta_step = (2*pi)/M;
	double r = diameter/2;
	double x, y;
	for (int i = 0; i < M; ++i)
		source_coords[i] = make_pair(r*cos(i*theta_step),r*sin(i*theta_step));
}

void generate_aux_gd(comp_matrix& Aux_GD, vector<cell* >& cells)
{
	//Set size as cells.size()*cells.size()
	double l1,l2;
	for (int i = 0; i < cells.size(); ++i)
	{
		for (int l = 0; l < cells.size(); ++l)
		{
			l1 = pow((cells[i]->coord.first - cells[l]->coord.first), 2);
			l2 = pow((cells[i]->coord.second - cells[l]->coord.second), 2);
			double p_il = sqrt(l1+l2);

			double a = sqrt(cells[i]->get_area())/(sqrt(pi));
			if(i == l)
				Aux_GD[i][l] = (k_b*k_b)*(jota/4.0)*(2/(k_b*k_b))*(pi*k_b*a*cyl_hankel_1(1,  k_b*a) + 2.0*jota);
			else
				Aux_GD[i][l] = (k_b*k_b)*(jota/4.0)*(2*pi*a/k_b)*(cyl_bessel_j(1, k_b*a))*cyl_hankel_1(0, k_b*p_il);
		}
	}
}

void generate_u_inc(comp_matrix& u_inc, vector<cell*>& domain, vector<double>& thetas)
{
	for (int i = 0; i < thetas.size(); ++i)
	{
		for (int j = 0; j < domain.size(); ++j)
		{
			u_inc[j][i] = exp(jota*k_b*( (domain[j]->coord.first)*cos(thetas[i]) + (domain[j]->coord.second)*sin(thetas[i])));
		}
	}
}

void generate_aux_gs(comp_matrix& Aux_GS, vector<cell* >& cells, vector<pair<double, double> >& source_coords, comp_matrix& u_inc)
{
	//Set size as M*cells.size()
	double l1,l2;
	for (int i = 0; i < M; ++i)
	{
		for (int l = 0; l < cells.size(); ++l)
		{

			l1 = pow(source_coords[i].first - cells[l]->coord.first, 2);
			l2 = pow(source_coords[i].second - cells[l]->coord.second, 2);
			double p_il = sqrt(l1+l2);
			double a = sqrt(cells[l]->get_area())/(sqrt(pi));
			Aux_GS[i][l] = (k_b*k_b)*(jota/4.0)*(2*pi*a/k_b)*(cyl_bessel_j(1, k_b*a)*cyl_hankel_1(0, k_b*p_il));
		}	

	}
}

void back_propagation(comp_matrix& Aux_GD, comp_matrix& Aux_GS, comp_matrix& f, comp_matrix& u_inc, vector<cell* >& cells)
{
	//Make Aux_GS* local vector of vector of dcomp
		//Set size as cells.size()*M
	int sz = cells.size();
	comp_matrix Aux_GS_star(sz, vector<dcomp>(M));
	for (int i = 0; i < Aux_GS.size(); ++i)
	{
		for (int j = 0; j < Aux_GS[i].size(); ++j)
			Aux_GS_star[j][i] = conj(Aux_GS[i][j]);	
	}

	comp_matrix numerator(sz, vector<dcomp>(M)); //N2xM
	comp_matrix denominator(M, vector<dcomp>(M)); //MxM
	multiply(Aux_GS_star, f, numerator);
	multiply(Aux_GS, numerator, denominator);

	comp_matrix w_0(sz, vector<dcomp>(M));

	for (int i = 0; i < M; ++i)
	{
		double norm_numerator = 0;
		double norm_denominator = 0;

		for (int l = 0; l < sz; ++l)
			norm_numerator += pow(abs(numerator[l][i]), 2);

		for (int l = 0; l < M; ++l)
			norm_denominator += pow(abs(denominator[l][i]), 2);

		//norm_numerator = sqrt(norm_numerator);
		//norm_denominator = sqrt(norm_denominator);

		for (int j = 0; j < sz; ++j)
		{
			w_0[j][i] = (norm_numerator/norm_denominator)*numerator[j][i];
		}

	}

	//vector<dcomp> w_sort;
	//for (int i = 0; i < sz; i++)
		//w_sort.push_back(w_0[i][1]);


	//sort(w_sort.begin(), w_sort.end(), sort_complex_abs);


	comp_matrix temp(sz, vector<dcomp>(M));
	comp_matrix u_0(sz, vector<dcomp>(M));
	multiply(Aux_GD, w_0, temp);

	for (int i = 0; i < sz; ++i)
		for (int j = 0; j < M; ++j)
			u_0[i][j] = u_inc[i][j] + temp[i][j];


	//vector<dcomp> chi_0(sz);
	for (int i = 0; i < sz; ++i){
		dcomp numerator = 0;
		double denominator = 0;
		for (int j = 0; j < M; ++j)
		{
			numerator += w_0[i][j]*conj(u_0[i][j]);
			denominator += pow(abs(u_0[i][j]), 2);
		}
		cells[i]->chi = numerator/denominator;
	}
}

void back_prop_master(vector<cell*>& domain, comp_matrix& f, vector<double>& thetas, vector<point>& source_coords)
{
	comp_matrix aux_GD(domain.size(), vector<dcomp>(domain.size()));
	comp_matrix aux_GS(M, vector<dcomp>(domain.size()));
	comp_matrix u_inc(domain.size(), vector<dcomp>(M));

	generate_aux_gd(aux_GD, domain);
	generate_u_inc(u_inc, domain, thetas);
	generate_aux_gs(aux_GS, domain, source_coords, u_inc);
	back_propagation(aux_GD, aux_GS, f, u_inc, domain);
}