#include <queue>
#include "math.h"
#include <algorithm>
#include <set>
#include <list>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/hankel.hpp>
//#include "init.cpp"

bool sort_complex_real(cell* a, cell* b){
	return real((a->chi)) < real((b->chi));
}


struct compareCellPtrs{
	bool operator()(const cell* a, const cell* b) const {
		double ax = a->coord.first;
		double ay = a->coord.second;
		double bx = b->coord.first;
		double by = b->coord.second;
		return (ax == bx ? ay - by : ax- bx);
	}
};


typedef set<cell*> set_cells;

void inline copy_vector(set<cell* >& a, vector<cell* >& b) {
	copy(a.begin(), a.end(), back_inserter(b));
}

/* scans through "cells" and pushes into populate, all the neighbors of *curr (including curr) */
void push_neigbors(cell* curr, set<cell* >& populate, vector<cell* >& cells) {
	int inspect = 0;
	for (int i = 0; i < cells.size(); i++) {
		inspect = populate.size();
		if (curr->neighbor(*cells[i]))
			populate.insert(cells[i]);
	}
}

void display_chi(vector<cell*>& domain)
{
	for (int i = 0; i < domain.size(); ++i)
	{
		cout<<domain[i]->chi<<"\n";
	}
}
/* inputs: the domain to be examined by multilevel sampling, as a vector<cell* > */
/* output: a vector<cell* > containing a collection of cells that need to be discretised further */
bool indices_to_flag(vector<cell* >& old_cells,vector<cell* >& new_flagged, double& cutoff) {
	vector<cell* > cells;
	int cutoff_index;
	double new_cutoff;
	bool stop = false;
	for (int i = 0; i < old_cells.size(); i++) {
		if (real(old_cells[i]->chi) >= cutoff)
			cells.push_back(old_cells[i]);
	}
	//pngwriter trial(width_image, width_image, 0, "trial.png");
	//display_chi(cells);
	//trial.close();

	/* cells contains all old_cells[i], that have a chi greater than cutoff */

	sort(cells.begin(), cells.end(), sort_complex_real);
	double smallest = -1; // smallest difference between consecutive cells
	for (int i = 0; i < cells.size()-1; ++i) {
		dcomp delta = (cells[i+1]->chi)- (cells[i]->chi);
		if (smallest == -1)
			smallest = real(delta);
		else
			smallest = (real(delta) < smallest ? real(delta) : smallest);
	}

	/* smallest = "abs of the smallest difference" between chi values of neighboring cells */

	/* decide the new cutoff value, and the cutoff index */
	for (int j = 1; j < cells.size() -1; ++j) {
		if(real((cells[j+1]->chi - cells[j]->chi)) >= M_mls*smallest){
			new_cutoff = real(cells[j+1]->chi);
			cutoff_index = j+1;
			break;
		}
	}

	set<cell* > unique_cells;
	int inspect = 0;
	for (int i = cutoff_index; i < cells.size(); i++) {
		inspect = unique_cells.size();
		push_neigbors(cells[i], unique_cells, cells);  
	}


	copy_vector(unique_cells, new_flagged);

	int inspect2 = new_flagged.size();
	double _eps = 0.0002; /* might want to change this */
	if (abs(cutoff- new_cutoff) < _eps)
		stop = true;

	cutoff = new_cutoff;
	return stop;

}


void divide(vector<cell*>& flagged, int discretization, int start, vector<double>& thetas, comp_matrix& f, vector<pair<double, double> >& source_coords) {
	vector<cell* > cofd;
	double cutoff = 0.0; 
	bool stop = false;
	for (int i = start; i < discretization && !stop; i++) {
		for (int j = 0; j < flagged.size(); j++) {
			cell* fr = flagged[j];
			fr->flag = true;
			typedef vector<cell* >::iterator itr;
			itr en = fr->children.end();
			for (itr it = fr->children.begin(); it != en; it++)
				cofd.push_back(*it);
		}

		flagged.clear();
		cout << cofd.size() << "size of current domain" << "\n";
		back_prop_master(cofd,f, thetas, source_coords);
		stop = indices_to_flag(cofd,flagged, cutoff); 
		cofd.clear();

	}
}

// vector "domain" consists of all the nodes that are finally part of the domain
void reconstruct_domain(vector<cell* >& final_domain, vector<cell* >& initial_domain) {
	queue<cell* > flagged_nd;
	for (int i = 0; i < initial_domain.size(); i++)
		flagged_nd.push(initial_domain[i]);


	for(;;) {
		if (flagged_nd.empty()) break;

		cell* fr = flagged_nd.front(); flagged_nd.pop();
		typedef vector<cell* >::iterator itr;
		itr en = fr->children.end();
		for (itr it = fr->children.begin(); it != en; it++) {
			if ((*it)->flag)
				flagged_nd.push(*it);
			else
				final_domain.push_back(*it);

		}
	}
}


// int main() {
// 	cell* root = new cell;
// 	double side_len = 20;
// 	root->set_coords(0,0);
// 	root->set_area(side_len*side_len);
// 	create_tree(root,side_len);
// 	pretty_print(root);
// 	return 0;
// 	cout << "HEEHHE";
// 	return 0;
// }