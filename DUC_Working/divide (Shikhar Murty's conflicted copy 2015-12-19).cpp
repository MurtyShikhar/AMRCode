#include <queue>
#include "math.h"
#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/hankel.hpp>

bool sort_complex_real(cell* a, cell* b){
	return real((a->chi)) < real((b->chi));
}

/* scans through "cells" and pushes into populate, all the neighbors of *curr (including curr) */
void push_neigbors(cell* curr, vector<cell* >& populate, vector<cell* > cells) {
	for (int i = 0; i < cells.size(); i++) {
		if (curr->neighbor(*cells[i]))
			populate.push_back(cells[i]);
	}
}


/* inputs: the domain to be examined by multilevel sampling, as a vector<cell* > */
/* output: a vector<cell* > containing a collection of cells that need to be discretised further */
bool indices_to_flag(vector<cell* >& old_cells,vector<cell* >& new_flagged, double& cutoff) {
	vector<cell* > cells;
	int cutoff_index;
	int new_cutoff;
	bool stop = false;
	for (int i = 0; i < old_cells.size(); i++) {
		if (abs(old_cells[i]->chi) >= cutoff)
			cells.push_back(old_cells[i]);
	}

	/* cells contains all old_cells[i], that have a chi greater than cutoff */

	sort(cells.begin(), cells.end(), sort_complex_real);
	double smallest = -1; // smallest difference between consecutive cells
	for (int i = 0; i < cells.size()-1; ++i) {
		dcomp delta = (cells[i+1]->chi)- (cells[i]->chi);
		if (smallest == -1)
			smallest = abs(delta);
		else
			smallest = (abs(delta) < smallest ? abs(delta) : smallest);
	}

	/* smallest = "abs of the smallest difference" between chi values of neighboring cells */


	/* decide the new cutoff value, and the cutoff index */
	for (int j = 1; j < cells.size() -1; ++j) {
		if(abs((cells[j+1]->chi - cells[j]->chi)) >= M_mls*smallest){
			new_cutoff = real(cells[j+1]->chi);
			cutoff_index = j+1;
			break;
		}
	}

	for (int i = cutoff_index; i < cells.size(); i++) {
		push_neigbors(cells[i], new_flagged, cells);  
	}

	double _eps = 0.0002; /* might want to change this */
	if (abs(cutoff- new_cutoff) < _eps)
		stop = true;

	cutoff = new_cutoff;
	return stop;

}


void divide(vector<cell*>& flagged, int discretization, int start, vector<double>& thetas, comp_matrix& f, vector<point>& source_coords) {
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
		back_prop_master(cofd,f, thetas, source_coords);
		stop = indices_to_flag(cofd,flagged, cutoff); 
		cofd.clear();

	}
}

// vector "domain" consists of all the nodes that are finally part of the domain
void reconstruct_domain(vector<cell* >& domain, cell* root) {
	queue<cell* > flagged_nd;
	flagged_nd.push(root);
	for(;;) {
		if (flagged_nd.empty()) break;

		cell* fr = flagged_nd.front(); flagged_nd.pop();
		typedef vector<cell* >::iterator itr;
		itr en = fr->children.end();
		for (itr it = fr->children.begin(); it != en; it++) {
			if ((*it)->flag)
				flagged_nd.push(*it);
			else
				domain.push_back(*it);

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