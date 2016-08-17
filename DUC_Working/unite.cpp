#include <cassert>
#include <algorithm>
#include <set>
using namespace std;

bool sort_by_y(const point& a, const point& b) {
	return a.second < b.second;
}

bool sort_by_x(const point& a, const point& b) {
	return a.first < b.first;
}

edge union_edge(const edge&e1, const edge& e2)
{
	vector<point> points(4);
	points[0] = e1.s;
	points[1] = e1.e;
	points[2] = e2.s;
	points[3] = e2.e;
	if(e1.v)
	{
		sort(points.begin(), points.end(), sort_by_y);
		edge e3(points[0], points[3], 1);
		return e3;
	}
	else
	{
		sort(points.begin(), points.end(), sort_by_x);
		edge e3(points[0], points[3], 0);
		return e3;
	}
}

void push_merged(set<edge>& edges, const edge& e1, const edge& e2) {
	vector<point> points(4);
	points[0] = e1.s;
	points[1] = e1.e;
	points[2] = e2.s;
	points[3] = e2.e;
	if (e1.v) {
		sort(points.begin(), points.end(), sort_by_y);
		edge e1(points[0], points[1], 1);
		edge e2(points[2], points[3], 1);
		edge e3(points[0], points[3], 1);

		if(points[1] == points[2])
		{
			edges.insert(e3);
			return;
		}
		else
		{
			if (e1.len() != 0) edges.insert(e1);
			if (e2.len() != 0) edges.insert(e2);
		return;
		}
	}
	else {
		sort(points.begin(), points.end(), sort_by_x);
		edge e1(points[0], points[1], 0);
		edge e2(points[2], points[3], 0);
		edge e3(points[0], points[3], 0);

		if(points[1] == points[2])
		{
			edges.insert(e3);
			return;
		}
		else
		{
			if (e1.len() != 0) edges.insert(e1);
			if (e2.len() != 0) edges.insert(e2);
		return;
		}
	}
}

bool common_point(const edge& a, const edge& b) {
	if (a.v != b.v)
		return false;
	else {
		if (a.v) {
			return (a.e.second == b.s.second || a.s.second == b.e.second) && (a.e.first == b.e.first) && !(a ==b);
		}

		else {
			return (a.e.first == b.s.first || a.s.first == b.e.first) && (a.e.second == b.e.second) && !(a == b);

		}
	}
}

/* resets the edges of cell c  */
void reset_edges(cell* c) {
	set<edge> edge_set;
	vector<edge> all_edges;

	int num_children = c->children.size();

	for (int i = 0; i < num_children; i++) {
		/* iterate over all children */
		int edge_sz = c->children[i]->edges.size();

		for (int j = 0; j < edge_sz; j++)  {
			/* push all its edges into both the edge_set and all_edges */
			edge_set.insert(c->children[i]->edges[j]);
			all_edges.push_back(c->children[i]->edges[j]);
		}
	}

	/* now we must compare all pairs of edges to see if they touch */
	int all_edges_sz = all_edges.size();
	for (int i = 0; i < all_edges_sz; i++ ) {
		for (int j = 0; j < all_edges_sz; j++) {
			if (i == j) continue;

			/* if 2 edges of children are the same, they cancel each other*/
			if (all_edges[i] == all_edges[j]) {
				set<edge>::iterator itr = edge_set.find(all_edges[i]);
				if (itr != edge_set.end()) {
					 edge_set.erase(itr);
				}
					
			}

			/* if they touch, then we delete both, and add the "merged edges" */
			else if (touch(all_edges[i], all_edges[j])) {
				edge_set.erase(all_edges[i]);
				edge_set.erase(all_edges[j]);
				push_merged(edge_set, all_edges[i], all_edges[j]);
			}

		}
	}

	set<edge>::iterator itr;
	set<edge>::iterator en = edge_set.end();
	c->edges.clear();
	for (itr = edge_set.begin(); itr != en; itr++) {
		c->add_edge(*itr);
	}

	for (int i = 0; i < c->edges.size(); ++i) {
		for (int j = 0; j < c->edges.size(); ++j) {
			if (i == j) continue;

			if (common_point(c->edges[i], c->edges[j])) {
				edge_set.erase(c->edges[i]);
				edge_set.erase(c->edges[j]);
				edge_set.insert(union_edge(c->edges[i], c->edges[j]));	
			}	
		}	
	}



	c->edges.clear();
	for (itr = edge_set.begin(); itr != en; itr++) {
		c->add_edge(*itr);
	}

}



/* BRUTE FORCE LCA ALGORITHM O(N) */
cell* lca(cell* a, cell* b) {
	int d1 = a->depth;
	int d2 = b->depth;
	cell* p1 = a;
	cell* p2 = b;
	while (d1 != d2) {
		if (d1 > d2)  {
			p1 = p1->parent;
			d1--;
		}
		else {
			p2 = p2->parent;
			d2--;
		}
	}
	while (p1 != p2) {
		p1 = p1->parent;
		p2 = p2->parent;
	}

	assert (p1 != NULL);
	return p1;
}

void recompute_depths(cell* a) {
	if (a == NULL) return;
	a->depth = (a->parent->depth)+ 1;
	int sz = a->children.size();
	for (int i = 0; i < sz; i++)
		recompute_depths(a->children[i]);
}

void inline set_coords(cell* united_cell, cell* a, cell* b) {
	double a1 = a->area, a2  = b->area;
	double x = (((a->coord).first)*a1 + ((b->coord).first)*a2)/(a1+a2);
	double y = (((a->coord).second)*a1 + ((b->coord).second)*a2)/(a1+a2);
	united_cell->coord = make_pair(x, y);
}


void fuse(cell* c) {
	cout << "HI\n";
	// if c has a parent
	if (c->depth > 0) {
		c->parent->add_child(c->children[0]);
		c->children[0]->parent = c->parent;
		c->parent->remove_child(c);
	}

	// c is the root. and it has only one child
	else {
		c->children[0]->parent = NULL;
	}

	recompute_depths(c->children[0]);
	delete c;
}



/* INPUT:  pointers to 2 neighboring cells a and b */
/* OUTPUT: a new_cell formed by merging the 2 cells */

/* 
   MODIFIES: splices a and b off its parents, and makes the united cell the new parent. 
   then adds the united cell to the LCA of a and b
*/

cell* unite(cell* a, cell* b) {
	cell* ancestor = lca(a, b); // lca of a and b returns the lowest common ancestor of a and b
	assert (a->neighbor(*b));
	cell* united_cell = new cell;
	united_cell->depth = (ancestor->depth) + 1;
	united_cell->add_child(a); // add a and b as children of united_cell
	united_cell->add_child(b);
	united_cell->area = a->area + b->area;
	set_coords(united_cell, a, b); 
	united_cell->parent = ancestor; // parent of united_cell is ancestor
	ancestor->add_child(united_cell); // add this as a child of ancestor

	reset_edges(united_cell); // since, we are uniting a and b, add new set of edges into the united cell

	a->parent->remove_child(a);
	b->parent->remove_child(b);
	
	cell* p1 = a->parent;
	cell* p2 = b->parent;
	int d1 = p1->depth;
	int d2 = p2->depth;

	while(d1 != d2) {

		if(d1 > d2) {
			reset_edges(p1);
			p1 = p1->parent;
			d1--;
		}

		else {
			reset_edges(p2);
			p2 = p2->parent;
			d2--;
		}
	}


	while(p1 != p2) {
		reset_edges(p1);
		reset_edges(p2);
		p1= p1->parent;
		p2= p2->parent;
	}

	reset_edges(ancestor);

	if (a->parent->children.size() == 1) {
		fuse(a->parent);
	}

	if (b->parent->children.size() == 1) {
		fuse(b->parent);
	}


	a->parent = united_cell;   	// splice a from tree, and 
	b->parent = united_cell;



	recompute_depths(a);
	recompute_depths(b);
	return united_cell;
}

