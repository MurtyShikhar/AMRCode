#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <complex>
#include <algorithm>
using namespace std;


typedef pair<double,double> point;
typedef complex<double> dcomp;
typedef vector<vector<dcomp> > comp_matrix;
const int discretization = 8;
//const int init_children = 4;

template <typename T>
struct pointer_values_equal
{
    const T* to_find;
    pointer_values_equal(T* a): to_find(a) {}

    bool operator()(const T* other) const
    {
        return *to_find == *other;
    }
};



// vectors are slowww
struct edge {
	bool v; // v == 1 if vertical edge, 0 if horizontal edge
	point s,e; // s == starting point, e == ending point

	edge(point s,point e, bool v): s(s), e(e), v(v) {}

	double len() const{
		if (v) 
			return abs(s.second - e.second);
		else return abs(s.first - e.first);
	}

	bool operator<(const edge& other) const{
		if (v == other.v)
			return (s < other.s) || (s == other.s && e < other.e);
		else
			return v < other.v;
	}

	bool operator==(const edge& other) const{
		return v == other.v && s == other.s && e == other.e;
	}
};

bool overlap(double a_f, double a_s, double b_f, double b_s) {
	if (b_s <= a_f || a_s <= b_f)
		return false;

	return true;
}



// return true if other and this edge touch each other
bool touch(const edge& first,const edge& other) {
	if (first.v == other.v) {
		if (first.v) 
			return (first.s.first == other.s.first) && overlap(first.s.second, first.e.second, other.s.second, other.e.second);
		else 
			return (first.s.second == other.s.second) && overlap(first.s.first, first.e.first, other.s.first, other.e.first);
	}
	return false;
}




struct cell {
	vector<cell*> children; // children of the cell
	cell* parent; // unique parent of the cell
	pair<double,double> coord; // coordinate of the centre of the cell.. not required
	double area; // area of the cell. useful for aux_GS, aux_GD calculation
	dcomp chi; // chi of the cell
	bool flag; 
	vector<edge> edges; // edges of the bounding polygon
	int depth;

	bool operator==(const cell& other) const {
		return (coord == other.coord) && (area == other.area);
	}
	cell(){
		flag = false;
	} 
	void set_coords(double x,double y){
		coord.first =x;
		coord.second = y;
	}

	// returns true if this cell and cell "other" are neighbors
	bool neighbor(const cell& other) {
		int sz1 = edges.size();
		int sz2 = other.edges.size();
		for (int i = 0; i < sz1; i++) {
			for (int j = 0; j < sz2; j++) {
				if (touch(edges[i], other.edges[j]))			
					return true;
			}
		}
		return false;
	}

	double get_x() const{
		return coord.first;
	}

	double get_y() const{
		return coord.second;
	}

	void set_area(double d){
		area = d;
	}

	double get_area() {
		return area;
	}

	void add_edge(const edge& ed) {
		edges.push_back(ed);
	}

	void add_child(cell* child){
		children.push_back(child);
	}

	void remove_child(cell* child) {
		pointer_values_equal<cell> p(child);

		vector<cell* >::iterator position =  find_if(children.begin(), children.end(), p);
		if (position != children.end())
			children.erase(position);
	}

	int num_children() {
		return (int) children.size();
	}
};

void create_cell(cell* curr, double side,int level) {
	int dx[2] = {-1,1};
	int dy[2] = {-1,1};
	double c = side/4;
	curr->depth = level; 
	//int k = 1 << (2*level+1);
	double area = side*side/4;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			cell* child = new cell;
			double xc = curr->coord.first + dx[i]*c;
			double yc = curr->coord.second + dy[j]*c;
			child->set_coords(xc,yc);
			child->set_area(area);
			/* add 4 edges since  we basically have a square cell of side length 2c X 2c */
			child->add_edge(edge( make_pair(xc-c, yc-c), make_pair(xc-c,yc+c) ,1));
			child->add_edge(edge( make_pair(xc+c, yc-c), make_pair(xc+c,yc+c) ,1));
			child->add_edge(edge( make_pair(xc-c,yc+c),  make_pair(xc+c,yc+c) ,0));
			child->add_edge(edge( make_pair(xc-c,yc-c),  make_pair(xc+c,yc-c) ,0));

			if (level < discretization)
				create_cell(child,side/2,level+1);
			curr->add_child(child);
			child->parent = curr;
		}

	}

}

void create_tree(cell* root, double dom_side) {
	create_cell(root,dom_side,0);
	return;
}

void pretty_print(cell* root,int indent = 0) {
	if (root != NULL) {
		for (int i = 0; i < root->num_children(); i++)
			pretty_print(root->children[i], indent+ 4);

		if (indent) {
			std:: cout << std::setw(indent) << ' ';
		}

		cout << "h" << "\n ";
	}
}

// void plot_mesh(vector<cell*>& domain, pngwriter& png)
// {
// 	//pngwriter png(300, 300, 0,"domainpicture.png");
// 	//png.line(5, 5, 5, 50, 1.0, 1.0, 1.0);

// 	for (int i = 0; i < domain.size(); i++)
// 		for (int j= 0; j< domain[i]->edges.size(); ++j)
// 			png.line(domain[i]->edges[j].s.first, domain[i]->edges[j].s.second, domain[i]->edges[j].e.first, domain[i]->edges[j].e.second, 1.0, 1.0, 1.0);
// }
