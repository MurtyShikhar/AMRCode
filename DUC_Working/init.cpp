#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <cassert>
#include "tree.cpp"
#include <algorithm>
#include "physics.cpp"
#include "divide.cpp"
#include <pngwriter.h>
#include "unite.cpp"
#include <cstring>

using namespace std;
typedef complex<double> dcomp;
typedef vector<vector<dcomp> > comp_matrix;
double r_meas, side, diameter, tr_x, tr_y, bl_x, bl_y, lambda;
double width_figure = 600;
double width = width_figure - 50;
/*
void multiply(comp_matrix& A, comp_matrix& B, comp_matrix& C);
void pop_source_coords(vector<pair<double, double> >& source_coords, double diameter, int M);
void generate_aux_gd(comp_matrix& Aux_GD, vector<cell* >& cells);
void generate_aux_gs(comp_matrix& Aux_GS, vector<cell* >& cells, vector<pair<double, double> >& source_coords, comp_matrix& u_inc);
void back_propagation(comp_matrix& Aux_GD, comp_matrix& Aux_GS, comp_matrix& f, comp_matrix& u_inc, vector<cell* >& cells);
*/



void print(cell* c) {

  cout << "COORDS: "<<c->get_x() << c->get_y() << " ";
  for (int i = 0; i < c->edges.size(); i++) {

    cout << c->edges[i].s.first << " ---> " << c->edges[i].s.second << " ";
    cout << c->edges[i].e.first << " ---> " << c->edges[i].e.second << " ";

  }
  cout << "\n";
}


void populate(vector<cell*>& domain, cell* curr, int level, int curr_level) {

  if (curr_level == level) {
      //print(curr);
      domain.push_back(curr);
  }
  else {
      for (int i = 0; i < curr->children.size(); i++) {
        //cout << "HERE: ";
        //print(curr->children[i]);
        //cout << "\n";
        populate(domain, curr->children[i], level, curr_level+1);
      }
  }
}


dcomp parse_to_complex(string s) {
	size_t sz;
	double real_part = stod(s, &sz);
	double img = stod(s.substr(sz));
	return dcomp(real_part, img);	
}

point modify(point e) {
  double e1 = (width_figure-width)/2 + ((e.first)*width/side) + width/2;
  double e2 = (width_figure-width)/2 + ((e.second)*width/side) + width/2; 
  return make_pair(e1, e2);
}

void plot(vector<cell* > domain, pngwriter& png) {
  for (int i = 0; i < domain.size(); i++)
    for (int j= 0; j< domain[i]->edges.size(); ++j) {
      point fir = modify(domain[i]->edges[j].s);
      point sec = modify(domain[i]->edges[j].e);
      //cout<<fir.first<<" "<<fir.second<<" "<<sec.first<<" "<<sec.second<<"\n";

      png.line(fir.first, fir.second, sec.first, sec.second, 1.0, 1.0, 1.0);
    }
} 

void read_data(comp_matrix& f, vector<double>& thetas) {
  ifstream file("square.txt");
  string lam, out_file, theta_list;

  vector<vector<string> > f_string;
  while(file)
  {
    string s;
    if(!getline(file, s)) break;

    stringstream ss(s);
    //cout<<s<<endl;
    if(s[0] == '#')
      continue;

    if(s.compare(0,6,"lambda") == 0)
    {
      //cout<<s<<endl;
      size_t pos = s.find("=");
      lam = s.substr(pos+2);
      cout<<lam<<endl;
      lambda = stod(lam);
      cout<<"ho\n";
    }
    if(s.compare(0,10,"outputfile") == 0)
    {
      size_t pos = s.find("=");
      //cout<<s<<endl;
      out_file = s.substr(pos+2);
    }

    if(s.compare(0,5,"theta") == 0)
    {
      size_t pos = s.find("=");
      theta_list = s.substr(pos+2);
    }
  }
  file.close();

  //cout<<lam<<endl<<out_file<<endl<<theta_list<<endl;

  ifstream infile(out_file.c_str());
  getline(infile, lam);
  while (infile)
  {
    string s;
    if (!getline( infile, s )) break;

    stringstream ss(s);
    vector<string> record;

    while (ss)
    {
      string s;
      if (!getline( ss, s, ',' )) break;
        record.push_back((s));
    //  cout<<s<< " ";
     // cout << s[0] << " " << s[2] << "  " << s.length() << " ";
         
    }
    //cout<<endl;
    f_string.push_back(record);
  } 

  f.resize(M, vector<dcomp>(M));
  for (int i = 0; i < M-1; ++i)
  {
    for (int j = 0; i < M; ++j)
    {
      f[i][j] = stod(f_string[j+M*i][2]) + stod(f_string[j+M*i][3])*jota;
    }
  }
  if (!infile.eof())
  {
    cerr << "Fooey!\n";
  }
  infile.close();


  stringstream sss(theta_list);
    cout<<sss<<endl;
    cout<<"Hello\n";

  while (sss)
      {
        string t;
        if (!getline(sss, t, ',' )) break;

        thetas.push_back(stod(t));
        //cout<<s<<endl;
      }

  ifstream infile3("extents.csv");
  string s;
  getline(infile3, s);
  cout<<s<<endl;
  getline(infile3, s);
  stringstream ss(s);
  getline(ss, s, ',');
  r_meas = stod(s);
  getline(ss, s, ',');
  tr_x = stod(s);
  getline(ss, s, ',');
  tr_y = stod(s);
  getline(ss, s, ',');
  bl_x = stod(s);
  getline(ss, s, ',');
  bl_y = stod(s);
  infile3.close();
}


int main()
{
	comp_matrix f;
	vector<double> thetas;

  read_data(f, thetas);

 	//double lambda = 1;
  	
    //cout<<"dod\n";
  	M = thetas.size();
    diameter = r_meas;
    side = 2*tr_x;

    vector<pair<double, double> > source_coords(M);

    pop_source_coords(source_coords, diameter, M);
    //cout<<"side : "<<side<<"\n";
    cell* root = new cell;
    root->set_coords(0,0);
    root->set_area(side*side);
    root->add_edge(edge(make_pair(-side/2,-side/2), make_pair(-side/2,side/2), 1));
    root->add_edge(edge(make_pair(side/2,-side/2), make_pair(side/2,side/2), 1));
    root->add_edge(edge(make_pair(-side/2,-side/2), make_pair(side/2,-side/2), 0));
    root->add_edge(edge(make_pair(-side/2,side/2), make_pair(side/2,side/2), 0));

    create_tree(root, side);
    //cout<<"whoo\n";

    vector<cell* > domain;

    unite(root->children[0], root->children[2]);

    populate(domain, root, 1, 0);



    pngwriter png(width_figure, width_figure, 0, "test.png");
    plot(domain,png);
    png.close();


    /* Domain plotted */

    comp_matrix aux_GD(domain.size(), vector<dcomp>(domain.size()));
    comp_matrix aux_GS(M, vector<dcomp>(domain.size()));
    comp_matrix u_inc(domain.size(), vector<dcomp>(M));

    k_b = 2*pi;



    generate_aux_gd(aux_GD, domain);
    //cout<<domain.size()<<"\n";

    generate_u_inc(u_inc, domain, thetas);
    //cout<<domain.size()<<"\n";
    generate_aux_gs(aux_GS, domain, source_coords, u_inc);

    back_propagation(aux_GD, aux_GS, f, u_inc, domain);



    /* divide algorithm */
    vector<cell* > initial_domain(domain);
    //display_chi(domain);
    divide(domain,7, 3, thetas, f, source_coords);
    reconstruct_domain(domain, initial_domain);



    vector<dcomp> chi_val;

    dcomp max = domain[0]->chi;
    int max_index = 0;
    for (int i = 1; i < domain.size(); ++i)
    {
      if(real(domain[i]->chi) > real(max))
      {
        max = domain[i]->chi;
        max_index = i;
      }
    }

    pngwriter grad(width_figure, width_figure, 0, "chi5.png");

    for (int i = 0; i < domain.size(); ++i){
    	point lowleft = modify(make_pair(domain[i]->get_x() - sqrt(domain[i]->area)/2, domain[i]->get_y() - sqrt(domain[i]->area)/2));
    	point topright = modify(make_pair(domain[i]->get_x() + sqrt(domain[i]->area)/2, domain[i]->get_y() + sqrt(domain[i]->area)/2));        
   		double color_red = (1+(real(domain[i]->chi))/(real(max)))/2;
    	grad.filledsquare(lowleft.first, lowleft.second, topright.first, topright.second, color_red, 0.0, 0.0);     
    }
    grad.close();




    return 0;
}





