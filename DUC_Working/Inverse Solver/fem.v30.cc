/* This is FEM v30, seeded by fem.v23, started 18/05/2015
 * This will be a total field formulation with adiabatic absorbers (not used)
 * Takes -- . Mesh file specifications are:
 * Edge lists;
	 * RBC: List of edges on outermost boundary, marked 'r'
	 * OUT: List of edges on inner border of AA above substrate, marked 'o'
	 * DIE: List of edges on air-substrate interface, marked 'd' (not used)
	 * TOP: List of edges from which far field is calculated, marked 't'
 * Element lists;
	 * TXX: List of elements belonging to the tissue number XX, marked 's'
	 * AAB: List of elements belonging to the AA layer, marked 'a' (not used)
 * The entire geometry must be centred at (0,0) for get_mesh_extents to work
 * FORWARD SOLVER
 */

#include <iostream>
//#include <fftw3.h>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <complex>
#include <time.h>
#include <queue>
#include <iomanip>
#include <algorithm>

//#define SELDON_WITH_BLAS
//#define SELDON_WITH_LAPACK
#define SELDON_WITH_MUMPS
#define SELDON_DEBUG_LEVEL_3

#include "Seldon.hxx"
#include "SeldonSolver.hxx"

using namespace std;
using namespace Seldon;

//#include "MersenneTwister.h"

typedef complex<double> cmp;

//Constants

const double PI = 3.14159265358979323846;
const double c_0 = 299792458.0;// m/s
const double eps_0 = 8.854187817e-12 ;// F/m
const double mu_0 = 12.566370614e-7;// kg m/s^2 / C/s
const double Z0 = 376.730313461; //Ohm
const cmp I = cmp(0,1.0);

#include "fem.v5.h"

int main(int argc, char* argv[])
{
	//TRY;	
	cout<<"* ========="<<endl<<"2D FEM v30 : synthetic data generation"<<endl;
	
	//Correct usage
	if(argc != 2){
		cout<<"** Correct way to invoke program is <exe> <config_file> **"<<endl;
		return -1;
	}
	clock_t time0=clock();
	
	//stock variables
	int i,j,k,ii,jj,kk,status;	
	//Create simulation input variables that will be passed to the wrapper:
	double lambda=0;
	//input mesh file
	ifstream inmesh; 
	//incidence angles measures from +x axis
	vector<double> thetas;
	//vector for tissue dielectric constants
	vector<cmp> tissues;
	// output file 
	ofstream  bircs;
	//Call the wrapper to read input variables
	status = input_configuration(argc, argv, &lambda, &inmesh, &thetas, &tissues,&bircs);
	if(status !=0 ){
		cout<<"** Error in reading input configuration **"<<endl;return -1;
	}
	
	//Now prepare to create the datastructures for the FEM
	vector<node> nodelist; vector<edge> edgelist; vector<element> elementlist;
	
	vector<int> rbc_ed,diel_ed,diel_el,top_ed,top_el,aa_el,out_ed;
	vector< vector<int> >tiss_els;
	
	status = create_fem_datastructures(&inmesh,lambda,&nodelist,&edgelist,&elementlist,
			&rbc_ed,&diel_ed,&diel_el,&top_ed,&top_el,&aa_el,&out_ed,&tiss_els);
			
	if(status !=0 ){
		cout<<"** Error in creating FEM datastructures **"<<endl;return -1;
	}
	else{
		print_mesh_stats(&elementlist, &edgelist,lambda);
	}
	int No=nodelist.size()-1, Ne=elementlist.size()-1, Ns=edgelist.size()-1;
	
	//Now create an edge adjacency list. Each node has listed the edges it is connected to.
	//vector< vector <int> > EAL(No+1, vector<int>(0,0)); 
	//for(i=1 ; i<=Ns ; i++){
		//j = edgelist[i].ni; k = edgelist[i].nj;
		//EAL[j].push_back(i); EAL[k].push_back(i);
	//}
	
	//Set the dielectric constants
	status = set_tissue_eps(&tissues,&tiss_els,&elementlist);
	
	if(status !=0 ){
		return -1;
	}
	
	//Determine mesh extents and properties
	//the int variables contain node numbers of mesh extents
	//int tr_aa,tr_va,bl_va,bl_aa;

	//get_mesh_extents(&out_ed,&rbc_ed,&tr_va,&tr_aa,&bl_va,&bl_aa,lambda,&EAL,&nodelist,&edgelist,&elementlist);
	get_mesh_extents(&diel_ed,&top_ed,&rbc_ed,lambda,&nodelist,&edgelist);
					
	//Now create contour of edge #s for RCS calculation
	if(top_ed.size()==0){
		cout<<"** Error: Contour for RCS calculation not specified **"<<endl; return -1;
	}
	
	//Sort top_ed to produce contour, necessary if field output is required
	status = sort_list('p',&top_ed,&nodelist,&edgelist);
	if( status != 0 ){
		return -1;
	}		
	vector<int> contour = top_ed;
	
	cout<<"RCS contour has "<<contour.size()<<" edges at radius "<<
	nodelist[edgelist[top_ed[0]].ni].node_dist(node(0,0))/lambda<<" from center."<<endl;
	
	//Set the AA region loss
   //set_aa_loss(lambda,&tr_va,&tr_aa,&bl_va,&bl_aa,&aa_el,&nodelist,&elementlist);	
	
	clock_t time1=clock(); 
			
	/* Start of computation, AX = B. For matrix assembly:	 */
	//	Matrix<cmp, General, ArrayRowSparse> A0tmp(Ns, Ns), A1tmp(Ns,Ns);
	Matrix<cmp, Symmetric, ArrayRowSymSparse> ATM(Ns, Ns), ATE(Ns,Ns);
	ATM.Zero();ATE.Zero();
	
	//C is the RHS, soln is in B.
	vector< Vector<cmp> > BTM(thetas.size(), Vector<cmp>(Ns)), BTE(thetas.size(), Vector<cmp>(Ns)); 
	
	cout<<"Building matrix: ";cout.flush();
	//Total field approach
	for(j=0; j<thetas.size(); j++){		
		//Dense RHS vectors need to filled to be initialized
		BTM[j].Fill(0.0); BTE[j].Fill(0.0);	
	}
	for(k=1; k<=Ne; k++){
		set_element_matrices(k,&nodelist,&edgelist,&elementlist);
		set_matrix_entries_tf(k,lambda,thetas,&ATM,&ATE,&BTM,&BTE,&nodelist,&edgelist,&elementlist);
		if(k%(Ne/10)==1){//print ten dots as progress
			cout<<".";cout.flush();
		}
	}
	
	clock_t time2=clock(); 
	// Use MUMPS direct solver
	cout<<" Solving for theta=";
	for(j=0; j<thetas.size(); j++){	 
		cout<<" "<<rtod(thetas[j]);
	}
	cout<<". TM";cout.flush();
	MatrixMumps<cmp> mat_lu;
	GetLU(ATM, mat_lu);
	int info = mat_lu.GetInfoFactorization();
	cout<<"("<<info<<")";
	
	if(info !=0){
		cout<<"**Error in matrix solution**, time = "<<int((clock()-time2)/CLOCKS_PER_SEC)<<"s."<<endl;
		return -1;
	}
	
	for(j=0; j<thetas.size(); j++){	 
		SolveLU(mat_lu, BTM[j]); cout<<"+";cout.flush();
	}
	mat_lu.Clear(); ATM.Clear();	
	
/*	cout<<", TE";cout.flush();
	GetLU(ATE, mat_lu);
	info = mat_lu.GetInfoFactorization();
	cout<<"("<<info<<")";
	if(info !=0){
		cout<<"**Error in matrix solution**, time = "<<int((clock()-time2)/CLOCKS_PER_SEC)<<"s."<<endl;
		return -1;
	}		
	
	for(j=0; j<thetas.size(); j++){	 
		SolveLU(mat_lu, BTE[j]); cout<<"+";cout.flush();
	}
	mat_lu.Clear(); ATE.Clear();
*/
	clock_t time4=clock();
	
	//Calculations done, now output the fields (either exact fields, or farfields)
	
	//Compute bistatic fields:
	//rec is the angle of the radar receiver from +x axis
	//going from rec_i to rec_f in steps of rec_s (all in degrees)
	double rec,rec_s,rec_i,rec_f,inc_ang;
	rec_i = dtor(0); rec_f = dtor(359); rec_s = dtor(1);		
	vector<cmp> fval;	
	
	//Print a header line
	ofstream bircsfull("bircsfull.csv");
	bircsfull<<"rec,Ez(abs),Ez(arg),Ht(abs),Ht(arg),Hz(abs),Hz(arg),Et(abs),Et(arg)"<<endl;
	bircs<<"inc,rec,Ez(real),Ez(imag),Ht(real),Ht(imag),Hz(real),Hz(imag),Et(real),Et(imag)"<<endl;
	//calculate bistatic fields or rcs
	
	for(j=0; j<thetas.size(); j++){
		//Put in one line to distinguish measurements of different incident fields
		inc_ang = thetas[j];
		bircsfull<<rtod(inc_ang)<<endl;
				
	/*	rec = rec_i;
		while( rec <= rec_f){			
			//Last argument to field_tf is 1 if we want to store scat filed, or 0 if we want total field
			fval = field_tf(&nodelist,&edgelist,&elementlist,&contour,&BTM[j],&BTE[j],rec,thetas[j],lambda,1);
			//flag 1 for scattered field, 0 for total field
			bircs<<rtod(rec)<<","<<abs(fval[0])<<","<<arg(fval[0])<<","<<abs(fval[1])<<","<<
			arg(fval[1])<<","<< abs(fval[2])<<","<<arg(fval[2])<<","<< abs(fval[3])<<","<<arg(fval[3])<<endl;

			rec += rec_s;
		}*/
		
		//Use this to store data for inverse solver, only a few rec angles needed
		//special case where every transmitter is also a receiver
		for(i=0; i<thetas.size(); i++){
			rec = thetas[i];
			fval = field_tf(&nodelist,&edgelist,&elementlist,&contour,&BTM[j],&BTE[j],rec,inc_ang,lambda,0);
			bircs<<rtod(inc_ang)<<","<<rtod(rec)<<","<<abs(fval[0])<<","<<arg(fval[0])<<","<<abs(fval[1])<<","<<
			arg(fval[1])<<","<< abs(fval[2])<<","<<arg(fval[2])<<","<< abs(fval[3])<<","<<arg(fval[3])<<endl;
		}

		//Store the full data here
		for(i=0; i<contour.size(); i++){
			node pmid = edgelist[contour[i]].midpoint(&nodelist), pi = nodelist[edgelist[contour[i]].ni];
			double ang = atan2(pmid.ny,pmid.nx);
			//if( ang < 0){ang += 2*PI;}
			//flag 1 for scattered field, 0 for total field
			fval = field_tf_at_ed(&nodelist,&edgelist,&elementlist,&contour,&BTM[j],&BTE[j],contour[i],thetas[j],lambda,0);
			bircsfull<<rtod(ang)<<","<<abs(fval[0])<<","<<arg(fval[0])<<","<<abs(fval[1])<<","<<
			arg(fval[1])<<","<< abs(fval[2])<<","<<arg(fval[2])<<","<< abs(fval[3])<<","<<arg(fval[3]);//<<endl;
			//For diagnostic purposes, also printing the scattered field
			fval = field_tf_at_ed(&nodelist,&edgelist,&elementlist,&contour,&BTM[j],&BTE[j],contour[i],thetas[j],lambda,1);
			bircsfull<<","<<abs(fval[0])<<","<<arg(fval[0])<<","<<abs(fval[1])<<","<<
			arg(fval[1])<<","<< abs(fval[2])<<","<<arg(fval[2])<<","<< abs(fval[3])<<","<<arg(fval[3])<<endl;
		}
	}
	bircs.close();
	bircsfull.close();
	
	//print_mesh_contours(&nodelist,&edgelist,&elementlist);
	
	//End of everything
	clock_t timeT=clock();
	k = int((timeT-time0)/CLOCKS_PER_SEC);
	j = int((timeT-time1)/(CLOCKS_PER_SEC));
	cout<<"Time(%): Setup "<<int((time1-time0)/CLOCKS_PER_SEC)<<
 	"s, execution "<<(j/60)%60<<"m "<<j%60<<"s "<<
	". Total: "<<k/3600<<"h "<<(k/60)%60<<"m "<<k%60<<"s."<<endl<<
	"========= *"<<endl;

	return 0;
}
