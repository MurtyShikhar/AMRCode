//Updated set_element_matrices (changed in orientation calculation) for fem.v22.c
//Functions
template <class T>
T Max(T a, T b){
	return (a>b?a:b);
}

template <class T>
T Min(T a, T b){
	return (a<b?a:b);
}

template <class T>
int sign(T a){
	return (a>=0?1:-1);
}

template <class T>
T sq(const T &a){
	return a*a;
}

template <class T> 
complex<T> iz(const complex<T>& z){ 
	return complex<T>(-z.imag(), z.real()); 
}

template<class T> 
complex<T> asin(const complex<T>& z)  { 
	return -I*(log(iz<T>(z) + sqrt(1.0 - z * z))); 
}

//this takes in any data type T, and returns it in string
template<class T>
string tos(T i)
{
    stringstream ss;
    string s;
    ss << i;
    s = ss.str();
    return s;
}
cmp polar(double r, cmp th)
{
	//normal polar fn, polar(r,th)=r e^{j th}=r(cos(th)+jsin(th)) where r,th are real.
	//s=r e^{j th}, th = a+jb, s= r e^{-b+ja} = r e^{-b}(cos(a)+jsin(a))=polar(r e^{-b},a)
	return polar(r*exp(-1.0*imag(th)),real(th));
}

double dtor(double d)
{
	return PI*d/180.0;
}

double rtod(double r)
{
	return 180.0*r/PI;
}

//Node class
class node{
	public:
	double nx,ny; //coordinates
	node(double,double);//constructor relevant for plain coordinate
	node(){};
	node& operator= (const node &);
	node operator- (const node &);
	bool operator== (const node &);
	double node_dist(const node &);
};

node::node(double x, double y){
	nx=x;ny=y;}

node& node::operator=(const node &rhs){ 
	if(this != &rhs)
	{
		nx=rhs.nx;ny=rhs.ny;
	}
	return *this;
}	

node node::operator-(const node &n){ 
	return node(nx-n.nx,ny-n.ny);
}	

bool node::operator==(const node &n1){
	if( (abs(nx-n1.nx) < 1e-5 ) && (abs(ny-n1.ny )< 1e-5) )
		return true;
	else
		return false;
}

double node::node_dist(const node &n1){
	return sqrt( sq(nx-n1.nx) + sq(ny-n1.ny));}

//Node related functions
node closest_node(node n0, const vector<node> *nlist)
{
	double min=n0.node_dist( (*nlist)[1]),tmp;
	int j=1;
	for(int i=2;i<(*nlist).size();i++)
	{
		tmp=n0.node_dist( (*nlist)[i]);
		if(tmp<min)
		{
			min=tmp;
			j=i;
		}
	}
	return (*nlist)[j];
}

//Return the signed area of a triangle of three nodes
double area_n(node * a, node * b, node * c)
{
	double ax=(*a).nx,ay=(*a).ny,
	bx=(*b).nx,by=(*b).ny,
	cx=(*c).nx,cy=(*c).ny;
	
	return 0.5*((bx*cy-by*cx) + ax*(by-cy) + ay*(cx-bx));
}

//Edge class
class edge{
public:
	int ni,nj; //node numbers of edge
	char edbtype;//p - PEC, r - RBC, m - metal, l - PML, d - dielectric
	double edlength;//edge length
	vector <int> els;//elements that share this edge
	edge(int,int);
	edge(int,int,const vector<node> *);
	edge(){};
	void setlength(const vector<node> *);
	bool operator== (const edge &);
	node midpoint(const vector<node> *);
};

edge::edge(int n1,int n2){
	ni=n1;nj=n2;}

edge::edge(int n1,int n2,const vector<node> * nlist){
	if(n1==0 || n2==0)
	{cout<<"**Node #0 requested -- problem.**";}
	ni=n1;nj=n2;
	
	double nix=(*nlist)[n1].nx,niy=(*nlist)[n1].ny,
	njx=(*nlist)[n2].nx,njy=(*nlist)[n2].ny;
	
	edlength=sqrt(sq(nix-njx)+sq(niy-njy));
}

void edge::setlength(const vector<node> * nlist){
	double nix=(*nlist)[ni].nx,niy=(*nlist)[ni].ny,
	njx=(*nlist)[nj].nx,njy=(*nlist)[nj].ny;
	edlength=sqrt(sq(nix-njx)+sq(niy-njy));
}

bool edge::operator== (const edge &e1){
	if(((ni==e1.ni)&&(nj==e1.nj))||((nj==e1.ni)&&(ni==e1.nj))) //safe because comparing integers
	{return true;}
	else
	{return false;}
}

node edge::midpoint(const vector<node> * nlist){
	node pi,pj,pk; pi = (*nlist)[ni]; pj = (*nlist)[nj];
	pk.nx = 0.5*(pi.nx+pj.nx); pk.ny = 0.5*(pi.ny+pj.ny);
	return pk;
}

//Edge related functions

// Return the edge number given an edge
int which_edge(edge *e,const vector<edge> *ed)
{
	for(int i=1;i<(*ed).size();i++){
		if((*e)==(*ed)[i])
		{return i;}
	}
	cerr<<"Error, returning edge #0"<<endl;
	return 0;
}

//Answers the question: Is this edge new?
bool new_edge(edge *e, vector<edge> *ed)
{
	for(int i=1;i<(*ed).size();i++){
		if((*e)==(*ed)[i])
		{return false;}
	}
	return true;
}

double node_edge_dist(node n, edge e, const vector<node> * nlist){
	node n1 = (*nlist)[e.ni], n2 = (*nlist)[e.nj];
	double ar = area_n(&n1,&n2,&n);
	return 2.0 * abs(ar) / e.edlength;
}

//Element class
class element{
public:
	//nodes and edges are stored in this order:
	//index:  0  1  2
	//nodes:  1  2  3
	//edges: 12 23 31 (with proper sign)
	vector<int> nodes;//node numbers from nodelist
	vector<int> edges;//edge numbers from edgelist. 
	char eltype;//boundary type: a - air, s - soil  
	double elarea;//element area
	cmp eleps;//element epsilon
	cmp elmu;//element mu
//	vector<cmp> elpml;//used only for PML
	vector<double> AA,BB,CC,LLen;//Computational details -- Save computation time at expense of memory
	//Functions
	element(){};
	element(int,int,int);
	element(int,int,int,const vector<node> *,const vector<edge> *);
	bool operator==(const element &);
	element& operator= (const element &);
	node elcenter(const vector<node> *);
	};

element::element(int n1,int n2,int n3){
	nodes.resize(3);edges.resize(3);
	nodes[0]=n1;nodes[1]=n2;nodes[2]=n3;
}

bool element::operator==(const element& e){
	int a=nodes[0],b=nodes[1],c=nodes[2];
	int aa=e.nodes[0],bb=e.nodes[1],cc=e.nodes[2];
	if( (a==aa && b==bb && c==cc)||
	   (a==bb && b==cc && c==aa)||
	   (a==cc && b==aa && c==bb)||
	   (a==aa && b==cc && c==bb)||
	   (a==bb && b==aa && c==cc)||
	   (a==cc && b==bb && c==aa) )
	{return true;}
	else 
	{return false;}
}

element& element::operator=(const element &rhs){ 
	if(this != &rhs)
	{
		nodes.resize(3);edges.resize(3); AA.resize(3); BB.resize(3); CC.resize(3);
		eleps = rhs.eleps; elmu = rhs.elmu; elarea = rhs.elarea; eltype = rhs.eltype;
		for(int i =0;i<3;i++)
		{
			nodes[i] = rhs.nodes[i]; edges[i] = rhs.edges[i]; 
			AA[i] = rhs.AA[i]; BB[i] = rhs.BB[i]; CC[i] = rhs.CC[i];
		}
	}
	return *this;
}	


node element::elcenter(const vector<node> * nodelist)
{
	node res;
	res.nx=((*nodelist)[nodes[0]].nx+(*nodelist)[nodes[1]].nx+(*nodelist)[nodes[2]].nx)/3.0;
	res.ny=((*nodelist)[nodes[0]].ny+(*nodelist)[nodes[1]].ny+(*nodelist)[nodes[2]].ny)/3.0;
	return res;
}
	
//Element related functions

//Return the element area
double area_e(element *e, vector<node> *nlist)
{
	double xi=(*nlist)[(*e).nodes[0]].nx,yi=(*nlist)[(*e).nodes[0]].ny,
	xj=(*nlist)[(*e).nodes[1]].nx,yj=(*nlist)[(*e).nodes[1]].ny,
	xk=(*nlist)[(*e).nodes[2]].nx,yk=(*nlist)[(*e).nodes[2]].ny,
	ar = 0.5 * ((xj*yk - xk*yj) - xi*(yk - yj) + yi*(xk - xj));
	return ar;
	//returning abs ar another waterloo
}

//Answers the question: Is this element new?
bool new_element(element *etmp,vector<element> * el){
	
	for(int i=1;i<(*el).size();i++){
		if((*etmp)==((*el)[i]))
		{return false;}
	}
	return true;
}

//Record class
class record{
public:
	int nn;
	node no;
	record(){};
	record(int,node);
};
record::record(int a, node b){
	nn = a;
	no = b;
}

//Layer class
class layer{
	public:
	double sa;
	double si;
	double cl;
	double th;
	double rb;
	layer(){};
	layer(double,double,double,double,double);
};
layer::layer(double a, double b, double c, double d, double e){
	sa = a; si = b; cl = c; th = d; rb = e;
}

//overloaded pow function for int, int, int. assuming e>0.
int pow(int b, int e)
{	
	if(e<0)
	{cout<<"Check, exponent less than 0"<<endl;}
	int i,p=1;
	for(i=0; i<e ; i++)
	{
		p *= b;
	}
	return p;
}

vector<double> integrate_lhs(element * el,vector <node> * nlist)
{
	vector<double> result(4,0.0);
	
	//Affine transform (x,y) -> (u,v) the given triangle to a rt. tr. and calculate the integral.
	//x = x1 + u*(x2-x1) + v*(x3-x1) = a + b*u +c*v
	//y = y1 + u*(y2-y1) + v*(y3-y1) = d + e*u +f*v
	// Functions in "result": {0 -> x, 1 -> y, 2 -> y^2, 3 -> x^2}
	// x^2=a^2 + b^2 u^2 + c^2 v^2 + 2ab u + 2ac v + 2bc uv
	// Over the unit triangle, Integral( x^q y^p) = q!r!/(q+r+2)!
	// So, Intg(x) becomes Intg(a+bu+cv) = a/2 + b/6 + c/6
	// Intg(y) = d/2 + e/6 + f/6
	// Intg(x^2) = a^2 /2 + b^2 /12 + c^2 /12 + ab/3 + ac/3 + bc /12
	// Intg(y^2) = d^2 /2 + e^2 /12 + f^2 /12 + ed/3 + fd/3 + ef /12
	double x1=(*nlist)[(*el).nodes[0]].nx,y1=(*nlist)[(*el).nodes[0]].ny,
	x2=(*nlist)[(*el).nodes[1]].nx, y2=(*nlist)[(*el).nodes[1]].ny,
	x3=(*nlist)[(*el).nodes[2]].nx, y3=(*nlist)[(*el).nodes[2]].ny;
	
	double a=x1, b=x2-x1, c=x3-x1, d=y1, e=y2-y1, f=y3-y1;
	double J = b * f - c * e; // This is the Jacobian
	result[0] = J * (a/2.0 + b/6.0 + c/6.0); //x
	result[1] = J * (d/2.0 + e/6.0 + f/6.0); //y 
	result[2] = J * (d*d /2.0 + e*e /12.0 + f*f /12.0 + d*e/3.0 + d*f/3.0 + e*f /12.0); //y^2
	result[3] = J * (a*a /2.0 + b*b /12.0 + c*c /12.0 + a*b/3.0 + a*c/3.0 + b*c /12.0); //x^2
	//above lines resulted in waterloo. almost.
	
	return result;
}



bool any_coupling(int ei, int ej,vector <edge> * edlist)	
{
	//Returns true if the two edges have any coupling
	int m,n;bool result;
	for( m=0;m<(*edlist)[ei].els.size();m++)
	{
		for( n=0;n<(*edlist)[ej].els.size();n++)
		{
			if( (*edlist)[ei].els[m] == (*edlist)[ej].els[n] )
			{return true;}
		}
	}
	return false;
}

double sec(double th)
{	
	return 1.0/cos(th);
}

double get_thorosG(double th, double L)
{
	//beam waist
	return L/10.0;//Max(L/4.00, lambda*6.0/pow(abs(sin(th)),1.5));
}

//for plane wave incident field: exp(-j (kx x + ky y))
cmp thoros(double x, double y, double th, double lambda){
	double K0 = 2.0* PI/lambda;
	//for plane wave
	return polar(1.0,-1.0*K0*(x*cos(th)+y*sin(th)));
}
//find partial x derivative: -j kx exp(-j (kx x + ky y))
cmp thoros_dx(double x, double y,double th, double lambda){	
	cmp T = thoros(x,y,th,lambda);
	double K0 = 2.0* PI/lambda;
	
	//for plane wave
	return -1.0*I*K0*cos(th) * T;
}
//find partial x derivative: -j ky exp(-j (kx x + ky y))
cmp thoros_dy(double x, double y, double th, double lambda){
	cmp T = thoros(x,y,th,lambda);
	double K0 = 2.0* PI/lambda;
	
	//for plane wave
	return -1.0*I*K0*sin(th) * T;
}

//For tapered incident field
cmp thoros(double x, double y, double g, double th, double y0,double lambda){
	
	//Center the beam about the soil center
	y -= y0;
	//change the theta to the one used in the formula
	th += PI/2.0;
	double K0 = 2.0* PI/lambda, w = (2.0*(sq((x+y*tan(th))/g) - 1.0) / sq(K0*g*cos(th)));
	
	//for plane wave
//	return polar(1.0,-1.0*K0*(x*sin(th)-y*cos(th)));
	
	//for tapered wave
	return exp(-1.0*sq((x+y*tan(th))/g)) * polar(1.0,-1.0*K0*(x*sin(th)-y*cos(th)) * (1.0 + w));
} 
cmp thoros_dx(double x, double y, double g, double th, double y0, double lambda){
	//find partial x derivative
	
	cmp T = thoros(x,y,g,th,y0,lambda);
	//Center the beam about the soil center
	y -= y0;
	//change the theta to the one used in the formula
	th += PI/2.0;
	double K0 = 2.0* PI/lambda;
	
	//for tapered wave
	return T*-1.0*
	(I*(2.0* sec(th)*(-2.0*x*y-(sq(g)-3.0*sq(x)+2.0*sq(y))*tan(th)+4.0*x*y*sq(tan(th))+sq(y*tan(th))*tan(th))
		+sq(g)*K0*(sq(g)*K0*sin(th)-2.0*I*(x+y*tan(th)))))/(g*g*g*g*K0);
} 
cmp thoros_dy(double x, double y, double g, double th, double y0, double lambda){
	//find partial x derivative
	
	cmp T = thoros(x,y,g,th,y0,lambda);
	//Center the beam about the soil center
	y -= y0;
	//change the theta to the one used in the formula
	th += PI/2.0;
	double K0 = 2.0* PI/lambda;
	
	//for tapered wave
 	return T*
	(I*(sq(g)*sq(g*K0)*cos(th) + 2.0*I*sq(g) *K0* tan(th)* (x + y* tan(th)) - 
		2.0* sec(th)* (sq(g) - sq(x) - 4.0* x* y* tan(th) + (2.0* sq(x) - 3.0*sq(y))*sq(tan(th)) + 
					   2.0* x* y* sq(tan(th))*tan(th))))/(g*g*g*g*K0);
}



//find the closest edge to a given point from a given list
int closest_edge(double x, double y, vector<int> * surf, vector<edge> * edlist, vector<node> * nlist)
{
	node p(x,y), ptmp = (*edlist)[(*surf)[0]].midpoint(nlist);
	double tmp = p.node_dist(ptmp);int i,j=0;
	
	for(i=1;i<(*surf).size();i++)
	{
		ptmp = (*edlist)[(*surf)[i]].midpoint(nlist);
		if( p.node_dist(ptmp) < tmp )
		{j=i;tmp = p.node_dist(ptmp);}
	}
	return j; // return vector index
}

//find the closest node to a given node from a given edge list
int closest_node(double x, double y, vector<int> * surf, vector<edge> * edlist, vector<node> * nlist)
{
	node p(x,y), ptmp = (*edlist)[(*surf)[0]].midpoint(nlist);
	double tmp = p.node_dist(ptmp);int i,j=0;
	
	for(i=1;i<(*surf).size();i++)
	{
		ptmp = (*edlist)[(*surf)[i]].midpoint(nlist);
		if( p.node_dist(ptmp) < tmp )
		{j=i;tmp = p.node_dist(ptmp);}
	}
	node ni = (*nlist)[(*edlist)[(*surf)[j]].ni], nj = (*nlist)[(*edlist)[(*surf)[j]].nj];

	//return node number
	if(p.node_dist(ni) < p.node_dist(nj) ){return (*edlist)[(*surf)[j]].ni;}
	else{ return (*edlist)[(*surf)[j]].nj; }
}

//find the closest edge to a given point
int closest_edge(double x, double y, vector<edge> * edlist, vector<node> * nlist)
{
	node p(x,y), ptmp = (*edlist)[1].midpoint(nlist);
	double tmp = p.node_dist(ptmp); int i,j=0;
	
	for(i=2;i<(*edlist).size();i++)
	{
		ptmp = (*edlist)[i].midpoint(nlist);
		if( p.node_dist(ptmp) < tmp )
		{j=i; tmp = p.node_dist(ptmp);}
	}
	return j; // return edge number
}

//find closest node to a given node
int closest_node(node n, vector<node> * nlist)
{
	double tmp = n.node_dist((*nlist)[1]); int i,j=1;
	for(i = 2; i < (*nlist).size(); i++)
	{
		if(n.node_dist((*nlist)[i]) < tmp)
		{
			tmp = n.node_dist((*nlist)[i]); j = i;
		}
	}
	return j;
}				  
//find the closest node to a given node from a given edge list

//find closest node to a given node	with a seed starting guess
//Not convinced of its optimality

int closest_node(node n, int s, vector<vector<int> > * val,vector<node> * nlist)
{
	double tmp, path; int i,j,k; node pi,pk;
	bool search = true; vector<int> npath(0);
	
	k = s;
	while(search)
	{
		npath.push_back(k);
		pk = (*nlist)[(*val)[k][0]];
		j = (*val)[k][1]; pi = (*nlist)[j];
		tmp = pk.node_dist(pi) + pi.node_dist(n);
		for(i=2;i<(*val)[k].size();i++)
		{
			pi = (*nlist)[(*val)[k][i]];
			path = pk.node_dist(pi) + pi.node_dist(n);
			if(path < tmp)
			{
				j = (*val)[k][i]; tmp = path;
			}
		}
		for(i=Max(0,int(npath.size()-2));i<npath.size();i++)//scan for oscillation
		{
			if(npath[i] == j)
				search = false;
		}
		k = j; pk = (*nlist)[(*val)[k][0]];

//		cout<<k<<","<<pk.nx<<","<<pk.ny<<","<<pk.node_dist(n)<<endl;
	}
	
	return npath[npath.size()-1];
}

vector<int> get_contour(vector<node> *nlist, vector<edge> *edlist, vector<vector<int> > * val,
						node pi, node pj,vector<int> * slist)
{
	/* This function produces a vector<int> list of edge numbers that represent a linear path from
	 node pi to node pj. */
	
	int nodei, edgei, nodej,edgej, i, j, k; node pk; double tmp=0,rms=0,mu=0; 
	vector<int> contour(0,0);
	
	/*	edgei = closest_edge(pi.nx,pi.ny,slist,edlist,nlist);
	 nodei = (*edlist)[(*slist)[edgei]].ni; i = (*edlist)[(*slist)[edgei]].nj;
	 if( (*nlist)[nodei].node_dist(pi) < (*nlist)[i].node_dist(pi))
	 {pi = (*nlist)[nodei];}
	 else
	 {pi = (*nlist)[i]; nodei = i;}
	 
	 edgej = closest_edge(pj.nx,pj.ny,slist,edlist,nlist);
	 nodej = (*edlist)[(*slist)[edgej]].ni; j = (*edlist)[(*slist)[edgej]].nj;
	 if( (*nlist)[nodej].node_dist(pj) < (*nlist)[j].node_dist(pj))
	 {pj = (*nlist)[nodej];}
	 else
	 {pj = (*nlist)[j]; nodej = j;} */
	
	//Brute force technique
	//	nodei = closest_node(pi,nlist); nodej = closest_node(pj,nlist);
	
	//clever technique
	i = closest_edge(pi.nx,pi.ny,slist,edlist,nlist); k = (*edlist)[i].ni;
	nodei = closest_node(pi,k,val,nlist); pi = (*nlist)[nodei];
	
	i = closest_edge(pj.nx,pj.ny,slist,edlist,nlist); k = (*edlist)[i].ni;
	nodej = closest_node(pj,k,val,nlist); pj = (*nlist)[nodej];
	
	edge e(nodei,nodej,nlist);
	cout<<"Contour: ("<<pi.nx<<","<<pi.ny<<") -> ("<<pj.nx<<","<<pj.ny<<"). ";
	
	k = nodei;
	while( k != nodej ) //traverse till pi reaches pj
	{
		//This algo maintains least distance from edge e
		/*		tmp = e.edlength; //initialize
		 for(i=1 ; i<(*val)[nodei].size(); i++) //run over nodes connected to nodei to find next nearest node
		 {
		 pi = (*nlist)[(*val)[nodei][i]];
		 if( pi.node_dist(pj) < (*nlist)[(*val)[nodei][0]].node_dist(pj) && node_edge_dist(pi,e,nlist) < tmp )
		 {
		 k = (*val)[nodei][i]; pk = (*nlist)[k]; tmp = node_edge_dist(pk,e,nlist);
		 }
		 } */
		
		//This algo is similar to A*
		k = (*val)[nodei][1]; pk = (*nlist)[k]; tmp = (*nlist)[nodei].node_dist(pk) + pk.node_dist(pj);
		
		for(i=2 ; i<(*val)[nodei].size(); i++) //run over nodes connected to nodei to find next nearest node
		{
			pi = (*nlist)[(*val)[nodei][i]];
			if( (*nlist)[nodei].node_dist(pi) + pi.node_dist(pj) < tmp )
			{
				k = (*val)[nodei][i]; pk = (*nlist)[k]; tmp = (*nlist)[nodei].node_dist(pi) + pi.node_dist(pj);
			}
		}
		//New node found
		edge ed(nodei,k); //make temp new edge
		edgei = which_edge(&ed,edlist); //get the edge number				
		nodei = k;  //update to the newly found node
		
		contour.push_back(edgei);
		mu += ed.midpoint(nlist).ny; rms += sq(ed.midpoint(nlist).ny);
		//		cout<<"("<<ed.midpoint(nlist).nx<<","<<ed.midpoint(nlist).ny<<") ";
		//		cout<<"("<<pk.nx<<","<<pk.ny<<") "<<pk.node_dist(pj)<<" "<<k<<endl; 
	}
	mu /= contour.size(); rms = sqrt(abs(rms/contour.size() - sq(mu)));
	cout<<contour.size()<<" edges. mu,rms = "<<mu<<","<<rms<<endl;
	//	cerr<<mu<<",";
	
	return contour;
}


//this returns the mu and stdev for a set of points
vector<double> stats(vector<double> * x){
	vector<double> R(2,0);
	
	for(int i=0; i<(*x).size(); i++)
	{
		R[0] += (*x)[i]; R[1] +=  sq((*x)[i]);
	}
	R[0] /= (*x).size(); R[1] = sqrt(R[1]/(*x).size() - sq(R[0]));
	
	return R;
}

//this returns the mu and stdev for a surface created by linearly interpolating
//a set of points, padded by one point on either side, each equal to the mean
vector<double> stats_surf(vector<double> * y){
	vector<double> R(2,0);
	int i, N  = (*y).size();
	
	for(i=0; i<N; i++){
		R[0] += (*y)[i];
	}
	R[0] /= N; 

	//create a zero mean set first
	vector<double> v(N,0);
	
	for(i=0; i<N; i++){
		v[i] = (*y)[i] - R[0];
	}
	
	//Now get the stdev, s^2 = 1/N( 2/3 * Sum_1^N vi^2 + 1/3 * Sum_1^N-1 vi * vi+1)
	for(i=0; i<N-1; i++){
		R[1] += 2.0 * sq(v[i]) + v[i]*v[i+1];
	}
	R[1] += 2.0 * sq(v[N-1]);
	
	R[1] = sqrt( R[1]/(3.0*N) );
	
	return R;
}


double get_correlation(vector<double> * z,int N, double dx, double l, int scorl)
{
	int i,j;
	vector<double> corr(int(5.0*(l/dx)),0), st = stats_surf(z);
	if( N> (*z).size() ){cerr<<"Asking for impossible correlation computation"<<endl; return 0;}
		
	for(i=0; i<corr.size(); i++)
	{
		for(j=0;j<N-i; j++)
		{
			corr[i] += ((*z)[j] - st[0]) * ((*z)[j+i] - st[0]);
		}
	}
	
	//estimate the correlation length
	double est_l;
	for(i=0;i<corr.size();i++)
	{
		if( abs(corr[i]/corr[0]) < exp(-1.0) )
			break;
	}
	//this means corr[i] = corr[0] exp(-(i*dx/l)^2)
	//i.e. 
/*	est_l = double(i-1) * dx /sqrt( log(corr[0]/corr[i-1]));
	cout<<"Internal variation in l: ("<<i-1<<") "<<est_l;
	est_l = double(i) * dx /sqrt( log(corr[0]/corr[i]));
	cout<<", ("<<i<<") "<<est_l<<". Correlation size: ("<<corr.size()<<")"<<endl; */
	
	//Found an unusual case where corr[i] is negative. Still don't understand why
	//maybe it is possible in small sample sizes
//	cout<<"["<<corr[0]<<","<<corr[i]<<"]";
	if(corr[i]<0 || isnan(corr[i]) ){
		return 0;
	}
	else{
		if(scorl == 0){
			return double(i) * dx /sqrt( log(corr[0]/corr[i]));//gaussian
		}
		else{
			return double(i) * dx /log(abs(corr[0]/corr[i])); //exponential
		}
	}
}

int other_node(int n, edge e){
	//Here n is the node, e is the edge to which it belongs
	//and we want to find the other node number
	if( e.ni == n ){
		return e.nj;
	}
	else{
		return e.ni;
	}
}

double det(node u, node v){
	//det(u,v) = ux vy - uy vx
	return u.nx * v.ny - u.ny *v.nx;
}

//returns true if node v is inside a triangle defined by n0,n1,n2.
bool in_triangle(node n0, node n1, node n2, node v){
	node v1 = n1-n0, v2 = n2-n0;
	//express tst = v0 + a* v1 + b* v2, where v1 and v2 are vectors from n0 to n1 and n2.
	//if a,b>0 and a+b <=1, tst is in interior
	double a = (det(v,v2)-det(n0,v2))/det(v1,v2), b = (det(n0,v1)-det(v,v1))/det(v1,v2), tol=1e-9;
	if( a<0 || b<0 || (a+b)>1.0){
		return false;
	}
	else{
		return true;
	}
}

bool interiorpt(vector <node> * nlist, node v, element e){
	node n0 = (*nlist)[e.nodes[0]], n1 = (*nlist)[e.nodes[1]], n2 = (*nlist)[e.nodes[2]];
	return in_triangle(n0,n1,n2,v);
}

//Brute force closest node find
int closest_node(node exct, vector<int> * itf, vector<node> * nlist,vector<element> * elist){ 
	int i,j,n,k=0; double dist = tan(dtor(89.9999)), tmp;
	for(i=0; i<(*itf).size(); i++){
		for(j=0; j<3; j++){
			n = (*elist)[(*itf)[i]].nodes[j];
			tmp = (*nlist)[n].node_dist(exct);
			if( tmp < dist){
				dist = tmp;
				k = n;
			}
		}
	}
	return k;
}			

//find closest node to a given node	with a seed starting guess
//this is a greedy best first search and thus fast.
int closest_node(node exct, int s, vector<vector<int> > * eal,vector<node> * nlist,vector<edge> * edlist){
	
	int i,j,k,newn, ed,oldn = s;
	bool searching = true;
	double sx = tan(dtor(89.9999)), npathl = sx, opathl = sx, tmpl = sx; //sx is some large number for init
	node pi, pj; vector<int> edpath;
	
	while(searching){
		//pi is starting node
		pi = (*nlist)[oldn]; opathl = pi.node_dist(exct);
//		cout<<"Starting with node "<<oldn<<" ("<<pi.nx<<","<<pi.ny<<"). "<<opathl<<endl;
		for(i=0; i<(*eal)[oldn].size(); i++){
			//pj is a node connected to pi
			j = other_node(oldn,(*edlist)[(*eal)[oldn][i]]); pj = (*nlist)[j];
			tmpl = pi.node_dist(pj) + pj.node_dist(exct);
			//within these nodes, update to a closest node
			if( tmpl < npathl ){
				newn = j; npathl = tmpl; 
				ed = (*eal)[oldn][i];
//				cout<<"Candidate node "<<j<<" might be closer ("<<pj.nx<<","<<pj.ny<<"). "
//				<<pj.node_dist(exct)<<endl;
			}
		}
		//add this node if it is closer
		if( (*nlist)[newn].node_dist(exct) < opathl){
			oldn = newn; npathl = sx;
			edpath.push_back(ed); 
//			cout<<"Candidate "<<newn<<" was closer, added edge "<<ed<<endl;
		}
		else{
			//searching is over, 
			searching = false;
		}
	}
	
	//See if the guess can be refined in an nbd of 3
	opathl = (*nlist)[oldn].node_dist(exct); 
	vector<int> ngbs; int ni,nj,nk,rguess = oldn;

	//Collect all nodes in this nbd
	for(i=0; i<(*eal)[oldn].size();i++){
		// ni are the ngbs of oldn
		ni = other_node(oldn,(*edlist)[(*eal)[oldn][i]]);
		for(j=0; j<(*eal)[ni].size(); j++){
			// nj are the ngbs of ni
			nj = other_node(ni,(*edlist)[(*eal)[ni][j]]);
			for(k=0; k<(*eal)[nj].size(); k++){
				//nk are the ngbs of nj
				nk = other_node(nj,(*edlist)[(*eal)[nj][k]]);
				ngbs.push_back(nk);
			}
			ngbs.push_back(nj);
		}
		ngbs.push_back(ni);
	}
	
	//See if anyone is closer
	for(i=0; i<ngbs.size(); i++){
		pi = (*nlist)[ngbs[i]];
		npathl = pi.node_dist(exct);
		if( npathl < opathl){
			rguess = ngbs[i]; npathl = opathl;
		}
	}
			
/*	if(rguess != oldn){
		cout<<"Closer node"<<endl;
	} */
	return rguess;
}

vector<int> greedy_path(int s, int end,  vector<vector<int> > * eal,vector<node> * nlist,vector<edge> * edlist){
	
	int i,j,newn, ed,oldn = s;
	bool searching = true;
	double sx = tan(dtor(89.9999)), npathl = sx, opathl = sx, tmpl = sx; //sx is some large number for init
	node pi, pj, exct = (*nlist)[end]; vector<int> edpath, npath;
	
	while(searching){
		//pi is starting node
		pi = (*nlist)[oldn]; opathl = pi.node_dist(exct);
		npath.push_back(oldn);
		//		cout<<"Starting with node "<<oldn<<" ("<<pi.nx<<","<<pi.ny<<"). "<<opathl<<endl;
		for(i=0; i<(*eal)[oldn].size(); i++){
			//pj is a node connected to pi
			j = other_node(oldn,(*edlist)[(*eal)[oldn][i]]); pj = (*nlist)[j];
			tmpl = pi.node_dist(pj) + pj.node_dist(exct);
			//within these nodes, update to a closest node
			if( tmpl < npathl ){
				newn = j; npathl = tmpl; 
				ed = (*eal)[oldn][i];
				//				cout<<"Candidate node "<<j<<" might be closer ("<<pj.nx<<","<<pj.ny<<"). "
				//				<<pj.node_dist(exct)<<endl;
			}
		}
		//add this node if it is closer
		if( (*nlist)[newn].node_dist(exct) < opathl){
			oldn = newn; npathl = sx;
			edpath.push_back(ed); 
			//			cout<<"Candidate "<<newn<<" was closer, added edge "<<ed<<endl;
		}
		else{
			//searching is over, 
			searching = false;
		}
	}
	/* //Print greedy path
	double len=0;
	cout<<"Greedy path is ";
	for(i=0;i<npath.size();i++){
		cout<<npath[npath.size()-1-i]<<" ";
		pi = (*nlist)[npath[npath.size()-1-i]];
		if(i!=0){
			len += (*nlist)[npath[i]].node_dist( (*nlist)[npath[i-1]] );
		}
 //		cerr<<pi.nx<<","<<pi.ny<<endl;
	}
	cout<<endl<<"Length is "<<len<<endl; */
	
	return edpath;
}

//This class stores the four currents for a given dl-segment along
//the rcs contour. A vector of this will contain the full field for
//a given monte carlo instance. Currents, M = E x n, J = n x H
//First two for TM, next two for TE
class currents{
	public:
	cmp Mt;//EzT for TM
	cmp Jz;//HtT for TM
	cmp Jt;//-HzT for TE
	cmp Mz;//-EtT for TE
	currents(){
		Mt = 0.0; Jz = 0.0; Jt = 0.0; Mz = 0.0;
	}
	currents(cmp a, cmp b, cmp c, cmp d){
		Mt = a; Jz = b; Jt = c; Mz = d;
	}
};

class card{
public:
	int nn;//the node number
	int pn;//the parent number
	int ed;//the edge from parent
	double g;//g cost
	double h;//heuristic cost
	card(){
		nn = 0; pn = 0; ed = 0; g = 0; h = 0;
	}
	card(int a, int b, int e, double c, double d){
		nn = a; pn = b; ed = e; g = c; h = d;
	}
	card(int a, int b, int e, double c, int goal, vector <node> * nlist){
		nn = a; pn = b; ed = e; g = c; h = (*nlist)[nn].node_dist( (*nlist)[goal] );
	}
	bool operator()(card& a, card& b){
		if( a.g+a.h > b.g+b.h ){return true;} // This is because lower h has higher priority
		else{return false;}
	}
	card& operator=(const card &rhs){
		if( this != &rhs){
			nn = rhs.nn; pn = rhs.pn; ed = rhs.ed; g = rhs.g; h = rhs.h;
		}
		return *this;
	}
};

//this is a shortest path algo based on A*
vector<int> shortest_path(int sta, int end, vector< vector <int> > *eal, vector<node> *nlist, vector<edge> *edlist){
	
	int i,j,ngb;
	vector<int> edpath(0);
	
	if(sta == end){return edpath;}
	
	//open and closed set creation
	vector<card> openset;
	vector<card> closedset;
	//starting node
	card firstnode(sta,sta,0,0,end,nlist); openset.push_back(firstnode);
	card current = firstnode;
	
	//find the shortest path
 	while( openset.size() != 0 ){
		current = openset.front(); //choose max priority for current
		openset.front() = openset.back(); //save the least priority card to front
		openset.pop_back(); //now delete the last entry
		closedset.push_back(current); //add current node to closed set
		
		if( current.nn == end ){
			break;
		}
		
		for(i = 0; i< (*eal)[current.nn].size(); i++){
//			double cost = current.g + (*edlist)[(*eal)[current.nn][i]].edlength;
			ngb = other_node(current.nn, (*edlist)[(*eal)[current.nn][i]]);
			double cost = current.g + (*nlist)[current.nn].node_dist( (*nlist)[ngb] ) ;

			bool in_openset = false, in_closedset = false;
			for(j=0; j<openset.size(); j++){
				if( ngb == openset[j].nn && cost < openset[j].g ){
					card update(ngb,current.nn,(*eal)[current.nn][i],cost,end,nlist);
					openset[j] = update;
					in_openset = true; break;
				}
			}
			for(j=0; j<closedset.size(); j++){
				if( ngb == closedset[j].nn ){
					in_closedset = true; break;
				}
			}
			if(!in_openset && !in_closedset){
				card nextone(ngb,current.nn,(*eal)[current.nn][i],cost,end,nlist);
				openset.push_back(nextone);
			}
		}
		make_heap(openset.begin(), openset.end(), card());// do the sorting only once
	}
	
	/* //Print optimal path
	card retrace = current;
	cout<<"Optimal path is "<<retrace.nn<<" ";
	node pi = (*nlist)[retrace.nn]; double len = 0;
 	cerr<<pi.nx<<","<<pi.ny<<endl;
	while(retrace.nn != sta){
		for(i=0; i< closedset.size(); i++){
			if(retrace.pn == closedset[i].nn){
				len += (*nlist)[closedset[i].nn].node_dist( (*nlist)[retrace.nn] );
				retrace = closedset[i]; cout<<retrace.nn<<" ";
				pi = (*nlist)[retrace.nn];
				cerr<<pi.nx<<","<<pi.ny<<endl;
				break;
			}
		}
	}
	cout<<endl<<"Length is "<<len<<endl; */
	
	// Save optimal path in terms of edges.
	card retrace = current; 
	while(retrace.nn != sta){
		for(i=0; i< closedset.size(); i++){
			if(retrace.pn == closedset[i].nn){
				edpath.push_back( retrace.ed );
				retrace = closedset[i];
				break;
			}
		}
	}
	return edpath;
}


vector<int> gen_smooth_segment(int oldn,int newn,vector< vector <int> > *eal, vector<node> *nlist, vector<edge> *edlist){
	//first generate a path, then see if there is more than 1 edge in it
	//if so, we smoothen it
	vector<int> segment = shortest_path(oldn,newn,eal,nlist,edlist);
	
	if(segment.size() == 1){
		return segment;
	}
}


void restore_record(vector<record> * rc, vector<node> *nlist){
	for(int i=0; i<(*rc).size(); i++){
		(*nlist)[(*rc)[i].nn] = (*rc)[i].no;
	}
}

//This function returns how many lines are there till
//the first delimiter is encountered
int linestill_delim(ifstream * in, string p){
	int i = 0;
	string line;
	if( (*in).is_open() ){
		(*in).seekg(0);(*in).clear();
		while( (*in).good() ){
			getline((*in),line);
			if( string(line) == p ){
				break;
			}
			else{
				i++;
			}
		}
		(*in).seekg(0);(*in).clear();
	}
	return i;
}

//Takes a list containing duplicates and creates a unique sorted list
//Note: original list is altered.
void sort_unq_list(vector<int> *inp){
	sort((*inp).begin(), (*inp).end());
	vector<int>::iterator it = unique((*inp).begin(),(*inp).end());
	(*inp).resize(it - (*inp).begin());
}

//Return sorted and non-duplicate list of nodes from a list of edges
vector<int> sort_edges(vector<int> * inp, vector<edge> * edgelist){
	vector<int> out(0);
	for(int i=0; i<(*inp).size(); i++){
		out.push_back((*edgelist)[(*inp)[i]].ni);
		out.push_back((*edgelist)[(*inp)[i]].nj);
	}
	sort_unq_list(&out);
	return out;
}

//Return sorted and non-duplicate list of nodes from a list of elements
vector<int> sort_elements(vector<int> * list, vector<element> * elist){
	vector<int> out(0);
	for(int i=0; i<(*list).size(); i++){
		for(int j=0; j<3; j++){
			out.push_back((*elist)[(*list)[i]].nodes[j]);
		}
	}
	sort_unq_list(&out);
	return out;
}

vector<cmp> print_anal(double theta,cmp eps)
{
	cmp eta, R_te, R_tm; double k;
	eta = sqrt(eps);
	k = sin(abs(theta)) / sin(acos(cos(abs(theta))/real(sqrt(eps))));
	R_te = (eta * k - 1.0)/(eta * k + 1.0);
	R_tm = (k - eta)/(k + eta);
	
//	cout<<"Analytic R_TM = ("<<abs(R_tm)<<","<<arg(R_tm)<<"), R_TE = ("<<abs(R_te)<<","<<arg(R_te)<<")"<<endl;
//	cout<<"Analytic R_TM = "<<sq(abs(R_tm))<<", R_TE = "<<sq(abs(R_te))<<endl;
	vector <cmp> R(2,0); R[0] = R_tm; R[1] = R_te; return R;
}

//This returns the bistatic RCS for a total field approach
vector<double> radar_tf(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, double soilw, double rec,double theta, double G,bool output_currents,
	double lambda)
{
	int i,j,k, N = (*ff).size(),ortn,sgn;
	node pi,pj,pm;	
	double psi,phi_i,phi_j,tphi,nphi,rp,dl,min_x=0,max_x=0,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda, g=K0*cos(theta),h=K0*sin(theta),ar,li,H_min=1,E_min=Z0; 
	vector <double> result(4,0);		
	cmp Mx,My,Mz,Jx,Jy,Jz,EtT,Ez,HtT,Hz,EzT,HzT,Hzs,Ets,Ezs,Hts,Mt,Jt,inc,er,mr; 
	vector <cmp> rval(2,cmp(0,0));
	
	//Assuming that the contour is always above the origin.
	//ortn stores match between contour orientation (c.ckwise) and edge orientation (i->j)
	
	//	ofstream cut; cut.open("cut.csv");
	ofstream opc;
	if(output_currents){
		opc.open("currents.v12.csv");
	}
	
	for(i=0; i<N; i++)
	{
		pm = (*edlist)[(*ff)[i]].midpoint(nlist);
		pi = (*nlist)[(*edlist)[(*ff)[i]].ni]; pj = (*nlist)[(*edlist)[(*ff)[i]].nj];
		//angles w.r.t center
		phi_i = atan2(pi.ny,pi.nx); phi_j = atan2(pj.ny,pj.nx);
		
		//match between contour and edge conv. 
		ortn = sign(phi_j - phi_i);
		if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
		{ortn *= -1;}
		
		//physical angle of tangent and outward normal (c.ckwise)
		tphi = atan2(double(ortn) * (pj.ny-pi.ny), double(ortn)*(pj.nx-pi.nx));
		nphi = tphi - PI/2.0; 
		tx = cos(tphi); ty = sin(tphi); nx = cos(nphi); ny = sin(nphi);
		
		//first get the total field 
		HtT = (*A)((*ff)[i]-1) * double(ortn); //H_min = Min(H_min,abs(HtT));
		EtT = (*B)((*ff)[i]-1) * double(ortn); //E_min = Min(E_min,abs(EtT));
		
		EzT = 0.0; HzT = 0.0;
		for(j=0; j<2; j++) //total field is avg over two elements
		{
			Ez = 0.0; Hz = 0.0;
			element ei = (*elist)[(*edlist)[(*ff)[i]].els[j]];
			ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
			
			for(k=0; k<3; k++)
			{
				sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
				li = (*edlist)[ei.edges[k]].edlength;
				Ez += double(sgn) * li * (*A)(ei.edges[k]-1);
				Hz += double(sgn) * li * (*B)(ei.edges[k]-1);
			}
			Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
			EzT += Ez; HzT += Hz;
		}
		EzT /= 2.0; HzT /= 2.0;
		
		//total currents, M = E x n, J = n x H
		Mt = EzT; 
		Jz = HtT;
		Jt = -1.0 * HzT; 
		Mz = -1.0 * EtT;
		
		//print if required
		if(output_currents){
			opc<<pm.nx<<","<<pm.ny<<","<<theta<<","<<abs(Mt)<<","<<arg(Mt)<<","<<abs(Jt)<<","<<arg(Jt)<<endl;
		}
		
		//rcs integrand needs these
		psi = atan2(pm.ny,pm.nx);
		rp = sqrt(sq(pm.ny)+sq(pm.nx));
		dl = sqrt(sq(pi.nx-pj.nx)+sq(pi.ny-pj.ny));
		
		//rcs
		rval[0] += (Mt * sin(tphi-rec) - Jz *Z0) * polar(dl,K0*rp*cos(psi-rec));
		rval[1] += (-1.0 * Jt * sin(tphi-rec) - Mz /Z0) * polar(dl,K0*rp*cos(psi-rec));
		
		//test
		//		cout<<phi_i<<" "<<phi_j<<" "<<ortn<<" "<<rval[0]<<" "<<rval[1]<<endl;
		//		cout<<abs(Jz)<<","<<abs(Mx)<<","<<abs(Mz)<<","<<abs(Jx)<<endl;
		//		cut<<pm.nx<<","<<pm.ny<<","<<abs(Mt)<<","<<abs(Jz)<<","<<abs(Jt)<<","<<abs(Mz)<<endl;
	}
	// 	cout<<"Min fields H*Z0, E/Z0 "<<H_min*Z0<<","<<E_min/Z0<<endl;
	if(output_currents){
		opc.close();
	}
	
	result[0] = 10.0*log10(0.25*K0*sq(abs(rval[0]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	result[1] = 10.0*log10(0.25*K0*sq(abs(rval[1]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	//return the phases
	result[2] = rtod(arg(rval[0]));
	result[3] = rtod(arg(rval[1]));
//	result[2] = rtod(arg(rval[0]/rval[1]));
	//	cut.close();
	return result;
}

//this is assuming PMC in case of TM and PEC in case of TM
vector<double> radar_pc(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, double soilw, double rec,double theta, double G,bool output_currents,double lambda)
{
	int i,j,k, N = (*ff).size(),ortn,sgn;
	node pi,pj,pm;	
	double psi,phi_i,phi_j,tphi,nphi,rp,dl,min_x=0,max_x=0,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda, g=K0*cos(theta),h=K0*sin(theta),ar,li,H_min=1,E_min=Z0; 
	vector <double> result(4,0);		
	cmp Mx,My,Mz,Jx,Jy,Jz,EtT,Ez,HtT,Hz,EzT,HzT,Hzs,Ets,Ezs,Hts,Mt,Jt,inc,er,mr; 
	vector <cmp> rval(2,cmp(0,0));
	
	//Assuming that the contour is always above the origin.
	//ortn stores match between contour orientation (c.ckwise) and edge orientation (i->j)
	
	//	ofstream cut; cut.open("cut.csv");
	ofstream opc;
	if(output_currents){
		opc.open("currents.v12.csv");
	}
	
	for(i=0; i<N; i++)
	{
		pm = (*edlist)[(*ff)[i]].midpoint(nlist);
		pi = (*nlist)[(*edlist)[(*ff)[i]].ni]; pj = (*nlist)[(*edlist)[(*ff)[i]].nj];
		//angles w.r.t center
		phi_i = atan2(pi.ny,pi.nx); phi_j = atan2(pj.ny,pj.nx);
		
		//match between contour and edge conv. 
		ortn = sign(phi_j - phi_i);
		if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
		{ortn *= -1;}
		
		//physical angle of tangent and outward normal (c.ckwise)
		tphi = atan2(double(ortn) * (pj.ny-pi.ny), double(ortn)*(pj.nx-pi.nx));
		nphi = tphi - PI/2.0; 
		tx = cos(tphi); ty = sin(tphi); nx = cos(nphi); ny = sin(nphi);
		
		//first get the total field 
		HtT = (*A)((*ff)[i]-1) * double(ortn); //H_min = Min(H_min,abs(HtT));
		EtT = (*B)((*ff)[i]-1) * double(ortn); //E_min = Min(E_min,abs(EtT));
		
		EzT = 0.0; HzT = 0.0;
		for(j=0; j<2; j++) //total field is avg over two elements
		{
			Ez = 0.0; Hz = 0.0;
			element ei = (*elist)[(*edlist)[(*ff)[i]].els[j]];
			ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
			
			for(k=0; k<3; k++)
			{
				sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
				li = (*edlist)[ei.edges[k]].edlength;
				Ez += double(sgn) * li * (*A)(ei.edges[k]-1);
				Hz += double(sgn) * li * (*B)(ei.edges[k]-1);
			}
			Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
			EzT += Ez; HzT += Hz;
		}
		EzT /= 2.0; HzT /= 2.0;
		
		//total currents, M = E x n, J = n x H
		//for TM, PMC => M=0, J=2xnxH. TM, PEC => J=0, M=2xExn
		Mt = 0.0 * EzT; 
		Jz = 2.0 * HtT;
		Jt = 0.0 * HzT; 
		Mz = -2.0 * EtT;
		
		//print if required
		if(output_currents){
			opc<<pm.nx<<","<<pm.ny<<","<<theta<<","<<abs(Mt)<<","<<arg(Mt)<<","<<abs(Jt)<<","<<arg(Jt)<<endl;
		}
		
		//rcs integrand needs these
		psi = atan2(pm.ny,pm.nx);
		rp = sqrt(sq(pm.ny)+sq(pm.nx));
		dl = sqrt(sq(pi.nx-pj.nx)+sq(pi.ny-pj.ny));
		
		//rcs
		rval[0] += (Mt * sin(tphi-rec) - Jz *Z0) * polar(dl,K0*rp*cos(psi-rec));
		rval[1] += (-1.0 * Jt * sin(tphi-rec) - Mz /Z0) * polar(dl,K0*rp*cos(psi-rec));
		
		//test
		//		cout<<phi_i<<" "<<phi_j<<" "<<ortn<<" "<<rval[0]<<" "<<rval[1]<<endl;
		//		cout<<abs(Jz)<<","<<abs(Mx)<<","<<abs(Mz)<<","<<abs(Jx)<<endl;
		//		cut<<pm.nx<<","<<pm.ny<<","<<abs(Mt)<<","<<abs(Jz)<<","<<abs(Jt)<<","<<abs(Mz)<<endl;
	}
	// 	cout<<"Min fields H*Z0, E/Z0 "<<H_min*Z0<<","<<E_min/Z0<<endl;
	if(output_currents){
		opc.close();
	}
	
	result[0] = 10.0*log10(0.25*K0*sq(abs(rval[0]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	result[1] = 10.0*log10(0.25*K0*sq(abs(rval[1]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	//return the phases
	result[2] = rtod(arg(rval[0]));
	result[3] = rtod(arg(rval[1]));
//	result[2] = rtod(arg(rval[0]/rval[1]));
	//	cut.close();
	return result;
}

//This function returns the component of the in-plane incident field
//along a given (vector) edge for both polarizations
vector<cmp> inc_field_along_edge(edge ed,double theta,double soilw,
	double lambda, double y0,vector<node> * nodelist){
	//In TM, Hi = j/(kZ) (dEi/dy, -dEi/dx), Ei = tapered thorsos
	//In TE, Ei = -jZ/k (dHi/dy, -dHi/dx), Hi = tapered thorsos
	//T has unit length, so T = (cos(phi), sin(phi))
	//we want T.Ui
	node pi = (*nodelist)[ed.ni], pj = (*nodelist)[ed.nj];
	double g = get_thorosG(theta,soilw), K0=2.0*PI/lambda,
	phi = atan2(pj.ny - pi.ny,pj.nx - pi.nx), u_phi = theta - PI/2.0;
	cmp douy = thoros_dy(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),g,theta,y0,lambda);
	cmp doux = thoros_dx(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),g,theta,y0,lambda);
	cmp u_d = I/(K0*Z0) * sqrt(douy*conj(douy)+doux*conj(doux));
	
	vector<cmp> result(2, 0.0);
	result[0] = I/(K0*Z0) * (cos(phi)*(douy) - sin(phi)*(doux));
	result[1] = -1.0 * Z0 * Z0 * result[0];
	return result;
}

//When the incident field is a plane wave
vector<cmp> inc_field_along_edge(edge ed,double theta,double lambda,vector<node> * nodelist){
	//In TM, Hi = j/(kZ) (dEi/dy, -dEi/dx), Ei =  thorsos
	//In TE, Ei = -jZ/k (dHi/dy, -dHi/dx), Hi =  thorsos
	//T has unit length, so T = (cos(phi), sin(phi))
	//we want T.Ui
	node pi = (*nodelist)[ed.ni], pj = (*nodelist)[ed.nj];
	double K0=2.0*PI/lambda,phi = atan2(pj.ny - pi.ny,pj.nx - pi.nx);
	
	cmp douy = thoros_dy(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),theta,lambda);
	cmp doux = thoros_dx(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),theta,lambda);
	
	vector<cmp> result(2, 0.0);
	result[0] = I/(K0*Z0) * (cos(phi)*(douy) - sin(phi)*(doux));
	result[1] = -1.0 * Z0 * Z0 * result[0];
	return result;
}

void store_currents(vector<vector<vector<currents> > > *all_currents,vector<int> *ff, 
vector<Vector<cmp> >*A, vector<Vector<cmp> >*B, double soilw, double lambda,double y0,int mcarlo, 
vector<double> * thetas, vector<node> *nlist, vector<edge> *edlist, vector<element> * elist){
	
	int th,i,j,k, N = (*ff).size(),ortn,sgn,ed;
	node pi,pj,pm;	
	double phi_i,phi_j,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda,ar,li,theta,G; 
	cmp inc,er,mr,HtT,EtT,HzT,EzT,Hz,Ez;
	
	for(th=0; th<(*thetas).size(); th++){
		theta = (*thetas)[th];
		G = get_thorosG(theta,soilw);
		for(i=0; i<N; i++){
			ed = (*ff)[i];
			pm = (*edlist)[ed].midpoint(nlist); pi = (*nlist)[(*edlist)[ed].ni]; 
			pj = (*nlist)[(*edlist)[ed].nj];
		
			//angles w.r.t center
			phi_i = atan2(pi.ny-y0,pi.nx); phi_j = atan2(pj.ny-y0,pj.nx);	
			//match between contour and edge conv. ortn stores match between
			//contour orientation (c.ckwise) and edge orientation (i->j)
			ortn = sign(phi_j - phi_i);
			if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
			{ortn *= -1;}
		
			//first get the total field = scat + inc
			vector<cmp> inc_comp = inc_field_along_edge((*edlist)[ed],theta,soilw,lambda,y0,nlist);
			HtT = ( (*A)[th](ed-1) + inc_comp[0] ) * double(ortn); 
			EtT = ( (*B)[th](ed-1) + inc_comp[1] ) * double(ortn); 
			
			EzT = 0.0; HzT = 0.0;
			for(j=0; j<2; j++){ //scat field is avg over two elements
				Ez = 0.0; Hz = 0.0;
				element ei = (*elist)[(*edlist)[(*ff)[i]].els[j]];
				ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
				for(k=0; k<3; k++){
					sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
					li = (*edlist)[ei.edges[k]].edlength;
					Ez += double(sgn) * li * (*A)[th](ei.edges[k]-1);
					Hz += double(sgn) * li * (*B)[th](ei.edges[k]-1);
				}
				Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
				EzT += Ez; HzT += Hz;
			}
			EzT /= 2.0; HzT /= 2.0;
			//add the inc fields now
			inc = thoros(pm.nx,pm.ny,G,theta,y0,lambda);
			EzT += inc; HzT += inc;
			
			//total currents, M = E x n, J = n x H
			(*all_currents)[th][mcarlo][i].Mt = EzT; 
			(*all_currents)[th][mcarlo][i].Jz = HtT;
			(*all_currents)[th][mcarlo][i].Jt = -1.0 * HzT; 
			(*all_currents)[th][mcarlo][i].Mz = -1.0 * EtT; 
		}
	}
}

void radar_sf_incoherent(double soilw, double y0, double lambda,
vector<int> * ff, vector<node> *nlist, vector<edge> *edlist,vector<double> *thetas,
vector<vector<vector< currents> > >* all_currents){
	double rec,rec_i,rec_f,rec_s,phi_i,psi,phi_j,tphi,K0=2.0*PI/lambda,G,rp,dl;
	int no_th = (*all_currents).size(), mcarlo = (*all_currents)[0].size(),
	segs = (*all_currents)[0][0].size(),i,j,k,ed,ortn;
	node pi,pj,pm; cmp Mt, Mz, Jt, Jz;
	
	//First generate the averages
	vector<vector<vector<currents> > >avg_currents(no_th,
	vector<vector<currents> >(1,
	vector<currents>(segs)));
	
	for(i=0; i<no_th; i++){
		for(j=0; j<mcarlo; j++){
			for(k=0; k<segs; k++){
				avg_currents[i][0][k].Jt += (*all_currents)[i][j][k].Jt/double(mcarlo);
				avg_currents[i][0][k].Jz += (*all_currents)[i][j][k].Jz/double(mcarlo);
				avg_currents[i][0][k].Mt += (*all_currents)[i][j][k].Mt/double(mcarlo);
				avg_currents[i][0][k].Mz += (*all_currents)[i][j][k].Mz/double(mcarlo);
			}
		}
	}
	
	ofstream ibircs("ibircs.v12.csv");
	//for each theta
	for(i=0; i<no_th; i++){
		rec_i = dtor(30); rec_f = dtor(150); rec_s = dtor(1); 
		G = get_thorosG((*thetas)[i],soilw);
		//for each instance 
		for(j=0; j<mcarlo; j++){
		ibircs<<rtod((*thetas)[i]);rec = rec_i;
		//calculate bistatic rcs
		while( rec <= rec_f){
			vector<cmp> rval(2,0);
			//from the far field line integral
			for(k=0; k<segs; k++){
				ed = (*ff)[k]; pm = (*edlist)[ed].midpoint(nlist); 
				pi = (*nlist)[(*edlist)[ed].ni]; pj = (*nlist)[(*edlist)[ed].nj];
				//angles w.r.t center
				phi_i = atan2(pi.ny-y0,pi.nx); phi_j = atan2(pj.ny-y0,pj.nx);
				
				//match between contour and edge conv. 
				ortn = sign(phi_j - phi_i);
				if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
				{ortn *= -1;}
				
				//physical angle of tangent and outward normal (c.ckwise)
				tphi = atan2(double(ortn) * (pj.ny-pi.ny), double(ortn)*(pj.nx-pi.nx));
				
				Mt = (*all_currents)[i][j][k].Mt - avg_currents[i][0][k].Mt; 
				Jz = (*all_currents)[i][j][k].Jz - avg_currents[i][0][k].Jz; 
				Jt = (*all_currents)[i][j][k].Jt - avg_currents[i][0][k].Jt; 
				Mz = (*all_currents)[i][j][k].Mz - avg_currents[i][0][k].Mz; 
				
				//rcs integrand needs these
				psi = atan2(pm.ny,pm.nx);
				rp = sqrt(sq(pm.ny)+sq(pm.nx));
				dl = abs((*edlist)[ed].edlength);
				
				//rcs
				rval[0] += (Mt * sin(tphi-rec) - Jz *Z0) * polar(dl,K0*rp*cos(psi-rec));
				rval[1] += (-1.0 * Jt * sin(tphi-rec) - Mz /Z0) * polar(dl,K0*rp*cos(psi-rec));
			}
			rval[0] = 10.0*log10(0.25*K0*sq(abs(rval[0]))/(G*abs(sin((*thetas)[i]))*sqrt(PI/2.0)));
			rval[1] = 10.0*log10(0.25*K0*sq(abs(rval[1]))/(G*abs(sin((*thetas)[i]))*sqrt(PI/2.0)));
			//since we are re-using rval, now as a double
			ibircs<<","<<real(rval[0])<<","<<real(rval[1]);
			rec += rec_s;
		}
		ibircs<<endl; 
		}
	} 
	ibircs.close();
	
}

//This returns the coherent bistatic fields for a scattered field approach
vector<cmp> farfield(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, double rec,double theta, double lambda){
	
	int i,j,k, N = (*ff).size(),ortn,sgn,ed;
	node pi,pj,pm;	
	double psi,phi_i,phi_j,tphi,nphi,rp,dl,min_x=0,max_x=0,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda,ar,li,H_min=1,E_min=Z0; 		
	cmp Mx,My,Mz,Jx,Jy,Jz,EtT,Ez,HtT,Hz,EzT,HzT,Hzs,Ets,Ezs,Hts,Mt,Jt,inc,er,mr; 
	vector<cmp> rval(2,cmp(0,0));
	//ortn stores match between contour orientation (c.ckwise) and edge orientation (i->j)
		
	for(i=0; i<N; i++){
		ed = (*ff)[i];
		pm = (*edlist)[ed].midpoint(nlist); pi = (*nlist)[(*edlist)[ed].ni]; 
		pj = (*nlist)[(*edlist)[ed].nj];
		//angles w.r.t center
		phi_i = atan2(pi.ny,pi.nx); phi_j = atan2(pj.ny,pj.nx);
		
		//match between contour and edge conv. 
		ortn = sign(phi_j - phi_i);
		if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
		{ortn *= -1;}
		
		//physical angle of tangent and outward normal (c.ckwise)
		tphi = atan2(double(ortn) * (pj.ny-pi.ny), double(ortn)*(pj.nx-pi.nx));
		nphi = tphi - PI/2.0; 
		tx = cos(tphi); ty = sin(tphi); nx = cos(nphi); ny = sin(nphi);
		
		//first get the total field = scat + inc
		//t is tangential, z is z-component, T is total
		vector<cmp> inc_comp = inc_field_along_edge((*edlist)[ed],theta,lambda,nlist);
		HtT = ( (*A)(ed-1) + inc_comp[0] ) * double(ortn); 
		EtT = ( (*B)(ed-1) + inc_comp[1] ) * double(ortn); 
		
		EzT = 0.0; HzT = 0.0;
		for(j=0; j<2; j++) //scat field is avg over two elements
		{
			Ez = 0.0; Hz = 0.0;
			element ei = (*elist)[(*edlist)[(*ff)[i]].els[j]];
			ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
			for(k=0; k<3; k++)
			{
				sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
				li = (*edlist)[ei.edges[k]].edlength;
				Ez += double(sgn) * li * (*A)(ei.edges[k]-1);
				Hz += double(sgn) * li * (*B)(ei.edges[k]-1);
			}
			Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
			EzT += Ez; HzT += Hz;
		}
		EzT /= 2.0; HzT /= 2.0;
		//add the inc fields now
		inc = thoros(pm.nx,pm.ny,theta,lambda);
		EzT += inc; HzT += inc;
		
		//total currents, M = E x n, J = n x H
		Mt = EzT; 
		Jz = HtT;
		Jt = -1.0 * HzT; 
		Mz = -1.0 * EtT;
		
		//rcs integrand needs these
		psi = atan2(pm.ny,pm.nx);
		rp = sqrt(sq(pm.ny)+sq(pm.nx));
		dl = sqrt(sq(pi.nx-pj.nx)+sq(pi.ny-pj.ny));
		
		//rcs
		rval[0] += (Mt * sin(tphi-rec) - Jz *Z0) * polar(dl,K0*rp*cos(psi-rec));
		rval[1] += (-1.0 * Jt * sin(tphi-rec) - Mz /Z0) * polar(dl,K0*rp*cos(psi-rec));
		
	}
	
	vector<cmp> result(2,0);
	result[0] = sqrt(0.25*K0)*rval[0];
	result[1] = sqrt(0.25*K0)*rval[1];
	
	return result;
}

//This returns the coherent bistatic fields for a total field approach
vector<cmp> farfield_tf(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, double rec,double theta, double lambda){
	
	int i,j,k, N = (*ff).size(),ortn,sgn,ed;
	node pi,pj,pm;	
	double psi,phi_i,phi_j,tphi,nphi,rp,dl,min_x=0,max_x=0,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda,ar,li,H_min=1,E_min=Z0; 		
	cmp Mx,My,Mz,Jx,Jy,Jz,EtT,Ez,HtT,Hz,EzT,HzT,Hzs,Ets,Ezs,Hts,Mt,Jt,inc,er,mr; 
	vector<cmp> rval(2,cmp(0,0));
	//ortn stores match between contour orientation (c.ckwise) and edge orientation (i->j)
		
	for(i=0; i<N; i++){
		ed = (*ff)[i];
		pm = (*edlist)[ed].midpoint(nlist); pi = (*nlist)[(*edlist)[ed].ni]; 
		pj = (*nlist)[(*edlist)[ed].nj];
		//angles w.r.t center
		phi_i = atan2(pi.ny,pi.nx); phi_j = atan2(pj.ny,pj.nx);
		
		//match between contour and edge conv. 
		ortn = sign(phi_j - phi_i);
		if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
		{ortn *= -1;}
		
		//physical angle of tangent and outward normal (c.ckwise)
		tphi = atan2(double(ortn) * (pj.ny-pi.ny), double(ortn)*(pj.nx-pi.nx));
		nphi = tphi - PI/2.0; 
		tx = cos(tphi); ty = sin(tphi); nx = cos(nphi); ny = sin(nphi);
		
		//first get the total field = scat + inc
		//t is tangential, z is z-component, T is total
	//	vector<cmp> inc_comp = inc_field_along_edge((*edlist)[ed],theta,lambda,nlist);
		HtT = ( (*A)(ed-1)  ) * double(ortn); 
		EtT = ( (*B)(ed-1)  ) * double(ortn); 
		
		EzT = 0.0; HzT = 0.0;
		for(j=0; j<2; j++) //scat field is avg over two elements
		{
			Ez = 0.0; Hz = 0.0;
			element ei = (*elist)[(*edlist)[(*ff)[i]].els[j]];
			ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
			for(k=0; k<3; k++)
			{
				sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
				li = (*edlist)[ei.edges[k]].edlength;
				Ez += double(sgn) * li * (*A)(ei.edges[k]-1);
				Hz += double(sgn) * li * (*B)(ei.edges[k]-1);
			}
			Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
			EzT += Ez; HzT += Hz;
		}
		EzT /= 2.0; HzT /= 2.0;
		//add the inc fields now
		//inc = thoros(pm.nx,pm.ny,theta,lambda);
		//EzT += inc; HzT += inc;
		
		//total currents, M = E x n, J = n x H
		Mt = EzT; 
		Jz = HtT;
		Jt = -1.0 * HzT; 
		Mz = -1.0 * EtT;
		
		//rcs integrand needs these
		psi = atan2(pm.ny,pm.nx);
		rp = sqrt(sq(pm.ny)+sq(pm.nx));
		dl = sqrt(sq(pi.nx-pj.nx)+sq(pi.ny-pj.ny));
		
		//rcs
		rval[0] += (Mt * sin(tphi-rec) - Jz *Z0) * polar(dl,K0*rp*cos(psi-rec));
		rval[1] += (-1.0 * Jt * sin(tphi-rec) - Mz /Z0) * polar(dl,K0*rp*cos(psi-rec));
		
	}
	
	vector<cmp> result(2,0);
	result[0] = sqrt(0.25*K0)*rval[0];
	result[1] = sqrt(0.25*K0)*rval[1];
	
	return result;
}

vector<cmp> field_sf(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, double rec,double theta, double lambda, bool flag){
	
	int i,j,k, N = (*ff).size(),ortn,sgn,ed;
	node pi,pj,pm;	
	double psi,phi_i,phi_j,tphi,nphi,rp,dl,min_x=0,max_x=0,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda,ar,li,H_min=1,E_min=Z0; 		
	cmp Mx,My,Mz,Jx,Jy,Jz,EtT,Ez,HtT,Hz,EzT,HzT,Hzs,Ets,Ezs,Hts,Mt,Jt,inc,er,mr; 
	vector<cmp> rval(4,cmp(0,0));
	//ortn stores match between contour orientation (c.ckwise) and edge orientation (i->j)
		
	for(i=0; i<N; i++){
		ed = (*ff)[i];
		pm = (*edlist)[ed].midpoint(nlist); 
		//angles w.r.t center
		psi = atan2(pm.ny,pm.nx);
		//psi=int(psi+360)% 360;
		if (psi<0)
			psi=psi+dtor(360);
		if (abs(psi-rec)<dtor(1)){
			//cout<<"  index : "<< i<<endl;
			break;
			}
	}
	//cout<<"ind"<< i<<endl;
		pi = (*nlist)[(*edlist)[ed].ni]; 
		pj = (*nlist)[(*edlist)[ed].nj];
		phi_i = atan2(pi.ny,pi.nx); phi_j = atan2(pj.ny,pj.nx);
		//match between contour and edge conv. 
		ortn = sign(phi_j - phi_i);
		if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
		{ortn *= -1;}
		//physical angle of tangent and outward normal (c.ckwise)
		//tphi = atan2(double(ortn) * (pj.ny-pi.ny), double(ortn)*(pj.nx-pi.nx));
		//nphi = tphi - PI/2.0; 
		//tx = cos(tphi); ty = sin(tphi); nx = cos(nphi); ny = sin(nphi);
		
		//first get the total field = scat + inc
		//t is tangential, z is z-component, T is total
		vector<cmp> inc_comp = inc_field_along_edge((*edlist)[ed],theta,lambda,nlist);
		HtT = ( (*A)(ed-1) + double(flag)*inc_comp[0] ) * double(ortn); 
		EtT = ( (*B)(ed-1) + double(flag)*inc_comp[1] ) * double(ortn); 
		EzT = 0.0; HzT = 0.0;
		for(j=0; j<2; j++) //scat field is avg over two elements
		{
			Ez = 0.0; Hz = 0.0;
			element ei = (*elist)[(*edlist)[(*ff)[i]].els[j]];
			ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
			for(k=0; k<3; k++)
			{
				sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
				li = (*edlist)[ei.edges[k]].edlength;
				Ez += double(sgn) * li * (*A)(ei.edges[k]-1);
				Hz += double(sgn) * li * (*B)(ei.edges[k]-1);
			}
			Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
			EzT += Ez; HzT += Hz;
		}
		EzT /= 2.0; HzT /= 2.0;
		//add the inc fields now
		if(flag==1){
			inc = thoros(pm.nx,pm.ny,theta,lambda);
			EzT += inc; HzT += inc;
		}
		
		//total currents, M = E x n, J = n x H
		Mt = EzT; 
		Jz = HtT;
		Jt = -1.0 * HzT; 
		Mz = -1.0 * EtT;
		
		
		
		//rcs
		rval[0] = Mt ;
		rval[1] = Jz;
		rval[2] =Jt;
		rval[3] =Mz;
	return rval;
}

//Store the scattered or total fields on the TOP contour in the case of the total-field (tf) formulation
//flag=0 => total field, flag=1 => scattered field (usually we are interested in this)
//A corresponds to TM, B corresponds to TE pol, and ed corresponds to the edge number where the field is desired
vector<cmp> field_tf_at_ed(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, int ed,double theta, double lambda, bool flag){
	
	int i,j,k, N = (*ff).size(),ortn,sgn;
	node pi,pj,pm;	
	double psi,phi_i,phi_j,tphi,nphi,rp,dl,min_x=0,max_x=0,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda,ar,li,H_min=1,E_min=Z0; 		
	cmp EtS,Ez,HtS,Hz,EzS,HzS,Hzs,Ets,Ezs,Hts,inc,er,mr; 
	vector<cmp> rval(4,cmp(0,0));
		
	pi = (*nlist)[(*edlist)[ed].ni]; 
	pj = (*nlist)[(*edlist)[ed].nj];
	phi_i = atan2(pi.ny,pi.nx); phi_j = atan2(pj.ny,pj.nx);
	
	//'ortn' stores match between contour orientation (c.ckwise)
	//and edge orientation. Convention is that node i points to node j (i->j)
	ortn = sign(phi_j - phi_i);
	if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
	{ortn *= -1;}
	
	//t is tangential, z is z-component, S is scattered (total-inc)
	//Get the tangential component first (TM - Ht, TE- Et)
	vector<cmp> inc_comp = inc_field_along_edge((*edlist)[ed],theta,lambda,nlist);
	HtS = ( (*A)(ed-1) - double(flag)*inc_comp[0] ) * double(ortn) ;
	EtS = ( (*B)(ed-1) - double(flag)*inc_comp[1] ) * double(ortn) ;
		
	//take the avg the z-component over two elements that share the edge (TM - Ez, TE - Hz)
	EzS = 0.0; HzS = 0.0;
	for(j=0; j<2; j++){
		Ez = 0.0; Hz = 0.0;
		element ei = (*elist)[(*edlist)[ed].els[j]];
		ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
		for(k=0; k<3; k++){
			//For an element with nodes i j k, the edges are defined 
			//as ij jk ki. Edge convention is smaller node # to larger node #.
			//sgn stores this convention
			sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
			li = (*edlist)[ei.edges[k]].edlength;
			Ez += double(sgn) * li * (*A)(ei.edges[k]-1);
			Hz += double(sgn) * li * (*B)(ei.edges[k]-1);
		}
		Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
		EzS += Ez; HzS += Hz;
	}
	EzS /= 2.0; HzS /= 2.0;
	//remove the inc field contribution if required (depending on flag)
	if(flag==1){
		inc = thoros(pm.nx,pm.ny,theta,lambda);
	//	inc = polar(1.0,-2.0*PI/lambda*(cos(theta)*pm.nx,sin(theta)*pm.ny));
		EzS -= inc; HzS -= inc;
	}

	//rcs
	rval[0] = EzS; //TM pol, Ez
	rval[1] = HtS; //TM pol, Ht
	rval[2] = HzS; //TE Pol, Hz
	rval[3] = EtS; //TE Pol, Et
	
	return rval;
}


//In the total field formulation, this function is asking for the field components at a receiver angle
//rec. This function finds out the closest edge number, ed, and calls field_tf_at_ed to get the result
vector<cmp> field_tf(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, double rec,double theta, double lambda, bool flag){
	
	int i,j,N = (*ff).size(),ed;
	node pm;	
	double psi;
	vector<cmp> rval(4,cmp(0,0));

/*	for(i=0; i<N; i++){
		ed = (*ff)[i];
		pm = (*edlist)[ed].midpoint(nlist); 
		//angles w.r.t center
		psi = atan2(pm.ny,pm.nx);
		//psi=int(psi+360)% 360;
		if (psi<0)
			psi=psi+dtor(360);
		if (abs(psi-rec)<dtor(1)){
			//cout<<"  index : "<< i<<endl;
			break;
		}
	} */
	//convert angles to [0,360]
	if(rec<0){rec += dtor(360.0);}
	//find the nearest edge
	double tmp = 2.0*PI; j=N;
	for(i=0; i<N; i++){
		ed = (*ff)[i];
		pm = (*edlist)[ed].midpoint(nlist); 
		psi = atan2(pm.ny,pm.nx);
		if( psi<0 ){ psi += dtor(360);}
		if( abs(psi-rec) < tmp ){ //update the nearest find
			tmp = abs(psi-rec);
			j = i;
		}
	}
	if( j==N ){cerr<<"** Something wrong in field_tf: no matching edge found! **"<<endl;}
	else{ed = (*ff)[j];}
	
	rval = field_tf_at_ed(nlist,edlist,elist,ff,A,B,ed,theta,lambda,flag);
	return rval;
}


//This returns the coherent bistatic RCS for a scattered field approach
vector<double> radar_sf(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, double soilw, double rec,double theta, double G,bool output_currents,
	double lambda,double y0){
	
	int i,j,k, N = (*ff).size(),ortn,sgn,ed;
	node pi,pj,pm;	
	double psi,phi_i,phi_j,tphi,nphi,rp,dl,min_x=0,max_x=0,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda,ar,li,H_min=1,E_min=Z0; 
	vector <double> result(4,0);		
	cmp Mx,My,Mz,Jx,Jy,Jz,EtT,Ez,HtT,Hz,EzT,HzT,Hzs,Ets,Ezs,Hts,Mt,Jt,inc,er,mr; 
	vector <cmp> rval(2,cmp(0,0));
	
	//Assuming that the contour is always above the origin.
	//ortn stores match between contour orientation (c.ckwise) and edge orientation (i->j)
	
	//	ofstream cut; cut.open("cut.csv");
	ofstream opc;
	if(output_currents){
		opc.open("currents.v12.csv");// cout<<"t,G,l,y "<<theta<<","<<G<<","<<lambda<<","<<y0<<endl;
	}
	
	for(i=0; i<N; i++){
		ed = (*ff)[i];
		pm = (*edlist)[ed].midpoint(nlist); pi = (*nlist)[(*edlist)[ed].ni]; 
		pj = (*nlist)[(*edlist)[ed].nj];
		//angles w.r.t center
		phi_i = atan2(pi.ny-y0,pi.nx); phi_j = atan2(pj.ny-y0,pj.nx);
		
		//match between contour and edge conv. 
		ortn = sign(phi_j - phi_i);
		if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
		{ortn *= -1;}
		
		//physical angle of tangent and outward normal (c.ckwise)
		tphi = atan2(double(ortn) * (pj.ny-pi.ny), double(ortn)*(pj.nx-pi.nx));
		nphi = tphi - PI/2.0; 
		tx = cos(tphi); ty = sin(tphi); nx = cos(nphi); ny = sin(nphi);
		
		//first get the total field = scat + inc
		vector<cmp> inc_comp = inc_field_along_edge((*edlist)[ed],theta,soilw,lambda,y0,nlist);
		HtT = ( (*A)(ed-1) + inc_comp[0] ) * double(ortn); 
		EtT = ( (*B)(ed-1) + inc_comp[1] ) * double(ortn); 
		
		EzT = 0.0; HzT = 0.0;
		for(j=0; j<2; j++) //scat field is avg over two elements
		{
			Ez = 0.0; Hz = 0.0;
			element ei = (*elist)[(*edlist)[(*ff)[i]].els[j]];
			ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
			for(k=0; k<3; k++)
			{
				sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
				li = (*edlist)[ei.edges[k]].edlength;
				Ez += double(sgn) * li * (*A)(ei.edges[k]-1);
				Hz += double(sgn) * li * (*B)(ei.edges[k]-1);
			}
			Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
			EzT += Ez; HzT += Hz;
		}
		EzT /= 2.0; HzT /= 2.0;
		//add the inc fields now
		inc = thoros(pm.nx,pm.ny,G,theta,y0,lambda);
		EzT += inc; HzT += inc;
		
		//total currents, M = E x n, J = n x H
		Mt = EzT; 
		Jz = HtT;
		Jt = -1.0 * HzT; 
		Mz = -1.0 * EtT;
		
		//print if required
		if(output_currents){
			opc<<pm.nx<<","<<pm.ny<<","<<theta<<","<<abs(Mt)<<","<<arg(Mt)<<","<<abs(Jz)<<
			","<<arg(Jz)<<","<<abs(Jt)<<","<<arg(Jt)<<","<<abs(Mz)<<","<<arg(Mz)<<endl;
		//	opc<<pm.nx<<","<<pm.ny<<","<<theta<<","<<abs((*A)(ed-1))<<","<<arg((*A)(ed-1))
		//	<<","<<abs((*B)(ed-1))<<","<<arg((*B)(ed-1))<<","<<abs(inc_comp[0])<<","
		//	<<abs(inc_comp[1])<<","<<abs(inc)<<endl;
		}
		
		//rcs integrand needs these
		psi = atan2(pm.ny,pm.nx);
		rp = sqrt(sq(pm.ny)+sq(pm.nx));
		dl = sqrt(sq(pi.nx-pj.nx)+sq(pi.ny-pj.ny));
		
		//rcs
		rval[0] += (Mt * sin(tphi-rec) - Jz *Z0) * polar(dl,K0*rp*cos(psi-rec));
		rval[1] += (-1.0 * Jt * sin(tphi-rec) - Mz /Z0) * polar(dl,K0*rp*cos(psi-rec));
		
		//test
		//		cout<<phi_i<<" "<<phi_j<<" "<<ortn<<" "<<rval[0]<<" "<<rval[1]<<endl;
		//		cout<<abs(Jz)<<","<<abs(Mx)<<","<<abs(Mz)<<","<<abs(Jx)<<endl;
		//		cut<<pm.nx<<","<<pm.ny<<","<<abs(Mt)<<","<<abs(Jz)<<","<<abs(Jt)<<","<<abs(Mz)<<endl;
	}
	// 	cout<<"Min fields H*Z0, E/Z0 "<<H_min*Z0<<","<<E_min/Z0<<endl;
	if(output_currents){
		opc.close();
	}
	
	result[0] = 10.0*log10(0.25*K0*sq(abs(rval[0]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	result[1] = 10.0*log10(0.25*K0*sq(abs(rval[1]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	//return the phases
	result[2] = rtod(arg(rval[0]));
	result[3] = rtod(arg(rval[1]));
//	result[2] = rtod(arg(rval[0]/rval[1]));
	//	cut.close();
	return result;
}

//This returns the coherent bistatic RCS for a scattered field approach
//this is for contour on the surface itself
vector<double> radar_sf_surf(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, double soilw, double rec,double theta, double G,bool output_currents,
	double lambda,double y0){
	
	int i,j,k, N = (*ff).size(),ortn,sgn,ed;
	node pi,pj,pm;	
	double psi,phi_i,phi_j,tphi,nphi,rp,dl,min_x=0,max_x=0,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda,ar,li,H_min=1,E_min=Z0; 
	vector <double> result(4,0);		
	cmp Mx,My,Mz,Jx,Jy,Jz,EtT,Ez,HtT,Hz,EzT,HzT,Hzs,Ets,Ezs,Hts,Mt,Jt,inc,er,mr; 
	vector <cmp> rval(2,cmp(0,0));
	
	//Assuming that the contour is always above the origin.
	//ortn stores match between contour orientation (c.ckwise) and edge orientation (i->j)
	
	//	ofstream cut; cut.open("cut.csv");
	ofstream opc;
	if(output_currents){
		opc.open("currents.v12.csv");// cout<<"t,G,l,y "<<theta<<","<<G<<","<<lambda<<","<<y0<<endl;
	}
	
	for(i=0; i<N; i++){
		ed = (*ff)[i];
		pm = (*edlist)[ed].midpoint(nlist); pi = (*nlist)[(*edlist)[ed].ni]; 
		pj = (*nlist)[(*edlist)[ed].nj];
		//angles w.r.t center
		phi_i = atan2(pi.ny-10*lambda,pi.nx); phi_j = atan2(pj.ny-10*lambda,pj.nx);
		
		//match between contour and edge conv. 
		ortn = sign(phi_j - phi_i);
		if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
		{ortn *= -1;}
		
		//physical angle of tangent and outward normal (c.ckwise)
		tphi = atan2(double(ortn) * (pj.ny-pi.ny), double(ortn)*(pj.nx-pi.nx));
		nphi = tphi - PI/2.0; 
		tx = cos(tphi); ty = sin(tphi); nx = cos(nphi); ny = sin(nphi);
		
		//first get the total field = scat + inc
		vector<cmp> inc_comp = inc_field_along_edge((*edlist)[ed],theta,soilw,lambda,y0,nlist);
		HtT = ( (*A)(ed-1) + inc_comp[0] ) * double(ortn); 
		EtT = ( (*B)(ed-1) + inc_comp[1] ) * double(ortn); 
		
		EzT = 0.0; HzT = 0.0;
		for(j=0; j<2; j++) //scat field is only vac field
		{
			Ez = 0.0; Hz = 0.0;
			element ei = (*elist)[(*edlist)[(*ff)[i]].els[j]];
			if( ei.eltype != 'v'){
				continue;
			}
			ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
			for(k=0; k<3; k++)
			{
				sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
				li = (*edlist)[ei.edges[k]].edlength;
				Ez += double(sgn) * li * (*A)(ei.edges[k]-1);
				Hz += double(sgn) * li * (*B)(ei.edges[k]-1);
			}
			Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
			EzT += Ez; HzT += Hz;
		}
		//EzT /= 2.0; HzT /= 2.0;
		//add the inc fields now
		inc = thoros(pm.nx,pm.ny,G,theta,y0,lambda);
		EzT += inc; HzT += inc;
		
		//total currents, M = E x n, J = n x H
		Mt = EzT; 
		Jz = HtT;
		Jt = -1.0 * HzT; 
		Mz = -1.0 * EtT;
		
		//print if required
		if(output_currents){
			opc<<pm.nx<<","<<pm.ny<<","<<theta<<","<<abs(Mt)<<","<<arg(Mt)<<","<<abs(Jz)<<
			","<<arg(Jz)<<","<<abs(Jt)<<","<<arg(Jt)<<","<<abs(Mz)<<","<<arg(Mz)<<endl;
		//	opc<<pm.nx<<","<<pm.ny<<","<<theta<<","<<abs((*A)(ed-1))<<","<<arg((*A)(ed-1))
		//	<<","<<abs((*B)(ed-1))<<","<<arg((*B)(ed-1))<<","<<abs(inc_comp[0])<<","
		//	<<abs(inc_comp[1])<<","<<abs(inc)<<endl;
		}
		
		//rcs integrand needs these
		psi = atan2(pm.ny,pm.nx);
		rp = sqrt(sq(pm.ny)+sq(pm.nx));
		dl = sqrt(sq(pi.nx-pj.nx)+sq(pi.ny-pj.ny));
		
		//rcs
		rval[0] += (Mt * sin(tphi-rec) - Jz *Z0) * polar(dl,K0*rp*cos(psi-rec));
		rval[1] += (-1.0 * Jt * sin(tphi-rec) - Mz /Z0) * polar(dl,K0*rp*cos(psi-rec));
		
		//test
		//		cout<<phi_i<<" "<<phi_j<<" "<<ortn<<" "<<rval[0]<<" "<<rval[1]<<endl;
		//		cout<<abs(Jz)<<","<<abs(Mx)<<","<<abs(Mz)<<","<<abs(Jx)<<endl;
		//		cut<<pm.nx<<","<<pm.ny<<","<<abs(Mt)<<","<<abs(Jz)<<","<<abs(Jt)<<","<<abs(Mz)<<endl;
	}
	// 	cout<<"Min fields H*Z0, E/Z0 "<<H_min*Z0<<","<<E_min/Z0<<endl;
	if(output_currents){
		opc.close();
	}
	
	result[0] = 10.0*log10(0.25*K0*sq(abs(rval[0]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	result[1] = 10.0*log10(0.25*K0*sq(abs(rval[1]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	//return the phases
	result[2] = rtod(arg(rval[0]));
	result[3] = rtod(arg(rval[1]));
//	result[2] = rtod(arg(rval[0]/rval[1]));
	//	cut.close();
	return result;
}


//This returns the bistatic RCS for a scattered field approach radar_pc version for PEC/PMC
vector<double> radar_pc(vector<node> *nlist, vector<edge> *edlist, vector<element> * elist, vector<int> * ff,
	Vector<cmp> * A, Vector<cmp> * B, double soilw, double rec,double theta, double G,bool output_currents,
	double lambda,double y0){
	
	int i,j,k, N = (*ff).size(),ortn,sgn,ed;
	node pi,pj,pm;	
	double psi,phi_i,phi_j,tphi,nphi,rp,dl,min_x=0,max_x=0,px,py,tmp,nx,ny,tx,ty,
	K0 = 2.0*PI/lambda,ar,li,H_min=1,E_min=Z0; 
	vector <double> result(4,0);		
	cmp Mx,My,Mz,Jx,Jy,Jz,EtT,Ez,HtT,Hz,EzT,HzT,Hzs,Ets,Ezs,Hts,Mt,Jt,inc,er,mr; 
	vector <cmp> rval(2,cmp(0,0));
	
	//Assuming that the contour is always above the origin.
	//ortn stores match between contour orientation (c.ckwise) and edge orientation (i->j)
	
	//	ofstream cut; cut.open("cut.csv");
	ofstream opc;
	if(output_currents){
		opc.open("currents.v12.csv");
	}
	
	for(i=0; i<N; i++){
		ed = (*ff)[i];
		pm = (*edlist)[ed].midpoint(nlist); pi = (*nlist)[(*edlist)[ed].ni]; pj = (*nlist)[(*edlist)[ed].nj];
		//angles w.r.t center
		phi_i = atan2(pi.ny,pi.nx); phi_j = atan2(pj.ny,pj.nx);
		
		//match between contour and edge conv. 
		ortn = sign(phi_j - phi_i);
		if( (phi_i > PI/2.0 && phi_j < PI/-2.0) || (phi_i < PI/-2.0 && phi_j > PI/2.0) )
		{ortn *= -1;}
		
		//physical angle of tangent and outward normal (c.ckwise)
		tphi = atan2(double(ortn) * (pj.ny-pi.ny), double(ortn)*(pj.nx-pi.nx));
		nphi = tphi - PI/2.0; 
		tx = cos(tphi); ty = sin(tphi); nx = cos(nphi); ny = sin(nphi);
		
		//first get the total field = scat + inc
		vector<cmp> inc_comp = inc_field_along_edge((*edlist)[ed],theta,soilw,lambda,y0,nlist);
		HtT = ( (*A)(ed-1) + inc_comp[0] ) * double(ortn); 
		EtT = ( (*B)(ed-1) + inc_comp[1] ) * double(ortn); 
		
		EzT = 0.0; HzT = 0.0;
		for(j=0; j<2; j++) //scat field is avg over two elements
		{
			Ez = 0.0; Hz = 0.0;
			element ei = (*elist)[(*edlist)[(*ff)[i]].els[j]];
			ar = ei.elarea; mr = ei.elmu; er = ei.eleps;
			for(k=0; k<3; k++)
			{
				sgn = sign( ei.nodes[(k+1)%3] - ei.nodes[k%3]);
				li = (*edlist)[ei.edges[k]].edlength;
				Ez += double(sgn) * li * (*A)(ei.edges[k]-1);
				Hz += double(sgn) * li * (*B)(ei.edges[k]-1);
			}
			Ez *= -1.0*I*Z0/(ar*K0*er); Hz *= I/(ar*Z0*K0*mr);
			EzT += Ez; HzT += Hz;
		}
		EzT /= 2.0; HzT /= 2.0;
		//add the inc fields now
		inc = thoros(pm.nx,pm.ny,G,theta,y0,lambda);
		EzT += inc; HzT += inc;
		
		//total currents, M = E x n, J = n x H
		//for TM, PMC => M=0, J=2xnxH. TM, PEC => J=0, M=2xExn
		Mt = 0.0 * EzT; 
		Jz = 2.0 * HtT;
		Jt = 0.0 * HzT; 
		Mz = -2.0 * EtT;
		
		//print if required
		if(output_currents){
			opc<<pm.nx<<","<<pm.ny<<","<<theta<<","<<abs(Mt)<<","<<arg(Mt)<<","<<abs(Jt)<<","<<arg(Jt)<<endl;
		}
		
		//rcs integrand needs these
		psi = atan2(pm.ny,pm.nx);
		rp = sqrt(sq(pm.ny)+sq(pm.nx));
		dl = sqrt(sq(pi.nx-pj.nx)+sq(pi.ny-pj.ny));
		
		//rcs
		rval[0] += (Mt * sin(tphi-rec) - Jz *Z0) * polar(dl,K0*rp*cos(psi-rec));
		rval[1] += (-1.0 * Jt * sin(tphi-rec) - Mz /Z0) * polar(dl,K0*rp*cos(psi-rec));
		
		//test
		//		cout<<phi_i<<" "<<phi_j<<" "<<ortn<<" "<<rval[0]<<" "<<rval[1]<<endl;
		//		cout<<abs(Jz)<<","<<abs(Mx)<<","<<abs(Mz)<<","<<abs(Jx)<<endl;
		//		cut<<pm.nx<<","<<pm.ny<<","<<abs(Mt)<<","<<abs(Jz)<<","<<abs(Jt)<<","<<abs(Mz)<<endl;
	}
	// 	cout<<"Min fields H*Z0, E/Z0 "<<H_min*Z0<<","<<E_min/Z0<<endl;
	if(output_currents){
		opc.close();
	}
	
	result[0] = 10.0*log10(0.25*K0*sq(abs(rval[0]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	result[1] = 10.0*log10(0.25*K0*sq(abs(rval[1]))/(G*abs(sin(theta))*sqrt(PI/2.0)));
	//return the phases
	result[2] = rtod(arg(rval[0]));
	result[3] = rtod(arg(rval[1]));
//	result[2] = rtod(arg(rval[0]/rval[1]));
	//	cut.close();
	return result;
}

vector<double> print_spm(double r,double corl,cmp eps,double theta, double rec,double lambda,int scorl){
	//theta and rec are measured from the +x axis, where as 
	//ti and ts are inc and scat angles resp. measured from +y axis and +ve.
	double ti = abs(PI/2.0 - abs(theta));
	double ts = abs(PI/2.0 - abs(rec));
	double ps; //0 in fwd direction, PI in bwd direction
	
	if(abs(theta) < PI/2.0){
		if(abs(rec) < PI/2.0){ps = 0;}
		else{ps = PI;}
	}
	else{
		if(abs(rec) < PI/2.0){ps = PI;}
		else{ps = 0;}
	}	
	
	double K0 = 2.0 * PI/lambda, h = r/K0; //rms height
	double l = corl*h; //surface correlation length
	//W(k) = l/(2sqrt(pi)) exp(-(kl)^2/4)
	//sigma \propto W(2ksin(th)) in b.s. or W(k(sin(ts)cos(ps)-sin(ti))) in general.
	double W;
	if(scorl == 0){
		W = l * sq(h) / (2.0*sqrt(PI)) * exp(-0.25*sq(K0*(sin(ts)*cos(ps)-sin(ti))*l)); //Gaussian
	}
	else{
		W = sq(h) * l / (PI*(sq(K0*(sin(ts)*cos(ps)-sin(ti))*l)+1.0)) ; //Exponential	
	}
	//big mistake above.jakob's book has for 2d,where as for 1d case it is sqrt(2) in l, and a 1/2 in W.
	
	cmp a_tm,a_te;
	//formula in Jakob's book only for back-scatter and 2D
	//	a_tm = (eps-1.0)/sq(cos(th)+sqrt(eps-sq(sin(th))));
	//	a_te = (eps-1.0)*((eps-1.0)*sq(sin(th))+eps)/sq(eps*cos(th)+sqrt(eps-sq(sin(th))));
	
	//formula in Ishimaru's book
	a_tm = (eps-1.0)*cos(ps)/((cos(ti)+sqrt(eps-sq(sin(ti))))*(cos(ts)+sqrt(eps-sq(sin(ts)))));
	a_te = (eps-1.0)*(eps*sin(ts)*sin(ti)-cos(ps)*sqrt(eps-sq(sin(ti)))*sqrt(eps-sq(sin(ts))))/
	((eps*cos(ti)+sqrt(eps-sq(sin(ti))))*(eps*cos(ts)+sqrt(eps-sq(sin(ts)))));
	
	double s_tm = 8.0*PI*K0* cos(ti) * sq(cos(ts) * K0 * abs(a_tm)) * W;
	double s_te = 8.0*PI*K0* cos(ti) * sq(cos(ts) * K0 * abs(a_te)) * W;
	
	vector <double> rval(2,0);
	rval[0] = 10.0*log10(s_tm);
	rval[1] = 10.0*log10(s_te);
	//	cout<<"SPM rcs for h,l = ("<<h<<","<<l<<") is "<<rval[0]<<","<<rval[1]<<endl;
	return rval;
}
	
//Return true if i is in the integest list
bool in_list(int i,vector<int> *lst){
	for(int j=0;j<(*lst).size();j++){
		if(i == (*lst)[j] ){
			return true;
		}
	}
	return false;
}

//Returns true if node p is inside an ellipse of center c, axis a,b 
bool in_ellipse(node p, node c, double a, double b){
	if( sq((p.nx-c.nx)/a) + sq((p.ny-c.ny)/b) <= 1.0 ){
		return true;
	}
	else{
		return false;
	}
}

//Brute force way to define an ellipse
void make_boulder(node p, double r, double bbya, vector<int> * itf, vector<element> * elist, vector<node> * nlist){
	int i;
	//Here, r is the semi major axis, 'a'
	//eccentricity is e = sqrt(1-(b/a)^2) => b = a sqrt(1-e^2)
	//so e=0 is a circle
	//if( e<-0.00001 || e> 0.9999){
	//	cout<<"**Unphysical ellipse being constructed."<<endl;
	//	return;
	//}
	if( bbya < -0.001){
		cout<<"**Unphysical ellipse being constructed."<<endl;
		return;
	}
		
	double a = r, b = a * bbya;
	for(i=0; i<(*itf).size(); i++){
		node pm = (*elist)[(*itf)[i]].elcenter(nlist);
		if( sq((p.nx-pm.nx)/a)+sq((p.ny-pm.ny)/b) < 1.0 ){
			(*elist)[(*itf)[i]].eltype = 's';
		}
	}
}

//Smarter way to define an ellipse centered at p
void make_boulder(node p, double r, double bbya, vector<node> * nlist,vector<edge> * edlist,
				  vector<element> * elist, vector<vector<int> > * eal){
	//Here, r is the semi major axis, 'a'
	//eccentricity is e = sqrt(1-(b/a)^2) => b = a sqrt(1-e^2)
	//so e=0 is a circle
	//if( e<-0.00001 || e> 0.9999){
	//	cout<<"**Unphysical ellipse being constructed."<<endl;
	//	return;
	//}
	if( bbya < -0.001){
		cout<<"**Unphysical ellipse being constructed."<<endl;
		return;
	}	
	double a = r, b = a * bbya;
	int current,i,j,k,ngb;
	vector<int> openset, closedset, elset;
	//find the center, any node (here "1") as a seed guess 
	current = closest_node(p,1,eal,nlist,edlist);
	//initialize the openset with the center
	openset.push_back(current);
	
	while( openset.size() !=0 ){
		//Move the back of the openset to the closedset
		current = openset.back();
		openset.pop_back();
		closedset.push_back(current);
		//Examine the neighbours
		for(i=0; i< (*eal)[current].size(); i++){
			ngb = other_node(current, (*edlist)[(*eal)[current][i]]);
			node pi = (*nlist)[ngb];
			//if this ngb is outside, skip it
			if( !in_ellipse(pi,p,a,b) ){
				continue;
			}
			//add the associated elements
			for(j=0; j<(*edlist)[(*eal)[current][i]].els.size(); j++){
				k = (*edlist)[(*eal)[current][i]].els[j];
				if(!in_list(k,&elset)){
					elset.push_back(k);
				}
			}
			//if it is either closed or open set, skip it
			if( in_list(ngb,&openset) || in_list(ngb,&closedset) ){
				continue;
			}
			//add it since it is a legitimate node
			openset.push_back(ngb);
		}
	}
	
	//Now set the properties
	for(i=0;i<elset.size();i++){
		//refine the elements to have nodes only from the closed set.
		int n0 = (*elist)[elset[i]].nodes[0], n1 = (*elist)[elset[i]].nodes[1], n2 = (*elist)[elset[i]].nodes[2];
		if( in_list(n0,&closedset) && in_list(n1,&closedset) && in_list(n2,&closedset) ){
			(*elist)[elset[i]].eltype = 's';
		}
	}
}				   

bool ensemble_converge( vector<double> * thetas, vector< vector< double> > * tm, vector< vector< double> > * te, int mc ){
	//This function returns true as long as ensemble has not converged to within 1 dB.
	int i,j,thetae = (*thetas).size(), inst = (*tm)[0].size();
	
	if( inst == 0 ){
		return true;
	}
	
	//A peculiar case happens when the iteration has filled tm/te with inf values
	//probably a memory issue
	bool junk = false;
	for(i=0; i<thetae; i++){
		if( abs((*tm)[i][inst-1]) > 1000 || abs((*te)[i][inst-1]) > 1000 ){
			junk = true; break;
		}
	}
	if(junk){
		for(i=0; i<thetae; i++){
			(*tm)[i].pop_back(); (*te)[i].pop_back();
		}
		cout<<"*Inst="<<inst<<" has produced junk. C'est la vie."<<endl;
		return true;
	}
	
	//Build up stats
	vector<double> mn_tm(thetae,0.0), mn_te(thetae,0.0), sd_tm(thetae,0.0), sd_te(thetae,0.0);
	for(i=0; i<thetae; i++){
		//first find the mean 
		for(j=0; j<inst; j++){
			mn_tm[i] += (*tm)[i][j];
			mn_te[i] += (*te)[i][j];
		}
		mn_tm[i] /= double(inst);
		mn_te[i] /= double(inst);
		//now the std dev
		for(j=0; j<inst; j++){
			sd_tm[i] += sq((*tm)[i][j] - mn_tm[i]);
			sd_te[i] += sq((*te)[i][j] - mn_te[i]);
		}
		sd_tm[i] = sqrt(sd_tm[i]/double(inst));
		sd_te[i] = sqrt(sd_te[i]/double(inst));
		//the std dev of eqv CLT process is
		sd_tm[i] /= sqrt(double(inst));
		sd_te[i] /= sqrt(double(inst));
	}
	cout<<"*Inst="<<inst<<",theta=(";
	for(i=0; i<thetae-1; i++){cout<<rtod((*thetas)[i])<<",";}
	cout<<rtod((*thetas)[thetae-1])<<"), TM_mu=(";
	for(i=0; i<thetae-1; i++){cout<<mn_tm[i]<<",";}
	cout<<mn_tm[thetae-1]<<"), TE_mu=(";
	for(i=0; i<thetae-1; i++){cout<<mn_te[i]<<",";}
	cout<<mn_te[thetae-1]<<"), TM_sd=(";
	for(i=0; i<thetae-1; i++){cout<<sd_tm[i]<<",";}
	cout<<sd_tm[thetae-1]<<"), TE_sd=(";
	for(i=0; i<thetae-1; i++){cout<<sd_te[i]<<",";}
	cout<<sd_te[thetae-1]<<")"<<endl;
	
	if(inst == mc){return false;}
	
	bool result = true; 
	
/*	double tol = 1.5;
	for(i=0; i<thetae; i++){
		if(sd_tm[i] > tol|| sd_te[i] > tol){
			result = true; break;
		}
	} */
	
	return result;
}
//This function returns a quadratic fit for soil moisture
//given the sm measurements at three different depths
//sm[0] + sm[1] z + sm[2] z^2, z=0 is soil top, z +ve.
vector<double> quadratic_sm(vector<double> s){
	double z1 = 0.04, z2 = 0.13, z3 = 0.30;
	double s1 = s[0], s2 = s[1], s3 = s[2];
	vector<double> sm(3,0);
//	cout<<s1<<" "<<s2<<" "<<s3;
	sm[0] = (s3*(-sq(z1)* z2 + z1* sq(z2)))/(-sq(z1)* z2 + z1* sq(z2) + sq(z1)* z3 - sq(z2)* z3 - 
    z1* sq(z3) + z2* sq(z3)) + (
   s2* (sq(z1)* z3 - z1* sq(z3)))/(-sq(z1)* z2 + z1* sq(z2) + sq(z1)* z3 - sq(z2)* z3 - 
    z1* sq(z3) + z2* sq(z3)) + (
   s1* (-sq(z2)* z3 + z2* sq(z3)))/(-sq(z1)* z2 + z1* sq(z2) + sq(z1)* z3 - sq(z2)* z3 - 
    z1* sq(z3) + z2* sq(z3));
    
    sm[1] = (s3* (sq(z1) - sq(z2)))/(-sq(z1)* z2 + z1* sq(z2) + sq(z1)* z3 - sq(z2)* z3 - 
    z1* sq(z3) + z2* sq(z3)) + (
   s1* (sq(z2) - sq(z3)))/(-sq(z1)* z2 + z1* sq(z2) + sq(z1)* z3 - sq(z2)* z3 - 
    z1* sq(z3) + z2* sq(z3)) + (
   s2* (-sq(z1) + sq(z3)))/(-sq(z1)* z2 + z1* sq(z2) + sq(z1)* z3 - sq(z2)* z3 - 
    z1* sq(z3) + z2* sq(z3));
    
    sm[2] = (s3* (-z1 + z2))/(-sq(z1)* z2 + z1* sq(z2) + sq(z1)* z3 - sq(z2)* z3 - z1* sq(z3) + 
    z2* sq(z3)) + (
   s2* (z1 - z3))/(-sq(z1)* z2 + z1* sq(z2) + sq(z1)* z3 - sq(z2)* z3 - z1* sq(z3) + 
    z2* sq(z3)) + (
   s1* (-z2 + z3))/(-sq(z1)* z2 + z1* sq(z2) + sq(z1)* z3 - sq(z2)* z3 - z1* sq(z3) + 
    z2* sq(z3));
    
    return sm;
}

//This fn performs Ep = Ep - E0 given an int list of locations of non-zero E0
void make_deltaE(vector<int> * diel_ed, Vector<cmp> * E0_TM, Vector<cmp> * E0_TE, 
	Vector<cmp> * Ep_TM, Vector<cmp> * Ep_TE){
	
	int i,j;
	for(i=0; i<(*diel_ed).size(); i++){
		j = (*diel_ed)[i];
		//j-1 below because edge numbers begin from 1
		(*Ep_TM)(j-1) -= (*E0_TM)(j-1);
		(*Ep_TE)(j-1) -= (*E0_TE)(j-1);
	}
}
//This fn performs C = A + B
void AddVec(Vector<cmp> * A, Vector<cmp> * B,Vector<cmp> * C){
	int i;
	for(i=0; i<(*A).GetSize(); i++){
		(*C)(i) = (*A)(i) + (*B)(i);
	}
}
	
		
//These functions return the mesh to pristine condition!
void refresh_mesh(vector<edge> * edlist, vector<int> *intf_eds){

	int i;
	for(i=0; i<(*intf_eds).size(); i++){
		//change 'd' to 'i'
		(*edlist)[(*intf_eds)[i]].edbtype = 'i';
	}
}

void refresh_mesh(vector<int> * itf_el,vector<int> * mat_el, 
	vector<edge> * edlist, vector<element> * elist, vector<int> *intf_eds){

	int i;
	//Reset all element types to 's'
	for(i=0; i<(*mat_el).size(); i++){
		(*elist)[(*mat_el)[i]].eltype = 's';
	}
	for(i=0; i<(*itf_el).size(); i++){
		(*elist)[(*itf_el)[i]].eltype = 's';
	}	
	refresh_mesh(edlist,intf_eds);
}

void refresh_mesh(vector<record> * changes,vector<int> * itf_el,vector<int> * mat_el, 
	vector<node> * nlist, vector<edge> * edlist, vector<element> * elist, vector<int> *intf_eds){
	restore_record(changes,nlist); 
	refresh_mesh(itf_el,mat_el,edlist,elist,intf_eds);
}

cmp icemvtoeps(double lambda, double mv){
	cmp eps;
	//p is from 0.1 to 0.5 g cm^-3, density.
	//mv = 0 to 20 % (just a guess for the upper layer)
	if(mv <0.0 || mv >1.0){
		cout<<"** Unphysical ice mv being requested. **"<<endl;
		return 0;
	}
	mv *= 100.0; //the program is called with mv in cm^3/cm^3, formula has it in %.
	double B = 0.073, f0 = 9.07e9, f = c_0/lambda, p = 0.3,A = 1 + 1.83*p + 0.02*pow(mv,1.015);
	
	eps.real() = A + B*(pow(mv,1.31)/(1.0 + sq(f/f0)));
	eps.imag() =-1.0* B*(f/f0)*(pow(mv,1.31)/(1.0+sq(f/f0)));
	
	return eps;
}

cmp soilmvtoeps_peplinsky(double lambda, double mv, double sand, double clay, double rb){
	//valid from 23cm to 1m wavelengths.
/*************************************************************************************
! THIS SUBROUTINE CALCULATES THE REAL AND IMAGINARY PARTS OF THE DIELCTRIC CONSTANT
! OF SOIL GIVEN THE SOIL PARAMETERS. THE FORMULATION IS BASED ON
! N. R. PEPLINSKY, F. T. ULABY, AND M. C. DOBSON, ''DIELECTRIC PROPERTIES OF SOILS IN
! 0.3-3.0 GHZ RANGE,'' IEEE TRANS. GEOSCIENCE REMOTE SENSING, VOL. 33, NO. 3, MAY 1995.
! THE PARAMETERS ARE
!	FREQ: FREQUENCY OF MEASUREMENT
!	RB: BULK DENSITY (G/CM3)
!	T: WATER TEMPERATURE (C)
!	SA: SALINITY OF WATER (4-35) (GRAMS SALT PER KG WATER)
!	S, C: MASS FRACTION OF SAND AND CLAY
!	MV: WATER VOLUME FRACTION

! CODE BY JOEL JOHNSON, THE OHIO STATE UNIVERSITY.*/
	double freq,K0,bp,bpp,rs,es,e0,t,sa,einf,relax,sigma1,efw,sigma,efwp,emp,empp,alpha;
	t = 10.0; sa = 4.0; rs = 2.66; alpha = 0.65;
	
	freq = c_0/lambda; K0 = 2.0*PI/lambda;
	bp = 1.2748-0.519*sand - 0.152*clay;
	bpp = 1.33797 - 0.603*sand - 0.166*clay;
	es = sq(1.01+0.44*rs) - 0.062;

	e0=(87.134 - (1.949e-1)*t - (1.276e-2)*sq(t) + (2.491e-4)*t*sq(t))* 
	(1.0 + (1.613e-5)*t*sa - (3.656e-3)*sa + (3.21e-5)*sq(sa) - (4.232e-7)*sa*sq(sa));

	einf = 4.9;

	relax = ((1.1109e-10) - (3.824e-12)*t + (6.398e-14)*sq(t) - (5.096e-16)*t*sq(t))* 
	(1.0 + (2.282e-5)*t*sa - (7.638e-4)*sa - (7.760e-6)*sq(sa) + (1.105e-8)*sa*sq(sa));

	sigma1 = 0.0467 + 0.2204*rb -0.4111*sand + 0.6614*clay;
	efw = (e0-einf)/(1+sq(relax*freq)) + einf;

	sigma = sigma1*(rs-rb)/rs/mv;
	efwp = (relax*freq*(e0-einf))/(1.0+sq(relax*freq)) + sigma/2.0/PI/freq/eps_0;
	
	emp = pow((1.0+rb/rs*(pow(es,alpha) - 1.0) + pow(mv,bp)*pow(efw,alpha) - mv),(1.0/alpha));
	empp = pow((pow(mv,bpp)*pow(efwp,alpha)),(1.0/alpha));
	
	emp = 1.15*emp-0.68;
	return cmp(emp,-1.0*empp);
}

cmp soilmvtoeps_hallikenen(int stype, double lambda,double mv, double sand, double clay){
	/*Data from Hallikenen's paper
	 s is sand, c is clay */
	 /*
	 soil types:
	 1: sandy loam, s = 0.515, c = 0.135
	 2: loam, s = 0.42, c = 0.085
	 3: silt loam, s = 0.306, c = 0.135
	 4: silt loam(2), s = 0.172, c = 0.19
	 5: silty clay, s = 0.05, c = 0.474
	 */
	double s,c;
	switch(stype)
	{
		case 0:
		s = sand; c = clay; break;
		case 1:
		s = 0.515; c = 0.135; break;
		case 2:
		s = 0.420; c = 0.085; break;
		case 3:
		s = 0.306; c = 0.135; break;
		case 4:
		s = 0.172; c = 0.190; break;
		case 5:
		s = 0.050; c = 0.474; break;
		case 6:
		//AG002 and it has 45.5% sand, 41.4% silt, 13.4% clay
		s = 0.455; c = 0.134; break;
		default:
		cout<<"** Unknown soil type **"<<endl; break;
	}
	double a0,a1,a2,b0,b1,b2,c0,c1,c2,//real coefficients
	aa0,aa1,aa2,bb0,bb1,bb2,cc0,cc1,cc2;//imag coefficients
	
	if( abs(100.0*lambda - 21.4) < 3.0 ){
	//This is at 1.4GHz -> 21.4cm use for 24cm also
	a0=2.862,a1=-0.012,a2=0.001,b0=3.803,b1=0.462,b2=-0.341,c0=119.006,c1=-0.5,c2=0.633;
	aa0=0.356,aa1=-0.003,aa2=-0.008,bb0=5.507,bb1=0.044,bb2=-0.002,cc0=17.753,cc1=-0.313,cc2=0.206;
	}
	else if( abs(100.0*lambda - 5.0) < 1.2 ){
	//This is at 6 GHz ->  5cm
	a0=1.993,a1=0.002,a2=0.015,b0=38.086,b1=-0.176,b2=-0.633,c0=10.720,c1=1.256,c2=1.522;
	aa0=-0.123,aa1=0.002,aa2=0.003,bb0=7.502,bb1=-0.058,bb2=-0.116,cc0=2.942,cc1=0.452,cc1=0.543;
	}
	else if( abs(100.0*lambda - 2.5) < 1.0 ){
	//This is at 12 GHz -> 2.5cm
	a0=2.2,a1=-0.001,a2=0.012,b0=26.473,b1=0.013,b2=-0.523,c0=34.333,c1=0.284,c2=1.062;
	aa0=-0.142,aa1=0.001,aa2=0.003,bb0=11.868,bb1=-0.059,bb2=-0.225,cc0=7.817,cc1=0.570,cc2=0.801;
	}
	else{
		cout<<"** Soil dielectric undefined at "<<lambda*100.0<<" cm wavelength"<<endl;
	}
	
	if( s>1||c>1||mv>1||s<-0.001||c<-0.001||mv<0 )
	{cout<<"** Unreal soil parameters being asked, mv "<<mv<<" **"<<endl;}
	
	//	cout<<(a0+a1*s+a2*c)<<" "<<(b0+b1*s+b2*c)<<" "<<(c0+c1*s+c2*c)<<endl;
	//	cout<<(aa0+aa1*s+aa2*c)<<" "<<(bb0+bb1*s+bb2*c)<<" "<<(cc0+cc1*s+cc2*c)<<endl;
	
	return ((a0+a1*s+a2*c)+(b0+b1*s+b2*c)*mv+(c0+c1*s+c2*c)*mv*mv) - 
	I*((aa0+aa1*s+aa2*c)+(bb0+bb1*s+bb2*c)*mv+(cc0+cc1*s+cc2*c)*mv*mv);
}

//This function checks if mv is unphysical and clipping is applied
void check_mv(int mtype, double ymax, double ymin, vector<vector<double> >* mprofiles,
double mvmin, double mvmax){
	
	double dy = (ymax-ymin)/100, y, mv; int i,j;
	
//	if(mtype == 1 && (*mprofile).size() > 1 ){
	if(mtype == 1 ){
		for(i=0; i<(*mprofiles).size(); i++){	
			bool minclip = false, maxclip = false; y = 0;
			while(y<ymax-ymin){
				mv = 0;
				for(j=0; j<(*mprofiles)[i].size(); j++){
					mv += (*mprofiles)[i][j] * pow(y,double(j));
				}
				if(mv < mvmin){
					minclip = true;
				}
				if(mv > mvmax){
					maxclip = true;
				}
				y += dy;
			//	cout<<y<<" "<<mv<<endl;
			}
			if(minclip){
				cout<<"*Caution: applying min mv "<<mvmin<<" * in profile "<<i+1<<endl;
			}
			if(maxclip){
				cout<<"*Caution: applying max mv "<<mvmax<<" * in profile "<<i+1<<endl;
			}
		}
	}	
}

//This function determines the eps at a given point
cmp get_eps(node p, double lambda, int mtype,cmp stat_eps,double refht,vector<vector<layer> >* Lyrs, 
vector<vector< double> >* mprofiles,double mvmin, double mvmax,double soilw,vector<int> seq){
		
	int i,j; double y,mv,ly,sa,cl,rb,x,dx; 
	//note that if refht is maxht, y is always +ve, starting from 0
	//else if refht is meanht, mv is taken to be const till meanht.
	y = refht - p.ny; 
	mv = 0; //initialize
	x = p.nx + soilw/2.0, dx = soilw/seq.size(); //offset to make x +ve to bin it properly
	i = int(x/dx);//this is the index of seq that the node belongs to
	if(i < 0 || x < 0){
	//	cout<<"** Caution, out of bounds x-region in substrate "<<p.nx<<","<<p.ny<<" *"<<endl;
		i = 0;
	}
	else if( i>= seq.size() ){
	//	cout<<"** Caution, out of bounds x+region in substrate "<<p.nx<<","<<p.ny<<" *"<<endl;
		i = seq.size() -1 ;
	}
	
	if( mtype == 2 ){ //This for ice
	//	stat_eps = cmp(3.1,-3.1/100000);
		return stat_eps;
	}
	else if( mtype == 3 ){ // This for wetice
		return icemvtoeps(lambda, (*mprofiles)[0][0]);
	}
	else if( mtype == 4 ){ // This is for misc
		return stat_eps;
	}
	else if( mtype == 1 ){ //This for soil
		if( y<0 ){
			mv = (*mprofiles)[seq[i]][0];
		}
		else{
			for(j=0; j<(*mprofiles)[seq[i]].size(); j++){
				mv += (*mprofiles)[seq[i]][j] * pow(y,double(j));
			}
		}
		//Check if mv values are physical
		//clip them to min/max values
		if( mv > mvmax ){
			mv = mvmax;
		}
		if( mv < mvmin ){
			mv = mvmin;
		}
		//end of physical check
		ly = 0;
		//only the last layer can have a "0" value for 
		//layer thickness, which implies infinity
		for(j=0; j<(*Lyrs)[seq[i]].size(); j++){
			ly += (*Lyrs)[seq[i]][j].th;
			if( y < ly || abs((*Lyrs)[seq[i]][j].th) < 1e-5){
				//This will also admit y<0 regions
				sa = (*Lyrs)[seq[i]][j].sa;
				cl = (*Lyrs)[seq[i]][j].cl;
				rb = (*Lyrs)[seq[i]][j].rb;
				break;
			}
		}
		//	cout<<seq[i]<<","<<x<<","<<y<<","<<mv<<","<<sa<<","<<cl<<" "<<rb<<" "<<endl;
		if( 100.0*lambda > 26.0 ){
			return soilmvtoeps_peplinsky(lambda,mv,sa,cl,rb);
		}
		else{
			return soilmvtoeps_hallikenen(0,lambda,mv,sa,cl);
		}
		//End of soil choices
	}
}

//This function sets the eps for the entire substrate, possibly heterogeneous
void set_eps(double lambda, int mtype, cmp stat_eps, double maxht, vector<vector< layer> >* Lyrs, vector<vector< double> >* mprofiles,
vector<int> * itf_el, vector<int> * mat_el, vector<node> *nlist, vector<element> *elist,
double mvmin,double mvmax,double soilw, vector<int> seq){	

	int i, elno, j; node pm;
	for(i=0; i<(*itf_el).size()+(*mat_el).size(); i++){
		if( i<(*itf_el).size() ){
			elno = (*itf_el)[i];
		}
		else{
			elno = (*mat_el)[i-(*itf_el).size()];
		}
		
		if( (*elist)[elno].eltype == 's' ){
			pm = (*elist)[elno].elcenter(nlist);
			(*elist)[elno].eleps = get_eps(pm,lambda,mtype,stat_eps,maxht,Lyrs,mprofiles,mvmin,mvmax,soilw,seq);
		}
		else{
			(*elist)[elno].eleps = cmp(1.0,0);
		}
	}
}

//This function sets the eps for only the itf layer
void set_eps(double lambda, int mtype, cmp stat_eps,double maxht, vector<vector< layer> >* Lyrs, vector<vector< double> >* mprofiles,
vector<int> * itf_el, vector<node> *nlist, vector<element> *elist,
double mvmin,double mvmax,double soilw, vector<int> seq){	

	vector<int> dummy(0);
	set_eps(lambda,mtype,stat_eps,maxht,Lyrs,mprofiles,itf_el,&dummy,nlist,elist,mvmin,mvmax,soilw,seq);
}


//This function removes white space in a string after = sign
string stripe_spaces_after_equals(string line){
	if( line.find("=") == string::npos ){
		cout<<"** Found no = sign in string, possible error in config file **"<<endl;
	}
	// extract string after = sign
	string tmp = string(line,1+line.find("="),string::npos); 
	// remove white spaces 
	tmp.erase(remove(tmp.begin(),tmp.end(),' '),tmp.end()); 
	return tmp;
}

//This function finds out how many delimiters of a given kind
//are in the give file. The delimiter must take up an entire line
int howmany_delims(ifstream * in, string p){
	string line;
	int i=0;
	if( (*in).is_open() ){
		while( (*in).good() ){
			getline((*in),line);
			if( string(line) == p ){
				i++;
			}
		}
	}
	else{
		cout<<"** Unable to open rsurface file to count delimiters **";
		return 0;
	}
	//clear flags and rewind the file
	(*in).clear(); (*in).seekg(0);
	return i;
}

/*
//This function is a wrapper that reads the config file for a substrate
//homogeneous in y, homogeneous in x
int input_configuration(int argc, char* argv[], double * lambda, ifstream * in, vector<double> * thetas, 
	int * mtype, cmp * stat_eps, int * mcarlo,double * rgh, double * corl,ifstream * rs,
	int * maxrs){
	//Initialize
	int i,stype=0; double sand=0, clay=0, mv_const = 0;
	(*mcarlo) = 0; (*lambda) = 0; (*rgh) = 0; (*corl) = 0;
	bool is_sand = false, is_clay = false; (*maxrs) = 0;
	//Read the config file
	ifstream config(argv[1]);
	if( config.is_open() ){
		cout<<"Config: '"<<argv[1]<<"'";
	}
	else{
		cout<<"** Error opening config file **"<<endl;
		return -1;
	}
	
	bool onetime = true;
	string line,tmp; 
	stringstream buffer, sline;
	buffer<<config.rdbuf();
	
	while(buffer){
		while(getline(buffer,line,'\n')){
			if(line.size()==0){
				continue; //if there is an empty line
			}
			if( string(line,0,1) == "#"){
				if( onetime ){
					cout<<", description: '"<<string(line,1,string::npos)<<"'"<<endl;
					onetime = false;
				}
				continue; // all comments begin with #
			}
			if( line.find("lambda") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				(*lambda) = atof(tmp.c_str());
				cout<<"Lambda="<<100.0*(*lambda)<<"cm, ";
			}
			if( line.find("meshfile") != string::npos ){
				
				tmp = stripe_spaces_after_equals(line);					
				cout<<"Mesh: '"<<tmp<<"'";
				(*in).open(tmp.c_str());					
			}
			if( line.find("theta") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				
				sline<<tmp;
				while(getline(sline,tmp,',')){ 
					// reusing tmp since sline has it
					(*thetas).push_back(dtor(atof(tmp.c_str())));
				}
				sline.clear();sline.str("");
				cout<<", Inc angle(s)=";
				for(i=0;i<(*thetas).size()-1;i++){
					cout<<rtod((*thetas)[i])<<",";
				}
				cout<<rtod((*thetas)[i]);	
			}
			if( line.find("substrate") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				if( tmp == "soil" ){
					(*mtype) = 1;
				}
				else if( tmp == "ice" ){
					(*mtype) = 2;
				}
				else if( tmp == "wetice" ){
					(*mtype) = 3;
				}
				else if( tmp == "misc" ){
					(*mtype) = 4;
				}
				else{
					(*mtype) = 0;
				}
				cout<<", Substrate: '"<<tmp<<"'";
			}
			if( (*mtype) == 1 ){
				if( line.find("soiltype") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					stype = atoi(tmp.c_str());
				}
				else if( line.find("sand") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					sand = atof(tmp.c_str()); 
					is_sand = true;
				}
				else if( line.find("clay") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					clay = atof(tmp.c_str());
					is_clay = true;
				}
			}
			if( (*mtype) == 1 || (*mtype) == 3 ){
				if( line.find("moisture") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					mv_const = atof(tmp.c_str());
					cout<<", moisture="<<mv_const<<"";
				}
			}
			if( line.find("mcarlo") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				(*mcarlo) = atoi(tmp.c_str());
			//	cout<<", instances="<<(*mcarlo);
			}
			if( line.find("roughness") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				(*rgh) = atof(tmp.c_str());
			}
			if( line.find("correlation") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				(*corl) = atof(tmp.c_str());
			}
			if( line.find("surfaces") != string::npos){
				tmp = stripe_spaces_after_equals(line);
				cout<<", surfaces: '"<<tmp<<"'";
				(*rs).open(tmp.c_str());
				(*maxrs) = howmany_delims(rs,"#");
			}
			//Add any more config file read-ins here
			//
		}
		line.clear();
	}
	
	//Set mcarlo to be the min of specified mcarlo, and available rsurfs
	if( (*maxrs) != 0 ){
		(*mcarlo) = Min( (*maxrs), (*mcarlo));
	}
	//Check if a required item is not specified.
	if( abs((*lambda)) <1e-3 ){
		cout<<endl<<"** Error wavelength not specified **"<<endl;
	}
	if( !(*in).is_open() ){
		cout<<endl<<"** Error in opening mesh file **"<<endl;
		return -1;
	}
	if( (*mtype) == 0 ){
		cout<<endl<<"** Error as substrate type not specified **"<<endl;
		return -1;
	}
	if( (*mtype) == 1 && stype == 0 && !is_sand && !is_clay ){
		cout<<endl<<"** Error as soil type or composition not specified **"<<endl;
		return -1;
	}
	if( (*mtype) == 1 && (abs(mv_const) < 1e-3 || mv_const>0.999) ){
		cout<<endl<<"** Error as soil moisture is unphysical **"<<endl;
		return -1;
	}
	if( (*mcarlo) == 0 && (*maxrs) == 0){
		cout<<endl<<"* Warning, number of Monte Carlo runs not specified, set to 1"<<endl;
		(*mcarlo) = 1;
	}
	cout<<", instances="<<(*mcarlo);
//	if( abs((*rgh))<1e-4 ||abs((*corl))<1e-3 ){
	if( abs((*corl))<1e-3 ){
		cout<<endl<<"** Error as roughness/correlation not specified **"<<endl;
		return -1;
	}
	//End of checks, now fix the dielctric constant
	(*stat_eps) = set_eps((*lambda),(*mtype),stype,mv_const,sand,clay);
	//test section
//	(*stat_eps) = cmp(1,0);
//	cout<<"**TESTING EPS IS VAC**"<<endl;
	//end section
	cout<<", const_eps="<<(*stat_eps)<<" "<<endl;
	return 0;
}

//This function is a wrapper that reads the config file for a substrate
//heterogeneous in y, homogeneous in x
int input_configuration(int argc, char* argv[], double * lambda, ifstream * in, vector<double> * thetas, 
	int * mtype, int * mcarlo,double * rgh, double * corl,ifstream * rs,
	int * maxrs, vector<layer> * lyrs, vector<double> * mprofile){
	//Initialize
	int i, lno = 0; double sand=0, clay=0, mv_const = 0;
	(*mcarlo) = 0; (*lambda) = 0; (*rgh) = 0; (*corl) = 0;
	(*maxrs) = 0;
	//Read the config file
	ifstream config(argv[1]);
	if( config.is_open() ){
		cout<<"Config: '"<<argv[1]<<"'";
	}
	else{
		cout<<"** Error opening config file **"<<endl;
		return -1;
	}
	
	bool onetime = true;
	string line,tmp; 
	stringstream buffer, sline;
	buffer<<config.rdbuf();
	
	while(buffer){
		while(getline(buffer,line,'\n')){
			if(line.size()==0){
				continue; //if there is an empty line
			}
			if( string(line,0,1) == "#"){
				if( onetime ){
					cout<<", description: '"<<string(line,1,string::npos)<<"'"<<endl;
					onetime = false;
				}
				continue; // all comments begin with #
			}
			if( line.find("lambda") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				(*lambda) = atof(tmp.c_str());
				cout<<"Lambda="<<100.0*(*lambda)<<"cm, ";
			}
			if( line.find("meshfile") != string::npos ){
				
				tmp = stripe_spaces_after_equals(line);					
				cout<<"Mesh: '"<<tmp<<"'";
				(*in).open(tmp.c_str());					
			}
			if( line.find("theta") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				
				sline<<tmp;
				while(getline(sline,tmp,',')){ 
					// reusing tmp since sline has it
					(*thetas).push_back(dtor(atof(tmp.c_str())));
				}
				sline.clear();sline.str("");
				cout<<", Inc angle(s)=";
				for(i=0;i<(*thetas).size()-1;i++){
					cout<<rtod((*thetas)[i])<<",";
				}
				cout<<rtod((*thetas)[i]);	
			}
			if( line.find("substrate") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				if( tmp == "soil" ){
					(*mtype) = 1;
				}
				else if( tmp == "ice" ){
					(*mtype) = 2;
				}
				else if( tmp == "wetice" ){
					(*mtype) = 3;
				}
				else if( tmp == "misc" ){
					(*mtype) = 4;
				}
				else{
					(*mtype) = 0;
				}
				cout<<", Substrate: '"<<tmp<<"'";
			}
			if( (*mtype) == 1 ){
				if( line.find("layers") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					(*lyrs).resize(atoi(tmp.c_str()));
				}
				else if( line.find("depth") != string::npos ){
					lno++;//update the layer number
					tmp = stripe_spaces_after_equals(line);
					(*lyrs)[lno-1].th = atof(tmp.c_str());
				}
				else if( line.find("sand") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					(*lyrs)[lno-1].sa = atof(tmp.c_str());
				}
				else if( line.find("clay") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					(*lyrs)[lno-1].cl = atof(tmp.c_str());
				}
				else if( line.find("density") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					(*lyrs)[lno-1].rb = atof(tmp.c_str());
				}
			}
			if( (*mtype) == 1 || (*mtype) == 3 ){
				if( line.find("moisture") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					istringstream word(tmp);
					while(getline(word,tmp,',')){
						(*mprofile).push_back( double(atof(tmp.c_str())) );
					}
					//Special case is if size == 3, I use this to signify
					//sm sensor readings at three depths, not the polyn coeff
					if( (*mprofile).size() == 3){
						vector<double> t = (*mprofile);
						(*mprofile) = quadratic_sm( t );
					}
					cout<<", moisture=(";
					for(i=0; i<(*mprofile).size(); i++){
						cout<<(*mprofile)[i]<<",";
					}
					cout<<'\b'<<")";
				}
			}
			if( line.find("mcarlo") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				(*mcarlo) = atoi(tmp.c_str());
			//	cout<<", instances="<<(*mcarlo);
			}
			if( line.find("roughness") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				(*rgh) = atof(tmp.c_str());
			}
			if( line.find("correlation") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				(*corl) = atof(tmp.c_str());
			}
			if( line.find("surfaces") != string::npos){
				tmp = stripe_spaces_after_equals(line);
				cout<<", surfaces: '"<<tmp<<"'";
				(*rs).open(tmp.c_str());
				(*maxrs) = howmany_delims(rs,"#");
			}
			//Add any more config file read-ins here
			//
		}
		line.clear();
	}
	
	//Set mcarlo to be the min of specified mcarlo, and available rsurfs
	if( (*maxrs) != 0 ){
		(*mcarlo) = Min( (*maxrs), (*mcarlo));
	}
	//Print layer info
	if( (*mtype) == 1 ){
		cout<<", layers (z,sa,si,cl,rb):";
		for(i=0; i<(*lyrs).size();i++){
			(*lyrs)[i].si = 1.0 - (*lyrs)[i].sa - (*lyrs)[i].cl	;
			cout<<"("<<(*lyrs)[i].th<<","<<(*lyrs)[i].sa<<","<<(*lyrs)[i].si
			<<","<<(*lyrs)[i].cl<<","<<(*lyrs)[i].rb<<"),";
		}
		cout<<'\b';
	}
	//Check if a required item is not specified.
	if( abs((*lambda)) <1e-3 ){
		cout<<endl<<"** Error wavelength not specified **"<<endl;
	}
	if( !(*in).is_open() ){
		cout<<endl<<"** Error in opening mesh file **"<<endl;
		return -1;
	}
	if( (*mtype) == 0 ){
		cout<<endl<<"** Error as substrate type not specified **"<<endl;
		return -1;
	}
	if( (*mtype) == 1 && (*lyrs).size() == 0 ){
		cout<<endl<<"** Error as soil layer composition not specified **"<<endl;
		return -1;
	}
	if( (*mcarlo) == 0 && (*maxrs) == 0){
		cout<<endl<<"* Warning, number of Monte Carlo runs not specified, set to 1"<<endl;
		(*mcarlo) = 1;
	}
	cout<<", instances="<<(*mcarlo)<<endl;
	if( abs((*corl))<1e-3 ){
		cout<<endl<<"** Error as roughness/correlation not specified **"<<endl;
		return -1;
	}
	//end section
	return 0;
}
* */

//This function is a wrapper that reads the config file for a substrate
//heterogeneous in y, heterogeneous in x
int input_configuration(int argc, char* argv[], double * lambda, ifstream * in, vector<double> * thetas, 
vector<cmp> * tissues, ofstream * outp
){
	//Initialize
	int i,j,lno; 
	//Read the config file
	ifstream config(argv[1]);
	if( config.is_open() ){
		cout<<"Config: '"<<argv[1]<<"'";
	}
	else{
		cout<<"** Error opening config file **"<<endl;
		return -1;
	}
	
	bool onetime = true, real_specified = false, imag_specified = false;
	cmp tis_eps;
	string line,tmp; 
	stringstream buffer, sline;
	buffer<<config.rdbuf();
	
	while(buffer){
		while(getline(buffer,line,'\n')){
			if(line.size()==0){
				continue; //if there is an empty line
			}
			if( string(line,0,1) == "#"){
				if( onetime ){
					cout<<", description: '"<<string(line,1,string::npos)<<"'"<<endl;
					onetime = false;
				}
				continue; // all comments begin with #
			}
			if( line.find("lambda") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				(*lambda) = atof(tmp.c_str());
				cout<<"Lambda="<<100.0*(*lambda)<<"cm, ";
			}
			if( line.find("meshfile") != string::npos ){
				tmp = stripe_spaces_after_equals(line);					
				cout<<"Mesh: '"<<tmp<<"'";
				(*in).open(tmp.c_str());					
			}
			if( line.find("theta") != string::npos ){
				tmp = stripe_spaces_after_equals(line);
				
				sline<<tmp;
				while(getline(sline,tmp,',')){ 
					// reusing tmp since sline has it
					(*thetas).push_back(dtor(atof(tmp.c_str())));
				}
				sline.clear();sline.str("");
				cout<<", Inc angle(s)=";
				for(i=0;i<(*thetas).size()-1;i++){
					cout<<rtod((*thetas)[i])<<",";
				}
				cout<<rtod((*thetas)[i]);	
			}
			//Read in tissue dielectric constants
			if( line.find("real") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					tis_eps.real() = double(atof(tmp.c_str()));
					cout<<", tissue"<<(*tissues).size()<<"=("<<tis_eps.real();
					real_specified = true;
				}
				else if( line.find("imag") != string::npos ){
					tmp = stripe_spaces_after_equals(line);
					tis_eps.imag() = double(atof(tmp.c_str()));
					cout<<","<<tis_eps.imag()<<")";
					imag_specified = true;
			}
			if(real_specified && imag_specified){
				(*tissues).push_back(tis_eps);
				real_specified = false; imag_specified = false;
			}
			
			//Add any more config file read-ins here
			//
			if( line.find("outputfile") != string::npos ){
				tmp = stripe_spaces_after_equals(line);					
				cout<<",  output: '"<<tmp<<"'";
				(*outp).open(tmp.c_str());	
			}
		}
		line.clear();
	}
	
	
	//Check if a required item is not specified.
	if( abs((*lambda)) <1e-3 ){
		cout<<endl<<"** Error wavelength not specified **"<<endl;
	}
	if( !(*in).is_open() ){
		cout<<endl<<"** Error in opening mesh file **"<<endl;
		return -1;
	}
	if( (*tissues).size() == 0 ){
		cout<<endl<<"** Tissues not specified **"<<endl;
		return -1;
	}
	else{
		for(i=0; i<(*tissues).size(); i++){
			if( (*tissues)[i].imag() > 0.00001 ){
				cout<<endl<<"*** Unphysical tissue eps ***"<<endl;
				return -1;
			}
		}
	}
	//end section
	return 0;
}

int create_fem_datastructures(ifstream * in, double lambda, vector<node> * nodelist, vector<edge> * edgelist, 
	vector<element> * elementlist, vector<int> * rbc_ed, vector<int> * diel_ed,vector<int> * diel_el,vector<int> * top_ed,vector<int> * top_el, vector<int> * aa_el,vector<int> * out_ed, vector< vector<int> > * tiss_els
	){
	
	int No=-1,Ne=-1,Ns=0,i,j,k,n, orig, trans,no_aa=0, tiss_no=0; double px,py;
	string line,tmp; 
	stringstream buffer, sline;
	bool read_no = false, read_el = false, read_rbc = false,
	read_diel = false, read_aa = false, str_len,read_top=false,read_out = false, read_tiss=false;
	vector<int> temp;
//	vector<bool> read_tiss( (*tiss_els).size(), false);
	
	//First pass through the file to figure out size of lists and allocate memory.
	buffer<<(*in).rdbuf();
	while(buffer)
	{
		while(getline(buffer,line,'\n'))
		{
			if( string(line,0,5) == "*NODE")
			{read_no=true;continue;}
			else if( string(line,0,8) == "*ELEMENT")
			{read_no=false;read_el=true;continue;}
			else if( string(line,0,6) == "*ELSET")
			{read_el=false;break;}
			
			if(read_no){No++;}
			else if(read_el){Ne++;}
		}
	}
	if(No == -1){
		cout<<"**Error in file/no nodes read**"<<endl; return -1;
	}
	//rewind the buffer after clearing error flags
	buffer.clear();	buffer.seekg(0);
	//Allocate memory
	(*nodelist).reserve(No); (*nodelist).push_back(node());
	(*elementlist).reserve(Ne); (*elementlist).push_back(element());
	(*edgelist).reserve(int(1.75*Ne)); (*edgelist).push_back(edge());
	
	while(buffer){
		while(getline(buffer,line,'\n')){
			if(line.size()>8) //to prevent string exception in ELSET line
			{str_len = true;}
			else
			{str_len = false;}
			
			if( string(line,0,5) == "*NODE")
			{
				read_no = true;continue;
			}
			else if( string(line,0,8) == "*ELEMENT")
			{
				read_el = true;continue;
			}
			else if(str_len)
			{
				string elset = string(line,8,9);
				if( elset == "ELSET=RBC" )
				{
					read_rbc = true; 
					orig = atoi( string(line.end()-1,line.end()).c_str() ) - 1;//Which edge
					trans = (orig+2)%3; continue;
				}
				else if( elset == "ELSET=DIE" )
				{
					read_diel = true; 
					orig = atoi( string(line.end()-1,line.end()).c_str() ) - 1;//Which edge
					trans = (orig+2)%3; continue;
				}
				else if( elset == "ELSET=TOP" )
				{
					read_top = true; 
					orig = atoi( string(line.end()-1,line.end()).c_str() ) - 1;//Which edge
					trans = (orig+2)%3; continue;
				}
				else if( elset == "ELSET=OUT" )
				{
					read_out = true; 
					orig = atoi( string(line.end()-1,line.end()).c_str() ) - 1;//Which edge
					trans = (orig+2)%3; continue;
				}
				else if( elset == "ELSET=AAB" )
				{
					read_aa = true;
					continue;
				}
				else if( string(line,8,6) == "ELSET=" )
				{
					read_tiss = true;
					continue;
				}
			
				else if( string(line,0,8) == "*SURFACE")
				{
					string sface = string(line.end()-3,line.end());
					if( sface == "RBC" )
					{read_rbc = false;}
					else if( sface == "DIE" )
					{read_diel = false;}
					else if( sface == "TOP" )
					{read_top = false;}
					else if( sface == "OUT" )
					{read_out = false;}
					else if( sface == "AAB" )
					{read_aa = false;}
					else
					{read_tiss = false; tiss_no++;
					(*tiss_els).push_back(temp);
					cout<< "tisse # (size) :" <<tiss_no<<"("<< temp.size()<<")"<<endl;
					temp.clear();
					}
					//The above assumes that there are NO other sidesets other than
					//RBC, DIE, TOP, OUT, AAB
					continue;
				}
				
			}
			else if( string(line,0,2) == "**")
			{
				if(read_no)
				{read_no = false; No=(*nodelist).size()-1; cout<<"Created: Nodes("<<No<<")";cout.flush();}
				else if(read_el)
				{read_el = false; Ne=(*elementlist).size()-1; cout<<", elements("<<Ne<<")";cout.flush();
					//node and element lists are created, now create edge list
					
					//* create unsorted indicator list
					vector< vector<int> >indicator(3*Ne,vector<int>(4,0));
					for(i=0; i<indicator.size(); i++)
					{
						indicator[i][1] = Min((*elementlist)[i/3+1].nodes[i%3],(*elementlist)[i/3+1].nodes[(i+1)%3]); //node n1
						indicator[i][2] = Max((*elementlist)[i/3+1].nodes[i%3],(*elementlist)[i/3+1].nodes[(i+1)%3]); //node n2
						indicator[i][3] = i/3+1; // element number
						indicator[i][0] = indicator[i][1] * indicator[i][2]; //the indicator
					}					
					//* sort the indicator list using shell sort
					vector<int> tmp_ind(4,0);
					vector<int> cols; int h=1; n = indicator.size();
					while(h<n)
					{
						cols.push_back(h); h = pow(2,1+cols.size())-1; 
					}
					for (k=cols.size()-1; k>=0; k--)
					{
						h=cols[k];
						for (i=h; i<n; i++)
						{
							tmp_ind = indicator[i];
							j=i;
							while (j>=h && indicator[j-h][0]>tmp_ind[0])
							{
								indicator[j]=indicator[j-h];
								j=j-h;
							}
							indicator[j]=tmp_ind;
						}
					}
					cols.clear();
					//* number the edges and create (*edgelist)
					bool add_edge=false; vector<int> same_ind(0,0);
					for(i=0; i<indicator.size(); i++)
					{
						for(j=i-1; j>=0; j--) //store all the entries with the same indicator
						{
							if(indicator[j][0] == indicator[i][0])
							{same_ind.push_back(j);}
							else
							{break;}
						}
						add_edge = true;
						for(j=0; j<same_ind.size(); j++)
						{
							if(indicator[i][1] == indicator[same_ind[j]][1])//see if one node is common
							{
								for(k=(*edgelist).size()-1;k>=0;k--)
								{
									if((*edgelist)[k] == edge(indicator[i][1],indicator[i][2]))
									{
										(*edgelist)[k].els.push_back(indicator[i][3]); break;
									}
								}
								add_edge = false; break;
							}
						}
						if( add_edge )
						{
							edge ed(indicator[i][1],indicator[i][2]);
							ed.els.reserve(2); ed.els.push_back(indicator[i][3]);
							ed.edbtype = 'i';// i for interior
							ed.setlength(&(*nodelist));
							(*edgelist).push_back(ed);
						}
						same_ind.clear();
					}
					indicator.clear();
					
					//* add edge information to element list
		 			for(i=1; i<(*edgelist).size();i++)
					{
						for(j=0 ; j<(*edgelist)[i].els.size(); j++) //run over elements of this edge
						{
							n = (*edgelist)[i].els[j];
							for(k=0; k<3; k++)//run over nodes of this element
							{
								edge ek((*elementlist)[n].nodes[k],(*elementlist)[n].nodes[(k+1)%3]);
								if( ek == (*edgelist)[i] ) //assign the edge to the element
								{
									(*elementlist)[n].edges[k] = i;
								}
							}
						}
					}  
					//* Done */
					Ns=(*edgelist).size()-1;cout<<", edges("<<Ns<<").";cout.flush();
				}
				continue;
			}
			sline.clear();sline.str (""); // <-- VERY imp
			sline << line;
			if(read_no)
			{
				//	Multiply by lambda to get real coordinates
				getline(sline,tmp,',');//this is the node number. pass
				getline(sline,tmp,',');//this is the x-coodrinate.
				px = lambda * (atof(tmp.c_str())); // float to double conversion.
				getline(sline,tmp,',');//this is the y-coodrinate.
				py = lambda * (atof(tmp.c_str()));
				(*nodelist).push_back(node(px,py));
			}
			if(read_el)
			{
				getline(sline,tmp,',');//this is element number. pass.
				element el(0,0,0);
				for(i=0;i<3;i++)
				{
					getline(sline,tmp,',');
					el.nodes[i] = atoi(tmp.c_str());
				}
				el.elarea = area_e(&el,&(*nodelist));
				el.eleps = cmp(1.0,0.0); el.elmu = cmp(1.0,0.0); el.eltype = 'v';//by default make all vacuum
				(*elementlist).push_back(el);					
			}
			if(read_rbc)
			{
				while(getline(sline,tmp,','))
				{
					i = atoi(tmp.c_str());
					(*rbc_ed).push_back((*elementlist)[i].edges[orig]);
	//				(*rbc_el).push_back(i);
					(*edgelist)[(*elementlist)[i].edges[orig]].edbtype = 'r';
				}
			}
			if(read_diel)
			{
				while(getline(sline,tmp,','))
				{
					i = atoi(tmp.c_str());
					(*diel_ed).push_back((*elementlist)[i].edges[orig]);
				//	(*diel_el).push_back(i);
					(*edgelist)[(*elementlist)[i].edges[orig]].edbtype = 'd';
				}
			}
			if(read_top)
			{	
				while(getline(sline,tmp,','))
				{
					i = atoi(tmp.c_str());
					(*top_ed).push_back((*elementlist)[i].edges[orig]);
					(*top_el).push_back(i);
					(*edgelist)[(*elementlist)[i].edges[orig]].edbtype = 't'; 
				}
			}
			if(read_out)
			{
				while(getline(sline,tmp,','))
				{
					i = atoi(tmp.c_str());
					(*out_ed).push_back((*elementlist)[i].edges[orig]);
					//(*inn_el).push_back(i);
					(*edgelist)[(*elementlist)[i].edges[orig]].edbtype = 'o';
				}
			}
			if(read_aa)
			{
				while(getline(sline,tmp,','))
				{
					i = atoi(tmp.c_str());
					(*elementlist)[i].eltype = 'a';
					(*aa_el).push_back(i);
					no_aa++;
				}
			}
			
			if(read_tiss)
			{
				//cout<<"tissue #  :"<<tiss_no<<endl;
				while(getline(sline,tmp,','))
				{
					i = atoi(tmp.c_str());
					(*elementlist)[i].eltype = 's';
					(temp).push_back(i);					
				}
			}	

			line.clear();
		}
	}
	buffer.clear();buffer.str("");
	(*in).close();
	if(No == 1){
		cout<<"** Error in file/no nodes read **"<<endl; return -1;
	}
	//Remove duplicates from top_ed and diel_ed
	sort_unq_list(top_ed);
	sort_unq_list(diel_ed);
	return 0;
}

//This function sorts an int list of edges using insertion sort
int sort_list(char p, vector<int> * list, vector<node> * nodelist, vector<edge> * edgelist){
	int i,j,k,ii,pos;
	node pi, pj;
	double ti, tj;
	
	if( (*list).size() == 0 ){
		cout<<"** Warning: empty list sent for sorting **"<<endl;
		return 0;
	}
	
	for(j=1 ; j<(*list).size(); j++){
		pos = j;//Find the right place for [j]
		pj = (*edgelist)[(*list)[j]].midpoint(nodelist);
		for(i=j-1 ; i>=0 ; i--){
			pi = (*edgelist)[(*list)[i]].midpoint(nodelist);
			if( p == 'x' ){
				//Sorting according to x-coordinates
				ti = pi.nx; tj = pj.nx;
			}
			else if( p == 'y'){
				//Sorting according to y-coordinates
				ti = pi.ny; tj = pj.ny;
			}
			else if( p == 'p'){
				//Sorting according to phi-coordinate of edge
				ti = atan2(pi.ny,pi.nx);
				tj = atan2(pj.ny,pj.nx); 
			}
			else{
				cout<<"** Error: Unknown sorting criteria **"<<endl;
				return -1;
			}
			// "<" for increasing
			if( tj < ti ){
				pos = i;
			}
		}
		//Insert [j] at the right place [pos]
		ii = (*list)[j];
		for(k=j-1 ; k >= pos ; k--){
			(*list)[k+1] = (*list)[k];
		}
		(*list)[pos] = ii;
	}
	return 0;
}

//This function sets the element matrices for a single element
void set_element_matrices(int i, vector<node> * nodelist, vector<edge> * edgelist, vector<element> * elementlist){
	int j,k,jj,kk,sgn;
	double px,py; 
	node pi,pj,pk;
	vector <double> a(3,0.0),b(3,0.0),c(3,0.0),A(3,0.0),B(3,0.0),C(3,0.0);
		
	//i is the element number
	element ei = (*elementlist)[i];
	for(j=0; j<3 ;j++){
		px = (*nodelist)[ei.nodes[j]].nx; py = (*nodelist)[ei.nodes[j]].ny;
		//cyclic convention
		jj = (j+1) % 3; kk = (j+2) % 3;
		pj = (*nodelist)[ei.nodes[jj]]; pk = (*nodelist)[ei.nodes[kk]];
		a[j] = pj.nx * pk.ny - pj.ny * pk.nx;
		b[j] = pj.ny - pk.ny;
		c[j] = pk.nx - pj.nx;
	}
	for(j=0; j<3 ;j++){
		(*elementlist)[i].AA.resize(3);(*elementlist)[i].BB.resize(3);(*elementlist)[i].CC.resize(3);
		//node j and j+1 compose edge j
		jj = j ; kk = (j+1) % 3; 
		sgn = sign((*elementlist)[i].nodes[kk]-(*elementlist)[i].nodes[jj]);
		A[j] = a[jj] * b[kk] - a[kk] * b[jj];
		B[j] = c[jj] * b[kk] - c[kk] * b[jj];
		C[j] = a[jj] * c[kk] - a[kk] * c[jj];
		A[j] *= sgn; (*elementlist)[i].AA = A;
		B[j] *= sgn; (*elementlist)[i].BB = B;
		C[j] *= sgn; (*elementlist)[i].CC = C;
	}
}

//Overloaded function to set matrices for a list of elements
void set_element_matrices(vector<int> * list, vector<node> * nodelist, vector<edge> * edgelist, 
	vector<element> * elementlist){
	for(int i=0; i<(*list).size(); i++){
		set_element_matrices( (*list)[i],nodelist,edgelist,elementlist );
	}
}

//This sets the ad ab losses for the entire adab region use in fem.v20
void set_aa_loss(double lambda,int* tr_va, int* tr_aa, int* bl_va, int* bl_aa,
	vector<int> * list, vector<node> * nodelist,vector<element> * elementlist){
	
	double K0 = 2.0 * PI/ lambda, px, py, adaS=5.0/K0;
	cmp eps,loss; int i,j=0,k=0;
	node trva = (*nodelist)[*tr_va], traa = (*nodelist)[*tr_aa],
	blva = (*nodelist)[*bl_va], blaa = (*nodelist)[*bl_aa], pm;
	
	for(i=0; i<(*list).size(); i++){
		pm = (*elementlist)[(*list)[i]].elcenter(nodelist);
		
		px = abs(pm.nx) - trva.nx; py = abs(pm.ny) - trva.ny;
		
		if( abs(pm.ny) < trva.ny ){// x adab
			loss = cmp(0.0,-adaS*pow(px/(traa.nx-trva.nx),2.0));
		}
		else if( abs(pm.nx) < trva.nx ){//y adab
			loss = cmp(0.0,-adaS*pow(py/(traa.ny-trva.ny),2.0));
		}
		else{// corner case adab
			loss = cmp(0.0, -adaS*pow(sq(px/(traa.nx-trva.nx))+sq(py/(traa.ny-trva.ny)),1.0));
		}
		//Set the loss
	//	cout<<(*list)[i]<<" a "<<sqrt((*elementlist)[(*list)[i]].eleps)<<" ";
		(*elementlist)[(*list)[i]].eleps = (cmp(1.0,0) + loss); 
		(*elementlist)[(*list)[i]].elmu = (cmp(1.0,0) + loss); 
		j++;
	//	cout<<sqrt((*elementlist)[(*list)[i]].eleps)<<endl;
		
	}
	cout<<"Set AdAb loss for "<<j<<"+"<<k<<" elements.";
}
/*
//This function sets the loss for the elements of the adiabatic absorber in substrate only, use in v11
void set_aa_loss(double lambda, double pmlxin, double pmlyin, double pmlxout, double pmlyout, cmp stat_eps,
	vector<int> * list, vector<node> * nodelist, vector<element> * elementlist){
	
	double K0 = 2.0 * PI/ lambda, px, py, adaS=5.0/K0;
	cmp loss;
	for( int i=0; i<(*list).size(); i++){		
		element ei = (*elementlist)[(*list)[i]];
		node pm = ei.elcenter(nodelist);
		px = abs(pm.nx)-pmlxin; py = abs(pm.ny)-pmlyin;
		
		if( abs(pm.ny) < pmlyin ){ // x pml
			loss = cmp(0.0,-adaS*pow(px/(pmlxout-pmlxin),2.0));
	//		cout<<sq(px/(pmlxout-pmlxin))<<endl;
		}
		else if( abs(pm.nx) < pmlxin){ //y pml
			loss = cmp(0.0,-adaS*pow(py/(pmlyout-pmlyin),2.0));
	//		cout<<sq(py/(pmlyout-pmlyin))<<endl;
		}
		else{// corner
			loss = cmp(0.0,-adaS*pow(sq(px/(pmlxout-pmlxin))+sq(py/(pmlyout-pmlyin)),1.0));
	//		cout<<sq(px/(pmlxout-pmlxin))+sq(py/(pmlyout-pmlyin))<<endl;
		}
		//Set the loss
		(*elementlist)[(*list)[i]].eleps = sq(sqrt(stat_eps) + loss);
		//Test, set no loss
	//	(*elementlist)[(*list)[i]].eleps = stat_eps;
	}
	cout<<"Set AdAb loss for "<<(*list).size()<<" elements.";
//	cout<<"** TESTING AA NO LOSS **"<<endl;
}
*/

		

//Returns the angle of the edge
double edge_angle(edge ed, vector<node> *nodelist){
	node pi = (*nodelist)[ed.ni], pj = (*nodelist)[ed.nj];
	return atan2(pj.ny - pi.ny, pj.nx - pi.nx);
}

//This function get and sets the extents of the mesh region with fem.v12
//assumes inn_ed is in the substrate region, and out_ed in the air
void get_mesh_extents(vector<int> * out_ed,vector<int> *rbc_ed,
 int *tr_va,int *tr_aa,int *bl_va, int*bl_aa, double lambda,vector<vector <int> > * EAL, 
 vector<node> * nodelist, vector<edge> * edgelist, vector<element> * elementlist){
	
	/* First make a list of the nodes in the data sets out_ed
	 * and rbc_ed, then sort and remove duplicate nodes.
	 * For each node in this list of 
	 * edges (each tagged "o"), examine EAL to see if any other edge of 
	 * that node is also tagged the same. If not, we are at the end.
	 * Determine the other nodes by seeing if the two similarly tagged
	 * edges are at right angles.
	*/
	int i,j,k,nd; vector<int> adj(0);
	vector<int> out_no = sort_edges(out_ed,edgelist), 
	rbc_no = sort_edges(rbc_ed,edgelist);
	
	for(i=0; i<out_no.size(); i++){
		//This holds the number of connected edges having of the correct type
		adj.clear();
		nd = out_no[i];//nd is the node number
		for(j=0; j<(*EAL)[nd].size(); j++){
			//k is the edge number
			k = (*EAL)[nd][j];
			if( (*edgelist)[k].edbtype == 'o' ){
				adj.push_back(k);
			}
		}
		if(adj.size()==2){
			if( 1-abs(sin( edge_angle((*edgelist)[adj[0]],nodelist)-edge_angle((*edgelist)[adj[1]],nodelist) )) < 1e-3 ){
				//This means we have found a corner
				if( (*nodelist)[nd].ny > 0 && (*nodelist)[nd].nx > 0 ){
					(*tr_va) = nd;
				}
				else if( (*nodelist)[nd].ny < 0 && (*nodelist)[nd].nx < 0){
					(*bl_va) = nd;
				}
			}
		}
		else{
			cout<<"** Error in getting mesh extents in out_no list **"<<endl;
		}
	}
	
	
	for(i=0; i<rbc_no.size(); i++){
		adj.clear();
		nd = rbc_no[i];//nd is the node number
		//investigat the EAL to find connected nodes whose edges are 'r'
		for(j=0; j<(*EAL)[nd].size(); j++){
			//k is the edge number
			k = (*EAL)[nd][j];
			if( (*edgelist)[k].edbtype == 'r' ){
				adj.push_back(k);
			}
		}
		if(adj.size() == 2){
			if( 1-abs(sin( edge_angle((*edgelist)[adj[0]],nodelist)-edge_angle((*edgelist)[adj[1]],nodelist) )) < 1e-3 ){
				//This means we have found a corner
				if( (*nodelist)[nd].ny > 0 && (*nodelist)[nd].nx > 0 ){
					(*tr_aa) = nd;
				}
				else if( (*nodelist)[nd].ny < 0 && (*nodelist)[nd].nx < 0){
					(*bl_aa) = nd;
				}
			}
		}
		else{
			cout<<"** Error in getting mesh extents in rbc_no list **"<<endl;
		}
	}
	
	//Output section
	if( abs((*nodelist)[(*tr_aa)].ny - abs((*nodelist)[(*bl_aa)].ny)) > lambda/20.0 ){
		cout<<"** Error: Mesh (AA) is not symmetric about x-axis **"<<endl;
	}
	if( abs((*nodelist)[(*tr_aa)].nx - abs((*nodelist)[(*bl_aa)].nx)) > lambda/20.0 ){
		cout<<"** Error: Mesh (AA) is not symmetric about y-axis **"<<endl;
	}
	if( abs((*nodelist)[(*tr_va)].ny - abs((*nodelist)[(*bl_va)].ny)) > lambda/20.0 ){
		cout<<"** Error: Mesh (VA) is not symmetric about x-axis **"<<endl;
	}
	if( abs((*nodelist)[(*tr_va)].nx - abs((*nodelist)[(*bl_va)].nx)) > lambda/20.0 ){
		cout<<"** Error: Mesh (VA) is not symmetric about y-axis **"<<endl;
	}
	
	cout<<setprecision(3)<<"Mesh(air_x,air_y)=("<<((*nodelist)[(*tr_va)].nx-(*nodelist)[(*bl_va)].nx)/lambda<<","
	<<((*nodelist)[(*tr_va)].ny-(*nodelist)[(*bl_va)].ny)/lambda<<"),"
	<<"AdAb(dx,dy)=("<<((*nodelist)[(*tr_aa)].nx-(*nodelist)[(*tr_va)].nx)/lambda<<","
	<<((*nodelist)[(*tr_aa)].ny-(*nodelist)[(*tr_va)].ny)/lambda<<")"<<endl;
}

//This function get and sets the extents of the mesh region with fem.v23
//assumes circular geometry 
void get_mesh_extents(vector<int> * out_ed,vector<int> *rbc_ed, double lambda,vector<node> * nodelist, vector<edge> * edgelist){

	double r_rbc, r_out;
	r_rbc = (*edgelist)[(*rbc_ed)[0]].midpoint(nodelist).node_dist(node(0,0));
	r_out = (*edgelist)[(*out_ed)[0]].midpoint(nodelist).node_dist(node(0,0));
	cout<<"Mesh extents: (rbc,out) "<<r_rbc/lambda<<", "<<r_out/lambda<<endl;
}

void get_mesh_extents(vector<int> *diel_ed, vector<int> *top_ed, vector<int> *rbc_ed, double lambda,vector<node> * nodelist, vector<edge> * edgelist){
	
	double r_rbc, r_top, tr_x, tr_y, bl_x, bl_y; node n;
	r_rbc = (*edgelist)[(*rbc_ed)[0]].midpoint(nodelist).node_dist(node(0,0));
	r_top = (*edgelist)[(*top_ed)[0]].midpoint(nodelist).node_dist(node(0,0));
	
	//Run through all the edges on diel
	for(int i=0; i<(*diel_ed).size(); i++){
		//choose each of the two nodes of the edge
		for(int j=0; j<2; j++){
			
			if(j==0){ n = (*nodelist)[(*edgelist)[(*diel_ed)[i]].ni];}
			else{ n = (*nodelist)[(*edgelist)[(*diel_ed)[i]].nj];}
				
			tr_x = Max(n.nx,tr_x); tr_y = Max(n.ny,tr_y);
			bl_x = Min(n.nx,bl_x); bl_y = Min(n.ny,bl_y);
		}
	}
	//Print on the screen in lambda-independent numbers
	cout<<"Mesh extents:(r_rbc,r_meas),(tr_x,tr_y),(bl_x,bl_y)=("
	<<r_rbc/lambda<<","<<r_top/lambda<<"),("<<tr_x/lambda<<","<<tr_y/lambda<<"),("<<bl_x/lambda<<","<<bl_y/lambda<<")"<<endl;
	
	//Hack for the inverse solver -- store them to file as well, but absolute numbers, not relative to wavelength
	ofstream extents("extents.csv");
	extents<<"r_meas,tr_x,tr_y,bl_x,bl_y"<<endl;
	extents<<r_top<<","<<tr_x<<","<<tr_y<<","<<bl_x<<","<<bl_y<<endl;
	extents.close();
}
		
		

//This function builds the LHS entries of the FEM matrix from element k
void set_element_LHS(int k, double lambda, Matrix<cmp,Symmetric,ArrayRowSymSparse> * A0tmp, 
Matrix<cmp,Symmetric,ArrayRowSymSparse> * A1tmp, vector<node> * nodelist, 
vector<edge> * edgelist, vector<element> * elementlist){	
	
	element el = (*elementlist)[k]; node pm = el.elcenter(nodelist),pi,pj;
	int i,j,m,n,tsign; cmp er = el.eleps, mr = el.elmu, zr = sqrt(mr/er);
	double ar = el.elarea,Ai,Bi,Ci,Li,Aj,Bj,Cj,Lj,K0=2.0*PI/lambda,I1,I2,I3,I4,I5,del,dot;
	vector<double> integral = integrate_lhs( &el, nodelist);
	
	for(i=0; i<3; i++){
		m = el.edges[i];
		Ai = el.AA[i]; Bi = el.BB[i]; Ci = el.CC[i]; Li = (*edgelist)[m].edlength;
		for(j=0; j<3; j++){
			n = el.edges[j];
			if( m>n ){continue;} //this is because of the symmetric nature only A(m,n) with m>=n is stored
			Aj = el.AA[j]; Bj = el.BB[j]; Cj = el.CC[j]; Lj = (*edgelist)[n].edlength;
			
			I1 = (Ai*Aj + Ci*Cj) * ar; //const
			I2 = -1.0 * (Ci*Bj + Cj*Bi) * integral[0]; //x
			I3 = (Ai*Bj + Aj*Bi) * integral[1]; //y
			I4 = (Bi*Bj) * integral[2]; //y^2
			I5 = (Bi*Bj) * integral[3]; //x^2 since Bi = - Di
			
			del = Li*Lj*Bi*Bj/(4.0*ar*ar*ar); dot = (I1+I2+I3+I4+I5) * Li*Lj/(16.0*ar*ar*ar*ar); 
			
			(*A0tmp)(m-1,n-1) += (del/er - K0*K0 * mr * dot);
			(*A1tmp)(m-1,n-1) += (del/mr - K0*K0 * er * dot);
			
			//RBC edge
			if( i==j && (*edgelist)[m].edbtype == 'r' && (*edgelist)[n].edbtype == 'r' ){
				(*A0tmp)(m-1,n-1) += I * K0 * Li * zr;
				(*A1tmp)(m-1,n-1) += I * K0 * Li / zr;
			}
			//PEC edge
			if( (*edgelist)[m].edbtype == 'p' || (*edgelist)[n].edbtype == 'p' ){
				if( m==n ) {
					(*A1tmp)(m-1,n-1) = 1.0;
				}
				else {
					(*A1tmp)(m-1,n-1) = 0.0;
				}
			}
		}//end j loop
	}//end i loop
	
}	

//Returns true if an edge belongs to the interior of the itf region
bool is_itf_int_edge(int edno, vector<node> * nlist,vector<edge> * edlist,vector<element> * elist, double ymin, double ymax){
	node p = (*edlist)[edno].midpoint(nlist);
	if( p.ny > ymax || p.ny < ymin){
	//	cout<<"y";	
		return false;
	}
	for(int i=0; i<(*edlist)[edno].els.size(); i++){
		if( (*elist)[(*edlist)[edno].els[i]].eltype == 'a' ){
	//		cout<<"x";
			return false;
		}
	}
	return true;
}

//This function sets LHS and RHS entries in the scattered field formulation for element k
void set_matrix_entries_sf(int k, double lambda,vector<double> thetas,double mn,double soilw, 
Matrix<cmp,Symmetric,ArrayRowSymSparse> * A0tmp, Matrix<cmp,Symmetric,ArrayRowSymSparse> * A1tmp, 
vector< Vector<cmp> > * BTM, vector< Vector<cmp> > * BTE,
vector<node> * nodelist, vector<edge> * edgelist, vector<element> * elementlist){
	
	element el = (*elementlist)[k]; node pi,pj,po(0,mn-10.0*lambda);
	int i,j,m,n,tsign,f; cmp er = el.eleps, mr = el.elmu, zr = sqrt(mr/er);
	double ar = el.elarea,Ai,Bi,Ci,Li,Aj,Bj,Cj,Lj,K0=2.0*PI/lambda,
	I1,I2,I3,I4,I5,del,dot,phi,phj;
	vector<double> integral = integrate_lhs( &el, nodelist);
	
	//Idea is to use the scattered field as the variable on the interface
	//Contour intergral contribution of interface is taken from vacuum element
	
	for(i=0; i<3; i++){
		m = el.edges[i];
		//Testing against edge m in element k
		Ai = el.AA[i]; Bi = el.BB[i]; Ci = el.CC[i]; Li = (*edgelist)[m].edlength;
		
		//Add RHS contour integral contri (from vacuum element) on interface edge
		//This is the first place where incidence field enters
		if( el.eltype == 'v' && (*edgelist)[m].edbtype == 'd' ){
			pi = (*nodelist)[(*edgelist)[m].ni]; pj = (*nodelist)[(*edgelist)[m].nj];
		/*	phi = atan2(pi.ny-po.ny,pi.nx-po.nx); 
			phj = atan2(pj.ny-po.ny,pj.nx-po.nx);
			//these angles are w.r.t po which is 10lambda below interface
			//edge points from pi to pj, tsign has the sign of dot product
			tsign = sign(phj-phi);
			if( (phi>PI/2.0 && phj<PI/(-2.0))||(phj>PI/2.0 && phi<PI/(-2.0)) ){
				tsign *= -1;
			} // due to discontinuity at +- PI */
			//instead, cos is +ve in +x quadrant for which the dot pr
			//is negative
			tsign = -1 * sign( cos(atan2(pj.ny-pi.ny,pj.nx-pi.nx)) );
			for(f=0; f<thetas.size(); f++){
				double g = get_thorosG(thetas[f],soilw);
				cmp Ui = thoros(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),g,thetas[f],mn,lambda);
				(*BTM)[f](m-1) += I*K0/Z0 * Ui * double(tsign) * Li;
				(*BTE)[f](m-1) -= I*K0*Z0 * Ui * double(tsign) * Li;
		//		cout<<" RHS1 "<<(*BTM)[f](m-1)<<endl;
			}
		}
		
		//Build LHS
		for(j=0; j<3; j++){
			n = el.edges[j];
			//Coupling to edge n
			
			Aj = el.AA[j]; Bj = el.BB[j]; Cj = el.CC[j]; Lj = (*edgelist)[n].edlength;
			
			I1 = (Ai*Aj + Ci*Cj) * ar; //const
			I2 = -1.0 * (Ci*Bj + Cj*Bi) * integral[0]; //x
			I3 = (Ai*Bj + Aj*Bi) * integral[1]; //y
			I4 = (Bi*Bj) * integral[2]; //y^2
			I5 = (Bi*Bj) * integral[3]; //x^2 since Bi = - Di
			
			del = Li*Lj*Bi*Bj/(4.0*ar*ar*ar); dot = (I1+I2+I3+I4+I5) * Li*Lj/(16.0*ar*ar*ar*ar); 
			
			//Add RHS contri from interface edge
			if( el.eltype == 's' && (*edgelist)[n].edbtype == 'd' ){
				for(f=0; f<thetas.size(); f++){
					vector<cmp> inc_comp = inc_field_along_edge((*edgelist)[n],thetas[f],soilw,lambda,mn,nodelist);
					(*BTM)[f](m-1) -= inc_comp[0] * (del/er - K0*K0 * mr * dot);
					(*BTE)[f](m-1) -= inc_comp[1] * (del/mr - K0*K0 * er * dot);
				//	cout<<" RHS2 "<<(*BTM)[f](m-1)<<endl;
				}
			}
			
			if( m>n ){continue;} 
			//this is because of the symmetric nature only A(m,n) 
			//with m<=n is stored, i.e. upper triangular
			(*A0tmp)(m-1,n-1) += (del/er - K0*K0 * mr * dot);
			(*A1tmp)(m-1,n-1) += (del/mr - K0*K0 * er * dot);
	//		cout<<" LHS "<<(*A0tmp)(m-1,n-1)<<" ";
						
			//RBC edge
			if( m==n && (*edgelist)[m].edbtype == 'r' ){
				(*A0tmp)(m-1,n-1) += I * K0 * Li * zr;
				(*A1tmp)(m-1,n-1) += I * K0 * Li / zr;
			}
			//PEC edge
			if( (*edgelist)[m].edbtype == 'p' || (*edgelist)[n].edbtype == 'p' ){
				if( m==n ) {
					(*A1tmp)(m-1,n-1) = 1.0;
				}
				else {
					(*A1tmp)(m-1,n-1) = 0.0;
				}
			}
		}//end j loop
	}//end i loop	
//	cout<<endl;
}
/*
//This function sets LHS and RHS entries in the total field formulation for element k
//soil formulation
void set_matrix_entries_tf(int k, double lambda,vector<double> thetas,double mn,double soilw, 
Matrix<cmp,Symmetric,ArrayRowSymSparse> * A0tmp, Matrix<cmp,Symmetric,ArrayRowSymSparse> * A1tmp, 
vector< Vector<cmp> > * BTM, vector< Vector<cmp> > * BTE,
vector<node> * nodelist, vector<edge> * edgelist, vector<element> * elementlist){
	
	element el = (*elementlist)[k]; node pi,pj,po(0,mn-10.0*lambda);
	int i,j,m,n,tsign,f; cmp er = el.eleps, mr = el.elmu, zr = sqrt(mr/er);
	double ar = el.elarea,Ai,Bi,Ci,Li,Aj,Bj,Cj,Lj,K0=2.0*PI/lambda,
	I1,I2,I3,I4,I5,del,dot,phi,phj,nphi;
	vector<double> integral = integrate_lhs( &el, nodelist);
	
	for(i=0; i<3; i++){
		m = el.edges[i];
		//Testing against edge m in element k
		Ai = el.AA[i]; Bi = el.BB[i]; Ci = el.CC[i]; Li = (*edgelist)[m].edlength;
		
		//Add RHS contour integral contri from RBC edges
		//This is the place where incidence field enters
		if( (*edgelist)[m].edbtype == 'r' ){
			pi = (*nodelist)[(*edgelist)[m].ni]; pj = (*nodelist)[(*edgelist)[m].nj];
			phi = atan2(pi.ny,pi.nx); 
			phj = atan2(pj.ny,pj.nx);
			//these angles are w.r.t po which is 10lambda below interface
			//edge points from pi to pj, tsign has the sign of dot product
			tsign = sign(phj-phi);
			if( (phi>PI/2.0 && phj<PI/(-2.0))||(phj>PI/2.0 && phi<PI/(-2.0)) ){
				tsign *= -1;
			} // due to discontinuity at +- PI 
			//instead, cos is +ve in +x quadrant for which the dot pr
			//is negative
			//tsign = -1 * sign( cos(atan2(pj.ny-pi.ny,pj.nx-pi.nx)) );
			nphi = atan2(tsign*(pj.ny-pi.ny),tsign*(pj.nx-pi.nx)) - PI/2.0;

			for(f=0; f<thetas.size(); f++){
				double G = get_thorosG(thetas[f],soilw);
				cmp taper = thoros(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),G,thetas[f],mn,lambda),
				taper_dx = thoros_dx(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),G,thetas[f],mn,lambda),
				taper_dy = thoros_dy(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),G,thetas[f],mn,lambda);
		
				//the incident field is represented by taper. assumption here is that contour integral
				//is equal to length * value at the center.
				(*BTM)[f](m-1) = (I*K0*taper + (cos(nphi)*taper_dx + sin(nphi)*taper_dy)) * double(tsign)* Li/Z0;
				(*BTE)[f](m-1) = (I*K0*taper + (cos(nphi)*taper_dx + sin(nphi)*taper_dy)) * double(tsign)* Li*Z0*-1.0;
			}
		}
		
		//Build LHS
		for(j=0; j<3; j++){
			n = el.edges[j];
			//Coupling to edge n
			
			Aj = el.AA[j]; Bj = el.BB[j]; Cj = el.CC[j]; Lj = (*edgelist)[n].edlength;
			
			I1 = (Ai*Aj + Ci*Cj) * ar; //const
			I2 = -1.0 * (Ci*Bj + Cj*Bi) * integral[0]; //x
			I3 = (Ai*Bj + Aj*Bi) * integral[1]; //y
			I4 = (Bi*Bj) * integral[2]; //y^2
			I5 = (Bi*Bj) * integral[3]; //x^2 since Bi = - Di
			
			del = Li*Lj*Bi*Bj/(4.0*ar*ar*ar); dot = (I1+I2+I3+I4+I5) * Li*Lj/(16.0*ar*ar*ar*ar); 
			
			if( m>n ){continue;} 
			//this is because of the symmetric nature only A(m,n) 
			//with m<=n is stored, i.e. upper triangular
			(*A0tmp)(m-1,n-1) += (del/er - K0*K0 * mr * dot);
			(*A1tmp)(m-1,n-1) += (del/mr - K0*K0 * er * dot);
	//		cout<<" LHS "<<(*A0tmp)(m-1,n-1)<<" ";
						
			//RBC edge
			if( m==n && (*edgelist)[m].edbtype == 'r' ){
				(*A0tmp)(m-1,n-1) += I * K0 * Li * zr;
				(*A1tmp)(m-1,n-1) += I * K0 * Li / zr;
			}
			//PEC edge
			if( (*edgelist)[m].edbtype == 'p' || (*edgelist)[n].edbtype == 'p' ){
				if( m==n ) {
					(*A1tmp)(m-1,n-1) = 1.0;
				}
				else {
					(*A1tmp)(m-1,n-1) = 0.0;
				}
			}
		}//end j loop
	}//end i loop	
//	cout<<endl;
} */

//This function sets LHS and RHS entries in the total field formulation for element k
//general formulation
void set_matrix_entries_tf(int k, double lambda,vector<double> thetas, 
Matrix<cmp,Symmetric,ArrayRowSymSparse> * A0tmp, Matrix<cmp,Symmetric,ArrayRowSymSparse> * A1tmp, 
vector< Vector<cmp> > * BTM, vector< Vector<cmp> > * BTE,
vector<node> * nodelist, vector<edge> * edgelist, vector<element> * elementlist){
	
	element el = (*elementlist)[k]; node pi,pj;
	int i,j,m,n,tsign,f; cmp er = el.eleps, mr = el.elmu, zr = sqrt(mr/er);
	double ar = el.elarea,Ai,Bi,Ci,Li,Aj,Bj,Cj,Lj,K0=2.0*PI/lambda,
	I1,I2,I3,I4,I5,del,dot,phi,phj,nphi;
	vector<double> integral = integrate_lhs( &el, nodelist);
	
	for(i=0; i<3; i++){
		m = el.edges[i];
		//Testing against edge m in element k
		Ai = el.AA[i]; Bi = el.BB[i]; Ci = el.CC[i]; Li = (*edgelist)[m].edlength;
		
		//Add RHS contour integral contri from RBC edges
		//This is the place where incidence field enters
		if( (*edgelist)[m].edbtype == 'r' ){
			pi = (*nodelist)[(*edgelist)[m].ni]; pj = (*nodelist)[(*edgelist)[m].nj];
			phi = atan2(pi.ny,pi.nx); 
			phj = atan2(pj.ny,pj.nx);
			//edge points from pi to pj, tsign has the sign of dot product
			tsign = sign(phj-phi);
			if( (phi>PI/2.0 && phj<PI/(-2.0))||(phj>PI/2.0 && phi<PI/(-2.0)) ){
				tsign *= -1;
			} // due to discontinuity at +- PI 
			//instead, cos is +ve in +x quadrant for which the dot pr
			//is negative
			//tsign = -1 * sign( cos(atan2(pj.ny-pi.ny,pj.nx-pi.nx)) );
			nphi = atan2(tsign*(pj.ny-pi.ny),tsign*(pj.nx-pi.nx)) - PI/2.0;

			for(f=0; f<thetas.size(); f++){
				cmp taper = thoros(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),thetas[f],lambda),
				taper_dx = thoros_dx(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),thetas[f],lambda),
				taper_dy = thoros_dy(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),thetas[f],lambda);
		
				//the incident field is represented by taper. assumption here is that contour integral
				//is equal to length * value at the center.
				(*BTM)[f](m-1) = (I*K0*taper + (cos(nphi)*taper_dx + sin(nphi)*taper_dy)) * double(tsign)* Li/Z0;
				(*BTE)[f](m-1) = (I*K0*taper + (cos(nphi)*taper_dx + sin(nphi)*taper_dy)) * double(tsign)* Li*Z0*-1.0;
			}
		}
		
		//Build LHS
		for(j=0; j<3; j++){
			n = el.edges[j];
			//Coupling to edge n
			
			Aj = el.AA[j]; Bj = el.BB[j]; Cj = el.CC[j]; Lj = (*edgelist)[n].edlength;
			
			I1 = (Ai*Aj + Ci*Cj) * ar; //const
			I2 = -1.0 * (Ci*Bj + Cj*Bi) * integral[0]; //x
			I3 = (Ai*Bj + Aj*Bi) * integral[1]; //y
			I4 = (Bi*Bj) * integral[2]; //y^2
			I5 = (Bi*Bj) * integral[3]; //x^2 since Bi = - Di
			
			del = Li*Lj*Bi*Bj/(4.0*ar*ar*ar); dot = (I1+I2+I3+I4+I5) * Li*Lj/(16.0*ar*ar*ar*ar); 
			
			if( m>n ){continue;} 
			//this is because of the symmetric nature only A(m,n) 
			//with m<=n is stored, i.e. upper triangular
			(*A0tmp)(m-1,n-1) += (del/er - K0*K0 * mr * dot);
			(*A1tmp)(m-1,n-1) += (del/mr - K0*K0 * er * dot);
	//		cout<<" LHS "<<(*A0tmp)(m-1,n-1)<<" ";
						
			//RBC edge
			if( m==n && (*edgelist)[m].edbtype == 'r' ){
				(*A0tmp)(m-1,n-1) += I * K0 * Li * zr;
				(*A1tmp)(m-1,n-1) += I * K0 * Li / zr;
			}
			//PEC edge
			if( (*edgelist)[m].edbtype == 'p' || (*edgelist)[n].edbtype == 'p' ){
				if( m==n ) {
					(*A1tmp)(m-1,n-1) = 1.0;
				}
				else {
					(*A1tmp)(m-1,n-1) = 0.0;
				}
			}
		}//end j loop
	}//end i loop	
//	cout<<endl;
}

//This function sets LHS and RHS entries in the scattered field formulation for element k
void set_matrix_entries(int k, double lambda,vector<double> thetas,
Matrix<cmp,Symmetric,ArrayRowSymSparse> * A0tmp, Matrix<cmp,Symmetric,ArrayRowSymSparse> * A1tmp, 
vector< Vector<cmp> > * BTM, vector< Vector<cmp> > * BTE,
vector<node> * nodelist, vector<edge> * edgelist, vector<element> * elementlist){
	
	element el = (*elementlist)[k]; node pi,pj;
	int i,j,m,n,tsign,f; cmp er = el.eleps, mr = el.elmu, zr = sqrt(mr/er);
	double ar = el.elarea,Ai,Bi,Ci,Li,Aj,Bj,Cj,Lj,K0=2.0*PI/lambda,
	I1,I2,I3,I4,I5,del,dot,phi,phj;
	vector<double> integral = integrate_lhs( &el, nodelist);
	
	//Idea is to use the scattered field as the variable on the interface
	//Contour intergral contribution of interface is taken from vacuum element
	
	for(i=0; i<3; i++){
		m = el.edges[i];
		//Testing against edge m in element k
		Ai = el.AA[i]; Bi = el.BB[i]; Ci = el.CC[i]; Li = (*edgelist)[m].edlength;
		
		//Add RHS contour integral contri (from vacuum element) on interface edge
		//This is the first place where incidence field enters
		if( el.eltype == 'v' && (*edgelist)[m].edbtype == 'd' ){
			pi = (*nodelist)[(*edgelist)[m].ni]; pj = (*nodelist)[(*edgelist)[m].nj];
			phi = atan2(pi.ny,pi.nx); 
			phj = atan2(pj.ny,pj.nx);
			//edge points from pi to pj, tsign has the sign of dot product
			tsign = sign(phj-phi);
			if( (phi>PI/2.0 && phj<PI/(-2.0))||(phj>PI/2.0 && phi<PI/(-2.0)) ){
				tsign *= -1;
			} // due to discontinuity at +- PI 
			//instead, cos is +ve in +x quadrant for which the dot pr
			//is negative
			//tsign = -1 * sign( cos(atan2(pj.ny-pi.ny,pj.nx-pi.nx)) );
			for(f=0; f<thetas.size(); f++){
				cmp Ui = thoros(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),thetas[f],lambda);
				(*BTM)[f](m-1) += I*K0/Z0 * Ui * double(tsign) * Li;
				(*BTE)[f](m-1) -= I*K0*Z0 * Ui * double(tsign) * Li;
			}
		}
		
		//Build LHS
		for(j=0; j<3; j++){
			n = el.edges[j];
			//Coupling to edge n
			
			Aj = el.AA[j]; Bj = el.BB[j]; Cj = el.CC[j]; Lj = (*edgelist)[n].edlength;
			
			I1 = (Ai*Aj + Ci*Cj) * ar; //const
			I2 = -1.0 * (Ci*Bj + Cj*Bi) * integral[0]; //x
			I3 = (Ai*Bj + Aj*Bi) * integral[1]; //y
			I4 = (Bi*Bj) * integral[2]; //y^2
			I5 = (Bi*Bj) * integral[3]; //x^2 since Bi = - Di
			
			del = Li*Lj*Bi*Bj/(4.0*ar*ar*ar); dot = (I1+I2+I3+I4+I5) * Li*Lj/(16.0*ar*ar*ar*ar); 
			
			//Add RHS contri from interface edge
			if( el.eltype == 's' && (*edgelist)[n].edbtype == 'd' ){
				for(f=0; f<thetas.size(); f++){
					vector<cmp> inc_comp = inc_field_along_edge((*edgelist)[n],thetas[f],lambda,nodelist);
					(*BTM)[f](m-1) -= inc_comp[0] * (del/er - K0*K0 * mr * dot);
					(*BTE)[f](m-1) -= inc_comp[1] * (del/mr - K0*K0 * er * dot);
				}
			}
			
			if( m>n ){continue;} 
			//this is because of the symmetric nature only A(m,n) 
			//with m<=n is stored, i.e. upper triangular
			(*A0tmp)(m-1,n-1) += (del/er - K0*K0 * mr * dot);
			(*A1tmp)(m-1,n-1) += (del/mr - K0*K0 * er * dot);
						
			//RBC edge
			if( m==n && (*edgelist)[m].edbtype == 'r' ){
				(*A0tmp)(m-1,n-1) += I * K0 * Li * zr;
				(*A1tmp)(m-1,n-1) += I * K0 * Li / zr;
			}
			//PEC edge
			if( (*edgelist)[m].edbtype == 'p' || (*edgelist)[n].edbtype == 'p' ){
				if( m==n ) {
					(*A1tmp)(m-1,n-1) = 1.0;
				}
				else {
					(*A1tmp)(m-1,n-1) = 0.0;
				}
			}

		}//end j loop
	}//end i loop	
}

void set_element_RHS(int ei,double lambda, double theta, double y0, double soilw,Vector<cmp> * BTM, Vector<cmp> * BTE, 
	vector<node> * nlist, vector<edge> * edlist, vector<element> *elist){
	//This function is for the total field formulation
	//This is the only place where the incident field appears. 
	//Convention for incident field is to zero the phase w.r.t y at the soil surface.
	int m,n,k,tsign; vector<cmp> result(2,0.0);
	double Li=(*edlist)[ei].edlength,phi,phj,nphi,K0 = 2.0*PI/lambda;
	cmp er,mr,zr,nz;
	
	node pi = (*nlist)[(*edlist)[ei].ni]; node pj = (*nlist)[(*edlist)[ei].nj];
	phi = atan2(pi.ny,pi.nx); phj = atan2(pj.ny,pj.nx);
	//edge points from pi to pj, tsign has the sign of dot product
	tsign = sign(phj-phi);
	if( (phi>PI/2.0 && phj<PI/(-2.0))||(phj>PI/2.0 && phi<PI/(-2.0)) )
	{tsign *= -1;} // due to discontinuity at +- PI
	
	nphi = atan2(tsign*(pj.ny-pi.ny),tsign*(pj.nx-pi.nx)) - PI/2.0;
//	tphi = atan2(tsign*(pj.ny-pi.ny),tsign*(pj.nx-pi.nx));
	
	n = (*edlist)[ei].els[0]; // on the boundary, so only one element
	er = (*elist)[n].eleps; mr = (*elist)[n].elmu; zr = sqrt(mr/er); nz = sqrt(er*mr);
	
	//Make a distinction between above and below ground
 	if( (*elist)[n].eltype == 'v'){ //only in air	
		double G = get_thorosG(theta,soilw);
		cmp taper = thoros(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),G,theta,y0,lambda),
		taper_dx = thoros_dx(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),G,theta,y0,lambda),
		taper_dy = thoros_dy(0.5*(pi.nx+pj.nx),0.5*(pi.ny+pj.ny),G,theta,y0,lambda);

		//the incident field is represented by taper. assumption here is that contour integral
		//is equal to length * value at the center.
		result[0] = (I*K0*taper + (cos(nphi)*taper_dx + sin(nphi)*taper_dy)) * double(tsign)* Li/Z0;
		result[1] = (I*K0*taper + (cos(nphi)*taper_dx + sin(nphi)*taper_dy)) * double(tsign)* Li*Z0*-1.0;
	}
	else { //incident field is dead by here
		result[0] = 0.0; // because the approximation is that there is no
		result[1] = 0.0; // more incident field left below the soil surface. needs to be refined on the sides
	}
	
	if( ((*edlist)[ei].edbtype == 'p')){
		result[0] = 0.0; // because curl H = E = 0
		result[1] = 0.0; // because curl E = 0
	}
	
	//the -1 is needed below because of numbering convention with MUMPS matrices 
	(*BTM)(ei-1) = result[0]; 
	(*BTE)(ei-1) = result[1]; 
}

//Print out mesh nodes for a given element list
void print_mesh_nodes(vector<int> * list, vector<node> * nlist, vector<element> * elist){
	//Get a non-duplicate list of nodes in these elements
	vector<int> pnodes = sort_elements(list,elist);
	//Now print them out
	ofstream mesho("mesh.csv");
	for(int i=0; i<pnodes.size(); i++){
		node p = (*nlist)[pnodes[i]];
		mesho<<p.nx<<","<<p.ny<<endl;
	}
	mesho.close();
}

//Print out mesh contours		 
void print_mesh_contours(vector<node> * nodelist, vector<edge> * edgelist, vector<element> * elementlist){
	ofstream mesho("contours.csv");

	for(int i=1; i<(*edgelist).size(); i++){
		node pi = (*edgelist)[i].midpoint(nodelist);
		if( (*edgelist)[i].els.size() == 2){
			if( (*elementlist)[(*edgelist)[i].els[0]].eltype != (*elementlist)[(*edgelist)[i].els[1]].eltype ){
				mesho<<pi.nx<<","<<pi.ny<<endl;
			}
		}
		else{
			mesho<<pi.nx<<","<<pi.ny<<endl;
		}
		if( (*edgelist)[i].edbtype == 't'){
			mesho<<pi.nx<<","<<pi.ny<<endl;
		}
	}
	mesho.close();
}

void print_mesh_stats(vector<element> * elist, vector<edge> * edlist, double lambda){
	double min_a=lambda, min_s=lambda, min_v=lambda, max_a=0, max_s=0, max_v=0,l;
	for(int i=1; i< (*edlist).size(); i++){
		l = (*edlist)[i].edlength;
		if( (*elist)[(*edlist)[i].els[0]].eltype == 'a'){
			min_a = Min(min_a,l); max_a = Max(max_a,l);
		}
		else if( (*elist)[(*edlist)[i].els[0]].eltype == 'v'){
			min_v = Min(min_v,l); max_v = Max(max_v,l);
		}
		else if( (*elist)[(*edlist)[i].els[0]].eltype == 's'){
			min_s = Min(min_s,l); max_s = Max(max_s,l);
		}
	}
	cout<<setprecision (3)<<"Min/max lengths/lambda, v=("
	<<min_v/lambda<<","<<max_v/lambda<<"), s="
	<<min_s/lambda<<","<<max_s/lambda<<"), a="
	<<min_a/lambda<<","<<max_a/lambda<<")."<<endl;
}

int set_tissue_eps(vector<cmp> * tissues, vector<vector< int> > * tiss_els, vector<element> * elist){
	
	cout<<"Tissue eps size "<<(*tissues).size()<<", elset number "<<(*tiss_els).size()<<endl;
	if( (*tissues).size() !=  (*tiss_els).size() ){
		cout<<"*** Mismatch between tissue dielectric specfications and elementset numbers ***"<<endl;
		return -1;
	}
	else{
		cout<<"Assigning "<<(*tiss_els).size()<<" tissues"<<endl;
	}
	
	for(int i=0; i< (*tiss_els).size(); i++){
		for(int j=0; j<(*tiss_els)[i].size(); j++){
			(*elist)[(*tiss_els)[i][j]].eleps = (*tissues)[i];
		}
	}
	return 0;
}
