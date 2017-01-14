/**************************************************************************************/
/* Fitzhugh-Nagumo reaction diffusion simulation with arbitrary vortex lines
   OPENMP VERSION
   Created by Carl Whitfield
   Last modified 03/01/17

   Operational order of the code:
   1) The code takes as input an stl file (defined in knot_filename) which defines an orientable surface with a boundary.
   2) This surface is scaled to fill a box of size xmax x ymax x zmax.
   3) A nunmerical integral is performed to calculate a phase field (phi_calc) on the 3D grid which winds around the boundary of the surface.
   4) This is then used to initialise the Fitzhugh-Nagumo set of partial differential equations such that on initialisation (uv_initialise):
   u = 2cos(phi) - 0.4        v = sin(phi) - 0.4
   The pde's used are
   dudt = (u - u^3/3 - v)/epsilon + Del^2 u
   dvdt = epsilon*(u + beta - gam v)
   5) The update method is Runge-Kutta fourth order (update_uv) unless RK4 is set to 0, otherwise Euler forward method is used.
   6) A parametric curve for the knot is found at each unit T


   The parameters epsilon, beta and gam are set to give rise to scroll waves (see Sutcliffe, Winfree, PRE 2003 and Maucher Sutcliffe PRL 2016) which eminate from a closed curve, on initialisation this curve corresponds to the boundary of the surface in the stl file.

   See below for various options to start the code from previous outputted data.*/
/**************************************************************************************/
#include "FN_knot_surface.h"    //contains some functions and all global variables
#include <omp.h>
#define RK4 1         //1 to use Runge Kutta 4th order method, 0 for euler forward method
//includes for the signal processing
#include <gsl/gsl_min.h>
#include <gsl/gsl_vector.h>
/*Available options:
FROM_PHI_FILE: Skip initialisation, input from previous run.
FROM_SURFACE_FILE: Initialise from input file(s) generated in surface evolver.
FROM_UV_FILE: Skip initialisation, run FN dynamics from uv file
FROM_KNOT_FILE: Initialise from parametric knot curve in .txt format (e.g. knotplot output)
 */

int option = FROM_UV_FILE;         //unknot default option
const bool periodic = false;                     //enable periodic boundaries in z

/**If FROM_SURFACE_FILE or FROM_KNOT_FILE chosen**/
string knot_filename = "zero1";      //if FROM_SURFACE_FILE assumed input filename format of "XXXXX.stl"
int ncomp = 1;                       //if FROM_KNOT_FILE assumed input filename format of "XXXXX.txt"
//if ncomp > 1 (no. of components) then component files should be separated to 'XXXXX.txt" "XXXXX2.txt", ....
/**IF FROM_PHI_FILE or FROM_UV_FILE chosen**/
string B_filename = "uv_plot20.vtk";    //filename for phi field or uv field

//Grid points
const int Nx = 241;   //No. points in x,y and z
const int Ny = 241;
const int Nz = 241;
const double TTime = 2000;         //total time of simulation (simulation units)
const double skiptime = 10;       //print out every # unit of time (simulation units)
const double starttime = 20;        //Time at start of simulation (non-zero if continuing from UV file)
const double dtime = 0.04;         //size of each time step

//System size parameters
const double lambda = 21.3;                //approx wavelength
const double size = 5*lambda;           //box size
const double h = size/(Nx-1);            //grid spacing
const double oneoverhsq = 1.0/(h*h);
const double epsilon = 0.3;                //parameters for F-N eqns
const double oneoverepsilon = 1.0/epsilon;
const double beta = 0.7;
const double gam = 0.5;

//Size boundaries of knot (now autoscaled)
double xmax = Nx*h/2.0;
double ymax = Ny*h/2.0;
double zmax = Nz*h/2.0;
int NK;   //number of surface points

//Unallocated matrices
vector<triangle> knotsurface;    //structure for storing knot surface coordinates
vector<knotpoint> knotcurve;
vector<knotpoint> knotcurveold; // a structure containing some number of knot curves, each curve a list of knotpoints
// vector<vector<knotpoint>> knotcurves; // a structure containing some number of knot curves, each curve a list of knotpoints
bool knotexists = true; //Always true at start
vector<double> X, Y, Z, dlx, dly,dlz;

double area;   //initial knot area
inline  int pt( int i,  int j,  int k)       //convert i,j,k to single index
{
	return (i*Ny*Nz+j*Nz+k);
}

int main (void)
{
	double *x, *y, *z, *phi, *u, *v, *ucvx, *ucvy, *ucvz;
	int i,j,k,n,l;
	int *missed;

	x = new double [Nx];
	y = new double [Ny];
	z = new double [Nz];
	phi = new double [Nx*Ny*Nz];  //scalar potential
	u = new double [Nx*Ny*Nz];
	v = new double [Nx*Ny*Nz];
	if(option==FROM_KNOT_FILE)
	{
		missed = new int [Nx*Ny*Nz];
		fill(missed,missed+Nx*Ny*Nz,0);
	}

	// output an info file on the run
	print_info(Nx, Ny, Nz, dtime, h, periodic, option, knot_filename, B_filename);

	// GSL initialization
	const gsl_min_fminimizer_type *Type;
	gsl_min_fminimizer *minimizerstate;
	Type = gsl_min_fminimizer_brent;
	minimizerstate = gsl_min_fminimizer_alloc (Type);
# pragma omp parallel shared ( x, y, z ) private ( i, j, k )
	{
#pragma omp for
		for(i=0;i<Nx;i++)           //initialise grid
		{
			x[i] = (i+0.5-Nx/2.0)*h;
		}
#pragma omp for
		for(j=0;j<Ny;j++)
		{
			y[j] = (j+0.5-Ny/2.0)*h;
		}
#pragma omp for
		for(k=0;k<Nz;k++)
		{
			z[k] = (k+0.5-Nz/2.0)*h;
		}
	}

	if (option == FROM_PHI_FILE)
	{
		cout << "Reading input file...\n";
		phi_file_read(phi);
	}
	else
	{
		if(option == FROM_UV_FILE)
		{
			cout << "Reading input file...\n";
			if(uvfile_read(u,v)) return 1;
		}
		else
		{
			//Initialise knot
			area = initialise_knot();
			if(area==0)
			{
				cout << "Error reading input option. Aborting...\n";
				return 1;
			}

			if(option == FROM_SURFACE_FILE) cout << "Total no. of surface points: ";
			else cout << "Total no. of knot points: ";
			cout << NK << '\n';

			//Calculate phi for initial conditions
			initial_cond(x,y,z,phi,missed);
		}

	}

	vector<triangle> ().swap(knotsurface);   //empty knotsurface memory
	vector<double> ().swap(X);   //empty initial knot curve memory
	vector<double> ().swap(Y);
	vector<double> ().swap(Z);
	vector<double> ().swap(dlx);
	vector<double> ().swap(dly);
	vector<double> ().swap(dlz);

	if(option!=FROM_UV_FILE)
	{
		cout << "Calculating u and v...\n";
		uv_initialise(phi,u,v,missed);
	}

	delete [] phi;
	if(option==FROM_KNOT_FILE) delete [] missed;

	ucvx = new double [Nx*Ny*Nz];
	ucvy = new double [Nx*Ny*Nz];
	ucvz = new double [Nx*Ny*Nz];
#if RK4
	double *ku, *kv, *kut, *kvt, *uold, *vold;
	ku = new double [Nx*Ny*Nz];
	kv = new double [Nx*Ny*Nz];
	kut = new double [Nx*Ny*Nz];
	kvt = new double [Nx*Ny*Nz];
	uold = new double [Nx*Ny*Nz];
	vold = new double [Nx*Ny*Nz];

#else
	double *D2u;

	D2u = new double [Nx*Ny*Nz];
#endif

	cout << "Updating u and v...\n";

	// initilialising counters
	int p=0;
	int q=0;
	n=0;

	// initialising timers
	time_t then = time(NULL);
	time_t rawtime;
	time (&rawtime);
	struct tm * timeinfo;
	ofstream wrout;
	wrout.open("writhe.txt");
	wrout << "Time\tWrithe\tTwist\tLength\n";
	wrout.close();

#if RK4
#pragma omp parallel default(none) shared ( x, y, z, u, v, uold, vold, n,ku,kv,kut,kvt,p,q,ucvx, ucvy, ucvz,cout, rawtime, timeinfo, knotcurve, knotcurveold, knotexists , minimizerstate)
#else
#pragma omp parallel default(none) shared ( x, y, z, u, v, n, D2u, p, q,ucvx, ucvy, ucvz, cout, rawtime, timeinfo, knotcurve, knotcurveold, knotexists ,minimizerstate)
#endif
	{
		while(n*dtime <= TTime)
		{
#pragma omp single
			{
				if(n*dtime >= q)  //Do this every unit T
				{
					if(knotexists) knotcurve.clear(); //empty vector with knot curve points
					crossgrad_calc(x,y,z,u,v,ucvx,ucvy,ucvz); //find Grad u cross Grad v
					cout << "T = " << n*dtime + starttime << endl;
					time (&rawtime);
					timeinfo = localtime (&rawtime);
					cout << "current time \t" << asctime(timeinfo) << "\n";
					if(n*dtime+starttime>=20 && knotexists)
					{
						find_knot_properties(x,y,z,ucvx,ucvy,ucvz,u,n*dtime+starttime,minimizerstate );      //find knot curve and twist and writhe
						print_knot(x,y,z,n*dtime+starttime,knotcurve);
					}
					q++;
				}

				if(n*dtime >= p*skiptime)
				{

					print_uv(x,y,z,u,v,ucvx,ucvy,ucvz,n*dtime+starttime);
					p++;
				}

				n++;
			}
#if RK4
			uv_update(u,v,ku,kv,kut,kvt,uold,vold);
#else
			uv_update_euler(u,v,D2u);
#endif
		}
	}
	time_t now = time(NULL);
	cout << "Time taken to complete uv part: " << now - then << " seconds.\n";

#if RK4
	delete [] uold;
	delete [] vold;
	delete [] ku;
	delete [] kv;
	delete [] kut;
	delete [] kvt;
#else
	delete [] D2u;
#endif
	delete [] x;
	delete [] y;
	delete [] z;
	delete [] u;
	delete [] v;
	delete [] ucvx;
	delete [] ucvy;
	delete [] ucvz;

	return 0;
}

/*************************Functions for knot initialisation*****************************/
double initialise_knot()
{
	double L;
	switch (option)
	{
		case FROM_SURFACE_FILE: L = init_from_surface_file();
					break;

		case FROM_KNOT_FILE: L = init_from_knot_file();
				     break;

		default: L=0;
			 break;
	}

	return L;
}

double init_from_surface_file(void)
{
	string filename, buff;
	stringstream ss;
	double A = 0;   //total area
	int i=0;
	int j;
	double r10,r20,r21,s,xcoord,ycoord,zcoord;
	string temp;
	ifstream knotin;
	/*  For recording max and min input values*/
	double maxxin = 0;
	double maxyin = 0;
	double maxzin = 0;
	double minxin = 0;
	double minyin = 0;
	double minzin = 0;

	ss.clear();
	ss.str("");
	ss << knot_filename << ".stl";

	filename = ss.str();
	knotin.open(filename.c_str());
	if(knotin.good())
	{
		if(getline(knotin,buff)) temp = buff;
	}
	else cout << "Error reading file\n";
	while(knotin.good())   //read in points for knot
	{
		if(getline(knotin,buff))  //read in surface normal
		{
			ss.clear();
			ss.str("");
			ss << buff;
			ss >> temp;
			if(temp.compare("endsolid") == 0) break;
			knotsurface.push_back(triangle());
			ss >> temp >> knotsurface[i].normal[0] >> knotsurface[i].normal[1] >> knotsurface[i].normal[2];
		}

		if(getline(knotin,buff)) temp = buff;   //read in "outer loop"
		knotsurface[i].centre[0] = 0;
		knotsurface[i].centre[1] = 0;
		knotsurface[i].centre[2] = 0;
		for(j=0;j<3;j++)
		{
			if(getline(knotin,buff))  //read in vertices
			{
				ss.clear();
				ss.str("");
				ss << buff;
				ss >> temp >> xcoord >> ycoord >> zcoord;

				if(xcoord>maxxin) maxxin = xcoord;
				if(ycoord>maxyin) maxyin = ycoord;
				if(zcoord>maxzin) maxzin = zcoord;
				if(xcoord<minxin) minxin = xcoord;
				if(ycoord<minyin) minyin = ycoord;
				if(zcoord<minzin) minzin = zcoord;

				knotsurface[i].xvertex[j] = xcoord;
				knotsurface[i].yvertex[j] = ycoord;
				knotsurface[i].zvertex[j] = zcoord;
				knotsurface[i].centre[0] += knotsurface[i].xvertex[j]/3.0;
				knotsurface[i].centre[1] += knotsurface[i].yvertex[j]/3.0;
				knotsurface[i].centre[2] += knotsurface[i].zvertex[j]/3.0;
			}
		}
		//cout << i << " (" << knotsurface[i].centre[0] << ',' << knotsurface[i].centre[1] << ',' << knotsurface[i].centre[2] << ") , (" << knotsurface[i].normal[0] << ',' << knotsurface[i].normal[1] << ',' << knotsurface[i].normal[2] << ") \n";

		if(getline(knotin,buff)) temp = buff;   //read in "outer loop"
		if(getline(knotin,buff)) temp = buff;   //read in "outer loop"

		i++;
	}

	NK = i;
	/* Work out space scaling for knot surface */
	double scale[3];
	if(maxxin-minxin>0) scale[0] = xmax/(maxxin-minxin);
	else scale[0] = 1;
	if(maxyin-minyin>0) scale[1] = ymax/(maxyin-minyin);
	else scale[1] = 1;
	if(maxzin-minzin>0) scale[2] = zmax/(maxzin-minzin);
	else scale[2] = 1;
	double midpoint[3];
	double norm;
	//double p1x,p1y,p1z,p2x,p2y,p2z,nx,ny,nz;
	midpoint[0] = 0.5*(maxxin+minxin);
	midpoint[1] = 0.5*(maxyin+minyin);
	midpoint[2] = 0.5*(maxzin+minzin);

	/*Rescale points and normals to fit grid properly*/
	for(i=0;i<NK;i++)
	{
		for(j=0;j<3;j++)
		{
			knotsurface[i].xvertex[j] = scale[0]*(knotsurface[i].xvertex[j] - midpoint[0]);
			knotsurface[i].yvertex[j] = scale[1]*(knotsurface[i].yvertex[j] - midpoint[1]);
			knotsurface[i].zvertex[j] = scale[2]*(knotsurface[i].zvertex[j] - midpoint[2]);
			knotsurface[i].centre[j] = scale[j]*(knotsurface[i].centre[j] - midpoint[j]);
		}

		norm = sqrt(scale[1]*scale[1]*scale[2]*scale[2]*knotsurface[i].normal[0]*knotsurface[i].normal[0] +
				scale[0]*scale[0]*scale[2]*scale[2]*knotsurface[i].normal[1]*knotsurface[i].normal[1] +
				scale[0]*scale[0]*scale[1]*scale[1]*knotsurface[i].normal[2]*knotsurface[i].normal[2]);

		knotsurface[i].normal[0] *= scale[1]*scale[2]/norm;
		knotsurface[i].normal[1] *= scale[0]*scale[2]/norm;
		knotsurface[i].normal[2] *= scale[0]*scale[1]/norm;

		/*Check surface normal is correct
		  p1x = knotsurface[i].xvertex[1] - knotsurface[i].xvertex[0];
		  p1y = knotsurface[i].yvertex[1] - knotsurface[i].yvertex[0];
		  p1z = knotsurface[i].zvertex[1] - knotsurface[i].zvertex[0];
		  p2x = knotsurface[i].xvertex[2] - knotsurface[i].xvertex[0];
		  p2y = knotsurface[i].yvertex[2] - knotsurface[i].yvertex[0];
		  p2z = knotsurface[i].zvertex[2] - knotsurface[i].zvertex[0];
		  nx = p1y*p2z - p2y*p1z;
		  ny = p1z*p2x - p2z*p1x;
		  nz = p1x*p2y - p2x*p1y;
		  norm = sqrt(nx*nx+ny*ny+nz*nz);
		  nx = nx/norm;
		  ny = ny/norm;
		  nz = nz/norm;
		  cout << nx*knotsurface[i].normal[0] + ny*knotsurface[i].normal[1] + nz*knotsurface[i].normal[2] << '\n';
		 */

		r10 = sqrt((knotsurface[i].xvertex[1]-knotsurface[i].xvertex[0])*(knotsurface[i].xvertex[1]-knotsurface[i].xvertex[0]) + (knotsurface[i].yvertex[1]-knotsurface[i].yvertex[0])*(knotsurface[i].yvertex[1]-knotsurface[i].yvertex[0]) + (knotsurface[i].zvertex[1]-knotsurface[i].zvertex[0])*(knotsurface[i].zvertex[1]-knotsurface[i].zvertex[0]));
		r20 = sqrt((knotsurface[i].xvertex[2]-knotsurface[i].xvertex[0])*(knotsurface[i].xvertex[2]-knotsurface[i].xvertex[0]) + (knotsurface[i].yvertex[2]-knotsurface[i].yvertex[0])*(knotsurface[i].yvertex[2]-knotsurface[i].yvertex[0]) + (knotsurface[i].zvertex[2]-knotsurface[i].zvertex[0])*(knotsurface[i].zvertex[2]-knotsurface[i].zvertex[0]));
		r21 = sqrt((knotsurface[i].xvertex[2]-knotsurface[i].xvertex[1])*(knotsurface[i].xvertex[2]-knotsurface[i].xvertex[1]) + (knotsurface[i].yvertex[2]-knotsurface[i].yvertex[1])*(knotsurface[i].yvertex[2]-knotsurface[i].yvertex[1]) + (knotsurface[i].zvertex[2]-knotsurface[i].zvertex[1])*(knotsurface[i].zvertex[2]-knotsurface[i].zvertex[1]));
		s = (r10+r20+r21)/2;
		knotsurface[i].area = sqrt(s*(s-r10)*(s-r20)*(s-r21));
		A += knotsurface[i].area;
	}

	cout << "Input scaled by: " << scale[0] << ' ' << scale[1] << ' ' << scale[2] << "in x,y and z\n";

	return A;
}

double init_from_knot_file(void)
{
	double xt,yt,zt,dx,dy,dz,dl,lseg,Lh;    //temporary variables
	vector<double> px,py,pz,dr,ntx,nty,ntz;  //points, distances and tangents
	int npts;  //counter
	string temp;
	double L=0;
	int n,m,t,s,NKh;
	ifstream knotin;   //knot file(s)
	string filename, buff;
	stringstream ss;
	double maxxin = 0;
	double maxyin = 0;
	double maxzin = 0;
	double minxin = 0;
	double minyin = 0;
	double minzin = 0;

	NK=0;    //count total no. of points

	for (m=1; m<=ncomp; m++)
	{
		npts=0;

		ss.clear();
		ss.str("");
		if (ncomp==1) ss << knot_filename << ".txt";
		else ss << knot_filename << m << ".txt";

		filename = ss.str();
		knotin.open(filename.c_str());

		while(knotin.good())   //read in points for knot
		{
			if(getline(knotin,buff))
			{
				ss.clear();
				ss.str("");
				ss << buff;
				ss >> xt >> yt >> zt;
			}
			else break;
			/*needs changing*/
			px.push_back(xt);
			py.push_back(yt);             //store points
			pz.push_back(zt);

			if(xt > maxxin) maxxin = xt;  //find max and min points
			if(yt > maxyin) maxyin = yt;
			if(zt > maxzin) maxzin = zt;
			if(xt < minxin) minxin = xt;
			if(yt < minyin) minyin = yt;
			if(zt < minzin) minzin = zt;

			npts++;
		}

		knotin.close();

		/*rescale knot*/
		double scale[3];
		if(maxxin-minxin>0) scale[0] = xmax/(maxxin-minxin);
		else scale[0] = 1;
		if(maxyin-minyin>0) scale[1] = ymax/(maxyin-minyin);
		else scale[1] = 1;
		if(maxzin-minzin>0) scale[2] = zmax/(maxzin-minzin);
		else scale[2] = 1;
		double midpoint[3];
		midpoint[0] = 0.5*(maxxin+minxin);
		midpoint[1] = 0.5*(maxyin+minyin);
		midpoint[2] = 0.5*(maxzin+minzin);

		Lh=0;
		for(t=0; t<npts; t++)
		{
			px[t] = scale[0]*(px[t]-midpoint[0]);
			py[t] = scale[1]*(py[t]-midpoint[1]);
			pz[t] = scale[2]*(pz[t]-midpoint[2]);
			if(t>0)
			{
				dx = px[t] - px[t-1];
				dy = py[t] - py[t-1];     //distance between points
				dz = pz[t] - pz[t-1];
				dr.push_back(sqrt(dx*dx + dy*dy + dz*dz));
				ntx.push_back(dx/dr[t-1]);    //tangent direction to next pt, goes from 0:npts-1
				nty.push_back(dy/dr[t-1]);
				ntz.push_back(dz/dr[t-1]);
				Lh += dr[t-1];                //total length of link component
			}
		}
		/*Do final point*/
		dx = px[0] - px[npts-1];
		dy = py[0] - py[npts-1];     //distance between points
		dz = pz[0] - pz[npts-1];
		dr.push_back(sqrt(dx*dx + dy*dy + dz*dz));
		ntx.push_back(dx/dr[npts-1]);    //tangent direction to next pt, goes from 0:npts-1
		nty.push_back(dy/dr[npts-1]);
		ntz.push_back(dz/dr[npts-1]);
		Lh += dr[npts-1];                //total length of link component

		NKh = ((int) (2*Lh/h));  //Number of points to define knot with ~h/2 spacing
		dl = Lh/NKh;       //Actual spacing ~h/2

		//Start at p0
		X.push_back(px[0]);
		Y.push_back(py[0]);
		Z.push_back(pz[0]);

		n = 0; //input pt counter

		for(t=1;t<NKh;t++)
		{
			s = NK+t;
			X.push_back(X[s-1] + dl*ntx[n]);     //interpolate between input points
			Y.push_back(Y[s-1] + dl*nty[n]);
			Z.push_back(Z[s-1] + dl*ntz[n]);
			lseg = sqrt((X[s] - px[n])*(X[s] - px[n]) + (Y[s] - py[n])*(Y[s] - py[n]) + (Z[s] -               pz[n])*(Z[s] - pz[n]));     //distance from last input point
			while(lseg>dr[n])   //if we have passed next input point
			{
				n=n+1;     //move to next input point
				X[s] = px[n] + (lseg-dr[n-1])*ntx[n];    //add extra bit in next direction
				Y[s] = py[n] + (lseg-dr[n-1])*nty[n];
				Z[s] = pz[n] + (lseg-dr[n-1])*ntz[n];
				lseg = sqrt((X[s] - px[n])*(X[s] - px[n]) + (Y[s] - py[n])*(Y[s] - py[n]) + (Z[s] - pz[n])*(Z[s] - pz[n]));    //recalculate segment length
			}
			//cout << X[s] << ' ' << Y[s] << ' ' << Z[s] << '\n';
		}

		px.clear();
		py.clear();
		pz.clear();
		ntx.clear();
		nty.clear();
		ntz.clear();
		dr.clear();

		//smooth curve
		/*double *Xnew,*Ynew,*Znew;
		  double d2x,d2y,d2z;

		  for(n=0;n<10000;n++)
		  {
		  for(t=0;t<NKh;t++)
		  {
		  d2x = X[NK+incp(t,1,NKh)] - 2*X[NK+t] + X[NK+incp(t,-1,NKh)];
		  d2y = Y[NK+incp(t,1,NKh)] - 2*Y[NK+t] + Y[NK+incp(t,-1,NKh)];
		  d2z = Z[NK+incp(t,1,NKh)] - 2*Z[NK+t] + Z[NK+incp(t,-1,NKh)];
		  X[NK+t] += 0.01*d2x;
		  Y[NK+t] += 0.01*d2y;
		  Z[NK+t] += 0.01*d2z;
		  }
		  }*/

		for(t=0;t<NKh;t++)
		{
			dlx.push_back(0.5*(X[NK+incp(t,1,NKh)] - X[NK+incp(t,-1,NKh)]));   //central diff for tangent
			dly.push_back(0.5*(Y[NK+incp(t,1,NKh)] - Y[NK+incp(t,-1,NKh)]));
			dlz.push_back(0.5*(Z[NK+incp(t,1,NKh)] - Z[NK+incp(t,-1,NKh)]));
		}
		NK += NKh;
		L += Lh;    //keep track of total length and npts to define link
	}

	ofstream oknotfile;

	oknotfile.open("knotfile.vtk");

	oknotfile << "# vtk DataFile Version 3.0\nKnotin\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	oknotfile << "POINTS " << NK << " float\n";

	for(t=0; t<NK; t++)
	{
		oknotfile << X[t] << ' ' << Y[t] << ' ' << Z[t] << '\n';
	}

	oknotfile.close();

	return L;
}

/*************************Functions for B and Phi calcs*****************************/

void initial_cond(double *x, double *y, double *z, double *phi, int* missed)
{
	if(option == FROM_KNOT_FILE)
	{
		int *ignore;  //Points to ignore
		int *ignore1;
		double *Bx;  //Mag field
		double *By;
		double *Bz;
		double *Bmag;

		ignore = new int [Nx*Ny*Nz];
		ignore1 = new int [Nx*Ny*Nz];
		Bx = new double [Nx*Ny*Nz];
		By = new double [Nx*Ny*Nz];
		Bz = new double [Nx*Ny*Nz];
		Bmag = new double [Nx*Ny*Nz];

		cout << "Calculating B field...\n";
		time_t then = time(NULL);
		B_field_calc(x,y,z,Bx, By, Bz, Bmag, ignore, ignore1, missed);
		time_t now = time(NULL);
		cout << "B field calc took " << now - then << " seconds.\n";
		cout << "Calculating scalar potential...\n";
		then = time(NULL);
		phi_calc_B(Bx, By, Bz, Bmag, ignore, ignore1, missed, phi);
		now = time(NULL);
		cout << "Phi field calc took " << now - then << " seconds.\n";
		cout << "Printing B and phi...\n";
		print_B_phi(x, y, z, phi, missed);

		delete [] ignore;
		delete [] ignore1;
		delete [] Bx;
		delete [] By;
		delete [] Bz;
		delete [] Bmag;
	}
	else
	{
		cout << "Calculating scalar potential...\n";
		time_t then = time(NULL);
		phi_calc(x,y,z,phi);
		time_t now = time(NULL);
		cout << "Initialisation took " << now - then << " seconds.\n";
		cout << "Printing B and phi...\n";
		print_B_phi(x, y, z, phi, missed);
	}
}

void phi_calc(double *x, double *y, double *z, double *phi)
{
	int i,j,k,n,s;
	double rx,ry,rz,r;


#pragma omp parallel default(none) shared ( x, y, z, knotsurface, phi, NK ) private ( i, j, k, n, s, rx, ry, rz , r)
	{
#pragma omp for
		for(i=0;i<Nx;i++)
		{
			for(j=0; j<Ny; j++)
			{
				for(k=0; k<Nz; k++)
				{
					n = pt(i,j,k);
					phi[n] = 0;
					for(s=0;s<NK;s++)
					{
						rx = knotsurface[s].centre[0]-x[i];
						ry = knotsurface[s].centre[1]-y[j];
						rz = knotsurface[s].centre[2]-z[k];
						r = sqrt(rx*rx+ry*ry+rz*rz);
						if(r>0) phi[n] += (rx*knotsurface[s].normal[0] + ry*knotsurface[s].normal[1] + rz*knotsurface[s].normal[2])*knotsurface[s].area/(2*r*r*r);
					}
					while(phi[n]>M_PI) phi[n] -= 2*M_PI;
					while(phi[n]<-M_PI) phi[n] += 2*M_PI;
				}
			}
		}
	}

}

void B_field_calc(double *x, double *y, double *z, double *Bx, double *By, double *Bz, double *Bmag, int *ignore, int *ignore1, int *missed)
{
	int i,j,k,n,t;
	double lx,ly,lz,lmag;
	double coresize = lambda/(2*M_PI);

#pragma omp parallel default(none) shared ( x, y, z, X, Y, Z, Bx, By, Bz, missed, Bmag, NK, coresize, ignore, ignore1, dlx, dly, dlz ) private ( i, j, k, n, t, lx, ly ,lz, lmag)
	{
#pragma omp for
		for(i=0;i<Nx;i++)
		{
			for(j=0;j<Ny;j++)
			{
				for(k=0;k<Nz;k++)
				{
					n = pt(i,j,k);    //3D counter
					Bx[n] = 0;
					By[n] = 0;
					Bz[n] = 0;
					missed[n] = 1;   //intialise
					for(t=0;t<NK;t++)  //integrate over line
					{
						lx = x[i]-X[t];    //distance to point on line
						ly = y[j]-Y[t];
						lz = z[k]-Z[t];
						lmag = sqrt(lx*lx + ly*ly + lz*lz);
						if (lmag < 2*coresize) ignore[n]=1;   //do not use these points first time
						if (lmag < 0.5*coresize) ignore1[n]=1; //do not use these at all
						Bx[n] += (ly*dlz[t] - lz*dly[t])/(2*lmag*lmag*lmag);
						By[n] += (lz*dlx[t] - lx*dlz[t])/(2*lmag*lmag*lmag);
						Bz[n] += (lx*dly[t] - ly*dlx[t])/(2*lmag*lmag*lmag);
					}
					Bmag[n] = sqrt(Bx[n]*Bx[n] + By[n]*By[n] + Bz[n]*Bz[n]);
				}
			}
		}
	}
}

void phi_calc_B(double *Bx, double *By, double *Bz, double *Bmag, int *ignore, int *ignore1, int *missed, double *phi)
{
	int i0=(Nx+1)/2;
	int j0=(Ny+1)/2;   //base point for path integral
	int k0=(Nz+1)/2;
	int i[2],j[2],k[2],id,jd,kd,c1,c2,c3,pathlength,t,nt,ntm;
	double Bxmid,Bymid,Bzmid;
	int *pi,*pj,*pk;
	int n = pt(i0,j0,k0);

	missed[n]=0;  //matrix to store points where phi is not calculated
	phi[n]=0;

	pi = new int [Nx+Ny+Nz];
	pj = new int [Nx+Ny+Nz];
	pk = new int [Nx+Ny+Nz];

	//#pragma omp parallel default(none) shared ( X, Y, Z, Bx, By, Bz, phi, ignore, ignore1, missed, Bmag, pi, pj, pk) private ( i, j, k, i0, j0, k0, id, jd, kd, n, c1, c2, c3, Bxmid, Bymid, Bzmid, t, nt, ntm, pathlength )
	{
		//#pragma omp for
		for(id=0; id<(Nx+1)/2; id++)    //from zero to half grid points
		{
			for(jd=0; jd<(Ny+1)/2; jd++)
			{
				for(kd=0; kd<(Nz+1)/2; kd++)
				{
					i[0] = id;
					i[1] = Nx-1-id;
					j[0] = jd;
					j[1] = Ny-1-jd;
					k[0] = kd;
					k[1] = Nz-1-kd;
					for(c1=0;c1<2;c1++)    //count inwards from corners
					{
						for(c2=0;c2<2;c2++)
						{
							for(c3=0;c3<2;c3++)
							{
								n = pt(i[c1],j[c2],k[c3]);

								if(missed[n]==1 && ignore[n]==0)
								{
									pathlength = pathfind(i0,j0,k0,i[c1],j[c2],k[c3],pi,pj,pk,ignore,Bx,By,Bz,Bmag);  //find path to current point
									for (t=1;t<=pathlength;t++)   //travel along path
									{
										nt = pt(pi[t],pj[t],pk[t]);     //this point
										ntm = pt(pi[t-1],pj[t-1],pk[t-1]); //prev pt
										Bxmid = 0.5*(Bx[nt]+Bx[ntm]);
										Bymid = 0.5*(By[nt]+By[ntm]);   //midpoint
										Bzmid = 0.5*(Bz[nt]+Bz[ntm]);
										phi[nt] = phi[ntm] + h*(Bxmid*(pi[t]-pi[t-1]) + Bymid*(pj[t]-pj[t-1]) + Bzmid*(pk[t]-pk[t-1]));    //integrate along
										missed[nt]=0;
										while(phi[nt]>M_PI) phi[nt] -= 2*M_PI;
										while(phi[nt]<-M_PI) phi[nt] += 2*M_PI;
									}
								}
							}
						}
					}
				}
			}
		}
		//#pragma omp for
		for(id=0; id<Nx; id++)    //fill in ignore points but not ignore1
		{
			for(jd=0; jd<Ny; jd++)
			{
				for(kd=0; kd<Nz; kd++)
				{
					n = pt(id,jd,kd);
					if(ignore1[n]==0 && missed[n]==1)
					{
						pathlength = pathfind(i0,j0,k0,id,jd,kd,pi,pj,pk,ignore1,Bx,By,Bz,Bmag);
						for (t=1;t<=pathlength;t++)
						{
							nt = pt(pi[t],pj[t],pk[t]);
							ntm = pt(pi[t-1],pj[t-1],pk[t-1]);
							Bxmid = 0.5*(Bx[nt]+Bx[ntm]);
							Bymid = 0.5*(By[nt]+By[ntm]);
							Bzmid = 0.5*(Bz[nt]+Bz[ntm]);
							phi[nt] = phi[ntm] + h*(Bxmid*(pi[t]-pi[t-1]) + Bymid*(pj[t]-pj[t-1]) + Bzmid*(pk[t]-pk[t-1]));
							missed[nt]=0;
							while(phi[nt]>M_PI) phi[nt] -= 2*M_PI;
							while(phi[nt]<-M_PI) phi[nt] += 2*M_PI;
						}
					}
				}
			}
		}
	}

	delete [] pi;
	delete [] pj;
	delete [] pk;
}
int pathfind(int i0, int j0, int k0, int ie, int je, int ke, int *pi, int *pj, int *pk, int *ignore, double *Bx, double *By, double *Bz, double *Bmag)
{
	int io,jo,ko,ip,jp,kp,n,np,nu,go,stop,t=0;
	int *track;
	double MAX,weight1,weight2;
	track = new int [Nx*Ny*Nz];

	fill(track,track+Nx*Ny*Nz,0);

	pi[0] = i0;    //starting point for path
	pj[0] = j0;
	pk[0] = k0;
	int di = ie - i0;
	int dj = je - j0;
	int dk = ke - k0;   //distance to go to final point

	while (t<Nx+Ny+Nz && (abs(di)>0 || abs(dj)>0 || abs(dk)>0))  //until reaches end of path or path too long
	{
		n = pt(pi[t],pj[t],pk[t]);
		nu = pt(pi[t]+sign(di),pj[t]+sign(dj),pk[t]+sign(dk));  //check direct route
		if(ignore[nu] + track[nu]==0)  //if next space is available
		{
			pi[t+1] = pi[t] + sign(di);
			pj[t+1] = pj[t] + sign(dj);
			pk[t+1] = pk[t] + sign(dk);
			t++;   //move to next point
			n = pt(pi[t],pj[t],pk[t]);
			track[n]=1;
		}
		else
		{
			MAX = -10;    //compare point values
			go = 0;
			for(ip=-1; ip<2; ip++)
			{
				for(jp=-1; jp<2; jp++)
				{
					for(kp=-1; kp<2; kp++)  //check all neighbours
					{
						np = pt(pi[t]+ip,pj[t]+jp,pk[t]+kp);
						if(pi[t]+ip<Nx && pi[t]+ip>0 && pj[t]+jp<Ny && pj[t]+jp>0 && pk[t]+kp<Nz && pk[t]+kp>0) //If it is in the simulation box
						{
							stop = ignore[np] + track[np];  //not allowed to visit ignore points or previously visited points
						}
						else stop = 1;
						if(stop==0)
						{
							go=1;
							//weigting for which point to favour
							//direction of final point weighting
							weight1 = (di*ip + dj*jp + dk*kp)/(sqrt(di*di + dj*dj + dk*dk)*sqrt(ip*ip + jp*jp + kp*kp));
							//direction of B field weighting (helps to choose a direction around a barrier)
							weight2 = (Bx[np]*ip + By[np]*jp + Bz[np]*kp)/(Bmag[np]*sqrt(ip*ip + jp*jp + kp*kp));
							if(weight1 + weight2 > MAX)
							{
								MAX = weight1+weight2;
								io = ip;
								jo = jp;    //store the most favourable point
								ko = kp;
								n = pt(pi[t],pj[t],pk[t]);
								track[n]=1;  //track points visited
							}
						}
					}
				}
			}
			if(go==1)   //found a point to move to
			{
				pi[t+1] = pi[t]+io;
				pj[t+1] = pj[t]+jo;
				pk[t+1] = pk[t]+ko;
				t++;   //move to next point
			}
			else
			{
				if(t==0)
				{
					cout << "Could not find path to" << ie << ' ' << je << ' ' << ke << endl;
					return 0;
				}
				else
				{
					t--;  //go back to refind previous point
				}
			}
		}
		di = ie - pi[t];
		dj = je - pj[t];
		dk = ke - pk[t];
	}

	if (t==Nx+Ny+Nz) t=0; //couldn't find path

	delete [] track;

	return t;
}


/*************************Functions for FN dynamics*****************************/

void uv_initialise(double *phi, double *u, double *v, int* missed)
{
	int n;

	for(n=0; n<Nx*Ny*Nz; n++)
	{
		u[n] = (2*cos(phi[n]) - 0.4);
		v[n] = (sin(phi[n]) - 0.4);
		if(option==FROM_KNOT_FILE && missed[n]==1)
		{
			u[n] = -0.4;
			v[n] = -0.4;
		}
	}
}

void crossgrad_calc(double *x, double *y, double *z, double *u, double *v, double *ucvx, double *ucvy, double *ucvz)
{
	int i,j,k,n,kup,kdwn;
	double dxu,dyu,dzu,dxv,dyv,dzv,ucvmag,ucvmax;
	knotcurve.push_back(knotpoint());
	ucvmax = -1.0; // should always be +ve, so setting it to an initially -ve # means it always gets written to once.  
	for(i=0;i<Nx;i++)
	{
		for(j=0; j<Ny; j++)
		{
			for(k=0; k<Nz; k++)   //Central difference
			{
				if(periodic)   //check for periodic boundaries
				{
					kup = incp(k,1,Nz);
					kdwn = incp(k,-1,Nz);
				}
				else
				{
					kup = incw(k,1,Nz);
					kdwn = incw(k,-1,Nz);
				}
				dxu = 0.5*(u[pt(incw(i,1,Nx),j,k)]-u[pt(incw(i,-1,Nx),j,k)])/h;
				dxv = 0.5*(v[pt(incw(i,1,Nx),j,k)]-v[pt(incw(i,-1,Nx),j,k)])/h;
				dyu = 0.5*(u[pt(i,incw(j,1,Ny),k)]-u[pt(i,incw(j,-1,Ny),k)])/h;
				dyv = 0.5*(v[pt(i,incw(j,1,Ny),k)]-v[pt(i,incw(j,-1,Ny),k)])/h;
				dzu = 0.5*(u[pt(i,j,kup)]-u[pt(i,j,kdwn)])/h;
				dzv = 0.5*(v[pt(i,j,kup)]-v[pt(i,j,kdwn)])/h;
				n = pt(i,j,k);
				ucvx[n] = dyu*dzv - dzu*dyv;
				ucvy[n] = dzu*dxv - dxu*dzv;    //Grad u cross Grad v
				ucvz[n] = dxu*dyv - dyu*dxv;
				ucvmag = sqrt(ucvx[n]*ucvx[n] + ucvy[n]*ucvy[n] + ucvz[n]*ucvz[n]);
				if(ucvmag > ucvmax)
				{
					ucvmax = ucvmag;
					knotcurve[0].xcoord=x[i];
					knotcurve[0].ycoord=y[j];
					knotcurve[0].zcoord=z[k];
				}
			}
		}
	}
	if(ucvmax<0.1) knotexists = false;
	else knotexists = true;
}

void find_knot_properties(double *x, double *y, double *z, double *ucvx, double *ucvy, double *ucvz, double *u, double t, gsl_min_fminimizer* minimizerstate)
{
	int i,j,k,idwn,jdwn,kdwn,m,pts,iinc,jinc,kinc,attempts;
	int s=1;
	bool finish=false;
	double ucvxs, ucvys, ucvzs, graducvx, graducvy, graducvz, prefactor, xd, yd ,zd, norm, fx, fy, fz, xdiff, ydiff, zdiff;

	/*calculate local direction of grad u x grad v (the tangent to the knot curve) at point s-1, then move to point s by moving along tangent + unit confinement force*/
	while (finish==false)
	{
		norm=0;
		/**Find nearest gridpoint**/
		idwn = (int) ((knotcurve[s-1].xcoord/h) - 0.5 + Nx/2.0);
		jdwn = (int) ((knotcurve[s-1].ycoord/h) - 0.5 + Ny/2.0);
		kdwn = (int) ((knotcurve[s-1].zcoord/h) - 0.5 + Nz/2.0);
		if(idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1) break;
		pts=0;
		ucvxs=0;
		ucvys=0;
		ucvzs=0;
		/*curve to gridpoint down distance*/
		xd = (knotcurve[s-1].xcoord - x[idwn])/h;
		yd = (knotcurve[s-1].ycoord - y[jdwn])/h;
		zd = (knotcurve[s-1].zcoord - z[kdwn])/h;
		for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
		{
			/* Work out increments*/
			iinc = m%2;
			jinc = (m/2)%2;
			kinc = (m/4)%2;
			/*Loop over nearest points*/
			i = incw(idwn, iinc, Nx);
			j = incw(jdwn, jinc, Ny);
			if(periodic) k = incp(kdwn,kinc, Nz);
			else k = incw(kdwn,kinc, Nz);
			prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
			/*interpolate grad u x grad v over nearest points*/
			ucvxs += prefactor*ucvx[pt(i,j,k)];
			ucvys += prefactor*ucvy[pt(i,j,k)];
			ucvzs += prefactor*ucvz[pt(i,j,k)];
		}
		norm = sqrt(ucvxs*ucvxs + ucvys*ucvys + ucvzs*ucvzs);
		ucvxs = ucvxs/norm; //normalise
		ucvys = ucvys/norm; //normalise
		ucvzs = ucvzs/norm; //normalise

		// okay we have our first guess, move forward in this direction
		double testx = knotcurve[s-1].xcoord + 0.5*ucvxs*lambda/(2*M_PI);
		double testy = knotcurve[s-1].ycoord + 0.5*ucvys*lambda/(2*M_PI);
		double testz = knotcurve[s-1].zcoord + 0.5*ucvzs*lambda/(2*M_PI);

		// now get the grad at this point
		idwn = (int) ((testx/h) - 0.5 + Nx/2.0);
		jdwn = (int) ((testy/h) - 0.5 + Ny/2.0);
		kdwn = (int) ((testz/h) - 0.5 + Nz/2.0);
		if(idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1) break;
		pts=0;
		graducvx=0;
		graducvy=0;
		graducvz=0;
		/*curve to gridpoint down distance*/
		xd = (testx - x[idwn])/h;
		yd = (testy - y[jdwn])/h;
		zd = (testz - z[kdwn])/h;
		for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
		{
			/* Work out increments*/
			iinc = m%2;
			jinc = (m/2)%2;
			kinc = (m/4)%2;
			/*Loop over nearest points*/
			i = incw(idwn, iinc, Nx);
			j = incw(jdwn, jinc, Ny);
			if(periodic) k = incp(kdwn,kinc, Nz);
			else k = incw(kdwn,kinc, Nz);
			prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
			/*interpolate gradients of |grad u x grad v|*/
			graducvx += prefactor*(sqrt(ucvx[pt(incw(i,1,Nx),j,k)]*ucvx[pt(incw(i,1,Nx),j,k)] + ucvy[pt(incw(i,1,Nx),j,k)]*ucvy[pt(incw(i,1,Nx),j,k)] + ucvz[pt(incw(i,1,Nx),j,k)]*ucvz[pt(incw(i,1,Nx),j,k)]) - sqrt(ucvx[pt(incw(i,-1,Nx),j,k)]*ucvx[pt(incw(i,-1,Nx),j,k)] + ucvy[pt(incw(i,-1,Nx),j,k)]*ucvy[pt(incw(i,-1,Nx),j,k)] + ucvz[pt(incw(i,-1,Nx),j,k)]*ucvz[pt(incw(i,-1,Nx),j,k)]))/(2*h);
			graducvy += prefactor*(sqrt(ucvx[pt(i,incw(j,1,Ny),k)]*ucvx[pt(i,incw(j,1,Ny),k)] + ucvy[pt(i,incw(j,1,Ny),k)]*ucvy[pt(i,incw(j,1,Ny),k)] + ucvz[pt(i,incw(j,1,Ny),k)]*ucvz[pt(i,incw(j,1,Ny),k)]) - sqrt(ucvx[pt(i,incw(j,-1,Ny),k)]*ucvx[pt(i,incw(j,-1,Ny),k)] + ucvy[pt(i,incw(j,-1,Ny),k)]*ucvy[pt(i,incw(j,-1,Ny),k)] + ucvz[pt(i,incw(j,-1,Ny),k)]*ucvz[pt(i,incw(j,-1,Ny),k)]))/(2*h);
			if(periodic) graducvz += prefactor*(sqrt(ucvx[pt(i,j,incp(k,1,Nz))]*ucvx[pt(i,j,incp(k,1,Nz))] + ucvy[pt(i,j,incp(k,1,Nz))]*ucvy[pt(i,j,incp(k,1,Nz))] + ucvz[pt(i,j,incp(k,1,Nz))]*ucvz[pt(i,j,incp(k,1,Nz))]) - sqrt(ucvx[pt(i,j,incp(k,-1,Nz))]*ucvx[pt(i,j,incp(k,-1,Nz))] + ucvy[pt(i,j,incp(k,-1,Nz))]*ucvy[pt(i,j,incp(k,-1,Nz))] + ucvz[pt(i,j,incp(k,-1,Nz))]*ucvz[pt(i,j,incp(k,-1,Nz))]))/(2*h);
			else graducvz += prefactor*(sqrt(ucvx[pt(i,j,incw(k,1,Nz))]*ucvx[pt(i,j,incw(k,1,Nz))] + ucvy[pt(i,j,incw(k,1,Nz))]*ucvy[pt(i,j,incw(k,1,Nz))] + ucvz[pt(i,j,incw(k,1,Nz))]*ucvz[pt(i,j,incw(k,1,Nz))]) - sqrt(ucvx[pt(i,j,incw(k,-1,Nz))]*ucvx[pt(i,j,incw(k,-1,Nz))] + ucvy[pt(i,j,incw(k,-1,Nz))]*ucvy[pt(i,j,incw(k,-1,Nz))] + ucvz[pt(i,j,incw(k,-1,Nz))]*ucvz[pt(i,j,incw(k,-1,Nz))]))/(2*h);
			
		}
		knotcurve.push_back(knotpoint());
		fx = (graducvx - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvxs);   //confining force perpendicular to curve direction
		fy = (graducvy - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvys);
		fz = (graducvz - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvzs);
		norm = sqrt(fx*fx + fy*fy + fz*fz);  //force is normalised, this could mean curve oscillates around the centre, but it does moke the confinement magnitude easier to control
		fx = fx/norm;
		fy = fy/norm;
		fz = fz/norm;

		// okay we have our direction to perfrom the line minimisation in

		gsl_vector* v = gsl_vector_alloc (3);
		gsl_vector_set (v, 0, testx);
		gsl_vector_set (v, 1, testy);
		gsl_vector_set (v, 2, testz);
		gsl_vector* f = gsl_vector_alloc (3);
		gsl_vector_set (f, 0, fx);
		gsl_vector_set (f, 1, fy);
		gsl_vector_set (f, 2, fz);

		struct parameters params;
		struct parameters* pparams = &params;
		pparams->x = x; pparams->y=y;pparams->z=z;
		pparams->ucvx=ucvx;pparams->ucvy=ucvy; pparams->ucvz = ucvz;
		pparams->v = v; pparams->f = f;
		// some initial values
		double minimum = 0;
		double a = -lambda/(8*M_PI);
		double b = lambda/(8*M_PI);
		gsl_function F;
		F.function = &my_f;
		F.params = (void*) pparams;
		gsl_min_fminimizer_set (minimizerstate, &F, minimum, a, b);

		int iter=0;
		int status =0;
		do
		{
			iter++;
			status = gsl_min_fminimizer_iterate (minimizerstate);

			minimum = gsl_min_fminimizer_x_minimum (minimizerstate);
			a = gsl_min_fminimizer_x_lower (minimizerstate);
			b = gsl_min_fminimizer_x_upper (minimizerstate);

			status = gsl_min_test_interval (a, b, 0.001, 0.0);

		}
		while (status == GSL_CONTINUE && iter < 500);


		gsl_vector_scale(f,minimum);
		gsl_vector_add(v,f);
		knotcurve[s].xcoord = gsl_vector_get(v, 0);
		knotcurve[s].ycoord= gsl_vector_get(v, 1);
		knotcurve[s].zcoord= gsl_vector_get(v, 2);
		
		gsl_vector_free(v);
		gsl_vector_free(f);
		
		xdiff = knotcurve[0].xcoord - knotcurve[s].xcoord;     //distance from start/end point
		ydiff = knotcurve[0].ycoord - knotcurve[s].ycoord;
		zdiff = knotcurve[0].zcoord - knotcurve[s].zcoord;
		if(sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff) < lambda/(2*M_PI) && s > 200) finish = true;
		if(s>50000) finish = true;
		s++;
	}
	/*Fill in remaining points*/
	double dx = xdiff/16.0;
	double dy = ydiff/16.0;
	double dz = zdiff/16.0;
	for(m=0; m<15; m++)
	{
		knotcurve.push_back(knotpoint());
		knotcurve[s+m].xcoord = knotcurve[s+m-1].xcoord + dx;
		knotcurve[s+m].ycoord = knotcurve[s+m-1].ycoord + dy;
		knotcurve[s+m].zcoord = knotcurve[s+m-1].zcoord + dz;
	}

	int NP = knotcurve.size();  //store number of points in knot curve

	/*******Erase some points********/

	//~ for(s=0; s<NP%8; s++) knotcurve.pop_back();    //delete last couple of elements
	//~ for(s=0; s<NP/8; s++)                          //delete 7 of every 8 elements
	//~ {
		//~ knotcurve.erase(knotcurve.end()-s-8,knotcurve.end()-s-1);
	//~ }

	/********************************/

	NP = knotcurve.size();  //update number of points in knot curve

	/*******Vertex averaging*********/

	double totlength, dl;
	//~ for(i=0;i<3;i++)   //repeat a couple of times because of end point
	//~ {
		//~ totlength=0;
		//~ for(s=0; s<NP; s++)   //Work out total length of curve
		//~ {
			//~ dx = knotcurve[incp(s,1,NP)].xcoord - knotcurve[s].xcoord;
			//~ dy = knotcurve[incp(s,1,NP)].ycoord - knotcurve[s].ycoord;
			//~ dz = knotcurve[incp(s,1,NP)].zcoord - knotcurve[s].zcoord;
			//~ totlength += sqrt(dx*dx + dy*dy + dz*dz);
		//~ }
		//~ dl = totlength/NP;
		//~ for(s=0; s<NP; s++)    //Move points to have spacing dl
		//~ {
			//~ dx = knotcurve[incp(s,1,NP)].xcoord - knotcurve[s].xcoord;
			//~ dy = knotcurve[incp(s,1,NP)].ycoord - knotcurve[s].ycoord;
			//~ dz = knotcurve[incp(s,1,NP)].zcoord - knotcurve[s].zcoord;
			//~ norm = sqrt(dx*dx + dy*dy + dz*dz);
			//~ knotcurve[incp(s,1,NP)].xcoord = knotcurve[s].xcoord + dl*dx/norm;
			//~ knotcurve[incp(s,1,NP)].ycoord = knotcurve[s].ycoord + dl*dy/norm;
			//~ knotcurve[incp(s,1,NP)].zcoord = knotcurve[s].zcoord + dl*dz/norm;
		//~ }
	//~ }

	/********************************/

	/*****Writhe and twist integrals******/
	NP = knotcurve.size();  //store number of points in knot curve
	//if(t==50) cout << "Number of points: " << NP << '\n';
	double totwrithe = 0;
	double tottwist = 0;
	double ds = 2*M_PI/NP;
	double dxds, dyds, dzds, dxdm, dydm, dzdm, dxu, dyu, dzu, dxup, dyup, dzup, bx, by, bz, check;
	totlength = 0;

	/******************Interpolate direction of grad u for twist calc*******/
	/**Find nearest gridpoint**/
	for(s=0; s<NP; s++)
	{
		idwn = (int) ((knotcurve[s].xcoord/h) - 0.5 + Nx/2.0);
		jdwn = (int) ((knotcurve[s].ycoord/h) - 0.5 + Ny/2.0);
		kdwn = (int) ((knotcurve[s].zcoord/h) - 0.5 + Nz/2.0);
		if(idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1) break;
		dxu=0;
		dyu=0;
		dzu=0;
		/*curve to gridpoint down distance*/
		xd = (knotcurve[s].xcoord - x[idwn])/h;
		yd = (knotcurve[s].ycoord - y[jdwn])/h;
		zd = (knotcurve[s].zcoord - z[kdwn])/h;
		for(m=0;m<8;m++)  //linear interpolation of 8 NNs
		{
			/* Work out increments*/
			iinc = m%2;
			jinc = (m/2)%2;
			kinc = (m/4)%2;
			/*Loop over nearest points*/
			i = incw(idwn, iinc, Nx);
			j = incw(jdwn, jinc, Ny);
			if(periodic) k = incp(kdwn,kinc, Nz);
			else k = incw(kdwn,kinc, Nz);
			prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);   //terms of the form (1-xd)(1-yd)zd etc. (interpolation coefficient)
			/*interpolate grad u over nearest points*/
			dxu += prefactor*0.5*(u[pt(incw(i,1,Nx),j,k)] -  u[pt(incw(i,-1,Nx),j,k)])/h;  //central diff
			dyu += prefactor*0.5*(u[pt(i,incw(j,1,Ny),k)] -  u[pt(i,incw(j,-1,Ny),k)])/h;
			if(periodic) prefactor*0.5*(u[pt(i,j,incp(k,1,Nz))] -  u[pt(i,j,incp(k,-1,Nz))])/h;
			else  dzu += prefactor*0.5*(u[pt(i,j,incw(k,1,Nz))] -  u[pt(i,j,incw(k,-1,Nz))])/h;
		}
		//project du onto perp of tangent direction first
		dx = 0.5*(knotcurve[incp(s,1,NP)].xcoord - knotcurve[incp(s,-1,NP)].xcoord);   //central diff as a is defined on the points
		dy = 0.5*(knotcurve[incp(s,1,NP)].ycoord - knotcurve[incp(s,-1,NP)].ycoord);
		dz = 0.5*(knotcurve[incp(s,1,NP)].zcoord - knotcurve[incp(s,-1,NP)].zcoord);
		dxup = dxu - (dxu*dx + dyu*dy + dzu*dz)*dx/(dx*dx+dy*dy+dz*dz);               //Grad u_j * (delta_ij - t_i t_j)
		dyup = dyu - (dxu*dx + dyu*dy + dzu*dz)*dy/(dx*dx+dy*dy+dz*dz);
		dzup = dzu - (dxu*dx + dyu*dy + dzu*dz)*dz/(dx*dx+dy*dy+dz*dz);
		/*Vector a is the normalised gradient of u, should point in direction of max u perp to t*/
		norm = sqrt(dxup*dxup+dyup*dyup+dzup*dzup);
		knotcurve[s].ax = dxup/norm;
		knotcurve[s].ay = dyup/norm;
		knotcurve[s].az = dzup/norm;
	}

	/***Do the integrals**/
	for(s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
	{
		dxds = (knotcurve[incp(s,1,NP)].xcoord - knotcurve[s].xcoord)/(ds);
		dyds = (knotcurve[incp(s,1,NP)].ycoord - knotcurve[s].ycoord)/(ds);
		dzds = (knotcurve[incp(s,1,NP)].zcoord - knotcurve[s].zcoord)/(ds);
		/*These quantities defined on line connecting points s and s+1*/
		knotcurve[s].writhe = 0;
		knotcurve[s].length = sqrt(dxds*dxds + dyds*dyds + dzds*dzds)*ds;  //actual length of thing
		bx = (knotcurve[incp(s,1,NP)].ax - knotcurve[s].ax)/ds;
		by = (knotcurve[incp(s,1,NP)].ay - knotcurve[s].ay)/ds;
		bz = (knotcurve[incp(s,1,NP)].az - knotcurve[s].az)/ds;
		knotcurve[s].twist = (dxds*(knotcurve[s].ay*bz - knotcurve[s].az*by) + dyds*(knotcurve[s].az*bx - knotcurve[s].ax*bz) + dzds*(knotcurve[s].ax*by - knotcurve[s].ay*bx))/(2*M_PI*sqrt(dxds*dxds + dyds*dyds + dzds*dzds));
		/*Check this is actually normal to tangent*/
		/*check = fabs(0.5*(knotcurve[s].ax + knotcurve[incp(s,1,NP)].ax)*dxds + 0.5*(knotcurve[s].ay + knotcurve[incp(s,1,NP)].ay)*dyds + 0.5*(knotcurve[s].az + knotcurve[incp(s,1,NP)].az)*dzds)/sqrt(dxds*dxds + dyds*dyds + dzds*dzds);
		  if(check>0.01) cout << s << ": (" << knotcurve[s].xcoord << ", " << knotcurve[s].ycoord << ", " << knotcurve[s].zcoord << "). Grad u . t = " << check << '\n';*/
		for(m=0; m<NP; m++)
		{
			if(s != m)
			{
				xdiff = 0.5*(knotcurve[incp(s,1,NP)].xcoord + knotcurve[s].xcoord - knotcurve[incp(m,1,NP)].xcoord - knotcurve[m].xcoord);   //interpolate, consistent with fwd diff
				ydiff = 0.5*(knotcurve[incp(s,1,NP)].ycoord + knotcurve[s].ycoord - knotcurve[incp(m,1,NP)].ycoord - knotcurve[m].ycoord);
				zdiff = 0.5*(knotcurve[incp(s,1,NP)].zcoord + knotcurve[s].zcoord - knotcurve[incp(m,1,NP)].zcoord - knotcurve[m].zcoord);
				dxdm = (knotcurve[incp(m,1,NP)].xcoord - knotcurve[m].xcoord)/(ds);
				dydm = (knotcurve[incp(m,1,NP)].ycoord - knotcurve[m].ycoord)/(ds);
				dzdm = (knotcurve[incp(m,1,NP)].zcoord - knotcurve[m].zcoord)/(ds);
				knotcurve[s].writhe += ds*(xdiff*(dyds*dzdm - dzds*dydm) + ydiff*(dzds*dxdm - dxds*dzdm) + zdiff*(dxds*dydm - dyds*dxdm))/(4*M_PI*(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)*sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff));
			}
		}
		/*Add on writhe, twist and length*/
		totwrithe += knotcurve[s].writhe*ds;
		totlength += knotcurve[s].length;
		tottwist  += knotcurve[s].twist*ds;
	}

	print_knot(x,y,z,t, knotcurve);

	/***Write values to file*******/
	ofstream wrout;
	wrout.open("writhe.txt",ios_base::app);
	wrout << t << '\t' << totwrithe << '\t' << tottwist << '\t' << totlength << '\n';
	wrout.close();
}

void uv_update(double *u, double *v, double *ku, double *kv, double *kut, double *kvt, double *uold, double *vold)
{
	int i,j,k,l,n,kup,kdwn;
	double D2u;

#pragma omp for
	for(i=0;i<Nx;i++)
	{
		for(j=0; j<Ny; j++)
		{
			for(k=0; k<Nz; k++)
			{
				n = pt(i,j,k);
				uold[n] = u[n];  //old value of u
				vold[n] = v[n];  //old value of v
				kut[n] = 0;
				kvt[n] = 0;
			}
		}
	}

	for(l=0;l<4;l++)  //u and v update for each fractional time step
	{
#pragma omp for
		for(i=0;i<Nx;i++)
		{
			for(j=0; j<Ny; j++)
			{
				for(k=0; k<Nz; k++)   //Central difference
				{
					n = pt(i,j,k);
					if(periodic)   //check for periodic boundaries
					{
						kup = incp(k,1,Nz);
						kdwn = incp(k,-1,Nz);
					}
					else
					{
						kup = incw(k,1,Nz);
						kdwn = incw(k,-1,Nz);
					}
					D2u = oneoverhsq*(u[pt(incw(i,1,Nx),j,k)] + u[pt(incw(i,-1,Nx),j,k)] + u[pt(i,incw(j,1,Ny),k)] + u[pt(i,incw(j,-1,Ny),k)] + u[pt(i,j,kup)] + u[pt(i,j,kdwn)] - 6.0*u[n]);
					ku[n] = oneoverepsilon*(u[n] - u[n]*u[n]*u[n]/3.0 - v[n]) + D2u;
					kv[n] = epsilon*(u[n] + beta - gam*v[n]);
				}
			}
		}

		switch (l)
		{
			case 0:
				{
					uv_add(u,v,uold,vold,ku,kv,kut,kvt,0.5,1.0);   //add k1 to uv and add to total k
				}
				break;

			case 1:
				{
					uv_add(u,v,uold,vold,ku,kv,kut,kvt,0.5,2.0);   //add k2 to uv and add to total k
				}
				break;

			case 2:
				{
					uv_add(u,v,uold,vold,ku,kv,kut,kvt,1.0,2.0);      //add k3 to uv and add to total k
				}
				break;

			case 3:
				{
#pragma omp for
					for(i=0;i<Nx;i++)
					{
						for(j=0; j<Ny; j++)
						{
							for(k=0; k<Nz; k++)  //update
							{
								n = pt(i,j,k);
								u[n] = uold[n] + dtime*sixth*(kut[n]+ku[n]);
								v[n] = vold[n] + dtime*sixth*(kvt[n]+kv[n]);
							}
						}
					}
				}
				break;

			default:
				break;
		}
	}
}

void uv_add(double *u, double *v, double* uold, double *vold, double *ku, double *kv, double *kut, double *kvt, double inc, double coeff)
{
	int i,j,k,n;

#pragma omp for
	for(i=0;i<Nx;i++)
	{
		for(j=0; j<Ny; j++)
		{
			for(k=0; k<Nz; k++)  //update
			{
				n = pt(i,j,k);
				u[n] = uold[n] + dtime*inc*ku[n];
				v[n] = vold[n] + dtime*inc*kv[n];
				kut[n] += coeff*ku[n];
				kvt[n] += coeff*kv[n];
			}
		}
	}

}

void uv_update_euler(double *u, double *v, double *D2u)
{
	int i,j,k,l,n,kup,kdwn;

#pragma omp for
	for(i=0;i<Nx;i++)
	{
		for(j=0; j<Ny; j++)
		{
			for(k=0; k<Nz; k++)
			{
				n = pt(i,j,k);
				if(periodic)   //check for periodic boundaries
				{
					kup = incp(k,1,Nz);
					kdwn = incp(k,-1,Nz);
				}
				else
				{
					kup = incw(k,1,Nz);
					kdwn = incw(k,-1,Nz);
				}
				D2u[n] = (u[pt(incw(i,1,Nx),j,k)] + u[pt(incw(i,-1,Nx),j,k)] + u[pt(i,incw(j,1,Ny),k)] + u[pt(i,incw(j,-1,Ny),k)] + u[pt(i,j,kup)] + u[pt(i,j,kdwn)] - 6*u[n])/(h*h);
			}
		}
	}

#pragma omp for
	for(i=0;i<Nx;i++)
	{
		for(j=0; j<Ny; j++)
		{
			for(k=0; k<Nz; k++)
			{
				n = pt(i,j,k);
				u[n] = u[n] + dtime*((u[n] - u[n]*u[n]*u[n]/3 - v[n])/epsilon + D2u[n]);
				v[n] = v[n] + dtime*(epsilon*(u[n] + beta - gam*v[n]));
			}
		}
	}
}

/*************************File reading and writing*****************************/

void print_uv(double *x, double *y, double *z, double *u, double *v, double *ucvx, double *ucvy, double *ucvz, double t)
{
	int i,j,k,n;
	stringstream ss;
	ss << "uv_plot" << t << ".vtk";
	ofstream uvout (ss.str().c_str());

	uvout << "# vtk DataFile Version 3.0\nUV fields\nASCII\nDATASET STRUCTURED_POINTS\n";
	uvout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
	uvout << "ORIGIN " << x[0] << ' ' << y[0] << ' ' << z[0] << '\n';
	uvout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
	uvout << "POINT_DATA " << Nx*Ny*Nz << '\n';
	uvout << "SCALARS u float\nLOOKUP_TABLE default\n";


	for(k=0; k<Nz; k++)
	{
		for(j=0; j<Ny; j++)
		{
			for(i=0; i<Nx; i++)
			{
				n = pt(i,j,k);
				uvout << u[n] << '\n';
			}
		}
	}

	uvout << "SCALARS v float\nLOOKUP_TABLE default\n";


	for(k=0; k<Nz; k++)
	{
		for(j=0; j<Ny; j++)
		{
			for(i=0; i<Nx; i++)
			{
				n = pt(i,j,k);
				uvout << v[n] << '\n';
			}
		}
	}

	uvout << "SCALARS ucrossv float\nLOOKUP_TABLE default\n";

	for(k=0; k<Nz; k++)
	{
		for(j=0; j<Ny; j++)
		{
			for(i=0; i<Nx; i++)
			{
				n = pt(i,j,k);
				uvout << sqrt(ucvx[n]*ucvx[n] + ucvy[n]*ucvy[n] + ucvz[n]*ucvz[n]) << '\n';
			}
		}
	}

	uvout.close();
}

void print_B_phi(double *x, double *y, double*z, double *phi, int* missed)
{
	int i,j,k,n;
	string fn = "phi.vtk";

	ofstream Bout (fn.c_str());

	Bout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
	Bout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
	Bout << "ORIGIN " << x[0] << ' ' << y[0] << ' ' << z[0] << '\n';
	Bout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
	Bout << "POINT_DATA " << Nx*Ny*Nz << '\n';
	Bout << "SCALARS Phi float\nLOOKUP_TABLE default\n";
	for(k=0; k<Nz; k++)
	{
		for(j=0; j<Ny; j++)
		{
			for(i=0; i<Nx; i++)
			{
				n = pt(i,j,k);
				Bout << phi[n] << '\n';
			}
		}
	}

	if(option==FROM_KNOT_FILE)
	{
		Bout << "\n\nSCALARS Missed int\nLOOKUP_TABLE default\n";
		for(k=0; k<Nz; k++)
		{
			for(j=0; j<Ny; j++)
			{
				for(i=0; i<Nx; i++)
				{
					n = pt(i,j,k);
					Bout << missed[n] << '\n';
				}
			}
		}
	}

	Bout.close();
}

void print_info(int Nx, int Ny, int Nz, double dtime, double h, const bool periodic,  int option, string knot_filename, string B_filename)
{
	string fn = "info.txt";

	time_t rawtime;
	struct tm * timeinfo;

	time (&rawtime);
	timeinfo = localtime (&rawtime);

	ofstream infoout (fn.c_str());

	infoout << "run started at\t" << asctime(timeinfo) << "\n";
	infoout << "Number of grid points\t" << Nx << '\t' << Ny << '\t' << Nz << '\n';
	infoout << "timestep\t" << dtime << '\n';
	infoout << "Spacing\t" << h << '\n';
	infoout << "Periodic\t" << periodic << '\n';
	infoout << "initoptions\t" << option << '\n';
	infoout << "knot filename\t" << knot_filename << '\n';
	infoout << "B or uv filename\t" << B_filename << '\n';
	infoout.close();
}

void print_knot(double *x, double *y, double *z, double t, std::vector<knotpoint>& knotcurve)
{
	stringstream ss;
	ss << "knotplot" << t << ".vtk";
	ofstream knotout (ss.str().c_str());

	int i;
	int n = knotcurve.size();

	knotout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	knotout << "POINTS " << n << " float\n";

	for(i=0; i<n; i++)
	{
		knotout << knotcurve[i].xcoord << ' ' << knotcurve[i].ycoord << ' ' << knotcurve[i].zcoord << '\n';
	}

	knotout << "\n\nCELLS " << n << ' ' << 3*n << '\n';

	for(i=0; i<n; i++)
	{
		knotout << 2 << ' ' << i << ' ' << incp(i,1,n) << '\n';
	}

	knotout << "\n\nCELL_TYPES " << n << '\n';

	for(i=0; i<n; i++)
	{
		knotout << "3\n";
	}

	knotout << "\n\nPOINT_DATA " << n << "\n\n";
	knotout << "\nVECTORS A float\n";
	for(i=0; i<n; i++)
	{
		knotout << knotcurve[i].ax << ' ' << knotcurve[i].ay << ' ' << knotcurve[i].az << '\n';
	}

	knotout << "\n\nCELL_DATA " << n << "\n\n";
	knotout << "\nSCALARS Writhe float\nLOOKUP_TABLE default\n";
	for(i=0; i<n; i++)
	{
		knotout << knotcurve[i].writhe << '\n';
	}

	knotout << "\nSCALARS Twist float\nLOOKUP_TABLE default\n";
	for(i=0; i<n; i++)
	{
		knotout << knotcurve[i].twist << '\n';
	}

	knotout << "\nSCALARS Length float\nLOOKUP_TABLE default\n";
	for(i=0; i<n; i++)
	{
		knotout << knotcurve[i].length << '\n';
	}

	knotout.close();
}

int phi_file_read(double *phi)
{
	string temp,buff;
	stringstream ss;
	ifstream fin (B_filename.c_str());
	int i,j,k,n;

	for(i=0;i<10;i++)
	{
		if(fin.good())
		{
			if(getline(fin,buff)) temp = buff;
		}
		else
		{
			cout << "Something went wrong!\n";
			return 1;
		}
	}

	for(k=0; k<Nz; k++)
	{
		for(j=0; j<Ny; j++)
		{
			for(i=0; i<Nx; i++)
			{
				n=pt(i,j,k);
				ss.clear();
				ss.str("");
				if(fin.good())
				{
					if(getline(fin,buff))
					{
						ss << buff;
						ss >> phi[n];
					}
				}
				else
				{
					cout << "Something went wrong!\n";
					return 1;
				}
			}
		}
	}

	fin.close();

	return 0;
}

int uvfile_read(double *u, double *v)
{
	string temp,buff;
	stringstream ss;
	ifstream fin (B_filename.c_str());
	int i,j,k,n;

	for(i=0;i<10;i++)
	{
		if(fin.good())
		{
			if(getline(fin,buff)) temp = buff;
		}
		else
		{
			cout << "Something went wrong!\n";
			return 1;
		}
	}

	for(k=0; k<Nz; k++)
	{
		for(j=0; j<Ny; j++)
		{
			for(i=0; i<Nx; i++)
			{
				n=pt(i,j,k);
				ss.clear();
				ss.str("");
				if(fin.good())
				{
					if(getline(fin,buff))
					{
						ss << buff;
						ss >> u[n];
					}
				}
				else
				{
					cout << "Something went wrong!\n";
					return 1;
				}
			}
		}
	}

	for(i=0;i<2;i++)
	{
		if(fin.good())
		{
			if(getline(fin,buff)) temp = buff;
		}
		else
		{
			cout << "Something went wrong!\n";
			return 1;
		}
	}

	for(k=0; k<Nz; k++)
	{
		for(j=0; j<Ny; j++)
		{
			for(i=0; i<Nx; i++)
			{
				n=pt(i,j,k);
				ss.clear();
				ss.str("");
				if(fin.good())
				{
					if(getline(fin,buff)) ss << buff;
					ss >> v[n];
				}
				else
				{
					cout << "Something went wrong!\n";
					return 1;
				}
			}
		}
	}

	fin.close();

	return 0;
}
int intersect3D_SegmentPlane( knotpoint SegmentStart, knotpoint SegmentEnd, knotpoint PlaneSegmentStart, knotpoint PlaneSegmentEnd, double& IntersectionFraction, std::vector<double>& IntersectionPoint )
{
	double ux = SegmentEnd.xcoord - SegmentStart.xcoord ;
	double uy = SegmentEnd.ycoord - SegmentStart.ycoord ;
	double uz = SegmentEnd.zcoord - SegmentStart.zcoord ;

	double wx= SegmentStart.xcoord - PlaneSegmentStart.xcoord ;
	double wy = SegmentStart.ycoord - PlaneSegmentStart.ycoord ;
	double wz = SegmentStart.zcoord - PlaneSegmentStart.zcoord ;

	double nx= PlaneSegmentEnd.xcoord  - PlaneSegmentStart.xcoord ;
	double ny = PlaneSegmentEnd.ycoord  - PlaneSegmentStart.ycoord ;
	double nz = PlaneSegmentEnd.zcoord  - PlaneSegmentStart.zcoord ;

	double D = nx*ux+ ny*uy + nz*uz;
	double N = - (nx*wx+ ny*wy + nz*wz);

	if (fabs(D) < 0.01)
	{           // segment is parallel to plane
		if (N == 0)                      // segment lies in plane
			return 2;
		else
			return 0;                    // no intersection
	}

	double sI = N / D;
	if (sI < 0 || sI > 1)
		return 0;                        // no intersection


	IntersectionFraction = sI;
	IntersectionPoint[0] = SegmentStart.xcoord + sI * ux;
	IntersectionPoint[1] = SegmentStart.ycoord + sI * uy;
	IntersectionPoint[2] = SegmentStart.zcoord + sI * uz;
	return 1;
}
// intersect3D_SegmentPlane(): find the 3D intersection of a segment and a plane
//    Input:  S = a segment, and Pn = a plane = {Point V0;  Vector n;}
//    Output: *I0 = the intersect point (when it exists)
//    Return: 0 = disjoint (no intersection)
//            1 =  intersection in the unique point *I0
//            2 = the  segment lies in the plane
//~ int
//~ intersect3D_SegmentPlane( Segment S, Plane Pn, Point* I )
//~ {
//~ Vector    u = S.P1 - S.P0;
//~ Vector    w = S.P0 - Pn.V0;
//~ float     D = dot(Pn.n, u);
//~ float     N = -dot(Pn.n, w);
//~ if (fabs(D) < SMALL_NUM) {           // segment is parallel to plane
//~ if (N == 0)                      // segment lies in plane
//~ return 2;
//~ else
//~ return 0;                    // no intersection
//~ }
//~ // they are not parallel
//~ // compute intersect param
//~ float sI = N / D;
//~ if (sI < 0 || sI > 1)
//~ return 0;                        // no intersection
//~ *I = S.P0 + sI * u;                  // compute segment intersect point
//~ return 1;
//~ }
/*********Naive FT doesn't work**********/
/*double *kxp, *kyp, *kzp, *kxpi, *kypi, *kzpi;
  int kmax;

  if(0.5*lambda/M_PI > totlength/NP)
  {
  kmax = ((int) (totlength*2*M_PI/lambda));  //Ignore oscillations of higher freq than 2pi/core diameter
  cout << "kmax: " << kmax << " NP: " << NP << '\n';
  kxp = new double [NP];
  kyp = new double [NP];
  kzp = new double [NP];
  kxpi = new double [NP];
  kypi = new double [NP];
  kzpi = new double [NP];
  for(k=0; k<kmax; k++)
  {
  kxp[k] = 0;
  kxpi[k] = 0;
  kyp[k] = 0;
  kypi[k] = 0;
  kzp[k] = 0;
  kzpi[k] = 0;
  for(s=0; s<NP; s++)
  {
  kxp[k] += knotcurve[s].xcoord*cos(2*M_PI*s*k/NP);
  kxpi[k] += knotcurve[s].xcoord*sin(2*M_PI*s*k/NP);
  kyp[k] += knotcurve[s].ycoord*cos(2*M_PI*s*k/NP);
  kypi[k] += knotcurve[s].ycoord*sin(2*M_PI*s*k/NP);
  kzp[k] += knotcurve[s].zcoord*cos(2*M_PI*s*k/NP);
  kzpi[k] += knotcurve[s].zcoord*sin(2*M_PI*s*k/NP);
  }
  }

  double px, py, pz;

  for(s=0; s<NP; s++)
  {
  px = 0;
  py = 0;
  pz = 0;
  for(k=0; k<kmax; k++)
  {
  px += (kxp[k]*cos(2*M_PI*s*k/NP) + kxpi[k]*sin(2*M_PI*s*k/NP))/kmax;
  py += (kyp[k]*cos(2*M_PI*s*k/NP) + kypi[k]*sin(2*M_PI*s*k/NP))/kmax;
  pz += (kzp[k]*cos(2*M_PI*s*k/NP) + kzpi[k]*sin(2*M_PI*s*k/NP))/kmax;
  }
  knotcurve[s].xcoord = px;
  knotcurve[s].ycoord = py;
  knotcurve[s].zcoord = pz;
  }

  cout << knotcurve.size() << '\n';

  delete [] kxp;
  delete [] kyp;
  delete [] kzp;
  delete [] kxpi;
  delete [] kypi;
  delete [] kzpi;
  }/*

/********Fourier transform******/   //may not be necessary
//~ /********Fourier transform******/
//~ NP = knotcurve.size();  //update number of points in knot curve
//~ vector<double> coord(NP);
//~ gsl_fft_real_wavetable * real;
//~ gsl_fft_halfcomplex_wavetable * hc;
//~ gsl_fft_real_workspace * work;
//~ work = gsl_fft_real_workspace_alloc (NP);
//~ real = gsl_fft_real_wavetable_alloc (NP);
//~ hc = gsl_fft_halfcomplex_wavetable_alloc (NP);
//~ for(int j=1; j<4; j++)
//~ {
//~ switch(j)
//~ {
//~ case 1 :
//~ for(int i=0; i<NP; i++) coord[i] =  knotcurve[i].xcoord ; break;
//~ case 2 :
//~ for(int i=0; i<NP; i++) coord[i] =  knotcurve[i].ycoord ; break;
//~ case 3 :
//~ for(int i=0; i<NP; i++) coord[i] =  knotcurve[i].zcoord ; break;
//~ }
//~ double* data = coord.data();
//~ // take the fft
//~ gsl_fft_real_transform (data, 1, NP, real, work);
//~ // 21/11/2016: make our low pass filter. To apply our filter. we should sample frequencies fn = n/Delta N , n = -N/2 ... N/2
//~ // this is discretizing the nyquist interval, with extreme frequency ~1/2Delta.
//~ // to cut out the frequencies of grid fluctuation size and larger we need a lengthscale Delta to
//~ // plug in above. im doing a rough length calc below, this might be overkill.
//~ // at the moment its just a hard filter, we can choose others though.
//~ // compute a rough length to set scale
//~ int roughlength =0;
//~ for (int  s =1 ; s<NP ; s++) roughlength += sqrt((knotcurve[s].xcoord - knotcurve[s-1].xcoord)*(knotcurve[s].xcoord - knotcurve[s-1].xcoord)+ (knotcurve[s].ycoord - knotcurve[s-1].ycoord)*(knotcurve[s].ycoord - knotcurve[s-1].ycoord) + (knotcurve[s].zcoord - knotcurve[s-1].zcoord)*(knotcurve[s].zcoord - knotcurve[s-1].zcoord));
//~ double filter;
//~ for (int i = 0; i < NP; ++i)
//~ {
//~ if (i < roughlength/h)
//~ {
//~ // passband
//~ filter = 1.0;
//~ }
//~ else
//~ {
//~ // stopband
//~ filter = 0.0;
//~ }
//~ data[i] *= filter;
//~ };
//~ // transform back
//~ gsl_fft_halfcomplex_inverse (data, 1, NP, hc, work);
//~ switch(j)
//~ {
//~ case 1 :
//~ for(int i=0; i<NP; i++)  knotcurve[i].xcoord = coord[i] ; break;
//~ case 2 :
//~ for(int i=0; i<NP; i++)  knotcurve[i].ycoord = coord[i] ; break;
//~ case 3 :
//~ for(int i=0; i<NP; i++)  knotcurve[i].zcoord = coord[i] ; break;
//~ }
//~ }
//~ gsl_fft_real_wavetable_free (real);
//~ gsl_fft_halfcomplex_wavetable_free (hc);
//~ gsl_fft_real_workspace_free (work);
/*******************************/
/* compute velocity vector, spin rate */
// Following Winfree (1990) review, we take the current arc, and look for where it punctures the local normal of the previous arc - this connects the two slightly
// displaced curves together, and tells us how to compute spin rate, and velocity of the filament with finite differences between these pairs.
// Practically, we have two vectors of possibly different lengths, representing an ordered bunch of segments forming the filaments.
// We know the curves are very similar - by aligning them initially (ie finding the same "start point" on each curve), then  travelling fractions of total arc length along them,
// we know we will be in roughly the same on each - thus we only need test a few segments for the desired interesection.
/*
   static bool first = true;

   if (!first)
   {

   int NPold = knotcurveold.size();

// align the two curves. minlocation will give the offset on the new curve.
double minlength =  (knotcurve[0].xcoord - knotcurveold[0].xcoord)*(knotcurve[0].xcoord - knotcurveold[0].xcoord) + (knotcurve[0].ycoord - knotcurveold[0].ycoord) * (knotcurve[0].ycoord - knotcurveold[0].ycoord) +  (knotcurve[0].zcoord - knotcurveold[0].zcoord) * (knotcurve[0].zcoord - knotcurveold[0].zcoord);
double templength = -1;
int offset = 0;
for(s = 1; s< NP; s++)
{
templength = (knotcurve[s].xcoord - knotcurveold[0].xcoord)*(knotcurve[s].xcoord - knotcurveold[0].xcoord) + (knotcurve[s].ycoord - knotcurveold[0].ycoord) * (knotcurve[s].ycoord - knotcurveold[0].ycoord) +  (knotcurve[s].zcoord - knotcurveold[0].zcoord) * (knotcurve[s].zcoord - knotcurveold[0].zcoord);
if (templength < minlength)
{
minlength = templength;
offset = s;
}
}

bool intersection = false;
double IntersectionFraction =-1;
std::vector<double> IntersectionPoint(3);
for(s = 0; s< knotcurveold.size(); s++)
{
intersection = false;
m = s + offset;
int stepnum = 0;
while(!intersection)
{
intersection = intersect3D_SegmentPlane( knotcurve[m%NP], knotcurve[(m+1)%NP], knotcurveold[s%NPold], knotcurveold[(s+1)%NPold], IntersectionFraction, IntersectionPoint );
if(intersection) break;
stepnum++;
stepnum%2? m = incp(m,-stepnum, NP): m = incp(m,stepnum, NP); // work outwards from our best guess

}
// linear interpolation of twist rate
double axinterpolated = knotcurve[(m+1)%NP].ax*IntersectionFraction + knotcurve[m%NP].ax*(1-IntersectionFraction);
double ayinterpolated = knotcurve[(m+1)%NP].ay*IntersectionFraction + knotcurve[m%NP].ay*(1-IntersectionFraction);
double azinterpolated = knotcurve[(m+1)%NP].az*IntersectionFraction + knotcurve[m%NP].az*(1-IntersectionFraction);

// remove tangential part of a

double ax = (axinterpolated - knotcurveold[s].ax);
double ay = (axinterpolated - knotcurveold[s].ay);
double az = (axinterpolated - knotcurveold[s].az);

double nx  =  knotcurveold[(s+1)%NPold].xcoord - knotcurveold[s%NPold].xcoord;
double ny  =  knotcurveold[(s+1)%NPold].ycoord - knotcurveold[s%NPold].ycoord;
double nz  =  knotcurveold[(s+1)%NPold].zcoord - knotcurveold[s%NPold].zcoord;

double proj = (ax*nx+ay*ny+az*nz)/(nx*nx+ny*ny+nz*nz);

ax = ax - proj*nx;
ay = ay - proj*ny;
az = ax - proj*nz;

norm = sqrt(ax*ax+ay*ay+az*az);

ax = ax/norm;
ay = ay/norm;
az = az/norm;

// work out velocity and twist rate
knotcurveold[s].vx = (IntersectionPoint[0] - knotcurveold[s].xcoord )/ dtime;
knotcurveold[s].vy = (IntersectionPoint[1] - knotcurveold[s].ycoord )/ dtime;
knotcurveold[s].vz = (IntersectionPoint[2] - knotcurveold[s].zcoord )/ dtime;

knotcurveold[s].spinrate = ((ax - knotcurveold[s].ax)/dtime)*((ax - knotcurveold[s].ax)/dtime) + ((ay - knotcurveold[s].ay)/dtime)*((ay - knotcurveold[s].ay)/dtime) + ((az - knotcurveold[s].az)/dtime)*((az - knotcurveold[s].az)/dtime);

}
}
first = false;
* */

double my_f(const double s, void* params)
{

	int i,j,k,idwn,jdwn,kdwn,m,pts,iinc,jinc,kinc;
	double ucvxs, ucvys, ucvzs,  xd, yd ,zd, xdiff, ydiff, zdiff, prefactor;
	struct parameters* myparameters = (struct parameters *) params;
	double* x= myparameters->x;
	double* y= myparameters->y;
	double* z= myparameters->z;
	double* ucvx= myparameters->ucvx;
	double* ucvy= myparameters->ucvy;
	double* ucvz= myparameters->ucvz;
	
	gsl_vector* tempf = gsl_vector_alloc (3);
	gsl_vector* tempv = gsl_vector_alloc (3);
	gsl_vector_memcpy (tempf,myparameters->f);
	gsl_vector_memcpy (tempv,myparameters->v);

	// s gives us how much of f to add to p
	gsl_vector_scale(tempf,s);
	gsl_vector_add(tempv,tempf);
	double px = gsl_vector_get(tempv, 0);
	double py = gsl_vector_get(tempv, 1);
	double pz = gsl_vector_get(tempv, 2);
	
	gsl_vector_free(tempf);
	gsl_vector_free(tempv);
	
	/**Find nearest gridpoint**/
	idwn = (int) ((px/h) - 0.5 + Nx/2.0);
	jdwn = (int) ((py/h) - 0.5 + Ny/2.0);
	kdwn = (int) ((pz/h) - 0.5 + Nz/2.0);
	pts=0;
	ucvxs=0;
	ucvys=0;
	ucvzs=0;
	/*curve to gridpoint down distance*/
	xd = (px - x[idwn])/h;
	yd = (py - y[jdwn])/h;
	zd = (pz - z[kdwn])/h;
	for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
	{
		/* Work out increments*/
		iinc = m%2;
		jinc = (m/2)%2;
		kinc = (m/4)%2;
		/*Loop over nearest points*/
		i = incw(idwn, iinc, Nx);
		j = incw(jdwn, jinc, Ny);
		if(periodic) k = incp(kdwn,kinc, Nz);
		else k = incw(kdwn,kinc, Nz);
		prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
		/*interpolate grad u x grad v over nearest points*/
		ucvxs += prefactor*ucvx[pt(i,j,k)];
		ucvys += prefactor*ucvy[pt(i,j,k)];
		ucvzs += prefactor*ucvz[pt(i,j,k)];
	}
	double ans = -1*sqrt(ucvxs*ucvxs + ucvys*ucvys + ucvzs*ucvzs);
	return  ans;
}
