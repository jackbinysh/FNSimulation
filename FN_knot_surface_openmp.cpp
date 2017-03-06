/**************************************************************************************/ /* Fitzhugh-Nagumo reaction diffusion simulation with arbitrary vortex lines OPENMP VERSION Created by Carl Whitfield Last modified 03/01/17 Operational order of the code: 1) The code takes as input an stl file (defined in knot_filename) which defines an orientable surface with a boundary.  2) This surface is scaled to fill a box of size xmax x ymax x zmax.
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
#include <math.h>
#include <string.h>
#define RK4 1         //1 to use Runge Kutta 4th order method, 0 for euler forward method
#define PRESERVE_RATIOS 1  //1 to scale input file preserving the aspect ratio
//includes for the signal processing
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

/*Available options:
FROM_PHI_FILE: Skip initialisation, input from previous run.
FROM_SURFACE_FILE: Initialise from input file(s) generated in surface evolver.
FROM_UV_FILE: Skip initialisation, run FN dynamics from uv file
FROM_KNOT_FILE: Initialise from parametric knot curve in .txt format (e.g. knotplot output)
FROM_FUNCTION: Initialise from some function which can be implemented by the user in phi_calc_manual. eg using theta(x) = artcan(y-y0/x-x0) to give a pole at x0,y0 etc..:wq
 */

const int option = FROM_SURFACE_FILE;         //unknot default option
const BoundaryType BoundaryType=ALLPERIODIC;

/** two rotation angles for the initial stl file, and a displacement vector for the file **/
const double initialthetarotation = 0;
const double initialphirotation = 0;
const double initialxdisplacement = 0;
const double initialydisplacement = 0;
const double initialzdisplacement = 0;

/**If FROM_SURFACE_FILE or FROM_KNOT_FILE chosen**/
string knot_filename = "five2";      //if FROM_SURFACE_FILE assumed input filename format of "XXXXX.stl"
int ncomp = 1;                       //if FROM_KNOT_FILE assumed input filename format of "XXXXX.txt"
//if ncomp > 1 (no. of components) then component files should be separated to 'XXXXX.txt" "XXXXX2.txt", ....
/**IF FROM_PHI_FILE or FROM_UV_FILE chosen**/
string B_filename = "uv_plot60.vtk";    //filename for phi field or uv field

//Grid points
const int Nx = 134;   //No. points in x,y and z
const int Ny = 134;
const int Nz = 134;
const double TTime = 700;       //total time of simulation (simulation units)
const double skiptime = 1;       //print out every # unit of time (simulation units)
const double starttime =0;        //Time at start of simulation (non-zero if continuing from UV file)
const double dtime = 0.02;         //size of each time step

//System size parameters
const double lambda = 21.3;                //approx wavelength
const double size = 4*lambda;   //box size
const double h = size/(Nx-1);            //grid spacing
const double oneoverhsq = 1.0/(h*h);
const double epsilon = 0.3;                //parameters for F-N eqns
const double oneoverepsilon = 1.0/epsilon;
const double beta = 0.7;
const double gam = 0.5;

//Size boundaries of knot (now autoscaled)
double xmax = 6*Nx*h/8.0;
double ymax = 6*Ny*h/8.0;
double zmax = 6*Nz*h/8.0;
int NK;   //number of surface points

//Unallocated matrices
vector<triangle> knotsurface;    //structure for storing knot surface coordinates
vector<knotcurve > knotcurves; // a structure containing some number of knot curves, each curve a list of knotpoints
vector<double> X, Y, Z, dlx, dly,dlz;

double area;   //initial knot area
inline  int pt( int i,  int j,  int k)       //convert i,j,k to single index
{
    return (i*Ny*Nz+j*Nz+k);
}

int main (void)
{
    double  *phi, *u, *v, *ucvx, *ucvy, *ucvz;
    int i,j,k,n,l;
    int *missed;

    phi = new double [Nx*Ny*Nz];  //scalar potential
    u = new double [Nx*Ny*Nz];
    v = new double [Nx*Ny*Nz];
    if(option==FROM_KNOT_FILE)
    {
        missed = new int [Nx*Ny*Nz];
        fill(missed,missed+Nx*Ny*Nz,0);
    }

    // GSL initialization
    const gsl_multimin_fminimizer_type *Type;
    gsl_multimin_fminimizer *minimizerstate;
    Type = gsl_multimin_fminimizer_nmsimplex2;
    minimizerstate = gsl_multimin_fminimizer_alloc (Type,2);

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
            if(option == FROM_FUNCTION)
            { 
                phi_calc_manual(phi);
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
                initial_cond(phi,missed);
            }
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
    double  *ku, *kv; 
    ku = new double [4*Nx*Ny*Nz];
    kv = new double [4*Nx*Ny*Nz];

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


    // let the wave system establish. we only want to do this if the system isnt already established! so don't bother if we are reading in from a uv file 
   /* if (option != FROM_UV_FILE)*/ //burnin(u,v,ucvx, ucvy,ucvz,ku, kv,minimizerstate);

    // okay, into the guts of the simulation
    #if RK4
    #pragma omp parallel default(none) shared (  u, v,  n,ku,kv,p,q,ucvx, ucvy, ucvz,cout, rawtime, timeinfo, knotcurves,   minimizerstate)
    #else
    #pragma omp parallel default(none) shared (  u, v, n, D2u, p, q,ucvx, ucvy, ucvz, cout, rawtime, timeinfo, knotcurves,  minimizerstate)
    #endif
    {
        while(n*dtime <= TTime)
        {
            #pragma omp single
            {
                if(n*dtime >= q)  //Do this every unit T
                {
                    crossgrad_calc(u,v,ucvx,ucvy,ucvz); //find Grad u cross Grad v
                    cout << "T = " << n*dtime + starttime << endl;
                    time (&rawtime);
                    timeinfo = localtime (&rawtime);
                    cout << "current time \t" << asctime(timeinfo) << "\n";
                    if(n*dtime+starttime>=10 )
                    {
                        find_knot_properties(ucvx,ucvy,ucvz,u,n*dtime+starttime,minimizerstate );      //find knot curve and twist and writhe
                    }
                    q++;
                }

                if(n*dtime >= p*skiptime)
                {

                    print_uv(u,v,ucvx,ucvy,ucvz,n*dtime+starttime);
                    p++;
                }

                n++;
            }
            #if RK4
            uv_update(u,v,ku,kv);
            #else
            uv_update_euler(u,v,D2u);
            #endif
        }
    }
    time_t now = time(NULL);
    cout << "Time taken to complete uv part: " << now - then << " seconds.\n";

#if RK4
    delete [] ku;
    delete [] kv;
#else
    delete [] D2u;
#endif
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
    double midpoint[3];
    double norm;
    scalefunction(scale,midpoint,maxxin,minxin,maxyin,minyin,maxzin,minzin);

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

        // apply any rotations and displacements  of the initial coniditions the user has specified
        for(j=0;j<3;j++) rotatedisplace(knotsurface[i].xvertex[j],knotsurface[i].yvertex[j],knotsurface[i].zvertex[j],initialthetarotation,initialphirotation,initialxdisplacement,initialydisplacement,initialzdisplacement);
        rotatedisplace(knotsurface[i].normal[0],knotsurface[i].normal[1],knotsurface[i].normal[2],initialthetarotation,initialphirotation,initialxdisplacement,initialydisplacement,initialzdisplacement);
        rotatedisplace(knotsurface[i].centre[0],knotsurface[i].centre[1],knotsurface[i].centre[2],initialthetarotation,initialphirotation,initialxdisplacement,initialydisplacement,initialzdisplacement);
    }

    cout << "Input scaled by: " << scale[0] << ' ' << scale[1] << ' ' << scale[2] << " in x,y and z\n";

    return A;
}

double init_from_knot_file(void)
{
    double xt,yt,zt,dx,dy,dz,dl,lseg,Lh;    //temporary variables
    vector<double> px,py,pz,dr,ntx,nty,ntz;  //points, distances and tangents
    int npts;  //counter
    string temp;
    double L=0;
    int i,n,m,t,s,NKh;
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
        double midpoint[3];
        scalefunction(scale,midpoint,maxxin,minxin,maxyin,minyin,maxzin,minzin); //function for calculating scale and midpoint

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
void scalefunction(double *scale, double *midpoint, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin)
{
    int i;
    bool nonzeroheight[3];  //marker: true if this dimension has non zero height in stl file
    if(maxxin-minxin>0) { scale[0] = xmax/(maxxin-minxin); nonzeroheight[0] = true; }
    else { scale[0] = 1;  nonzeroheight[0] = false; }
    if(maxyin-minyin>0) { scale[1] = ymax/(maxyin-minyin); nonzeroheight[1] = true; }
    else { scale[1] = 1;  nonzeroheight[1] = false; }
    if(maxzin-minzin>0) { scale[2] = zmax/(maxzin-minzin); nonzeroheight[2] = true; }
    else { scale[2] = 1;  nonzeroheight[2] = false; }
    //double p1x,p1y,p1z,p2x,p2y,p2z,nx,ny,nz;
    midpoint[0] = 0.5*(maxxin+minxin);
    midpoint[1] = 0.5*(maxyin+minyin);
    midpoint[2] = 0.5*(maxzin+minzin);
#if PRESERVE_RATIOS
    double minscale=1000000000;
    int imin=3;
    for(i = 0;i<3;i++)   //find minimum scale factor
    {
        if(scale[i] < minscale && nonzeroheight[i])
        {
            imin = i;
            minscale = scale[i];
        }
    }
    if(imin < 3)      //scale x,y, and z directions by same scale factor
    {
        for(i = 0;i<3;i++) scale[i] = scale[imin];
    }
#endif
}

/*************************Functions for B and Phi calcs*****************************/

void initial_cond(double *phi, int* missed)
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
        B_field_calc(Bx, By, Bz, Bmag, ignore, ignore1, missed);
        time_t now = time(NULL);
        cout << "B field calc took " << now - then << " seconds.\n";
        cout << "Calculating scalar potential...\n";
        then = time(NULL);
        phi_calc_B(Bx, By, Bz, Bmag, ignore, ignore1, missed, phi);
        now = time(NULL);
        cout << "Phi field calc took " << now - then << " seconds.\n";
        cout << "Printing B and phi...\n";
        print_B_phi(phi, missed);

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
        phi_calc(phi);
        time_t now = time(NULL);
        cout << "Initialisation took " << now - then << " seconds.\n";
        cout << "Printing B and phi...\n";
        print_B_phi(phi, missed);
    }
}

void phi_calc(double *phi)
{
    int i,j,k,n,s;
    double rx,ry,rz,r;


#pragma omp parallel default(none) shared ( knotsurface, phi, NK ) private ( i, j, k, n, s, rx, ry, rz , r)
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
                        rx = knotsurface[s].centre[0]-x(i);
                        ry = knotsurface[s].centre[1]-y(j);
                        rz = knotsurface[s].centre[2]-z(k);
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
void phi_calc_manual(double *phi)
{
    int i,j,k,n;
    for(i=0;i<Nx;i++)
    {
        for(j=0; j<Ny; j++)
        {
            for(k=0; k<Nz; k++)
            {
                n = pt(i,j,k);
                phi[n] = 0;
                double theta = 0.5;
                phi[n] = 0 ;
                while(phi[n]>M_PI) phi[n] -= 2*M_PI;
                while(phi[n]<-M_PI) phi[n] += 2*M_PI;
            }
        }
    }
}
void B_field_calc( double *Bx, double *By, double *Bz, double *Bmag, int *ignore, int *ignore1, int *missed)
{
    int i,j,k,n,t;
    double lx,ly,lz,lmag;
    double coresize = lambda/(2*M_PI);

#pragma omp parallel default(none) shared (  X, Y, Z, Bx, By, Bz, missed, Bmag, NK, coresize, ignore, ignore1, dlx, dly, dlz ) private ( i, j, k, n, t, lx, ly ,lz, lmag)
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
                        lx = x(i)-X[t];    //distance to point on line
                        ly = y(j)-Y[t];
                        lz = z(k)-Z[t];
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

void crossgrad_calc( double *u, double *v, double *ucvx, double *ucvy, double *ucvz)
{
    int i,j,k,n,kup,kdown;
    double dxu,dyu,dzu,dxv,dyv,dzv;
    for(i=0;i<Nx;i++)
    {
        for(j=0; j<Ny; j++)
        {
            for(k=0; k<Nz; k++)   //Central difference
            {
                kup = gridinc(k,1,Nz,2);
                kdown = gridinc(k,-1,Nz,2);
                dxu = 0.5*(u[pt(gridinc(i,1,Nx,0),j,k)]-u[pt(gridinc(i,-1,Nx,0),j,k)])/h;
                dxv = 0.5*(v[pt(gridinc(i,1,Nx,0),j,k)]-v[pt(gridinc(i,-1,Nx,0),j,k)])/h;
                dyu = 0.5*(u[pt(i,gridinc(j,1,Ny,1),k)]-u[pt(i,gridinc(j,-1,Ny,1),k)])/h;
                dyv = 0.5*(v[pt(i,gridinc(j,1,Ny,1),k)]-v[pt(i,gridinc(j,-1,Ny,1),k)])/h;
                dzu = 0.5*(u[pt(i,j,kup)]-u[pt(i,j,kdown)])/h;
                dzv = 0.5*(v[pt(i,j,kup)]-v[pt(i,j,kdown)])/h;
                n = pt(i,j,k);
                ucvx[n] = dyu*dzv - dzu*dyv;
                ucvy[n] = dzu*dxv - dxu*dzv;    //Grad u cross Grad v
                ucvz[n] = dxu*dyv - dyu*dxv;
            }
        }
    }
}
void burnin(double* u,double* v,double *ucvx,double* ucvy,double* ucvz,double* ku, double* kv,gsl_multimin_fminimizer* minimizerstate)
{
    int n = 0;
    int p = 0;
    bool markedpoints[Nx*Ny*Nz] = {false};
    #pragma omp parallel default(none) shared (X,Y,Z,u, v,p,n,ku,kv,ucvx, ucvy, ucvz, minimizerstate,markedpoints)
    {
        // evolve for timesteps to smooth things out
        //while(n*dtime <= 1)
        //{
        //    #if RK4
        //    uv_update(u,v,ku,kv);
        //    #else
        //    uv_update_euler(u,v,D2u);
        //    #endif
        //    #pragma omp single
        //    {
        //        n++;
        //    } 
        //}
        // get the region to exclude
        #pragma omp single
        {
            double status = init_from_knot_file();
            crossgrad_calc(u,v,ucvx,ucvy,ucvz); //find Grad u cross Grad v
     //       find_knot_properties_burnin(markedpoints,ucvx,ucvy,ucvz,minimizerstate);
        //    n = 0;
// for 1: length X, read off point, find nearest gridpoint mark region
            for(int c = 0; c<X.size() ;c++)
            {
                int idwn = (int) ((X[c]/h) - 0.5 + Nx/2.0);
                int jdwn = (int) ((Y[c]/h) - 0.5 + Ny/2.0);
                int kdwn = (int) ((Z[c]/h) - 0.5 + Nz/2.0);
                // idwn etc can be off the actual grid , into "ghost" grids around the real one. this is useful for knotcurve tracing over periodic boundaries
                // but we also need the corresponding real grid positions!
                  int modidwn = circularmod(idwn,Nx);
                int modjdwn = circularmod(jdwn,Ny);
                  int modkdwn = circularmod(kdwn,Nz);
               if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
               if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
                // mark these points , up to roughly a core radius in all directions, in the "marked" array
               // int delta = ceil((lambda/(8*M_PI))/h);
                int delta = 1;
                for(int i=gridinc(modidwn,-delta,Nx,0);i<=gridinc(modidwn,delta,Nx,0);i++)
                {
                    for(int j=gridinc(modjdwn,-delta,Ny,1); j<=gridinc(modjdwn,delta,Ny,1) ;j++)
                    {
                        for(int k=gridinc(modkdwn,-delta,Nz,2); k<=gridinc(modkdwn,delta,Nz,2); k++)   //Central difference
                        {
                            markedpoints[pt(i,j,k)] = true;
                            u[pt(i,j,k)] = -1.29;
                            v[pt(i,j,k)]  = -0.61 ;
                        }
                    }
                }
            }
        }
            print_uv(u,v,ucvx,ucvy,ucvz,-2,markedpoints);
        // okay we have a region to exclude. now lets run the dynamics , excluding these points
        while(n*dtime <= 40)
        {
            uv_update_burnin(u,v,ku,kv,markedpoints);
            #pragma omp single
            {
                if(n*dtime >= p*skiptime)
                {

                    print_uv(u,v,ucvx,ucvy,ucvz,(n*dtime)+starttime);
                    p++;
                }
                n++;
            } 
        }
    }
}
// a stripped down version of find_knot_properties, used to get a marked region for the burn in
void find_knot_properties_burnin(bool *markedpoints, double *ucvx, double *ucvy, double *ucvz, gsl_multimin_fminimizer* minimizerstate)
{
    int c =0;
    bool first = false;
    bool knotexists = true;
    bool cleanupneeded = false;
    //while(knotexists)
    {
        double  ucvmag, norm;
        double   ucvmax = -1.0; // should always be +ve, so setting it to an initially -ve # means it always gets written to once.
        int n,i,j,k,imax,jmax,kmax;
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k);
                    ucvmag = sqrt(ucvx[n]*ucvx[n] + ucvy[n]*ucvy[n] + ucvz[n]*ucvz[n]);
                    if(!markedpoints[n] &&(ucvmag > ucvmax))
                    {
                        ucvmax = ucvmag;
                        imax = i;
                        jmax = j;
                        kmax=k;
                    }
                }
            }
        }
        if(ucvmax<0.7) knotexists = false;
        else
        {
            knotexists = true; 
            cleanupneeded = true;
            // its the first time we enter here we know we are setting up the knot
            static bool firsttimehere = true;
            if(firsttimehere) first = true; firsttimehere = false;
        }
        if(knotexists)
        {
            knotcurves.push_back(knotcurve() );
            knotcurves[c].knotcurve.push_back(knotpoint());
            knotcurves[c].knotcurve[0].xcoord=x(imax);
            knotcurves[c].knotcurve[0].ycoord=y(jmax);
            knotcurves[c].knotcurve[0].zcoord=z(kmax);

            int s=1;
            bool finish=false;

            int idwn,jdwn,kdwn, modidwn, modjdwn, modkdwn,m,pts,iinc,jinc,kinc,attempts;
            double ucvxs, ucvys, ucvzs, graducvx, graducvy, graducvz, prefactor, xd, yd ,zd, fx, fy, fz, xdiff, ydiff, zdiff;
            /*calculate local direction of grad u x grad v (the tangent to the knot curve) at point s-1, then move to point s by moving along tangent + unit confinement force*/
            while (finish==false)
            {
                norm=0;
                /**Find nearest gridpoint**/
                idwn = (int) ((knotcurves[c].knotcurve[s-1].xcoord/h) - 0.5 + Nx/2.0);
                jdwn = (int) ((knotcurves[c].knotcurve[s-1].ycoord/h) - 0.5 + Ny/2.0);
                kdwn = (int) ((knotcurves[c].knotcurve[s-1].zcoord/h) - 0.5 + Nz/2.0);
                // idwn etc can be off the actual grid , into "ghost" grids around the real one. this is useful for knotcurve tracing over periodic boundaries
                // but we also need the corresponding real grid positions!
                modidwn = circularmod(idwn,Nx);
                modjdwn = circularmod(jdwn,Ny);
                modkdwn = circularmod(kdwn,Nz);
               if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
               if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
                // mark these points , up to roughly a core radius in all directions, in the "marked" array
               // int delta = ceil((lambda/(8*M_PI))/h);
                int delta = 0;
                for(i=gridinc(modidwn,-delta,Nx,0);i<=gridinc(modidwn,delta,Nx,0);i++)
                {
                    for(j=gridinc(modjdwn,-delta,Ny,1); j<=gridinc(modjdwn,delta,Ny,1) ;j++)
                    {
                        for(k=gridinc(modkdwn,-delta,Nz,2); k<=gridinc(modkdwn,delta,Nz,2); k++)   //Central difference
                        {
                            n = pt(i,j,k);
                            markedpoints[n] = true;
                        }
                    }
                }
                pts=0;
                ucvxs=0;
                ucvys=0;
                ucvzs=0;
                /*curve to gridpoint down distance*/
                xd = (knotcurves[c].knotcurve[s-1].xcoord - x(idwn))/h;
                yd = (knotcurves[c].knotcurve[s-1].ycoord - y(jdwn))/h;
                zd = (knotcurves[c].knotcurve[s-1].zcoord - z(kdwn))/h;
                for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
                {
                    /* Work out increments*/
                    iinc = m%2;
                    jinc = (m/2)%2;
                    kinc = (m/4)%2;
                    /*Loop over nearest points*/
                    i = gridinc(modidwn, iinc, Nx,0);
                    j = gridinc(modjdwn, jinc, Ny,1);
                    k = gridinc(modkdwn,kinc, Nz,2);
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
                double testx = knotcurves[c].knotcurve[s-1].xcoord + 0.1*ucvxs*lambda/(2*M_PI);
                double testy = knotcurves[c].knotcurve[s-1].ycoord + 0.1*ucvys*lambda/(2*M_PI);
                double testz = knotcurves[c].knotcurve[s-1].zcoord + 0.1*ucvzs*lambda/(2*M_PI);

                // now get the grad at this point
                idwn = (int) ((testx/h) - 0.5 + Nx/2.0);
                jdwn = (int) ((testy/h) - 0.5 + Ny/2.0);
                kdwn = (int) ((testz/h) - 0.5 + Nz/2.0);
                modidwn = circularmod(idwn,Nx);
                modjdwn = circularmod(jdwn,Ny);
                modkdwn = circularmod(kdwn,Nz);
                // again, bear in mind these numbers can be into the "ghost" grids
               if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
               if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
                pts=0;
                graducvx=0;
                graducvy=0;
                graducvz=0;
                /*curve to gridpoint down distance*/
                xd = (testx - x(idwn))/h;
                yd = (testy - y(jdwn))/h;
                zd = (testz - z(kdwn))/h;
                for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
                {
                    /* Work out increments*/
                    iinc = m%2;
                    jinc = (m/2)%2;
                    kinc = (m/4)%2;
                    /*Loop over nearest points*/
                    i = gridinc(modidwn, iinc, Nx,0);
                    j = gridinc(modjdwn, jinc, Ny,1);
                    k = gridinc(modkdwn,kinc, Nz,2);
                    prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
                    /*interpolate gradients of |grad u x grad v|*/
                    graducvx += prefactor*(sqrt(ucvx[pt(gridinc(i,1,Nx,0),j,k)]*ucvx[pt(gridinc(i,1,Nx,0),j,k)] + ucvy[pt(gridinc(i,1,Nx,0),j,k)]*ucvy[pt(gridinc(i,1,Nx,0),j,k)] + ucvz[pt(gridinc(i,1,Nx,0),j,k)]*ucvz[pt(gridinc(i,1,Nx,0),j,k)]) - sqrt(ucvx[pt(gridinc(i,-1,Nx,0),j,k)]*ucvx[pt(gridinc(i,-1,Nx,0),j,k)] + ucvy[pt(gridinc(i,-1,Nx,0),j,k)]*ucvy[pt(gridinc(i,-1,Nx,0),j,k)] + ucvz[pt(gridinc(i,-1,Nx,0),j,k)]*ucvz[pt(gridinc(i,-1,Nx,0),j,k)]))/(2*h);
                    graducvy += prefactor*(sqrt(ucvx[pt(i,gridinc(j,1,Ny,1),k)]*ucvx[pt(i,gridinc(j,1,Ny,1),k)] + ucvy[pt(i,gridinc(j,1,Ny,1),k)]*ucvy[pt(i,gridinc(j,1,Ny,1),k)] + ucvz[pt(i,gridinc(j,1,Ny,1),k)]*ucvz[pt(i,gridinc(j,1,Ny,1),k)]) - sqrt(ucvx[pt(i,gridinc(j,-1,Ny,1),k)]*ucvx[pt(i,gridinc(j,-1,Ny,1),k)] + ucvy[pt(i,gridinc(j,-1,Ny,1),k)]*ucvy[pt(i,gridinc(j,-1,Ny,1),k)] + ucvz[pt(i,gridinc(j,-1,Ny,1),k)]*ucvz[pt(i,gridinc(j,-1,Ny,1),k)]))/(2*h);
                    graducvz += prefactor*(sqrt(ucvx[pt(i,j,gridinc(k,1,Nz,2))]*ucvx[pt(i,j,gridinc(k,1,Nz,2))] + ucvy[pt(i,j,gridinc(k,1,Nz,2))]*ucvy[pt(i,j,gridinc(k,1,Nz,2))] + ucvz[pt(i,j,gridinc(k,1,Nz,2))]*ucvz[pt(i,j,gridinc(k,1,Nz,2))]) - sqrt(ucvx[pt(i,j,gridinc(k,-1,Nz,2))]*ucvx[pt(i,j,gridinc(k,-1,Nz,2))] + ucvy[pt(i,j,gridinc(k,-1,Nz,2))]*ucvy[pt(i,j,gridinc(k,-1,Nz,2))] + ucvz[pt(i,j,gridinc(k,-1,Nz,2))]*ucvz[pt(i,j,gridinc(k,-1,Nz,2))]))/(2*h);

                }
                knotcurves[c].knotcurve.push_back(knotpoint());
                // one of the vectors in the plane we wish to perfrom our minimisation in
                fx = (graducvx - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvxs); 
                fy = (graducvy - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvys);
                fz = (graducvz - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvzs);
                norm = sqrt(fx*fx + fy*fy + fz*fz);
                fx = fx/norm;
                fy = fy/norm;
                fz = fz/norm;

                // okay we have our direction to perfrom the line minimisation in
                // the point
                gsl_vector* v = gsl_vector_alloc (3);
                gsl_vector_set (v, 0, testx);
                gsl_vector_set (v, 1, testy);
                gsl_vector_set (v, 2, testz);
                // one vector in the plane we with to minimize in
                gsl_vector* f = gsl_vector_alloc (3);
                gsl_vector_set (f, 0, fx);
                gsl_vector_set (f, 1, fy);
                gsl_vector_set (f, 2, fz);
                // the ucv vector
                gsl_vector* ucv = gsl_vector_alloc (3);
                gsl_vector_set (ucv, 0, ucvxs);
                gsl_vector_set (ucv, 1, ucvys);
                gsl_vector_set (ucv, 2, ucvzs);
                // take a cross product to get the other vector in the plane 
                gsl_vector* b = gsl_vector_alloc (3);
                cross_product(f,ucv,b); 
                // initial conditions
                gsl_vector* minimum = gsl_vector_alloc (2);
                gsl_vector_set (minimum, 0, 0);
                gsl_vector_set (minimum, 1, 0);
                struct parameters params; struct parameters* pparams = &params;
                pparams->ucvx=ucvx;pparams->ucvy=ucvy; pparams->ucvz = ucvz;
                pparams->v = v; pparams->f = f;pparams->b=b;
                // some initial values
                gsl_multimin_function F;
                F.n=2;
                F.f = &my_f;
                F.params = (void*) pparams;
                gsl_vector* stepsize = gsl_vector_alloc (2);
                gsl_vector_set (stepsize, 0, lambda/(8*M_PI));
                gsl_vector_set (stepsize, 1, lambda/(8*M_PI));
                gsl_multimin_fminimizer_set (minimizerstate, &F, minimum, stepsize);

                int iter=0;
                int status =0;
                double minimizersize=0;
                do
                {
                    iter++;
                    status = gsl_multimin_fminimizer_iterate(minimizerstate);

                    if (status) 
                        break;

                    minimizersize = gsl_multimin_fminimizer_size (minimizerstate);
                    status = gsl_multimin_test_size (size, 1e-2);

                }
                while (status == GSL_CONTINUE && iter < 500);


                gsl_vector_scale(f,gsl_vector_get(minimizerstate->x, 0));
                gsl_vector_scale(b,gsl_vector_get(minimizerstate->x, 1));
                gsl_vector_add(f,b);
                gsl_vector_add(v,f);
                knotcurves[c].knotcurve[s].xcoord = gsl_vector_get(v, 0);
                knotcurves[c].knotcurve[s].ycoord= gsl_vector_get(v, 1);
                knotcurves[c].knotcurve[s].zcoord= gsl_vector_get(v, 2);

                gsl_vector_free(v);
                gsl_vector_free(f);
                gsl_vector_free(b);
                gsl_vector_free(ucv);
                gsl_vector_free(stepsize);

                xdiff = knotcurves[c].knotcurve[0].xcoord - knotcurves[c].knotcurve[s].xcoord;     //distance from start/end point
                ydiff = knotcurves[c].knotcurve[0].ycoord - knotcurves[c].knotcurve[s].ycoord;
                zdiff = knotcurves[c].knotcurve[0].zcoord - knotcurves[c].knotcurve[s].zcoord;
                if(sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff) <3*h  && s > 10) finish = true;
                if(s>50000) finish = true;
                s++;
            }
        }
        c++;
    }
    std::vector<int> permutation (1,0);
    print_knot(-1, knotcurves, permutation);
    knotcurves.clear(); //empty vector with knot curve points
}
void find_knot_properties( double *ucvx, double *ucvy, double *ucvz, double *u, double t, gsl_multimin_fminimizer* minimizerstate)
{
    int c =0;
    static bool xmarked[Nx] = {false};
    static bool ymarked[Ny] = {false};
    static bool zmarked[Nz] = {false};
    bool first = false;
    bool knotexists = true;
    bool cleanupneeded = false;
    while(knotexists)
    {
        double  ucvmag, norm;
        double   ucvmax = -1.0; // should always be +ve, so setting it to an initially -ve # means it always gets written to once.
        int n,i,j,k,imax,jmax,kmax;
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k);
                    ucvmag = sqrt(ucvx[n]*ucvx[n] + ucvy[n]*ucvy[n] + ucvz[n]*ucvz[n]);
                    if( ( !xmarked[i]||!ymarked[j]||!zmarked[k] )&& (ucvmag > ucvmax))
                    {
                        ucvmax = ucvmag;
                        imax = i;
                        jmax = j;
                        kmax=k;
                    }
                }
            }
        }
        if(ucvmax<0.1) knotexists = false;
        else
        {
            knotexists = true; 
            cleanupneeded = true;
            // its the first time we enter here we know we are setting up the knot
            static bool firsttimehere = true;
            if(firsttimehere) first = true; firsttimehere = false;
        }
        if(knotexists)
        {
            knotcurves.push_back(knotcurve() );
            knotcurves[c].knotcurve.push_back(knotpoint());
            knotcurves[c].knotcurve[0].xcoord=x(imax);
            knotcurves[c].knotcurve[0].ycoord=y(jmax);
            knotcurves[c].knotcurve[0].zcoord=z(kmax);

            int s=1;
            bool finish=false;

            int idwn,jdwn,kdwn, modidwn, modjdwn, modkdwn,m,pts,iinc,jinc,kinc,attempts;
            double ucvxs, ucvys, ucvzs, graducvx, graducvy, graducvz, prefactor, xd, yd ,zd, fx, fy, fz, xdiff, ydiff, zdiff;
            /*calculate local direction of grad u x grad v (the tangent to the knot curve) at point s-1, then move to point s by moving along tangent + unit confinement force*/
            while (finish==false)
            {
                norm=0;
                /**Find nearest gridpoint**/
                idwn = (int) ((knotcurves[c].knotcurve[s-1].xcoord/h) - 0.5 + Nx/2.0);
                jdwn = (int) ((knotcurves[c].knotcurve[s-1].ycoord/h) - 0.5 + Ny/2.0);
                kdwn = (int) ((knotcurves[c].knotcurve[s-1].zcoord/h) - 0.5 + Nz/2.0);
                // idwn etc can be off the actual grid , into "ghost" grids around the real one. this is useful for knotcurve tracing over periodic boundaries
                // but we also need the corresponding real grid positions!
                modidwn = circularmod(idwn,Nx);
                modjdwn = circularmod(jdwn,Ny);
                modkdwn = circularmod(kdwn,Nz);
               if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
               if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
                // mark these points , up to roughly a core radius in all directions, in the "marked" array
                int delta = ceil((lambda/(2*M_PI))/h);
                for(int q = - delta; q <= delta; q++)
                {

                    xmarked[gridinc(modidwn,q,Nx,0)] = true;
                    ymarked[gridinc(modjdwn,q,Ny,1)] = true;
                    zmarked[gridinc(modkdwn,q,Nz,2)] = true;

                }
                pts=0;
                ucvxs=0;
                ucvys=0;
                ucvzs=0;
                /*curve to gridpoint down distance*/
                xd = (knotcurves[c].knotcurve[s-1].xcoord - x(idwn))/h;
                yd = (knotcurves[c].knotcurve[s-1].ycoord - y(jdwn))/h;
                zd = (knotcurves[c].knotcurve[s-1].zcoord - z(kdwn))/h;
                for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
                {
                    /* Work out increments*/
                    iinc = m%2;
                    jinc = (m/2)%2;
                    kinc = (m/4)%2;
                    /*Loop over nearest points*/
                    i = gridinc(modidwn, iinc, Nx,0);
                    j = gridinc(modjdwn, jinc, Ny,1);
                    k = gridinc(modkdwn,kinc, Nz,2);
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
                double testx = knotcurves[c].knotcurve[s-1].xcoord + 0.1*ucvxs*lambda/(2*M_PI);
                double testy = knotcurves[c].knotcurve[s-1].ycoord + 0.1*ucvys*lambda/(2*M_PI);
                double testz = knotcurves[c].knotcurve[s-1].zcoord + 0.1*ucvzs*lambda/(2*M_PI);

                // now get the grad at this point
                idwn = (int) ((testx/h) - 0.5 + Nx/2.0);
                jdwn = (int) ((testy/h) - 0.5 + Ny/2.0);
                kdwn = (int) ((testz/h) - 0.5 + Nz/2.0);
                modidwn = circularmod(idwn,Nx);
                modjdwn = circularmod(jdwn,Ny);
                modkdwn = circularmod(kdwn,Nz);
                // again, bear in mind these numbers can be into the "ghost" grids
               if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
               if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
                pts=0;
                graducvx=0;
                graducvy=0;
                graducvz=0;
                /*curve to gridpoint down distance*/
                xd = (testx - x(idwn))/h;
                yd = (testy - y(jdwn))/h;
                zd = (testz - z(kdwn))/h;
                for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
                {
                    /* Work out increments*/
                    iinc = m%2;
                    jinc = (m/2)%2;
                    kinc = (m/4)%2;
                    /*Loop over nearest points*/
                    i = gridinc(modidwn, iinc, Nx,0);
                    j = gridinc(modjdwn, jinc, Ny,1);
                    k = gridinc(modkdwn,kinc, Nz,2);
                    prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
                    /*interpolate gradients of |grad u x grad v|*/
                    graducvx += prefactor*(sqrt(ucvx[pt(gridinc(i,1,Nx,0),j,k)]*ucvx[pt(gridinc(i,1,Nx,0),j,k)] + ucvy[pt(gridinc(i,1,Nx,0),j,k)]*ucvy[pt(gridinc(i,1,Nx,0),j,k)] + ucvz[pt(gridinc(i,1,Nx,0),j,k)]*ucvz[pt(gridinc(i,1,Nx,0),j,k)]) - sqrt(ucvx[pt(gridinc(i,-1,Nx,0),j,k)]*ucvx[pt(gridinc(i,-1,Nx,0),j,k)] + ucvy[pt(gridinc(i,-1,Nx,0),j,k)]*ucvy[pt(gridinc(i,-1,Nx,0),j,k)] + ucvz[pt(gridinc(i,-1,Nx,0),j,k)]*ucvz[pt(gridinc(i,-1,Nx,0),j,k)]))/(2*h);
                    graducvy += prefactor*(sqrt(ucvx[pt(i,gridinc(j,1,Ny,1),k)]*ucvx[pt(i,gridinc(j,1,Ny,1),k)] + ucvy[pt(i,gridinc(j,1,Ny,1),k)]*ucvy[pt(i,gridinc(j,1,Ny,1),k)] + ucvz[pt(i,gridinc(j,1,Ny,1),k)]*ucvz[pt(i,gridinc(j,1,Ny,1),k)]) - sqrt(ucvx[pt(i,gridinc(j,-1,Ny,1),k)]*ucvx[pt(i,gridinc(j,-1,Ny,1),k)] + ucvy[pt(i,gridinc(j,-1,Ny,1),k)]*ucvy[pt(i,gridinc(j,-1,Ny,1),k)] + ucvz[pt(i,gridinc(j,-1,Ny,1),k)]*ucvz[pt(i,gridinc(j,-1,Ny,1),k)]))/(2*h);
                    graducvz += prefactor*(sqrt(ucvx[pt(i,j,gridinc(k,1,Nz,2))]*ucvx[pt(i,j,gridinc(k,1,Nz,2))] + ucvy[pt(i,j,gridinc(k,1,Nz,2))]*ucvy[pt(i,j,gridinc(k,1,Nz,2))] + ucvz[pt(i,j,gridinc(k,1,Nz,2))]*ucvz[pt(i,j,gridinc(k,1,Nz,2))]) - sqrt(ucvx[pt(i,j,gridinc(k,-1,Nz,2))]*ucvx[pt(i,j,gridinc(k,-1,Nz,2))] + ucvy[pt(i,j,gridinc(k,-1,Nz,2))]*ucvy[pt(i,j,gridinc(k,-1,Nz,2))] + ucvz[pt(i,j,gridinc(k,-1,Nz,2))]*ucvz[pt(i,j,gridinc(k,-1,Nz,2))]))/(2*h);

                }
                knotcurves[c].knotcurve.push_back(knotpoint());
                // one of the vectors in the plane we wish to perfrom our minimisation in
                fx = (graducvx - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvxs); 
                fy = (graducvy - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvys);
                fz = (graducvz - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvzs);
                norm = sqrt(fx*fx + fy*fy + fz*fz);
                fx = fx/norm;
                fy = fy/norm;
                fz = fz/norm;

                // okay we have our direction to perfrom the line minimisation in
                // the point
                gsl_vector* v = gsl_vector_alloc (3);
                gsl_vector_set (v, 0, testx);
                gsl_vector_set (v, 1, testy);
                gsl_vector_set (v, 2, testz);
                // one vector in the plane we with to minimize in
                gsl_vector* f = gsl_vector_alloc (3);
                gsl_vector_set (f, 0, fx);
                gsl_vector_set (f, 1, fy);
                gsl_vector_set (f, 2, fz);
                // the ucv vector
                gsl_vector* ucv = gsl_vector_alloc (3);
                gsl_vector_set (ucv, 0, ucvxs);
                gsl_vector_set (ucv, 1, ucvys);
                gsl_vector_set (ucv, 2, ucvzs);
                // take a cross product to get the other vector in the plane 
                gsl_vector* b = gsl_vector_alloc (3);
                cross_product(f,ucv,b); 
                // initial conditions
                gsl_vector* minimum = gsl_vector_alloc (2);
                gsl_vector_set (minimum, 0, 0);
                gsl_vector_set (minimum, 1, 0);
                struct parameters params; struct parameters* pparams = &params;
                pparams->ucvx=ucvx;pparams->ucvy=ucvy; pparams->ucvz = ucvz;
                pparams->v = v; pparams->f = f;pparams->b=b;
                // some initial values
                gsl_multimin_function F;
                F.n=2;
                F.f = &my_f;
                F.params = (void*) pparams;
                gsl_vector* stepsize = gsl_vector_alloc (2);
                gsl_vector_set (stepsize, 0, lambda/(8*M_PI));
                gsl_vector_set (stepsize, 1, lambda/(8*M_PI));
                gsl_multimin_fminimizer_set (minimizerstate, &F, minimum, stepsize);

                int iter=0;
                int status =0;
                double minimizersize=0;
                do
                {
                    iter++;
                    status = gsl_multimin_fminimizer_iterate(minimizerstate);

                    if (status) 
                        break;

                    minimizersize = gsl_multimin_fminimizer_size (minimizerstate);
                    status = gsl_multimin_test_size (size, 1e-2);

                }
                while (status == GSL_CONTINUE && iter < 500);


                gsl_vector_scale(f,gsl_vector_get(minimizerstate->x, 0));
                gsl_vector_scale(b,gsl_vector_get(minimizerstate->x, 1));
                gsl_vector_add(f,b);
                gsl_vector_add(v,f);
                knotcurves[c].knotcurve[s].xcoord = gsl_vector_get(v, 0);
                knotcurves[c].knotcurve[s].ycoord= gsl_vector_get(v, 1);
                knotcurves[c].knotcurve[s].zcoord= gsl_vector_get(v, 2);

                gsl_vector_free(v);
                gsl_vector_free(f);
                gsl_vector_free(b);
                gsl_vector_free(ucv);
                gsl_vector_free(stepsize);

                xdiff = knotcurves[c].knotcurve[0].xcoord - knotcurves[c].knotcurve[s].xcoord;     //distance from start/end point
                ydiff = knotcurves[c].knotcurve[0].ycoord - knotcurves[c].knotcurve[s].ycoord;
                zdiff = knotcurves[c].knotcurve[0].zcoord - knotcurves[c].knotcurve[s].zcoord;
                if(sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff) <3*h  && s > 10) finish = true;
                if(s>50000) finish = true;
                s++;
            }

            int NP = knotcurves[c].knotcurve.size();  //store number of points in knot curve


            /*******Vertex averaging*********/

            double totlength, dl, dx,dy,dz;
            for(i=0;i<3;i++)   //repeat a couple of times because of end point
            {
                totlength=0;
                for(s=0; s<NP; s++)   //Work out total length of curve
                {
                    dx = knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[s].xcoord;
                    dy = knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[s].ycoord;
                    dz = knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[s].zcoord;
                    totlength += sqrt(dx*dx + dy*dy + dz*dz);
                }
                dl = totlength/NP;
                for(s=0; s<NP; s++)    //Move points to have spacing dl
                {
                    dx = knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[s].xcoord;
                    dy = knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[s].ycoord;
                    dz = knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[s].zcoord;
                    norm = sqrt(dx*dx + dy*dy + dz*dz);
                    knotcurves[c].knotcurve[incp(s,1,NP)].xcoord = knotcurves[c].knotcurve[s].xcoord + dl*dx/norm;
                    knotcurves[c].knotcurve[incp(s,1,NP)].ycoord = knotcurves[c].knotcurve[s].ycoord + dl*dy/norm;
                    knotcurves[c].knotcurve[incp(s,1,NP)].zcoord = knotcurves[c].knotcurve[s].zcoord + dl*dz/norm;
                }
            }

            /*************Curve Smoothing*******************/
            vector<double> coord(NP);
            gsl_fft_real_wavetable * real;
            gsl_fft_halfcomplex_wavetable * hc;
            gsl_fft_real_workspace * work;
            work = gsl_fft_real_workspace_alloc (NP);
            real = gsl_fft_real_wavetable_alloc (NP);
            hc = gsl_fft_halfcomplex_wavetable_alloc (NP);
            for(j=1; j<4; j++)
            {
                switch(j)
                {
                    case 1 :
                        for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].xcoord ; break;
                    case 2 :
                        for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ycoord ; break;
                    case 3 :
                        for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].zcoord ; break;
                }
                double* data = coord.data();
                // take the fft
                gsl_fft_real_transform (data, 1, NP, real, work);
                // 21/11/2016: make our low pass filter. To apply our filter. we should sample frequencies fn = n/Delta N , n = -N/2 ... N/2
                // this is discretizing the nyquist interval, with extreme frequency ~1/2Delta.
                // to cut out the frequencies of grid fluctuation size and larger we need a lengthscale Delta to
                // plug in above. im doing a rough length calc below, this might be overkill.
                // at the moment its just a hard filter, we can choose others though.
                // compute a rough length to set scale
                double filter;
                const double cutoff = 2*M_PI*(totlength/(6*lambda));
                for (i = 0; i < NP; ++i)
                {
                    filter = 1/sqrt(1+pow((i/cutoff),8));
                    data[i] *= filter;
                };
                // transform back
                gsl_fft_halfcomplex_inverse (data, 1, NP, hc, work);
                switch(j)
                {
                    case 1 :
                        for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].xcoord = coord[i] ; break;
                    case 2 :
                        for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ycoord = coord[i] ; break;
                    case 3 :
                        for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].zcoord = coord[i] ; break;
                }
            }

            /******************Interpolate direction of grad u for twist calc*******/
            /**Find nearest gridpoint**/
            double dxu, dyu, dzu, dxup, dyup, dzup;
            for(s=0; s<NP; s++)
            {
                idwn = (int) ((knotcurves[c].knotcurve[s].xcoord/h) - 0.5 + Nx/2.0);
                jdwn = (int) ((knotcurves[c].knotcurve[s].ycoord/h) - 0.5 + Ny/2.0);
                kdwn = (int) ((knotcurves[c].knotcurve[s].zcoord/h) - 0.5 + Nz/2.0);
                modidwn = circularmod(idwn,Nx);
                modjdwn = circularmod(jdwn,Ny);
                modkdwn = circularmod(kdwn,Nz);
               if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
               if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
                dxu=0;
                dyu=0;
                dzu=0;
                /*curve to gridpoint down distance*/
                xd = (knotcurves[c].knotcurve[s].xcoord - x(idwn))/h;
                yd = (knotcurves[c].knotcurve[s].ycoord - y(jdwn))/h;
                zd = (knotcurves[c].knotcurve[s].zcoord - z(kdwn))/h;
                for(m=0;m<8;m++)  //linear interpolation of 8 NNs
                {
                    /* Work out increments*/
                    iinc = m%2;
                    jinc = (m/2)%2;
                    kinc = (m/4)%2;
                    /*Loop over nearest points*/
                    i = gridinc(modidwn, iinc, Nx,0);
                    j = gridinc(modjdwn, jinc, Ny,1);
                    k = gridinc(modkdwn,kinc, Nz,2);
                    prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);   //terms of the form (1-xd)(1-yd)zd etc. (interpolation coefficient)
                    /*interpolate grad u over nearest points*/
                    dxu += prefactor*0.5*(u[pt(gridinc(i,1,Nx,0),j,k)] -  u[pt(gridinc(i,-1,Nx,0),j,k)])/h;  //central diff
                    dyu += prefactor*0.5*(u[pt(i,gridinc(j,1,Ny,1),k)] -  u[pt(i,gridinc(j,-1,Ny,1),k)])/h;
                    dzu += prefactor*0.5*(u[pt(i,j,gridinc(k,1,Nz,2))] -  u[pt(i,j,gridinc(k,-1,Nz,2))])/h;
                }
                //project du onto perp of tangent direction first
                dx = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[incp(s,-1,NP)].xcoord);   //central diff as a is defined on the points
                dy = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[incp(s,-1,NP)].ycoord);
                dz = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[incp(s,-1,NP)].zcoord);
                dxup = dxu - (dxu*dx + dyu*dy + dzu*dz)*dx/(dx*dx+dy*dy+dz*dz);               //Grad u_j * (delta_ij - t_i t_j)
                dyup = dyu - (dxu*dx + dyu*dy + dzu*dz)*dy/(dx*dx+dy*dy+dz*dz);
                dzup = dzu - (dxu*dx + dyu*dy + dzu*dz)*dz/(dx*dx+dy*dy+dz*dz);
                /*Vector a is the normalised gradient of u, should point in direction of max u perp to t*/
                norm = sqrt(dxup*dxup+dyup*dyup+dzup*dzup);
                knotcurves[c].knotcurve[s].ax = dxup/norm;
                knotcurves[c].knotcurve[s].ay = dyup/norm;
                knotcurves[c].knotcurve[s].az = dzup/norm;
            }

            for(j=1; j<4; j++)
            {
                switch(j)
                {
                    case 1 :
                        for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ax ; break;
                    case 2 :
                        for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ay ; break;
                    case 3 :
                        for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].az ; break;
                }
                double* data = coord.data();
                // take the fft
                gsl_fft_real_transform (data, 1, NP, real, work);
                // 21/11/2016: make our low pass filter. To apply our filter. we should sample frequencies fn = n/Delta N , n = -N/2 ... N/2
                // this is discretizing the nyquist interval, with extreme frequency ~1/2Delta.
                // to cut out the frequencies of grid fluctuation size and larger we need a lengthscale Delta to
                // plug in above. im doing a rough length calc below, this might be overkill.
                // at the moment its just a hard filter, we can choose others though.
                // compute a rough length to set scale
                double filter;
                const double cutoff = 2*M_PI*(totlength/(6*lambda));
                for (i = 0; i < NP; ++i)
                {
                    filter = 1/sqrt(1+pow((i/cutoff),8));
                    data[i] *= filter;
                };
                // transform back
                gsl_fft_halfcomplex_inverse (data, 1, NP, hc, work);
                switch(j)
                {
                    case 1 :
                        for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ax= coord[i] ; break;
                    case 2 :
                        for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ay= coord[i] ; break;
                    case 3 :
                        for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].az = coord[i] ; break;
                }
            }
            gsl_fft_real_wavetable_free (real);
            gsl_fft_halfcomplex_wavetable_free (hc);
            gsl_fft_real_workspace_free (work);

            /*****Writhe and twist integrals******/
            NP = knotcurves[c].knotcurve.size();  //store number of points in knot curve
            double totwrithe = 0;
            double tottwist = 0;
            double dxds, dyds, dzds, dxdm, dydm, dzdm,bx, by, bz;
            totlength = 0;
            /***Do the integrals**/
            double T[3][3];
            double N[2][3];
            double deltas[3]; double ds;
            double curvature[2];double torsion;
            for(s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
            {
                for(i=0;i<3;i++)
                {

                    dx = (knotcurves[c].knotcurve[incp(s,i+1,NP)].xcoord - knotcurves[c].knotcurve[incp(s,i,NP)].xcoord);   //central diff as a is defined on the points
                    dy = (knotcurves[c].knotcurve[incp(s,i+1,NP)].ycoord - knotcurves[c].knotcurve[incp(s,i,NP)].ycoord);
                    dz = (knotcurves[c].knotcurve[incp(s,i+1,NP)].zcoord - knotcurves[c].knotcurve[incp(s,i,NP)].zcoord);
                    deltas[i] = sqrt(dx*dx+dy*dy+dz*dz);
                    T[i][0] = dx/(deltas[i]);
                    T[i][1] = dy/(deltas[i]);
                    T[i][2] = dz/(deltas[i]);
                    if(i==0)
                    {
                        knotcurves[c].knotcurve[s].length = deltas[0];
                        dxds = T[0][0];
                        dyds = T[0][1];
                        dzds = T[0][2];
                        ds = deltas[0];
                    }
                }
                for(i=0;i<2;i++)
                {
                    N[i][0] = (T[i+1][0]-T[i][0])/deltas[i];
                    N[i][1] = (T[i+1][1]-T[i][1])/deltas[i];
                    N[i][2] = (T[i+1][2]-T[i][2])/deltas[i];
                    curvature[i] = sqrt(N[i][0]*N[i][0]+N[i][1]*N[i][1]+N[i][2]*N[i][2]);
                    N[i][0] /=curvature[i];
                    N[i][1] /=curvature[i];
                    N[i][2] /=curvature[i];
                }
                // lets get the torsion by computing the x component of the binomral, and looking at the scale factor between it and dn/ds+kn
                torsion= ((N[1][0]-N[0][0])/deltas[0] + curvature[0]*T[0][0])/(T[0][1]*N[0][2] - N[0][1]*T[0][2]);

                knotcurves[c].knotcurve[s].curvature = curvature[0];
                knotcurves[c].knotcurve[s].torsion = torsion;
                /*These quantities defined on line connecting points s and s+1*/
                knotcurves[c].knotcurve[s].writhe = 0;
                bx = (knotcurves[c].knotcurve[incp(s,1,NP)].ax - knotcurves[c].knotcurve[s].ax)/ds;
                by = (knotcurves[c].knotcurve[incp(s,1,NP)].ay - knotcurves[c].knotcurve[s].ay)/ds;
                bz = (knotcurves[c].knotcurve[incp(s,1,NP)].az - knotcurves[c].knotcurve[s].az)/ds;
                knotcurves[c].knotcurve[s].twist = (dxds*(knotcurves[c].knotcurve[s].ay*bz - knotcurves[c].knotcurve[s].az*by) + dyds*(knotcurves[c].knotcurve[s].az*bx - knotcurves[c].knotcurve[s].ax*bz) + dzds*(knotcurves[c].knotcurve[s].ax*by - knotcurves[c].knotcurve[s].ay*bx))/(2*M_PI*sqrt(dxds*dxds + dyds*dyds + dzds*dzds));
                /*Check this is actually normal to tangent*/
                /*check = fabs(0.5*(knotcurves[c].knotcurve[s].ax + knotcurves[c].knotcurve[incp(s,1,NP)].ax)*dxds + 0.5*(knotcurves[c].knotcurve[s].ay + knotcurves[c].knotcurve[incp(s,1,NP)].ay)*dyds + 0.5*(knotcurves[c].knotcurve[s].az + knotcurves[c].knotcurve[incp(s,1,NP)].az)*dzds)/sqrt(dxds*dxds + dyds*dyds + dzds*dzds);
                  if(check>0.01) cout << s << ": (" << knotcurves[c].knotcurve[s].xcoord << ", " << knotcurves[c].knotcurve[s].ycoord << ", " << knotcurves[c].knotcurve[s].zcoord << "). Grad u . t = " << check << '\n';*/
                for(m=0; m<NP; m++)
                {
                    if(s != m)
                    {
                        xdiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].xcoord + knotcurves[c].knotcurve[s].xcoord - knotcurves[c].knotcurve[incp(m,1,NP)].xcoord - knotcurves[c].knotcurve[m].xcoord);   //interpolate, consistent with fwd diff
                        ydiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].ycoord + knotcurves[c].knotcurve[s].ycoord - knotcurves[c].knotcurve[incp(m,1,NP)].ycoord - knotcurves[c].knotcurve[m].ycoord);
                        zdiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].zcoord + knotcurves[c].knotcurve[s].zcoord - knotcurves[c].knotcurve[incp(m,1,NP)].zcoord - knotcurves[c].knotcurve[m].zcoord);
                        dxdm = (knotcurves[c].knotcurve[incp(m,1,NP)].xcoord - knotcurves[c].knotcurve[m].xcoord)/(ds);
                        dydm = (knotcurves[c].knotcurve[incp(m,1,NP)].ycoord - knotcurves[c].knotcurve[m].ycoord)/(ds);
                        dzdm = (knotcurves[c].knotcurve[incp(m,1,NP)].zcoord - knotcurves[c].knotcurve[m].zcoord)/(ds);
                        knotcurves[c].knotcurve[s].writhe += ds*(xdiff*(dyds*dzdm - dzds*dydm) + ydiff*(dzds*dxdm - dxds*dzdm) + zdiff*(dxds*dydm - dyds*dxdm))/(4*M_PI*(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)*sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff));
                    }
                }
                /*Add on writhe, twist and length*/
                knotcurves[c].writhe += knotcurves[c].knotcurve[s].writhe*ds;
                knotcurves[c].length += knotcurves[c].knotcurve[s].length;
                knotcurves[c].twist  += knotcurves[c].knotcurve[s].twist*ds;
            }

            c++;
        }

    }

    // the order of the components within the knotcurves vector is not guaranteed to remain fixed from timestep to timestep. thus, componenet 0 at one timtestep could be 
    // components 1 at the next. the code needs a way of tracking which componenet is which.
    // at the moment, im doing this by fuzzily comparing summary stats on the components - at this point, the length twist and writhe. 

    // these variables have the summary stats from the last timestep

    static vector<double> oldwrithe(knotcurves.size()); 
    static vector<double> oldtwist(knotcurves.size());
    static vector<double> oldlength(knotcurves.size());

    vector<int> permutation(knotcurves.size());
    if(first)
    {
        for(int i = 0; i<knotcurves.size();i++)
        {
            oldwrithe[i] = knotcurves[i].writhe;
            oldlength[i] = knotcurves[i].length;
            oldtwist[i] = knotcurves[i].twist;
            permutation[i] = i;
        }

    }
    else
    {
        for(int i = 0; i<knotcurves.size();i++)
        {
			// do j=0 manually to ensure minscore gets written to.
            double minscore = ((knotcurves[0].length - oldlength[i])/oldlength[i])*((knotcurves[0].length - oldlength[i])/oldlength[i]) +((knotcurves[0].writhe - oldwrithe[i])/oldwrithe[i])*((knotcurves[0].writhe - oldwrithe[i])/oldwrithe[i]) +((knotcurves[0].twist - oldtwist[i])/oldtwist[i])*((knotcurves[0].twist - oldtwist[i])/oldtwist[i]); 
            permutation[i] = 0;
            for(int j = 1; j<knotcurves.size();j++)
            {    
                double score = ((knotcurves[j].length - oldlength[i])/oldlength[i])*((knotcurves[j].length - oldlength[i])/oldlength[i]) +((knotcurves[j].writhe - oldwrithe[i])/oldwrithe[i])*((knotcurves[j].writhe - oldwrithe[i])/oldwrithe[i]) +((knotcurves[j].twist - oldtwist[i])/oldtwist[i])*((knotcurves[j].twist - oldtwist[i])/oldtwist[i]); 
                if(score<minscore) permutation[i] = j; minscore = score; 
            }
        }
        // apply the permutation to the list of lengths etc. It is now "correct" in the sense that index [0] really is component 0 etc. these labellings are arbitrarlly set at the simulations start and
        // must be consistently carried forward
        for(int i = 0; i<knotcurves.size();i++)
        {
            oldwrithe[i] = knotcurves[permutation[i]].writhe;
            oldlength[i] = knotcurves[permutation[i]].length;
            oldtwist[i] = knotcurves[permutation[i]].twist;

        }
    }
    first = false;
    print_knot(t, knotcurves, permutation);

    if(cleanupneeded)
    {
        knotcurves.clear(); //empty vector with knot curve points
        memset(xmarked, 0, sizeof(xmarked));
        memset(ymarked, 0, sizeof(ymarked));
        memset(zmarked, 0, sizeof(zmarked));
    }
}

void uv_update(double *u, double *v,  double *ku, double *kv)
{
    int i,j,k,l,n,kup,kdown,iup,idown,jup,jdown;
    double D2u;
    const int arraysize = Nx*Ny*Nz;
// first loop. get k1, store (in testun] and testv[n], the value u[n]+h/2k1)
    #pragma omp for 
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k);
                    kup = gridinc(k,1,Nz,2);
                    kdown = gridinc(k,-1,Nz,2);
                    D2u = oneoverhsq*(u[pt(gridinc(i,1,Nx,0),j,k)] + u[pt(gridinc(i,-1,Nx,0),j,k)] + u[pt(i,gridinc(j,1,Ny,1),k)] + u[pt(i,gridinc(j,-1,Ny,1),k)] + u[pt(i,j,kup)] + u[pt(i,j,kdown)] - 6.0*u[n]);
                    ku[n] = oneoverepsilon*(u[n] - (ONETHIRD*u[n])*(u[n]*u[n]) - v[n]) + D2u;
                    kv[n] = epsilon*(u[n] + beta - gam*v[n]);
                }
            }
        }
// 2nd and 3rd loops
    double inc ;   
    for(l=1;l<=3;l++)  //u and v update for each fractional time step
    {
        switch (l)
        {
            case 1:
                {
                    inc=0.5;   //add k1 to uv and add to total k
                }
                break;

            case 2:
                {
                    inc=0.5 ;   //add k1 to uv and add to total k
                }
                break;
            case 3:
                {
                    inc=1 ;   //add k1 to uv and add to total k
                }
                break;
        }   
        #pragma omp for 
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k);

                     iup = pt(gridinc(i,1,Nx,0),j,k);
                    idown =pt(gridinc(i,-1,Nx,0),j,k);
                    jup = pt(i,gridinc(j,1,Ny,1),k);
                    jdown =pt(i,gridinc(j,-1,Ny,1),k);
                    kup = pt(i,j,gridinc(k,1,Nz,2));
                    kdown = pt(i,j,gridinc(k,-1,Nz,2));
                    double currentu = u[n] + dtime*inc*ku[(l-1)*arraysize+n];
                    double currentv = v[n] + dtime*inc*kv[(l-1)*arraysize+n];

                    D2u = oneoverhsq*((u[iup]+dtime*inc*ku[(l-1)*arraysize+iup]) + (u[idown]+dtime*inc*ku[(l-1)*arraysize+idown]) +(u[jup]+dtime*inc*ku[(l-1)*arraysize+jup]) +(u[jdown]+dtime*inc*ku[(l-1)*arraysize+jdown]) + (u[kup]+dtime*inc*ku[(l-1)*arraysize+kup]) + (u[kdown]+dtime*inc*ku[(l-1)*arraysize+kdown])- 6.0*(currentu));


                    ku[arraysize*l+n] = oneoverepsilon*(currentu - (ONETHIRD*currentu)*(currentu*currentu) - currentv) + D2u;
                    kv[arraysize*l+n] = epsilon*(currentu + beta - gam*currentv);
                }
            }
        }
    }
    #pragma omp for 
    for(n=0;n<Nx*Ny*Nz;n++)
    {

            u[n] = u[n] + dtime*sixth*(ku[n]+2*ku[arraysize+n]+2*ku[2*arraysize+n]+ku[3*arraysize+n]);
            v[n] = v[n] + dtime*sixth*(kv[n]+2*kv[arraysize+n]+2*kv[2*arraysize+n]+kv[3*arraysize+n]);
    }

}

void uv_update_burnin(double *u, double *v,  double *ku, double *kv,bool* markedpoints)
{
    int i,j,k,l,n,kup,kdown,iup,idown,jup,jdown;
    double D2u;
    const int arraysize = Nx*Ny*Nz;
// first loop. get k1, store (in testun] and testv[n], the value u[n]+h/2k1)
    #pragma omp for 
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k);
                    kup = gridinc(k,1,Nz,2);
                    kdown = gridinc(k,-1,Nz,2);
                    D2u = oneoverhsq*(u[pt(gridinc(i,1,Nx,0),j,k)] + u[pt(gridinc(i,-1,Nx,0),j,k)] + u[pt(i,gridinc(j,1,Ny,1),k)] + u[pt(i,gridinc(j,-1,Ny,1),k)] + u[pt(i,j,kup)] + u[pt(i,j,kdown)] - 6.0*u[n]);
                    ku[n] = oneoverepsilon*(u[n] - (ONETHIRD*u[n])*(u[n]*u[n]) - v[n]) + D2u;
                    kv[n] = epsilon*(u[n] + beta - gam*v[n]);
                }
            }
        }
// 2nd and 3rd loops
    double inc ;   
    for(l=1;l<=3;l++)  //u and v update for each fractional time step
    {
        switch (l)
        {
            case 1:
                {
                    inc=0.5;   //add k1 to uv and add to total k
                }
                break;

            case 2:
                {
                    inc=0.5 ;   //add k1 to uv and add to total k
                }
                break;
            case 3:
                {
                    inc=1 ;   //add k1 to uv and add to total k
                }
                break;
        }   
        #pragma omp for 
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k);

                     iup = pt(gridinc(i,1,Nx,0),j,k);
                    idown =pt(gridinc(i,-1,Nx,0),j,k);
                    jup = pt(i,gridinc(j,1,Ny,1),k);
                    jdown =pt(i,gridinc(j,-1,Ny,1),k);
                    kup = pt(i,j,gridinc(k,1,Nz,2));
                    kdown = pt(i,j,gridinc(k,-1,Nz,2));
                    double currentu = u[n] + dtime*inc*ku[(l-1)*arraysize+n];
                    double currentv = v[n] + dtime*inc*kv[(l-1)*arraysize+n];

                    D2u = oneoverhsq*((u[iup]+dtime*inc*ku[(l-1)*arraysize+iup]) + (u[idown]+dtime*inc*ku[(l-1)*arraysize+idown]) +(u[jup]+dtime*inc*ku[(l-1)*arraysize+jup]) +(u[jdown]+dtime*inc*ku[(l-1)*arraysize+jdown]) + (u[kup]+dtime*inc*ku[(l-1)*arraysize+kup]) + (u[kdown]+dtime*inc*ku[(l-1)*arraysize+kdown])- 6.0*(currentu));


                    ku[arraysize*l+n] = oneoverepsilon*(currentu - (ONETHIRD*currentu)*(currentu*currentu) - currentv) + D2u;
                    kv[arraysize*l+n] = epsilon*(currentu + beta - gam*currentv);
                }
            }
        }
    }
    #pragma omp for 
    for(n=0;n<Nx*Ny*Nz;n++)
    {
        if(!markedpoints[n])
        {
            u[n] = u[n] + dtime*sixth*(ku[n]+2*ku[arraysize+n]+2*ku[2*arraysize+n]+ku[3*arraysize+n]);
            v[n] = v[n] + dtime*sixth*(kv[n]+2*kv[arraysize+n]+2*kv[2*arraysize+n]+kv[3*arraysize+n]);
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
    int i,j,k,l,n,kup,kdown;

#pragma omp for
    for(i=0;i<Nx;i++)
    {
        for(j=0; j<Ny; j++)
        {
            for(k=0; k<Nz; k++)
            {
                n = pt(i,j,k);
                kup = gridinc(k,1,Nz,2);
                kdown = gridinc(k,-1,Nz,2);
                D2u[n] = (u[pt(gridinc(i,1,Nx,0),j,k)] + u[pt(gridinc(i,-1,Nx,0),j,k)] + u[pt(i,gridinc(j,1,Ny,1),k)] + u[pt(i,gridinc(j,-1,Ny,1),k)] + u[pt(i,j,kup)] + u[pt(i,j,kdown)] - 6*u[n])/(h*h);
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

void print_uv( double *u, double *v, double *ucvx, double *ucvy, double *ucvz, double t, bool *markedpoints )
{
    int i,j,k,n;
    stringstream ss;
    ss << "uv_plot" << t << ".vtk";
    ofstream uvout (ss.str().c_str());

    uvout << "# vtk DataFile Version 3.0\nUV fields\nASCII\nDATASET STRUCTURED_POINTS\n";
    uvout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    uvout << "ORIGIN " << x(0) << ' ' << y(0) << ' ' << z(0) << '\n';
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

if(markedpoints != NULL)
{
    uvout << "SCALARS markedpoints float\nLOOKUP_TABLE default\n";
    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k);
                uvout << markedpoints[n] << '\n';
            }
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

void print_B_phi( double *phi, int* missed)
{
    int i,j,k,n;
    string fn = "phi.vtk";

    ofstream Bout (fn.c_str());

    Bout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
    Bout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    Bout << "ORIGIN " << x(0) << ' ' << y(0) << ' ' << z(0) << '\n';
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


void print_knot( double t, vector<knotcurve>& knotcurves, vector<int>& permutation)
{
    for( int c=0; c <= (knotcurves.size()-1) ; c++)
    {

        /***Write values to file*******/
        stringstream ss;
        ss << "globaldata" << "_" << c <<  ".txt";
        ofstream wrout (ss.str().c_str(), std::ofstream::app);
        wrout << t << '\t' << knotcurves[permutation[c]].writhe << '\t' << knotcurves[permutation[c]].twist << '\t' << knotcurves[permutation[c]].length << '\n';
        wrout.close();

        ss.str("");
        ss.clear();       

        ss << "knotplot" << c << "_" << t <<  ".vtk";
        ofstream knotout (ss.str().c_str());

        int i;
        int n = knotcurves[permutation[c]].knotcurve.size();

        knotout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        knotout << "POINTS " << n << " float\n";

        for(i=0; i<n; i++)
        {
            knotout << knotcurves[permutation[c]].knotcurve[i].xcoord << ' ' << knotcurves[permutation[c]].knotcurve[i].ycoord << ' ' << knotcurves[permutation[c]].knotcurve[i].zcoord << '\n';
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

        knotout << "\nSCALARS Curvature float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[permutation[c]].knotcurve[i].curvature << '\n'; }

        knotout << "\nSCALARS Torsion float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[permutation[c]].knotcurve[i].torsion << '\n';
        }

        knotout << "\nVECTORS A float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[permutation[c]].knotcurve[i].ax << ' ' << knotcurves[permutation[c]].knotcurve[i].ay << ' ' << knotcurves[permutation[c]].knotcurve[i].az << '\n';
        }

        knotout << "\n\nCELL_DATA " << n << "\n\n";
        knotout << "\nSCALARS Writhe float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[permutation[c]].knotcurve[i].writhe << '\n';
        }

        knotout << "\nSCALARS Twist float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[permutation[c]].knotcurve[i].twist << '\n';
        }

        knotout << "\nSCALARS Length float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[permutation[c]].knotcurve[i].length << '\n';
        }
        knotout.close();
    }
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

double my_f(const gsl_vector* minimum, void* params)
{

    int i,j,k,idwn,jdwn,kdwn,modidwn,modjdwn,modkdwn,m,pts,iinc,jinc,kinc;
    double ucvxs, ucvys, ucvzs,  xd, yd ,zd, xdiff, ydiff, zdiff, prefactor;
    struct parameters* myparameters = (struct parameters *) params;
    double* ucvx= myparameters->ucvx;
    double* ucvy= myparameters->ucvy;
    double* ucvz= myparameters->ucvz;

    gsl_vector* tempf = gsl_vector_alloc (3);
    gsl_vector* tempv = gsl_vector_alloc (3);
    gsl_vector* tempb = gsl_vector_alloc (3);
    gsl_vector_memcpy (tempf,myparameters->f);
    gsl_vector_memcpy (tempv,myparameters->v);
    gsl_vector_memcpy (tempb,myparameters->b);

    // s gives us how much of f to add to p

    gsl_vector_scale(tempf,gsl_vector_get (minimum, 0));
    gsl_vector_scale(tempb,gsl_vector_get (minimum, 1));
    gsl_vector_add(tempf,tempb);
    gsl_vector_add(tempv,tempf);
    double px = gsl_vector_get(tempv, 0);
    double py = gsl_vector_get(tempv, 1);
    double pz = gsl_vector_get(tempv, 2);
    gsl_vector_free(tempf);
    gsl_vector_free(tempv);
    gsl_vector_free(tempb);

    /**Find nearest gridpoint**/
    idwn = (int) ((px/h) - 0.5 + Nx/2.0);
    jdwn = (int) ((py/h) - 0.5 + Ny/2.0);
    kdwn = (int) ((pz/h) - 0.5 + Nz/2.0);
    modidwn = circularmod(idwn,Nx);
    modjdwn = circularmod(jdwn,Ny);
    modkdwn = circularmod(kdwn,Nz);
    pts=0;
    ucvxs=0;
    ucvys=0;
    ucvzs=0;
    /*curve to gridpoint down distance*/
    xd = (px - x(idwn))/h;
    yd = (py - y(jdwn))/h;
    zd = (pz - z(kdwn))/h;
    for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
    {
        /* Work out increments*/
        iinc = m%2;
        jinc = (m/2)%2;
        kinc = (m/4)%2;
        /*Loop over nearest points*/
        i = gridinc(modidwn, iinc, Nx,0);
        j = gridinc(modjdwn, jinc, Ny,1);
        k = gridinc(modkdwn,kinc, Nz,1);
        prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
        /*interpolate grad u x grad v over nearest points*/
        ucvxs += prefactor*ucvx[pt(i,j,k)];
        ucvys += prefactor*ucvy[pt(i,j,k)];
        ucvzs += prefactor*ucvz[pt(i,j,k)];
    }
    double ans = -1*sqrt(ucvxs*ucvxs + ucvys*ucvys + ucvzs*ucvzs);
    return  ans;
}
void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product)
{
    double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
        - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

    double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
        - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

    double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
        - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

    gsl_vector_set(product, 0, p1);
    gsl_vector_set(product, 1, p2);
    gsl_vector_set(product, 2, p3);
}
void rotatedisplace(double& xcoord, double& ycoord, double& zcoord, const double theta, const double phi, const double dispx,const double dispy,const double dispz)
{

    double xprime = cos(phi)*cos(theta)*xcoord -sin(phi)*ycoord + cos(phi)*sin(theta)* zcoord;
    double yprime = sin(phi)*cos(theta)*xcoord +cos(phi)*ycoord + sin(phi)*sin(theta)*zcoord ;
    double zprime = -sin(theta)*xcoord  + cos(theta)*zcoord;
    xcoord = xprime;
    ycoord = yprime;
    zcoord = zprime; 

}
inline int circularmod(int i, int N)    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
{
    if(i<0) return N - ((-i)%N);
    else return i%N;
}
// inlined functions for incrementing things respecting boundaries
inline int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}

inline int incw(int i, int p, int N)    //increment with reflecting boundary between -1 and 0 and N-1 and N
{
    if(i+p<0) return (-(i+p+1));
    if(i+p>N-1) return (2*N-(i+p+1));
    return (i+p);
}

// this function is specifically designed to incremenet, in the direction specified, respecting the boundary conditions, which are global enums
inline int gridinc(int i, int p, int N, int direction )    //increment with reflecting boundary between -1 and 0 and N-1 and N
{

    if(BoundaryType == ALLREFLECTING)
    {
        return incw(i,p,N);
    }

    if(BoundaryType == ALLPERIODIC)
    {
        return incp(i,p,N);
    }

    if(BoundaryType == ZPERIODIC)
    {
        if(direction ==2) return incp(i,p,N);
        else return incw(i,p,N);
    }
return 0;
}
inline double x(int i)
{
    return (i+0.5-Nx/2.0)*h;
}
inline double y(int i)
{
    return (i+0.5-Ny/2.0)*h;
}
inline double z(int i)
{
    return (i+0.5-Nz/2.0)*h;
}
/*******Erase some points********/

// for(s=0; s<NP%8; s++) knotcurve.pop_back();    //delete last couple of elements
// for(s=0; s<NP/8; s++)                          //delete 7 of every 8 elements
// {
//     knotcurve.erase(knotcurve.end()-s-8,knotcurve.end()-s-1);
// }

/********************************/

//NP = knotcurves[c].size();  //update number of points in knot curve
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
