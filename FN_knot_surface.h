//
//  FN_knot_code.h
//
//
//  Created by Carl Whitfield on 17/05/2016.
//
//  Last modified 3/11/16

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
using namespace std;

#define FROM_PHI_FILE 0
#define FROM_SURFACE_FILE 1
#define FROM_UV_FILE 2
#define FROM_FUNCTION 3
const double sixth = 1.0/6.0;
const double ONETHIRD = 1.0/3.0;
enum BoundaryType {ALLREFLECTING, ZPERIODIC, ALLPERIODIC};

struct parameters
{
	gsl_vector *v,*f,*b;
	double *ucvx,*ucvy,*ucvz;	
};

struct triangle
{
    double xvertex[3];   //x components of vertices
    double yvertex[3];   //y components of vertices
    double zvertex[3];   //z components of vertices
    double normal[3];   //normal vector
    double area;      //surface area
    double centre[3];  //centre position vector
};

struct knotpoint
{
    double xcoord;   //position vector x coord
    double ycoord;   //position vector y coord
    double zcoord;   //position vector z coord
    double ax;       //grad vector x coord
    double ay;       //grad vector y coord
    double az;       //grad vector z coord
	double vx;       // x velocity
    double vy;       //y velocity
    double vz;       //z velocity
    double curvature;        // curvature
    double torsion;        // torsion
    double twist;    //local twist value
    double writhe;   //local writhe value
    double length;   //length of line
};

struct knotcurve
{
    std::vector<knotpoint> knotcurve; // the actual data of the curve
    // global data for the knot component
    double twist;    //total twist value
    double writhe;   //total  writhe value
    double length;   //total lengthh of line
};
/*************************General maths and integer functions*****************************/

// little inline guys
inline double x(int i);
inline double y(int i);
inline double z(int i);
inline int sign(int i);
inline  int pt( int i,  int j,  int k);       //convert i,j,k to single index
inline int circularmod(int i, int N);    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
inline int incp(int i, int p, int N);    //increment i with p for periodic boundary
inline int incw(int i, int p, int N);    //increment with reflecting boundary between -1 and 0 and N-1 and N
inline int gridinc(int i, int p, int N, int direction );    //increment with reflecting boundary between -1 and 0 and N-1 and N

void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product);
double my_f(const gsl_vector* minimum, void* params);
void rotatedisplace(double& xcoord, double& ycoord, double& zcoord, const double theta, const double phi, const double dispx,const double dispy,const double dispz);
/*************************Functions for knot initialisation*****************************/

double initialise_knot(std::vector<triangle>& knotsurface);

double init_from_surface_file(std::vector<triangle>& knotsurface);

void scalefunction(double *scale, double *midpoint, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin);

/*************************Functions for B and Phi calcs*****************************/

void phi_calc( double *phi,std::vector<triangle>& knotsurface);

void phi_calc_manual( double *phi);

//FitzHugh Nagumo functions
void uv_initialise(double *phi, double *u, double *v);
void crossgrad_calc( double *u, double *v, double *ucvx, double *ucvy, double *ucvz);
void find_knot_properties( double *ucvx, double *ucvy, double *ucvz, double* u,std::vector<knotcurve>& knotcurves,double t, gsl_multimin_fminimizer *minimizerstate);
void uv_update(double *u, double *v,  double *ku, double *kv);
// 3d geometry functions
int intersect3D_SegmentPlane( knotpoint SegmentStart, knotpoint SegmentEnd, knotpoint PlaneSegmentStart, knotpoint PlaneSegmentEnd, double& IntersectionFraction, std::vector<double>& IntersectionPoint );

/*************************File reading and writing*****************************/

void print_B_phi( double *phi);
void print_uv( double *u, double *v, double *ucvx, double *ucvy, double *ucvz, double t);
int phi_file_read(double *phi);
void print_knot( double t, vector<knotcurve >& knotcurves,vector<int>& permutation);
int uvfile_read(double *u,double *v);
