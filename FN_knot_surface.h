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
#include <gsl/gsl_min.h>
#include <gsl/gsl_vector.h>
using namespace std;

#define FROM_PHI_FILE 0
#define FROM_SURFACE_FILE 1
#define FROM_UV_FILE 2
#define FROM_KNOT_FILE 3

const double sixth = 1.0/6.0;

struct parameters
{
	gsl_vector *v,*f;
	double *x,*y,*z,*ucvx,*ucvy,*ucvz;	
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
    double twist;    //local twist value
    double writhe;   //local writhe value
    double length;   //length of line
};

/*************************General maths and integer functions*****************************/

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

inline int sign(int i)
{
    if(i==0) return 0;
    else return i/abs(i);
}

/*************************Functions for knot initialisation*****************************/

double initialise_knot();

double init_from_surface_file(void);

double init_from_knot_file(void);

/*************************Functions for B and Phi calcs*****************************/

void initial_cond(double *x, double *y, double *z, double *phi, int *missed);

void phi_calc(double *x, double *y, double *z, double *phi);

void B_field_calc(double *x, double *y, double *z, double *Bx, double *By, double *Bz, double *Bmag, int *ignore, int *ignore1, int *missed);

void phi_calc_B(double *Bx, double *By, double *Bz, double *Bmag, int *ignore, int *ignore1, int *missed, double *phi);

int pathfind(int i0, int j0, int k0, int ie, int je, int ke, int *pi, int *pj, int *pk, int *ignore, double *Bx, double *By, double *Bz, double *Bmag);

//FitzHugh Nagumo functions
void uv_initialise(double *phi, double *u, double *v, int* missed);
void crossgrad_calc(double *x, double *y, double *z, double *u, double *v, double *ucvx, double *ucvy, double *ucvz);
void find_knot_properties(double *x, double *y, double *z, double *ucvx, double *ucvy, double *ucvz, double* u,double t, gsl_min_fminimizer *minimizerstate);
void uv_update(double *u, double *v, double *ku, double *kv, double *kut, double *kvt, double *uold, double *vold);
void uv_add(double *u, double *v, double* uold, double *vold, double *ku, double *kv, double *kut, double *kvt, double inc, double coeff);
void uv_update_euler(double *u, double *v, double *D2u);
double my_f(const double s, void* params);
// 3d geometry functions
int intersect3D_SegmentPlane( knotpoint SegmentStart, knotpoint SegmentEnd, knotpoint PlaneSegmentStart, knotpoint PlaneSegmentEnd, double& IntersectionFraction, std::vector<double>& IntersectionPoint );

/*************************File reading and writing*****************************/

void print_B_phi(double *x, double *y, double*z, double *phi, int *missed);
void print_uv(double *x, double *y, double *z, double *u, double *v, double *ucvx, double *ucvy, double *ucvz, double t);
int phi_file_read(double *phi);
void print_knot(double *x, double *y, double *z, double t, std::vector<knotpoint>& knotcurve);
int uvfile_read(double *u,double *v);
void print_info(int Nx, int Ny, int Nz, double dtime, double h, const bool periodic,  int option, string knot_filename, string B_filename);
