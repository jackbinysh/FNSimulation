#include "FN_Constants.h"
#include "TriCubicInterpolator.h"
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
#include <complex.h>
#include <fftw3.h>
using namespace std;

#ifndef FNKNOT_H
#define FNKNOT_H
struct Plans
{

    // the data
    fftw_complex* uhat;
    fftw_complex* vhat;
    fftw_complex* uhatnext;
    fftw_complex* vhatnext;
    fftw_complex* uhattemp;
    double* utemp;

    double complex* L;
    double complex* Lhalf;

    // the plans
    fftw_plan utemp_to_uhattemp;
    fftw_plan uhattemp_to_utemp;
    fftw_plan uext_to_uhat;
    fftw_plan vext_to_vhat;
    fftw_plan uhattemp_to_uext;
    fftw_plan uhattemp_to_vext;
};

struct Griddata
{
    int Nx,Ny,Nz;
    double h;
};
struct parameters
{
	gsl_vector *v,*f,*b;
    likely::TriCubicInterpolator* ucvmag;
    Griddata mygriddata;
};

struct viewpoint
{
    double xcoord;
    double ycoord;
    double zcoord;
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
    double modxcoord;   //position vector x coord, modded out by the lattice
    double modycoord;   //position vector y coord, ""
    double modzcoord;   //position vector z coord, ""
    double ax;       //grad vector x coord
    double ay;       //grad vector y coord
    double az;       //grad vector z coord
    double tx;       //grad vector x coord
    double ty;       //grad vector y coord
    double tz;       //grad vector z coord
    double nx;       //grad vector x coord
    double ny;       //grad vector y coord
    double nz;       //grad vector z coord
    double bx;       //grad vector x coord
    double by;       //grad vector y coord
    double bz;       //grad vector z coord
    double vx;       //grad vector x coord
    double vy;       //grad vector y coord
    double vz;       //grad vector z coord
    double kappaNx;  //curvature vector x component
    double kappaNy;  //curvature vector x component
    double kappaNz;  //curvature vector x component
    double vdotnx;       //grad vector x coord
    double vdotny;       //grad vector x coord
    double vdotnz;       //grad vector y coord
    double vdotbx;       //grad vector x coord
    double vdotby;       //grad vector x coord
    double vdotbz;       //grad vector y coord
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
    double xavgpos;  // average position of the knot
    double yavgpos;
    double zavgpos;
};

struct Link
{
    std::vector<knotcurve> Components;
    int NumComponents;
    int NumPoints;
    double length;   //total length of line
    double writhe;
    double twist;
};

/*************************General maths and integer functions*****************************/

// little inline guys
double x(int i, const Griddata &griddata);
double y(int i, const Griddata &griddata);
double z(int i, const Griddata &griddata);
int sign(int i);
int pt(int i,  int j,  int k, const Griddata &griddata);       //convert i,j,k to single index
int coordstopt(double x, double y, double z, Griddata&griddata);
int circularmod(int i, int N);    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
int incp(int i, int p, int N);    //increment i with p for periodic boundary
int incw(int i, int p, int N);    //increment with reflecting boundary between -1 and 0 and N-1 and N
int gridinc(int i, int p, int N, int direction );    //increment with reflecting boundary between -1 and 0 and N-1 and N

void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product);
double my_f(const gsl_vector* minimum, void* params);
void rotatedisplace(double& xcoord, double& ycoord, double& zcoord, const double theta, const double dispx,const double dispy,const double dispz);


/*************************Functions for B and Phi calcs*****************************/


//FitzHugh Nagumo functions
void uv_initialise(vector<double>&phi, vector<double>&u, vector<double>&v,const Griddata& griddata);
void crossgrad_calc(vector<double>&u, vector<double>&v, vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>&ucvmag, const Griddata &griddata);
void find_knot_properties(vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>& ucvmag, vector<double>&u, vector<knotcurve>& knotcurves, double t, gsl_multimin_fminimizer* minimizerstate, const Griddata &griddata);
void find_knot_velocity(const vector<knotcurve>& knotcurves, vector<knotcurve>& knotcurvesold, const Griddata &griddata, const double deltatime);
void uhatvhat_initialise(const Plans& plans, const Griddata& Griddata);
void uv_update(const Plans &plans, const Griddata& Griddata);
void uv_update_external(vector<double>&u, vector<double>&v, const Plans &plans, const Griddata& Griddata);
void Initialise(vector<double>&u, vector<double>&v,Plans &plans, const Griddata& Griddata);
// 3d geometry functions
int intersect3D_SegmentPlane( knotpoint SegmentStart, knotpoint SegmentEnd, knotpoint PlaneSegmentStart, knotpoint PlaneSegmentEnd, double& IntersectionFraction, std::vector<double>& IntersectionPoint );

// things for the grown function

inline int incabsorb(int i, int p, int N);

#endif //FNKNOT_H
