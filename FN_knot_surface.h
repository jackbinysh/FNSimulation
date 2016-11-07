
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

using namespace std;

#define FROM_PHI_FILE 0
#define FROM_SURFACE_FILE 1
#define FROM_UV_FILE 2

struct triangle
{
    double xvertex[3];   //x components of vertices
    double yvertex[3];   //y components of vertices
    double zvertex[3];   //z components of vertices
    double normal[3];   //normal vector
    double area;      //surface area
    double centre[3];  //centre position vector
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

/*************************Functions for B and Phi calcs*****************************/

void initial_cond(double *x, double *y, double *z, double *phi, unsigned int *missed);

void phi_calc(double *x, double *y, double *z, unsigned int *missed, double *phi);

//FitzHugh Nagumo functions
void uv_initialise(double *phi, double *u, double *v, double *ucv, unsigned int *missed);
void crossgrad_calc(double *u, double *v, double *ucv);
void uv_update(double *u, double *v, double **ku, double **kv, double *uold, double *vold);

/*************************File reading and writing*****************************/

void print_B_phi(double *x, double *y, double*z, unsigned int *missed, double *phi);
void print_uv(double *u, double *v, double *ucv, double t);
int phi_file_read(double *phi, unsigned int *missed);
int uvfile_read(double *u,double *v);
void print_info(int Nx, int Ny, int Nz, double dtime, double h, const bool periodic, unsigned int option, string knot_filename, string B_filename);
