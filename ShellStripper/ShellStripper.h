#include "ShellStripper_Constants.h"
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

struct griddata 
{
	double Nx,Ny,Nz;	
};

/*************************General maths and integer functions*****************************/

// little inline guys
inline double x(int i, const griddata& griddata);
inline double y(int i,const griddata& griddata);
inline double z(int i,const griddata& griddata);
inline int sign(int i);
inline  int pt( int i,  int j,  int k,const griddata& griddata);       //convert i,j,k to single index
inline int circularmod(int i, int N);    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
inline int incp(int i, int p, int N);    //increment i with p for periodic boundary
inline int gridinc(int i, int p, int N, int direction );    //increment with reflecting boundary between -1 and 0 and N-1 and N

void resizebox(vector<double>&u,vector<double>&v,vector<double>&ucvx,vector<double>&ucvy,vector<double>&ucvz,vector<knotcurve>&knotcurves,vector<double>&ku,vector<double>&kv,griddata& oldgriddata);

void print_marked( vector<int>&marked,int shelllabel, const griddata& griddata);

void print_uv( vector<double>&u, vector<double>&v,vector<double>&ucvmag, double t, const griddata& griddata);
int uvfile_read_BINARY(vector<double>&u, vector<double>&v,vector<double>&ucvmag,const griddata& griddata);
float FloatSwap( float f );
void ByteSwap(const char* TobeSwapped, char* swapped );

// things for the grown function

inline int incabsorb(int i, int p, int N);
void growshell(vector<double>&u,vector<int>& marked,double ucrit, const griddata& griddata);
void grow(const vector<double>&u,vector<int>&marked,double ucrit,const griddata& griddata);

