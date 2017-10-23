#include "FN_Constants.h"
#include "FN_Knot.h"
using namespace std;

#ifndef READINGWRITING_H
#define READINGWRITING_H

void print_B_phi(vector<double>&phi, const Griddata &griddata);
void print_uv(vector<double>&u, vector<double>&v, vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>&ucvmag, double t, const Griddata &griddata);
void print_knot(double t, vector<knotcurve>& knotcurves, const Griddata &griddata);
int uvfile_read(vector<double>&u, vector<double>&v, vector<double>& ku, vector<double>& kv, vector<double>& ucvx, vector<double>& ucvy, vector<double>& ucvz, vector<double> &ucvmag, Griddata &griddata);
int uvfile_read_ASCII(vector<double>&u, vector<double>&v, const Griddata &griddata); // for legacy purposes
int uvfile_read_BINARY(vector<double>&u, vector<double>&v, const Griddata &griddata);
float FloatSwap( float f );
void ByteSwap(const char* TobeSwapped, char* swapped );


#endif //READINGWRITING_H
