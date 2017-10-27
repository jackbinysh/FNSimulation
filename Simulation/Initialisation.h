
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
#include "FN_Knot.h"
using namespace std;

#ifndef SOLIDANGLECALCULATION_H
#define SOLIDANGLECALCULATION_H

struct viewpoint
{
    double xcoord;
    double ycoord;
    double zcoord;
};

/*************************Functions for knot initialisation*****************************/

void InitialiseFromFile(struct Link& Curve);
void RefineCurve(struct Link& Curve);

/**********************Functions for curve geometry************************/

void ComputeLengths(struct Link& Curve);
void ComputeTangent(struct Link& Curve);
void ComputeKappaN(struct Link& Curve);
void ComputeWrithe(struct Link& Curve);

/***********************Functions for outputting the solid angle*************************/

double SolidAngleCalc(const Link& Curve, const viewpoint& View);
void phi_calc_curve(vector<double> &phi, const struct Link& Curve, const Griddata &griddata);
double init_from_surface_file(std::vector<triangle>& knotsurface);

#endif //SOLID_ANGLECALCULATION_H
