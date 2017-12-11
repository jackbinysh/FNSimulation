#include <string>

#ifndef FNCONSTANTS_H
#define FNCONSTANTS_H

#define FROM_SURFACE_FILE 0
#define FROM_CURVE_FILE 1
#define FROM_UV_FILE 2
#define FROM_FUNCTION 3
// the different boundary conditions
enum BoundaryType {ALLREFLECTING, ZPERIODIC, ALLPERIODIC};

/* CHANGE THESE OPTIONS */

// OPTION - // what kind of initialisation
/* the different initialisation options
Available options:
FROM_PHI_FILE: Skip initialisation, input from previous run.
FROM_SURFACE_FILE: Initialise from input file(s) generated in surface evolver.
FROM_UV_FILE: Skip initialisation, run FN dynamics from uv file
FROM_FUNCTION: Initialise from some function which can be implemented by the user in phi_calc_manual. eg using theta(x) = artcan(y-y0/x-x0) to give a pole at x0,y0 etc..:wq
 */
//if ncomp > 1 (no. of components) then component files should be separated to 'XXXXX.txt" "XXXXX2.txt", ....
const int option = INSERT_INITIALISATION_TYPE;         //unknot default option
const std::string knot_filename = "INSERT_SURFACE_FILENAME";      //if FROM_SURFACE_FILE assumed input filename format of "XXXXX.stl"
const std::string B_filename = "INSERT_UV_FILENAME";    //filename for phi field or uv field
const int NumComponents = INSERT_NUM_COMPONENTS;   //No. points in x,y and z

// OPTION - what kind of boundary condition
const BoundaryType BoundaryType=INSERT_BOUNDARY_TYPE;

//OPTION - do you want the geometry of the input file to be exactly preserved, or can it be scaled to fit the box better
#define PRESERVE_RATIOS INSERT_PRESERVE_RATIOS  //1 to scale input file preserving the aspect ratio

// OPTION - how long should it run, when do you want data printed, what time value should it start at
const double TTime = INSERT_RUNTIME;       //total time of simulation (simulation units)
const double UVPrintTime = INSERT_UVPRINTTIME;       //print out UV every # unit of time (simulation units)
const double VelocityKnotplotPrintTime = INSERT_VELOCITYPRINTTIME;       //print out the velocity every # unit of time (simulation units)
const double FrequentKnotplotPrintTime = INSERT_FREQUENTPRINTTIME; // print out the knot , without the velocity
const double InitialSkipTime = INSERT_SKIPTIME;       // amout to skip before beginning the curve tracing

// OPTION - what grid values do you want/ timestep
//Grid points
const double initialh = INSERT_GRIDSPACING;            //grid spacing
const int initialNx = INSERT_NX;   //No. points in x,y and z
const int initialNy = INSERT_NY;
const int initialNz = INSERT_NZ;

// OPTION - do you want to read in a coarse uv file , and interpolate onto a finer grid? If so,
// first, set the flag to 1 if you want, 0 if you dont.
// give the # points in each dimension, which should be > initialNx - the spacing will be set by (initialNx-1)*h/(interpolatedNx-1)
const int interpolationflag = INSERT_INTERPOLATION_FLAG ;
const int interpolatedNx = INSERT_INTERPOLATED_NX;   //No. points in x,y and z
const int interpolatedNy = INSERT_INTERPOLATED_NY;
const int interpolatedNz = INSERT_INTERPOLATED_NZ;


// timestep
const double dtime = INSERT_TIMESTEP;         //size of each time step

// OPTION - do you want to resize the box? if so, when?
const bool BoxResizeFlag = 0;
const double BoxResizeTime = 1000;


// OPTION - how big should the knot be in the box, do you want it tilted or displaced?
const double sensorxcoord = 35;
const double sensorycoord = 0;
const double sensorzcoord = 0;

// OPTION - how big should the knot be in the box, do you want it tilted or displaced?
//Size boundaries of knot (now autoscaled)
const double xmax = 8*initialNx*initialh/10.0;
const double ymax = 8*initialNy*initialh/10.0;
const double zmax = 8*initialNz*initialh/10.0;
/** two rotation angles for the initial stl file, and a displacement vector for the file **/
const double initialthetarotation = 0;
const double initialxdisplacement = 0;
const double initialydisplacement = 0;
const double initialzdisplacement = 0;
const double xscalefactortweak = 0;
const double yscalefactortweak = 0;
const double zscalefactortweak = 0;

// OPTION - what system params do you want . Don't touch these usually
//System size parameters
const double lambda = 21.3;                //approx wavelength
const double epsilon = 0.3;                //parameters for F-N eqns
const double beta = 0.7;
const double gam = 0.5;


#endif //FNCONSTANTS_H
