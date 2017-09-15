//Constants.h
#include <string>

// OPTION - // what kind of initialisation
/* the different initialisation options
Available options:
FROM_PHI_FILE: Skip initialisation, input from previous run.
FROM_SURFACE_FILE: Initialise from input file(s) generated in surface evolver.
FROM_UV_FILE: Skip initialisation, run FN dynamics from uv file
FROM_FUNCTION: Initialise from some function which can be implemented by the user in phi_calc_manual. eg using theta(x) = artcan(y-y0/x-x0) to give a pole at x0,y0 etc..:wq
 */
//if ncomp > 1 (no. of components) then component files should be separated to 'XXXXX.txt" "XXXXX2.txt", ....
std::string B_filename = "test.vtk";    //filename for phi field or uv field

// OPTION - what grid values do you want/ timestep
//Grid points
const double h = 0.641566;            //grid spacing
const int initialNx = 238;   //No. points in x,y and z
const int initialNy = 238;
const int initialNz = 238;

// OPTION - do you want to resize the box? if so, when?
const bool BoxResizeFlag = 0;
const double BoxResizeTime = 1000;


