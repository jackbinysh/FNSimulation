//Constants.h
/*Available options:
FROM_PHI_FILE: Skip initialisation, input from previous run.
FROM_SURFACE_FILE: Initialise from input file(s) generated in surface evolver.
FROM_UV_FILE: Skip initialisation, run FN dynamics from uv file
FROM_FUNCTION: Initialise from some function which can be implemented by the user in phi_calc_manual. eg using theta(x) = artcan(y-y0/x-x0) to give a pole at x0,y0 etc..:wq
 */

#include <string>

#define FROM_PHI_FILE 0
#define FROM_SURFACE_FILE 1
#define FROM_UV_FILE 2
#define FROM_FUNCTION 3
enum BoundaryType {ALLREFLECTING, ZPERIODIC, ALLPERIODIC};
const int option = FROM_SURFACE_FILE;         //unknot default option
const BoundaryType BoundaryType=ALLREFLECTING;

/** two rotation angles for the initial stl file, and a displacement vector for the file **/
const double initialthetarotation = 0;
const double initialphirotation = 0;
const double initialxdisplacement = 0;
const double initialydisplacement = 0;
const double initialzdisplacement = 0;

/**If FROM_SURFACE_FILE chosen**/
std::string knot_filename = "five1";      //if FROM_SURFACE_FILE assumed input filename format of "XXXXX.stl"
//if ncomp > 1 (no. of components) then component files should be separated to 'XXXXX.txt" "XXXXX2.txt", ....
/**IF FROM_PHI_FILE or FROM_UV_FILE chosen**/
std::string B_filename = "uv_plot0.vtk";    //filename for phi field or uv field

//Grid points
const int initialNx = 201;   //No. points in x,y and z
const int initialNy = 201;
const int initialNz = 91;
const double TTime = 20000;       //total time of simulation (simulation units)
const double skiptime = 500;       //print out every # unit of time (simulation units)
const double starttime =0;        //Time at start of simulation (non-zero if continuing from UV file)
const double dtime = 0.02;         //size of each time step

//System size parameters
const double lambda = 21.3;                //approx wavelength
const double size = 6*lambda;   //box size
const double h = size/(initialNx-1);            //grid spacing
const double oneoverhsq = 1.0/(h*h);
const double epsilon = 0.3;                //parameters for F-N eqns
const double oneoverepsilon = 1.0/epsilon;
const double beta = 0.7;
const double gam = 0.5;

//Size boundaries of knot (now autoscaled)
double xmax = 3*initialNx*h/4.0;
double ymax = 3*initialNy*h/4.0;
double zmax = 3*initialNz*h/4.0;
