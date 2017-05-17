#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FN_Configuration.h" 

//structure to hold size of simulation and arrays on the host
typedef struct DataArray
{
	int xmax;
	int ymax;
	int zmax;
	gridprecision* u;
	gridprecision* v;
} DataArray;

//Initial data for the simulation
typedef struct InitialData
{
	float da; //Diffusion coefficient for species A
	float dt; //Timestep
	float dx; //Grid spacing
	float epsilon; //Grid spacing
	float beta; //Grid spacing
	float gam; //Grid spacing
	float oneoverepsilon; //Grid spacing
	float onethird; //Grid spacing
} InitialData;


//Kernel configuration parameter for boundary conditions. All dimensions must be multiple of this.
#define BOUNDARYBLOCKSIZE 16	


# ifdef MOVINGTILES
//3D Laplace-computing kernel configuration parameters. Dimensions of a kernel is how many threads it contains
#define CELLW	16	//width of thread block
#define CELLH	5	//height of thred block
#define CELLD	5	//depth of thread block. At least 3, but value is arbitrary since CELLW*CELLH is already multiple of the warp size (32)
#endif //MOVINGTILES

# ifdef SHARED
//3D Laplace-computing kernel configuration parameters. Dimensions of a kernel is how many threads it contains
#define CELLW	8	//width of thread block
#define CELLH	8	//height of thred block
#define CELLD	8	//depth of thread block. At least 3, but value is arbitrary since CELLW*CELLH is already multiple of the warp size (32)
#endif //SHARED

// FOR SHARED
// all should by 8
//because of the overlap, reads in the laplacian computation kernel will never be coalesced.
//block dimensions are chosen to minimize redundancy in reading data.
//      On 2nd generation devices all global memory operations are coalesced, block dimensions are chosen to minimize read redundancy
// FOR MOVING TILES
//NOTE: On devices with compute capability 1.0 and 1.1 (1st generation devices) CELLW should be 16, CELLH and CELLD: 5.
//      this causes reading more data in every step, but the read will be coalesced!

//Preprocessor functions to compute 2nd power using multiplication instead of pow().
#define P2(x) ((x)*(x))
#define p2(x) ((x)*(x))

//functions in main.cu
void CudaCheck();
void Error(const char* text);

//functions in device_functions.cu
void gpuClose(DataArray& device1, DataArray& device2);
void gpuStep(DataArray& device1, DataArray& device2);
void gpuInitialize(InitialData& id, DataArray& host, DataArray& device1, DataArray& device2);

#endif //declarations_h
