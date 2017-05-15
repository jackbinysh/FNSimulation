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
} InitialData;


//Kernel configuration parameter for boundary conditions. All dimensions must be multiple of this.
#define BOUNDARYBLOCKSIZE 16	

//3D Laplace-computing kernel configuration parameters. Dimensions of a kernel is how many threads it contains
#define CELLW	8	//width of thread block
#define CELLH	8	//height of thred block
#define CELLD	8	//depth of thread block. At least 3, but value is arbitrary since CELLW*CELLH is already multiple of the warp size (32)
//because of the overlap, reads in the laplacian computation kernel will never be coalesced.
//block dimensions are chosen to minimize redundancy in reading data.

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
