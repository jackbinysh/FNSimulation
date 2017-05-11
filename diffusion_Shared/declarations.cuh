#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//structure to hold size of simulation and arrays on the host
typedef struct DataArray
{
	int xmax;
	int ymax;
	int zmax;
	float* a;
	//float* b;
	//float* c;
} DataArray;

//Initial data for the simulation
typedef struct InitialData
{
	float da; //Diffusion coefficient for species A
	float dt; //Timestep
	float dx; //Grid spacing
	float anull; //Initial concentration of species A
	float radius; //Radius of the sphere in which A is located initially 
	int itmax; //Max number of iteration steps
} InitialData;


//Kernel configuration parameter for boundary conditions. All dimensions must be multiple of this.
#define BLOCKSIZE 16	

//3D Laplace-computing kernel configuration parameters. Dimensions of a kernel is how many threads it contains
#define CELLW	8	//width of thread block
#define CELLH	8	//height of thred block
#define CELLD	8	//depth of thread block. At least 3, but value is arbitrary since CELLW*CELLH is already multiple of the warp size (32)
//because of the overlap, reads in the laplacian computation kernel will never be coalesced.
//block dimensions are chosen to minimize redundancy in reading data.

//Preprocessor functions to compute 2nd power using multiplication instead of pow().
#define P2(x) ((x)*(x))
#define p2(x) ((x)*(x))

//fundamental variables for addressing 3D coordinates in 1D arrays
extern int zsize, ysize, total;

//functions in main.cu
void CudaCheck();
void Error(const char* text);

//functions in host_functions.cu
void hostInitialize(DataArray& host);
void hostClose(DataArray& host);

//functions in device_functions.cu
void gpuClose(DataArray& device1, DataArray& device2);
void gpuStep(DataArray& device1, DataArray& device2);
void gpuInitialize(InitialData& id, DataArray& host, DataArray& device1, DataArray& device2);
void ExportCheckpoint(const char* name, DataArray& host, DataArray& device1, int l);

#endif //declarations_h
