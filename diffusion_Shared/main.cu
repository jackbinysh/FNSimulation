#include "declarations.cuh"

//data arrays and simulation parameters
DataArray host, device1, device2;
InitialData id;

//variables for memory allocation and addressing
int zsize, ysize, total;

//Function to check for error after a CUDA call.
void CudaCheck()
{
	cudaError_t er = cudaGetLastError();
	if (er!=0)
	{
		printf("CUDA Error: %s\n", cudaGetErrorString(er));
#ifdef _WIN32
		system("PAUSE");
#endif
		exit(-1);
	};	
};

//function to halt program with an error message
void Error(const char* text)
{
	printf("%s\n", text);
#ifdef _WIN32
		system("PAUSE");
#endif
	exit(-1);
}

//main program
int main()
{	
	//define dimensions of simulated space
	host.xmax = 192;
	host.ymax = 192;
	host.zmax = 192;
	
	//compute 3D addressing variables
	total = host.xmax * host.ymax * host.zmax;
	zsize = host.xmax * host.ymax;
	ysize = host.xmax;

	//define parameters
	id.da = 1.0f;
	id.dt = 0.02f;
	id.dx = 1.0f;
	id.anull = 10.0f;
	id.itmax = 10000;
	id.radius = 20.0f;

	//INITIALIZE
	hostInitialize(host);
	gpuInitialize(id, host, device1, device2);
	
	int reportInterval = 100;
	int exportInterval = 1000; //must be multiple of reportInterval

	//events for time measurement
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//COMPUTE
	cudaEventRecord(start, 0);
	for (int i = 1; i <= id.itmax; i++)
	{
		gpuStep(device1, device2);

		if (i % reportInterval == 0)
		{
			cudaEventRecord(stop, 0);
			cudaEventSynchronize(stop);
			float elapsed;
			cudaEventElapsedTime(&elapsed, start, stop); //gives in milliseconds
			elapsed /= 1000.0f; //convert to seconds
			float mpoints = (float)reportInterval * host.xmax * host.ymax * host.zmax / 1.0e6f;			
			printf("%d: %f Mpoints/s \n", i, mpoints / elapsed);

			if (i % exportInterval == 0) ExportCheckpoint("test", host, device1, i);
			cudaEventRecord(start, 0);
		}		
	}

	//FINISH
	gpuClose(device1, device2);
	hostClose(host);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
	return 0;
}