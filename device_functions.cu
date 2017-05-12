#include "declarations.cuh"
#include "kernels.cu"

//kernel configuration variables
// main compute blocks and threads
dim3 blocks_laplace;	
dim3 threads_laplace;
// face blocks and threads
dim3 blocks_boundary_X;
dim3 blocks_boundary_Y;
dim3 blocks_boundary_Z;
dim3 threads_boundary_X;
dim3 threads_boundary_Y;
dim3 threads_boundary_Z;
// edge blocks and threads
dim3 blocks_boundary_edge_X;
dim3 blocks_boundary__edge_Y;
dim3 blocks_boundary_edge_Z;
dim3 threads_boundary_edge_X;
dim3 threads_boundary_edge_Y;
dim3 threads_boundary_edge_Z;

//function to allocate device memory and initialize data
void gpuInitialize(InitialData& id, DataArray& host, DataArray& device1, DataArray& device2)
{
	if (host.xmax % BLOCKSIZE != 0 || host.ymax % BLOCKSIZE != 0 || host.zmax % BLOCKSIZE != 0)
	{
		char buf[1024];
		sprintf(buf, "All dimensions must be multiple of %d", BLOCKSIZE);
		Error(buf);
	}

	device1.xmax = host.xmax;
	device1.ymax = host.ymax;
	device1.zmax = host.zmax;

	device2.xmax = host.xmax;
	device2.ymax = host.ymax;
	device2.zmax = host.zmax;

    // Block and threads for the main compute laplacian
	blocks_laplace.x = (int)ceil((float)host.xmax/(CELLW-2));	
	blocks_laplace.y = (int)ceil((float)host.ymax/(CELLH-2));	
	threads_laplace.x = CELLW;
	threads_laplace.y = CELLH;
	threads_laplace.z = CELLD;

    // Block and threads for the face periodic bc's
	blocks_boundary_X.x = (int)ceil((float)host.ymax/(BOUNDARYBLOCKSIZE-2));	
	blocks_boundary_X.y = (int)ceil((float)host.zmax/(BOUNDARYBLOCKSIZE-2));	

	threads_boundary_X.x = 4;
	threads_boundary_X.y = BOUNDARYBLOCKSIZE;
	threads_boundary_X.z=BOUNDARYBLOCKSIZE;

	blocks_boundary_Y.x = (int)ceil((float)host.xmax/(BOUNDARYBLOCKSIZE-2));	
	blocks_boundary_Y.y = (int)ceil((float)host.zmax/(BOUNDARYBLOCKSIZE-2));	

	threads_boundary_Y.x = BOUNDARYBLOCKSIZE;
	threads_boundary_Y.y = 4;
	threads_boundary_Y.z=BOUNDARYBLOCKSIZE;

	blocks_boundary_Z.x = (int)ceil((float)host.xmax/(BOUNDARYBLOCKSIZE-2));	
	blocks_boundary_Z.y = (int)ceil((float)host.ymax/(BOUNDARYBLOCKSIZE-2));	

	threads_boundary_Z.x = BOUNDARYBLOCKSIZE;
	threads_boundary_Z.y = BOUNDARYBLOCKSIZE;
	threads_boundary_Z.z= 4;

    // Block and threads for the edge periodic bc's
	blocks_boundary_edge_X.x = host.xmax/BOUNDARYBLOCKSIZE;	
	blocks_boundary_edge_Y.x = host.ymax/BOUNDARYBLOCKSIZE;	
	blocks_boundary_edge_Z.x = host.zmax/BOUNDARYBLOCKSIZE;	

	threads_boundary_edge_X.x = BOUNDARYBLOCKSIZE;
	threads_boundary_edge_Y.x = BOUNDARYBLOCKSIZE;
	threads_boundary_edge_Z.x = BOUNDARYBLOCKSIZE;

	//allocate device arrays
    int total = host.xmax * host.ymax *host.zmax;
	int size = total * sizeof(float);	
	cudaMalloc((void**)&device1.u, size);
	cudaMalloc((void**)&device1.v, size);
	cudaMalloc((void**)&device2.u, size);
	cudaMalloc((void**)&device2.v, size);
	CudaCheck();
    // copy cpu memory to gpu memory
    cudaMemcpy(device1.u,host.u,size,cudaMemcpyHostToDevice);
    cudaMemcpy(device1.v,host.v,size,cudaMemcpyHostToDevice);

	//copy data to constant memory on device
	cudaMemcpyToSymbol(DT, &id.dt, sizeof(float));
	cudaMemcpyToSymbol(DX, &id.dx, sizeof(float));
	cudaMemcpyToSymbol(DA, &id.da, sizeof(float));

	cudaMemcpyToSymbol(XMAX, &host.xmax, sizeof(int));
	cudaMemcpyToSymbol(YMAX, &host.ymax, sizeof(int));
	cudaMemcpyToSymbol(ZMAX, &host.zmax, sizeof(int));

	cudaMemcpyToSymbol(Zsize, &zsize, sizeof(int));
	cudaMemcpyToSymbol(Ysize, &ysize, sizeof(int));

	CudaCheck();

	float t;
	t = id.da * id.dt/(6.0f*id.dx*id.dx);
	cudaMemcpyToSymbol(DaDt_6DxDx, &t, sizeof(float));
	CudaCheck();

	cudaThreadSynchronize();
	CudaCheck();
	
}

//compute one time step
void gpuStep(DataArray& device1, DataArray& device2)
{
    kernelPeriodicBCX<<<blocks_boundary_X, threads_boundary_X>>>(device1.u,device1.v,device2.u,device2.v);
    kernelPeriodicBCY<<<blocks_boundary_Y, threads_boundary_Y>>>(device1.u,device1.v,device2.u,device2.v);
    kernelPeriodicBCZ<<<blocks_boundary_Z, threads_boundary_Z>>>(device1.u,device1.v,device2.u,device2.v);
    kernelPeriodicBCEdgeX<blocks_boundary_edge_X,threads_boundary_edge_X>>>(device1.u,device1.v,device2.u,device2.v);
    kernelPeriodicBCEdgeY<blocks_boundary_edge_Y,threads_boundary_edge_Y>>>(device1.u,device1.v,device2.u,device2.v);
    kernelPeriodicBCEdgeZ<blocks_boundary_edge_Z,threads_boundary_edge_Z>>>(device1.u,device1.v,device2.u,device2.v);

	//swap
	float *temp;
    temp=device2.u; device2.u=device1.u; device1.u=temp;
	temp=device2.v; device2.v=device1.v; device1.v=temp;

	CudaCheck(); //because of asynchronous execution, there will be several steps before it can return an error, if any
}

//free device memory
void gpuClose(DataArray& device1, DataArray& device2)
{
	cudaFree(device1.u);
	cudaFree(device1.v);
	cudaFree(device2.u);
	cudaFree(device2.v);
}

