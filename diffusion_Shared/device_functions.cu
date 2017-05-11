#include "declarations.cuh"
#include "kernels.cu"

//kernel configuration variables
dim3 blocks_init;
dim3 threads_init;

dim3 blocks_laplace;	
dim3 threads_laplace;

dim3 blocks_boundary_x;
dim3 blocks_boundary_y;
dim3 blocks_boundary_z;
dim3 threads_boundary;

dim3 blocks_edge_x;
dim3 blocks_edge_y;
dim3 blocks_edge_z;
dim3 threads_edge;

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

	//Calculate kernel configurations
	blocks_init.x = host.xmax/BLOCKSIZE;
	blocks_init.y = host.ymax/BLOCKSIZE;
	threads_init.x = BLOCKSIZE;
	threads_init.y = BLOCKSIZE;

	blocks_laplace.x = (int)ceil((float)host.xmax/(CELLW-2));	
	blocks_laplace.y = (int)ceil((float)host.ymax/(CELLH-2));	
	threads_laplace.x = CELLW;
	threads_laplace.y = CELLH;
	threads_laplace.z = CELLD;

	blocks_boundary_x.x = host.ymax/BLOCKSIZE;
	blocks_boundary_x.y = host.zmax/BLOCKSIZE;

	blocks_boundary_y.x = host.xmax/BLOCKSIZE;
	blocks_boundary_y.y = host.zmax/BLOCKSIZE;

	blocks_boundary_z.x = host.xmax/BLOCKSIZE;
	blocks_boundary_z.y = host.ymax/BLOCKSIZE;

	threads_boundary.x = BLOCKSIZE;
	threads_boundary.y = BLOCKSIZE;

	blocks_edge_x.x = host.xmax/BLOCKSIZE;
	blocks_edge_y.x = host.ymax/BLOCKSIZE;
	blocks_edge_z.x = host.zmax/BLOCKSIZE;
	threads_edge.x = BLOCKSIZE;

	//allocate device arrays
	int size = total * sizeof(float);	
	cudaMalloc((void**)&device1.a, size);
	cudaMalloc((void**)&device2.a, size);
	CudaCheck();

	//copy data to constant memory on device
	cudaMemcpyToSymbol(DT, &id.dt, sizeof(float));
	cudaMemcpyToSymbol(DX, &id.dx, sizeof(float));
	cudaMemcpyToSymbol(DA, &id.da, sizeof(float));
	cudaMemcpyToSymbol(A_NULL, &id.anull, sizeof(float));
	cudaMemcpyToSymbol(ITMAX, &id.itmax, sizeof(int));

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

	kernelInitData<<<blocks_init, threads_init>>>(device1.a, id.radius);
	cudaThreadSynchronize();
	CudaCheck();
	
}

//compute one time step
void gpuStep(DataArray& device1, DataArray& device2)
{
	kernelBoundaryX<<<blocks_boundary_x, threads_boundary>>>(device1.a);
	kernelBoundaryY<<<blocks_boundary_y, threads_boundary>>>(device1.a);
	kernelBoundaryZ<<<blocks_boundary_z, threads_boundary>>>(device1.a);		
	kernelBoundaryEdgeX<<<blocks_edge_x, threads_edge>>>(device1.a);
	kernelBoundaryEdgeY<<<blocks_edge_y, threads_edge>>>(device1.a);
	kernelBoundaryEdgeZ<<<blocks_edge_z, threads_edge>>>(device1.a);

	kernelDiffusion<<<blocks_laplace, threads_laplace>>>(device1.a, device2.a);	

	//swap
	float *temp=device2.a; device2.a=device1.a; device1.a=temp;

	CudaCheck(); //because of asynchronous execution, there will be several steps before it can return an error, if any
}

//free device memory
void gpuClose(DataArray& device1, DataArray& device2)
{
	cudaFree(device1.a);
	cudaFree(device2.a);
}

//copy data from device to host and write it to a binary file
void ExportCheckpoint(const char* name, DataArray& host, DataArray& device1, int l)
{	
	cudaMemcpy(host.a, device1.a, total*sizeof(float), cudaMemcpyDeviceToHost);
	CudaCheck();
	
	//create a file
	char filename[4096];
	FILE *f;	
	sprintf(filename, "%s_%d.bin", name, l);
	f = fopen(filename, "wb");
	if (f==0) { printf("  error creating %s\n", filename); return; }

	//create 12-byte header(3x 4byte-integers): write dimensions
	fwrite(&host.xmax, sizeof(host.xmax), 1, f);
	fwrite(&host.ymax, sizeof(host.ymax), 1, f);
	fwrite(&host.zmax, sizeof(host.zmax), 1, f);

	//write data
	fwrite(host.a, total * sizeof(float), 1, f);

	fclose(f);
}
