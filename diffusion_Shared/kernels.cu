//Constants for the simulation defined in device Constant memory
__constant__ float DT, DX, DA, A_NULL, RADIUS; //same as in InitialData structure
__constant__ int XMAX, YMAX, ZMAX, ITMAX; //Dimensions of simulated space and time interval
__constant__ int Zsize, Ysize; //Helps to address memory.
__constant__ float DaDt_6DxDx; //Precomputed: DA * DT / (6 * DX * DX)

//Data initialization. This 2D kernel loops over Z and sets initial values
__global__ void kernelInitData(float *a_old, int radius)
{
	int x = blockIdx.x * BLOCKSIZE + threadIdx.x;
	int y = blockIdx.y * BLOCKSIZE + threadIdx.y;
	int z;

	for (z = 0; z < ZMAX; z++)
	{
		if (P2(x-XMAX/2) + P2(y-YMAX/2) + P2(z-ZMAX/2) < radius * radius)
		{
			a_old[z*Zsize+y*Ysize+x]=A_NULL;
		} else
		{
			a_old[z*Zsize+y*Ysize+x]=0.01f;
		};
	}
}

//Noflux boundary, perpendicular to X axis (YZ sides of space)
__global__ void kernelBoundaryX(float *a_old)
{
	int y=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int z=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int left_dest, left_src, right_dest, right_src;

	left_dest=z*Zsize+y*Ysize;
	left_src=z*Zsize+y*Ysize+1;	
	right_dest=z*Zsize+y*Ysize+(XMAX-1);
	right_src=z*Zsize+y*Ysize+(XMAX-2);

	a_old[left_dest]=a_old[left_src]; 
	a_old[right_dest]=a_old[right_src];
}

//Noflux boundary, perpendicular to Y axis (XZ sides of space)
__global__ void kernelBoundaryY(float *a_old)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int z=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int left_dest, left_src, right_dest, right_src;

	left_dest=z*Zsize+x;
	left_src=z*Zsize+Ysize+x;
	right_dest=z*Zsize+(YMAX-1)*Ysize+x;
	right_src=z*Zsize+(YMAX-2)*Ysize+x;

	a_old[left_dest]=a_old[left_src]; 
	a_old[right_dest]=a_old[right_src];
}

//Noflux boundary, perpendicular to Z axis (XY sides of space)
__global__ void kernelBoundaryZ(float *a_old)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int y=blockIdx.y*BLOCKSIZE+threadIdx.y;
	int left_dest, left_src, right_dest, right_src;

	left_dest=y*Ysize+x;
	left_src=Zsize+y*Ysize+x;
	right_dest=(ZMAX-1)*Zsize+y*Ysize+x;
	right_src=(ZMAX-2)*Zsize+y*Ysize+x;

	a_old[left_dest]=a_old[left_src]; 
	a_old[right_dest]=a_old[right_src];
}

//Noflux boundary on edges of space, parallel to X axis
__global__ void kernelBoundaryEdgeX(float *a_old)
{
	int x=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int dest_1, dest_2, dest_3, dest_4;
	int src_1, src_2, src_3, src_4;

	dest_1=x;
	src_1=Zsize+Ysize+x;
	dest_2=(YMAX-1)*Ysize+x;
	src_2=Zsize+(YMAX-2)*Ysize+x;
	dest_3=(ZMAX-1)*Zsize+(YMAX-1)*Ysize+x;
	src_3=(ZMAX-2)*Zsize+(YMAX-2)*Ysize+x;
	dest_4=(ZMAX-1)*Zsize+x;
	src_4=(ZMAX-2)*Zsize+Ysize+x;

	a_old[dest_1]=a_old[src_1];
	a_old[dest_2]=a_old[src_2];
	a_old[dest_3]=a_old[src_3];
	a_old[dest_4]=a_old[src_4];
}

//Noflux boundary on edges of space, parallel to Y axis
__global__ void kernelBoundaryEdgeY(float *a_old)
{
	int y=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int dest_1, dest_2, dest_3, dest_4;
	int src_1, src_2, src_3, src_4;

	dest_1=y*Ysize;
	src_1=Zsize+y*Ysize+1;
	dest_2=y*Ysize+(XMAX-1);
	src_2=Zsize+y*Ysize+(XMAX-2);
	dest_3=(ZMAX-1)*Zsize+y*Ysize+(XMAX-1);
	src_3=(ZMAX-2)*Zsize+y*Ysize+XMAX-2;
	dest_4=(ZMAX-1)*Zsize+y*Ysize;
	src_4=(ZMAX-2)*Zsize+y*Ysize+1;

	a_old[dest_1]=a_old[src_1];
	a_old[dest_2]=a_old[src_2];
	a_old[dest_3]=a_old[src_3];
	a_old[dest_4]=a_old[src_4];
}

//Noflux boundary on edges of space, parallel to Z axis
__global__ void kernelBoundaryEdgeZ(float *a_old)
{
	int z=blockIdx.x*BLOCKSIZE+threadIdx.x;
	int dest_1, dest_2, dest_3, dest_4;
	int src_1, src_2, src_3, src_4;

	dest_1=z*Zsize;
	src_1=z*Zsize+1*Ysize+1;
	dest_2=z*Zsize+(XMAX-1);
	src_2=z*Zsize+Ysize+(XMAX-2);
	dest_3=z*Zsize+(YMAX-1)*Ysize+(XMAX-1);
	src_3=z*Zsize+(YMAX-2)*Ysize+(XMAX-2);
	dest_4=z*Zsize+(YMAX-1)*Ysize;
	src_4=z*Zsize+(YMAX-2)*Ysize+1;

	a_old[dest_1]=a_old[src_1];
	a_old[dest_2]=a_old[src_2];
	a_old[dest_3]=a_old[src_3];
	a_old[dest_4]=a_old[src_4];
}

//main compute kernel to update concentration values on the grid
__global__ void kernelDiffusion(float *a_old, float *a_new)
{
	__shared__ float ao[CELLD][CELLH][CELLW];
	
	unsigned int p, k, x, y, z, kmax;
	bool ok_read, ok_compute;
	float laplace_a;

	//location and memory address. Note: neighbouring blocks overlap.
	x = blockIdx.x*(CELLW-2)+threadIdx.x;
	y = blockIdx.y*(CELLH-2)+threadIdx.y;
	z = threadIdx.z;
	p = z * Zsize + y * Ysize + x;	

	//precompute some conditions
	ok_read = (y<YMAX) & (x<XMAX);
	ok_compute = (threadIdx.y>0) & (threadIdx.y<CELLH-1) & (threadIdx.x>0) & (threadIdx.x<CELLW-1) & (threadIdx.z > 0) & (threadIdx.z < CELLD-1) & (y<YMAX-1) & (x<XMAX-1);

	//read first block-layer of data
	if (ok_read)
	{
		ao[threadIdx.z][threadIdx.y][threadIdx.x] = a_old[p];
	};

	__syncthreads();

	//loop down along z
	kmax = (int)ceilf((float)(ZMAX-2)/(CELLD-2));
	for (k=0; k < kmax; k++)
	{
		z = k * (CELLD-2) + threadIdx.z;
		p = z * Zsize + y * Ysize + x;		

		//calculate laplacian with 19-point stencil
		if (ok_compute & (z < ZMAX-1))
		{
			laplace_a =
			     ao[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
				 ao[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
				 ao[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
				 ao[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
				 ao[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
				 ao[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
				 ao[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
				 ao[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
				 ao[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
				 ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] +
				 ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
				 ao[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
			2.0f * ( ao[threadIdx.z-1][threadIdx.y][threadIdx.x] +
			     ao[threadIdx.z+1][threadIdx.y][threadIdx.x] +
				 ao[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
				 ao[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
				 ao[threadIdx.z][threadIdx.y][threadIdx.x-1] +
				 ao[threadIdx.z][threadIdx.y][threadIdx.x+1] 
			) -24.0f * ao[threadIdx.z][threadIdx.y][threadIdx.x]; 
			
			//write updated values to global memory
			a_new[p] = ao[threadIdx.z][threadIdx.y][threadIdx.x] + DaDt_6DxDx * laplace_a;
		};

		__syncthreads();

		//copy last two z-layers in the block to the first position
		if (threadIdx.z >= CELLD-2)
		{
			ao[threadIdx.z-CELLD+2][threadIdx.y][threadIdx.x] = ao[threadIdx.z][threadIdx.y][threadIdx.x];			
		};

		//No need to syncthreads() here because those warps that write shared mem above
		//are exactly the same warps that will load data to the same place afterwards.
		//Warps always execute code sequentially.
		//If block width*height is multiple of 32 then warps involved are not divergent.

		//read new z block-layer
		z = (k+1) * (CELLD-2) + threadIdx.z;
		p = z * Zsize + y * Ysize + x;
		if ((z<ZMAX) & ok_read & (threadIdx.z >= 2))
		{
			ao[threadIdx.z][threadIdx.y][threadIdx.x] = a_old[p];			
		}

		__syncthreads();
	}
}
