//Constants for the simulation defined in device Constant memory
__constant__ float DT, DX, DA, A_NULL, RADIUS; //same as in InitialData structure
__constant__ int XMAX, YMAX, ZMAX, ITMAX; //Dimensions of simulated space and time interval
__constant__ int Zsize, Ysize; //Helps to address memory.
__constant__ float DaDt_6DxDx; //Precomputed: DA * DT / (6 * DX * DX)

// BCs 

// perpendicular to the X axis (YZ sides of space)
__global__ void kernelPeriodicBCX(float *a_old, float *a_new)
{
    __shared__ float ao[BOUNDARYBLOCKSIZE][BOUNDARYBLOCKSIZE][4];

    unsigned int p, k, x, y, z, kmax;
    bool ok_read, ok_compute;
    float laplace_a;

    // recall that, irritatingly, block id's are always 2d, in x and y. since we are doing the yz boundary, but we have to label are blocks xy,
    // we initially define "real space" block id's here, mapping "block space" x -> real space y , "block space" y -> real space z

    int realspaceblockIdy = blockIdx.x;
    int realspaceblockIdz = blockIdx.y;

    // beyond this, our convention is as follows: LOCK THE THREAD ID'S TO REAL SPACE

    //location and memory address. Note: neighbouring blocks overlap.
    y=realspaceblockIdy*(BOUNDARYBLOCKSIZE-2)+threadIdx.y;
    z=realspaceblockIdz*(BOUNDARYBLOCKSIZE-2)+threadIdx.z;
    // ok, x runs as follows: if threadidx.z = 0, x=N-2, idx.z = 1. x = N-1, idx.z = 2, x = 0, idx.z = 3, x = 1
    // ie we move right from x = N-2, go over the boundary, and continue to x = 1
    x=(threadIdx.x+(XMAX-2))%XMAX;

    p = z * Zsize + y * Ysize + x;	

    //precompute some conditions
    ok_read = (y<YMAX) & (z<ZMAX);
    ok_compute = (threadIdx.y>0) & (threadIdx.y<BOUNDARYBLOCKSIZE-1) & (threadIdx.x>0) & (threadIdx.x<4-1) & (threadIdx.z > 0) & (threadIdx.z < BOUNDARYBLOCKSIZE-1) & (y<YMAX-1) & (z<ZMAX-1);

    //read first block-layer of data
    if (ok_read)
    {
        ao[threadIdx.z][threadIdx.y][threadIdx.x] = a_old[p];
    };

    __syncthreads();


    //calculate laplacian with 19-point stencil
    if (ok_compute )
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
            ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
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

}

// perpendicular to the Y axis (XZ sides of space)
__global__ void kernelPeriodicBCY(float *a_old, float *a_new)
{
    __shared__ float ao[BOUNDARYBLOCKSIZE][4][BOUNDARYBLOCKSIZE];

    unsigned int p, k, x, y, z, kmax;
    bool ok_read, ok_compute;
    float laplace_a;

    // recall that, irritatingly, block id's are always 2d, in x and y. since we are doing the yz boundary, but we have to label are blocks xy,
    // we initially define "real space" block id's here, mapping "block space" x -> real space y , "block space" y -> real space z

    int realspaceblockIdx = blockIdx.x;
    int realspaceblockIdz = blockIdx.y;

    // beyond this, our convention is as follows: LOCK THE THREAD ID'S TO REAL SPACE

    //location and memory address. Note: neighbouring blocks overlap.
    x=realspaceblockIdx*(BOUNDARYBLOCKSIZE-2)+threadIdx.x;
    z=realspaceblockIdz*(BOUNDARYBLOCKSIZE-2)+threadIdx.z;
    // ok, y runs as follows: if threadidx.y = 0, y=N-2, idx.y = 1. x = N-1, idx.y = 2, x = 0, idx.y = 3, y = 1
    // ie we move right from y = N-2, go over the boundary, and continue to y = 1
    y=(threadIdx.y+(YMAX-2))%YMAX;

    p = z * Zsize + y * Ysize + x;	

    //precompute some conditions
    ok_read = (x<XMAX) & (z<ZMAX);
    ok_compute = (threadIdx.y>0) & (threadIdx.y<4-1) & (threadIdx.x>0) & (threadIdx.x<BOUNDARYBLOCKSIZE -1) & (threadIdx.z > 0) & (threadIdx.z < BOUNDARYBLOCKSIZE-1) & (x<XMAX-1) & (z<ZMAX-1);

    //read first block-layer of data
    if (ok_read)
    {
        ao[threadIdx.z][threadIdx.y][threadIdx.x] = a_old[p];
    };

    __syncthreads();


    //calculate laplacian with 19-point stencil
    if (ok_compute )
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
            ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
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

}

// perpendicular to the Z axis (XY sides of space)
__global__ void kernelPeriodicBCZ(float *a_old, float *a_new)
{
    __shared__ float ao[4][BOUNDARYBLOCKSIZE][BOUNDARYBLOCKSIZE];

    unsigned int p, k, x, y, z, kmax;
    bool ok_read, ok_compute;
    float laplace_a;

    // recall that, irritatingly, block id's are always 2d, in x and y. since we are doing the yz boundary, but we have to label are blocks xy,
    // we initially define "real space" block id's here, mapping "block space" x -> real space y , "block space" y -> real space z

    int realspaceblockIdx = blockIdx.x;
    int realspaceblockIdy = blockIdx.y;

    // beyond this, our convention is as follows: LOCK THE THREAD ID'S TO REAL SPACE

    //location and memory address. Note: neighbouring blocks overlap.
    x=realspaceblockIdx*(BOUNDARYBLOCKSIZE-2)+threadIdx.x;
    y=realspaceblockIdy*(BOUNDARYBLOCKSIZE-2)+threadIdx.y;
    // ok, y runs as follows: if threadidx.y = 0, y=N-2, idx.y = 1. x = N-1, idx.y = 2, x = 0, idx.y = 3, y = 1
    // ie we move right from y = N-2, go over the boundary, and continue to y = 1
    z=(threadIdx.z+(ZMAX-2))%ZMAX;

    p = z * Zsize + y * Ysize + x;	

    //precompute some conditions
    ok_read = (x<XMAX) & (y<YMAX);
    ok_compute = (threadIdx.y>0) & (threadIdx.y<BOUNDARYBLOCKSIZE -1) & (threadIdx.x>0) & (threadIdx.x<BOUNDARYBLOCKSIZE -1) & (threadIdx.z > 0) & (threadIdx.z < 4-1) & (x<XMAX-1) & (y<YMAX-1);

    //read first block-layer of data
    if (ok_read)
    {
        ao[threadIdx.z][threadIdx.y][threadIdx.x] = a_old[p];
    };

    __syncthreads();


    //calculate laplacian with 19-point stencil
    if (ok_compute )
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
            ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
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

}

// Edge BC's

//periodic boundary of the edges of space, parallel to X axis
// enough smarts, lets just do this the dumb way
__global__ void kernelPeriodicBCEdgeX(float *a_old,float *a_new)
{
    int x=blockIdx.x*BLOCKSIZE+threadIdx.x;

    for(int y=0, y<YMAX, y += YMAX)
    {
        for(int z=0, z<ZMAX, z += ZMAX)
        {

            // the 7 point stencil
            x_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            xup_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xdown_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            x_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	

            // get the 4 extras around xup_y_z
            xup_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xup_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xup_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xup_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	

            // get the 4 extras around xdown_y_z
            xdown_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            xdown_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            xdown_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            xdown_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	

            // get the remaining guys
            x_yup_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_yup_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_ydown_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_ydown_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	

            //calculate laplacian with 19-point stencil
            laplace_a =
                aold[x_ydown_zup]+
                //       ao[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
                aold[xdown_y_zup]+
                //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
                aold[xup_y_zup]+
                //      ao[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
                aold[x_yup_zup]+
                //       ao[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
                aold[xdown_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
                aold[xup_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
                aold[xdown_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
                aold[xup_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
                aold[x_ydown_zdown]+
                //       ao[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
                aold[xdown_y_zdown] + aold[xup_y_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
                aold[x_yup_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
                2.0f* (aold[x_y_zdown]+
                        //  2.0f * ( ao[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                aold[x_y_zup]+
                //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                aold[x_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                aold[x_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                aold[xdown_y_z]+
                //       ao[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                aold[xup_y_z]
                //       ao[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                ) - 24.0f* aold[x_y_z];
                        //  ) -24.0f * ao[threadIdx.z][threadIdx.y][threadIdx.x]; 

            //write updated values to global memory
            a_new[x_y_z] = aold[x_y_z] + DaDt_6DxDx * laplace_a;
        }
    }
}

//periodic boundary of the edges of space, parallel to Y axis
// enough smarts, lets just do this the dumb way
__global__ void kernelPeriodicBCBoundaryEdgeY(float *a_old,float *a_new)
{
    int y=blockIdx.x*BLOCKSIZE+threadIdx.x;

    for(int x=0, x<XMAX, x += XMAX)
    {
        for(int z=0, z<ZMAX, z += ZMAX)
        {

            // the 7 point stencil
            x_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            xup_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xdown_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            x_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	

            // get the 4 extras around xup_y_z
            xup_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xup_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xup_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xup_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	

            // get the 4 extras around xdown_y_z
            xdown_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            xdown_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            xdown_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            xdown_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	

            // get the remaining guys
            x_yup_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_yup_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_ydown_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_ydown_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	

            //calculate laplacian with 19-point stencil
            laplace_a =
                aold[x_ydown_zup]+
                //       ao[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
                aold[xdown_y_zup]+
                //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
                aold[xup_y_zup]+
                //      ao[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
                aold[x_yup_zup]+
                //       ao[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
                aold[xdown_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
                aold[xup_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
                aold[xdown_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
                aold[xup_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
                aold[x_ydown_zdown]+
                //       ao[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
                aold[xdown_y_zdown] + aold[xup_y_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
                aold[x_yup_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
                2.0f* (aold[x_y_zdown]+
                        //  2.0f * ( ao[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                aold[x_y_zup]+
                //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                aold[x_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                aold[x_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                aold[xdown_y_z]+
                //       ao[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                aold[xup_y_z]
                //       ao[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                ) - 24.0f* aold[x_y_z];
                        //  ) -24.0f * ao[threadIdx.z][threadIdx.y][threadIdx.x]; 

            //write updated values to global memory
            a_new[x_y_z] = aold[x_y_z] + DaDt_6DxDx * laplace_a;
        }
    }
}

//periodic boundary of the edges of space, parallel to Z axis
// enough smarts, lets just do this the dumb way
__global__ void kernelPeriodicBCBoundaryEdgeZ(float *a_old,float *a_new)
{
    int z=blockIdx.x*BLOCKSIZE+threadIdx.x;

    for(int x=0, x<XMAX, x += XMAX)
    {
        for(int z=0, z<ZMAX, z += ZMAX)
        {
            // the 7 point stencil
            x_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            xup_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xdown_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            x_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	

            // get the 4 extras around xup_y_z
            xup_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xup_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xup_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            xup_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	

            // get the 4 extras around xdown_y_z
            xdown_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            xdown_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            xdown_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            xdown_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	

            // get the remaining guys
            x_yup_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_yup_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_ydown_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            x_ydown_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	

            //calculate laplacian with 19-point stencil
            laplace_a =
                aold[x_ydown_zup]+
                //       ao[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
                aold[xdown_y_zup]+
                //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
                aold[xup_y_zup]+
                //      ao[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
                aold[x_yup_zup]+
                //       ao[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
                aold[xdown_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
                aold[xup_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
                aold[xdown_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
                aold[xup_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
                aold[x_ydown_zdown]+
                //       ao[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
                aold[xdown_y_zdown] + aold[xup_y_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
                aold[x_yup_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
                2.0f* (aold[x_y_zdown]+
                        //  2.0f * ( ao[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                aold[x_y_zup]+
                //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                aold[x_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                aold[x_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                aold[xdown_y_z]+
                //       ao[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                aold[xup_y_z]
                //       ao[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                ) - 24.0f* aold[x_y_z];
                        //  ) -24.0f * ao[threadIdx.z][threadIdx.y][threadIdx.x]; 

            //write updated values to global memory
            a_new[x_y_z] = aold[x_y_z] + DaDt_6DxDx * laplace_a;
        }
    }
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

// Two little helper functions to do the periodc bc's
__device__ int Gridinc(int p,int inc, int N)
{
    return(GridMod(p+inc,N));
}


// a little helper function to do the periodc bc's
__device__ int GridMod(int p, int N)
{
    int LessThanZero = (p<0);
    int ZeroOrAbove = 1 - LessThanZero;

    return ZeroOrAbove*(p%N) + LessThanZero*(N-p);
}
