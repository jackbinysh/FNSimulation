//Constants for the simulation defined in device Constant memory
#include "declarations.cuh"
__constant__ gridprecision DT, DX, DA, EPSILON,BETA,GAMMA,ONEOVEREPSILON,ONETHIRD ; //same as in InitialData structure
__constant__ int XMAX, YMAX, ZMAX; //Dimensions of simulated space and time interval
__constant__ int Zsize, Ysize; //Helps to address memory.
__constant__ gridprecision DaDt_6DxDx; //Precomputed: DA * DT / (6 * DX * DX)

// Two little helper functions to do the periodc bc's

__device__ int GridMod(int p, int N)
{
    int LessThanZero = (p<0);
    int ZeroOrAbove = 1 - LessThanZero;

    return ZeroOrAbove*(p%N) + LessThanZero*(N-p);
}
__device__ int GridInc(int p,int inc, int N)
{
    return(GridMod(p+inc,N));
}

// BCs 

// perpendicular to the X axis (YZ sides of space)
__global__ void kernelPeriodicBCX(gridprecision *u_old,gridprecision *v_old, gridprecision *u_new, gridprecision *v_new)
{
    __shared__ gridprecision u_shared[BOUNDARYBLOCKSIZE][BOUNDARYBLOCKSIZE][4];
    __shared__ gridprecision v_shared[BOUNDARYBLOCKSIZE][BOUNDARYBLOCKSIZE][4];

    unsigned int p, x, y, z;
    bool ok_read, ok_compute;
    gridprecision laplace_u,u,v;

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
        u_shared[threadIdx.z][threadIdx.y][threadIdx.x] = u_old[p];
        v_shared[threadIdx.z][threadIdx.y][threadIdx.x] = v_old[p];
    };

    __syncthreads();


    //calculate laplacian with 19-point stencil
    if (ok_compute )
    {
        laplace_u =
            u_shared[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
            u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
            u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
            u_shared[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
            u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
            u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
            u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
            u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
            u_shared[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
            u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
            u_shared[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
            2.0f * ( u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                    u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                    u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                    u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                    u_shared[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                    u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                   ) -24.0f * u_shared[threadIdx.z][threadIdx.y][threadIdx.x]; 

        //write updated values to global memory
        u =  u_shared[threadIdx.z][threadIdx.y][threadIdx.x] ;
        v =  v_shared[threadIdx.z][threadIdx.y][threadIdx.x] ;
        u_new[p] = u + (DaDt_6DxDx * laplace_u)+ ONEOVEREPSILON*(u - ONETHIRD*u*u*u - v)*DT;
        v_new[p] = v +EPSILON*(u + BETA - GAMMA*v)*DT;
    };

}

// perpendicular to the Y axis (XZ sides of space)
__global__ void kernelPeriodicBCY(gridprecision *u_old,gridprecision *v_old, gridprecision *u_new, gridprecision *v_new)
{
    __shared__ gridprecision u_shared[BOUNDARYBLOCKSIZE][4][BOUNDARYBLOCKSIZE];
    __shared__ gridprecision v_shared[BOUNDARYBLOCKSIZE][4][BOUNDARYBLOCKSIZE];

    unsigned int p,  x, y, z;
    bool ok_read, ok_compute;
    gridprecision laplace_u,u,v;

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
        u_shared[threadIdx.z][threadIdx.y][threadIdx.x] = u_old[p];
        v_shared[threadIdx.z][threadIdx.y][threadIdx.x] = v_old[p];
    };

    __syncthreads();


    //calculate laplacian with 19-point stencil
    if (ok_compute )
    {
        laplace_u =
            u_shared[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
            u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
            u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
            u_shared[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
            u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
            u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
            u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
            u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
            u_shared[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
            u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
            u_shared[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
            2.0f * ( u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                    u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                    u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                    u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                    u_shared[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                    u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                   ) -24.0f * u_shared[threadIdx.z][threadIdx.y][threadIdx.x]; 

        //write updated values to global memory
        u =  u_shared[threadIdx.z][threadIdx.y][threadIdx.x] ;
        v =  v_shared[threadIdx.z][threadIdx.y][threadIdx.x] ;
        u_new[p] = u + (DaDt_6DxDx * laplace_u)+ ONEOVEREPSILON*(u - ONETHIRD*u*u*u - v)*DT;
        v_new[p] = v +EPSILON*(u + BETA - GAMMA*v)*DT;
    };

}

// perpendicular to the Z axis (XY sides of space)
__global__ void kernelPeriodicBCZ(gridprecision *u_old,gridprecision *v_old, gridprecision *u_new, gridprecision *v_new)
{
    __shared__ gridprecision u_shared[4][BOUNDARYBLOCKSIZE][BOUNDARYBLOCKSIZE];
    __shared__ gridprecision v_shared[4][BOUNDARYBLOCKSIZE][BOUNDARYBLOCKSIZE];

    unsigned int p, x, y, z;
    bool ok_read, ok_compute;
    gridprecision laplace_u,u,v;

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
        u_shared[threadIdx.z][threadIdx.y][threadIdx.x] = u_old[p];
        v_shared[threadIdx.z][threadIdx.y][threadIdx.x] = v_old[p];
    };

    __syncthreads();


    //calculate laplacian with 19-point stencil
    if (ok_compute )
    {
        laplace_u =
            u_shared[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
            u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
            u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
            u_shared[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
            u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
            u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
            u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
            u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
            u_shared[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
            u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
            u_shared[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
            2.0f * ( u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                    u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                    u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                    u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                    u_shared[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                    u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                   ) -24.0f * u_shared[threadIdx.z][threadIdx.y][threadIdx.x]; 

        //write updated values to global memory
        u =  u_shared[threadIdx.z][threadIdx.y][threadIdx.x] ;
        v =  v_shared[threadIdx.z][threadIdx.y][threadIdx.x] ;
        u_new[p] = u + (DaDt_6DxDx * laplace_u)+ ONEOVEREPSILON*(u - ONETHIRD*u*u*u - v)*DT;
        v_new[p] = v +EPSILON*(u + BETA - GAMMA*v)*DT;
    };

}

// Edge BC's

//periodic boundary of the edges of space, parallel to X axis
// enough smarts, lets just do this the dumb way
__global__ void kernelPeriodicBCEdgeX(gridprecision *u_old,gridprecision *v_old, gridprecision *u_new, gridprecision *v_new)
{
    int x=blockIdx.x*BOUNDARYBLOCKSIZE+threadIdx.x;
    gridprecision laplace_u,u,v;

    for(int y=0; y<YMAX; y += YMAX)
    {
        for(int z=0; z<ZMAX; z += ZMAX)
        {

            // the 7 point stencil
            int x_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int xup_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xdown_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int x_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	

            // get the 4 extras around xup_y_z
            int xup_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xup_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xup_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xup_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	

            // get the 4 extras around xdown_y_z
            int xdown_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int xdown_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int xdown_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int xdown_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	

            // get the remaining guys
            int x_yup_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_yup_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_ydown_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_ydown_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	

            //calculate laplacian with 19-point stencil
            laplace_u =
                u_old[x_ydown_zup]+
                //       ao[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
                u_old[xdown_y_zup]+
                //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
                u_old[xup_y_zup]+
                //      ao[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
                u_old[x_yup_zup]+
                //       ao[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
                u_old[xdown_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
                u_old[xup_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
                u_old[xdown_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
                u_old[xup_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
                u_old[x_ydown_zdown]+
                //       ao[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
                u_old[xdown_y_zdown] + u_old[xup_y_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
                u_old[x_yup_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
                2.0f* (u_old[x_y_zdown]+
                        //  2.0f * ( ao[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                        u_old[x_y_zup]+
                        //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                        u_old[x_ydown_z]+
                        //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                        u_old[x_yup_z]+
                        //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                        u_old[xdown_y_z]+
                        //       ao[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                        u_old[xup_y_z]
                        //       ao[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                      ) - 24.0f* u_old[x_y_z];
                      //  ) -24.0f * ao[threadIdx.z][threadIdx.y][threadIdx.x]; 

                      //write updated values to global memory
            u =  u_old[x_y_z] ;
            v= v_old[x_y_z] ;
            u_new[x_y_z] = u + (DaDt_6DxDx * laplace_u)+ ONEOVEREPSILON*(u - ONETHIRD*u*u*u - v)*DT;
            v_new[x_y_z] = v +EPSILON*(u + BETA - GAMMA*v)*DT;
        }
    }
}

//periodic boundary of the edges of space, parallel to Y axis
// enough smarts, lets just do this the dumb way
__global__ void kernelPeriodicBCEdgeY(gridprecision *u_old,gridprecision *v_old, gridprecision *u_new, gridprecision *v_new)
{
    int y=blockIdx.x*BOUNDARYBLOCKSIZE+threadIdx.x;
    gridprecision laplace_u,u,v;

    for(int x=0; x<XMAX; x += XMAX)
    {
        for(int z=0; z<ZMAX; z += ZMAX)
        {

            // the 7 point stencil
            int x_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int xup_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xdown_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int x_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	

            // get the 4 extras around xup_y_z
            int xup_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xup_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xup_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xup_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	

            // get the 4 extras around xdown_y_z
            int xdown_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int xdown_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int xdown_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int xdown_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	

            // get the remaining guys
            int x_yup_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_yup_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_ydown_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_ydown_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	

            //calculate laplacian with 19-point stencil
            laplace_u =
                u_old[x_ydown_zup]+
                //       ao[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
                u_old[xdown_y_zup]+
                //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
                u_old[xup_y_zup]+
                //      ao[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
                u_old[x_yup_zup]+
                //       ao[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
                u_old[xdown_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
                u_old[xup_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
                u_old[xdown_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
                u_old[xup_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
                u_old[x_ydown_zdown]+
                //       ao[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
                u_old[xdown_y_zdown] + u_old[xup_y_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
                u_old[x_yup_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
                2.0f* (u_old[x_y_zdown]+
                        //  2.0f * ( ao[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                        u_old[x_y_zup]+
                        //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                        u_old[x_ydown_z]+
                        //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                        u_old[x_yup_z]+
                        //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                        u_old[xdown_y_z]+
                        //       ao[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                        u_old[xup_y_z]
                        //       ao[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                      ) - 24.0f* u_old[x_y_z];
                      //  ) -24.0f * ao[threadIdx.z][threadIdx.y][threadIdx.x]; 

                      //write updated values to global memory
            u =  u_old[x_y_z] ;
            v= v_old[x_y_z] ;
            u_new[x_y_z] = u + (DaDt_6DxDx * laplace_u)+ ONEOVEREPSILON*(u - ONETHIRD*u*u*u - v)*DT;
            v_new[x_y_z] = v +EPSILON*(u + BETA - GAMMA*v)*DT;
        }
    }
}

//periodic boundary of the edges of space, parallel to Z axis
// enough smarts, lets just do this the dumb way
__global__ void kernelPeriodicBCEdgeZ(gridprecision *u_old,gridprecision *v_old, gridprecision *u_new, gridprecision *v_new)
{
    int z=blockIdx.x*BOUNDARYBLOCKSIZE+threadIdx.x;
    gridprecision laplace_u,u,v;

    for(int x=0; x<XMAX; x += XMAX)
    {
        for(int y=0; y<YMAX; y += YMAX)
        {
            // the 7 point stencil
            int x_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int xup_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xdown_y_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int x_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,0,XMAX);	

            // get the 4 extras around xup_y_z
            int xup_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xup_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xup_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	
            int xup_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,1,XMAX);	

            // get the 4 extras around xdown_y_z
            int xdown_yup_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int xdown_ydown_z = GridInc(z,0,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int xdown_y_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	
            int xdown_y_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,0,YMAX) * Ysize + GridInc(x,-1,XMAX);	

            // get the remaining guys
            int x_yup_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_yup_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_ydown_zup = GridInc(z,1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	
            int x_ydown_zdown = GridInc(z,-1,ZMAX) * Zsize + GridInc(y,-1,YMAX) * Ysize + GridInc(x,0,XMAX);	

            //calculate laplacian with 19-point stencil
            laplace_u =
                u_old[x_ydown_zup]+
                //       ao[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
                u_old[xdown_y_zup]+
                //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
                u_old[xup_y_zup]+
                //      ao[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
                u_old[x_yup_zup]+
                //       ao[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
                u_old[xdown_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
                u_old[xup_ydown_z]+
                //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
                u_old[xdown_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
                u_old[xup_yup_z]+
                //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
                u_old[x_ydown_zdown]+
                //       ao[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
                u_old[xdown_y_zdown] + u_old[xup_y_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y][threadIdx.x-1] + ao[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
                u_old[x_yup_zdown]+
                //      ao[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
                2.0f* (u_old[x_y_zdown]+
                        //  2.0f * ( ao[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                        u_old[x_y_zup]+
                        //       ao[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                        u_old[x_ydown_z]+
                        //       ao[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                        u_old[x_yup_z]+
                        //       ao[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                        u_old[xdown_y_z]+
                        //       ao[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                        u_old[xup_y_z]
                        //       ao[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                      ) - 24.0f* u_old[x_y_z];
                      //  ) -24.0f * ao[threadIdx.z][threadIdx.y][threadIdx.x]; 

                      //write updated values to global memory
            u =  u_old[x_y_z] ;
            v= v_old[x_y_z] ;
            u_new[x_y_z] = u + (DaDt_6DxDx * laplace_u)+ ONEOVEREPSILON*(u - ONETHIRD*u*u*u - v)*DT;
            v_new[x_y_z] = v +EPSILON*(u + BETA - GAMMA*v)*DT;
        }
    }
}

//main compute kernel to update concentration values on the grid
__global__ void kernelDiffusion_Shared(gridprecision *u_old,gridprecision *v_old, gridprecision *u_new, gridprecision *v_new)
{
    __shared__ gridprecision u_shared[CELLD][CELLH][CELLW];
    __shared__ gridprecision v_shared[CELLD][CELLH][CELLW];

    unsigned int p, k, x, y, z, kmax;
    bool ok_read, ok_compute;
    gridprecision laplace_u,u,v;

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
        u_shared[threadIdx.z][threadIdx.y][threadIdx.x] = u_old[p];
        v_shared[threadIdx.z][threadIdx.y][threadIdx.x] = v_old[p];
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
            laplace_u =
                u_shared[threadIdx.z+1][threadIdx.y-1][threadIdx.x] + 
                u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x-1] +
                u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
                u_shared[threadIdx.z+1][threadIdx.y+1][threadIdx.x] + 
                u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x-1] +
                u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
                u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x-1] +
                u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
                u_shared[threadIdx.z-1][threadIdx.y-1][threadIdx.x] + 
                u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x-1] +
                u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
                u_shared[threadIdx.z-1][threadIdx.y+1][threadIdx.x] + 
                2.0f * ( u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                        u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                        u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x] + 
                        u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x] + 
                        u_shared[threadIdx.z][threadIdx.y][threadIdx.x-1] +
                        u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1] 
                       ) -24.0f * u_shared[threadIdx.z][threadIdx.y][threadIdx.x]; 

            //write updated values to global memory
            u =  u_shared[threadIdx.z][threadIdx.y][threadIdx.x] ;
            v =  v_shared[threadIdx.z][threadIdx.y][threadIdx.x] ;
            u_new[p] = u + (DaDt_6DxDx * laplace_u)+ ONEOVEREPSILON*(u - ONETHIRD*u*u*u - v)*DT;
            v_new[p] = v +EPSILON*(u + BETA - GAMMA*v)*DT;
        };

        __syncthreads();

        //copy last two z-layers in the block to the first position
        if (threadIdx.z >= CELLD-2)
        {
            u_shared[threadIdx.z-CELLD+2][threadIdx.y][threadIdx.x] = u_shared[threadIdx.z][threadIdx.y][threadIdx.x];			
            v_shared[threadIdx.z-CELLD+2][threadIdx.y][threadIdx.x] = v_shared[threadIdx.z][threadIdx.y][threadIdx.x];			
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
            u_shared[threadIdx.z][threadIdx.y][threadIdx.x] = u_old[p];			
            v_shared[threadIdx.z][threadIdx.y][threadIdx.x] = v_old[p];			
        }

        __syncthreads();
    }
}
__global__ void kernelDiffusion_MovingTiles(gridprecision *u_old,gridprecision *v_old, gridprecision *u_new, gridprecision *v_new)
{
    __shared__ float u_shared[CELLD][CELLH][CELLW*2+1];

    unsigned int p, k, x, y, z;
    bool ok_read, ok_compute;
    float laplace_u;

    //3D position and memory address
    z = blockIdx.y*(CELLD-2)+threadIdx.z; 
    y = blockIdx.x*(CELLH-2)+threadIdx.y; 
    x = threadIdx.x;
    p = z * Zsize + y * Ysize + x;

    //precompute conditions
    ok_read = (z<ZMAX) & (y<YMAX);
    ok_compute = (threadIdx.y>0) & (threadIdx.y<CELLH-1) & (threadIdx.z>0) & (threadIdx.z<CELLD-1) & (y<YMAX-1) & (z<ZMAX-1);

    //read first two tiles of data
    if (ok_read)
    {
        u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1] = u_old[p];
        v_shared[threadIdx.z][threadIdx.y][threadIdx.x] = v_old[p];
        p+=CELLW;       
        u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1+CELLW] = u_old[p];
    };
    __syncthreads();

    //Move Tiles to the right: computing, writing result and reading new data in each iteration step.
    for (k=0; k < XMAX/CELLW; k++)
    {
        x = k*CELLW + threadIdx.x; 
        p = z * Zsize + y * Ysize + x;

        //calculate
        if (ok_compute & (x>0) & (x<XMAX-1))
        {
            laplace_u =

                u_shared[threadIdx.z+1][threadIdx.y-1][threadIdx.x+1] + 
                u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x] +
                u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x+2] +
                u_shared[threadIdx.z+1][threadIdx.y+1][threadIdx.x+1] + 
                u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x] +
                u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x+2] + 
                u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x] +
                u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x+2] + 
                u_shared[threadIdx.z-1][threadIdx.y-1][threadIdx.x+1] + 
                u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x] +
                u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x+2] +
                u_shared[threadIdx.z-1][threadIdx.y+1][threadIdx.x+1] + 
                2.0f*( u_shared[threadIdx.z-1][threadIdx.y][threadIdx.x+1] +
                        u_shared[threadIdx.z+1][threadIdx.y][threadIdx.x+1] +
                        u_shared[threadIdx.z][threadIdx.y-1][threadIdx.x+1] + 
                        u_shared[threadIdx.z][threadIdx.y+1][threadIdx.x+1] + 
                        u_shared[threadIdx.z][threadIdx.y][threadIdx.x] +
                        u_shared[threadIdx.z][threadIdx.y][threadIdx.x+2] 
                     ) -24.0f*u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1]; 

            u =  u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1] ;
            v =  v_shared[threadIdx.z][threadIdx.y][threadIdx.x] ;
            u_new[p] = u + (DaDt_6DxDx * laplace_u)+ ONEOVEREPSILON*(u - ONETHIRD*u*u*u - v)*DT;
            v_new[p] = v +EPSILON*(u + BETA - GAMMA*v)*DT;
        };

        __syncthreads();

        //copy last column of first tile to the extra column before the first tile (no bank conflict) 
        if (threadIdx.x==CELLW-1)
        {
            u_shared[threadIdx.z][threadIdx.y][0] = u_shared[threadIdx.z][threadIdx.y][CELLW];
        };

        //no need to syncthreads() here because threads (warps) that read in the memcopy above
        //are exactly the ones that will write to the same address in the following memcopy

        //moving the tile: copy the second tile onto the first
        //no bank conflict -> this is as fast as setting a new value to a register (in every thread)
        u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1] = u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1 + CELLW];

        //read new data into the second tile
        x = (k+2) * CELLW + threadIdx.x;
        p = z * Zsize + y * Ysize + x;
        if (ok_read & (x<XMAX)) 
        {
            v_shared[threadIdx.z][threadIdx.y][threadIdx.x] = v_old[p];

            if (k < XMAX/CELLW -2) //don't read in last two iterations
            {
                u_shared[threadIdx.z][threadIdx.y][threadIdx.x+1 + CELLW] = u_old[p];
            }
        }

        __syncthreads();
    }
}
