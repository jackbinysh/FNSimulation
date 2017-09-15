#include "ShellStripper.h"    //contains user defined variables for the simulation, and the parameters used 
#include <omp.h>
#include <math.h>
#include <string.h>

int main (int argc, char *argv[])
{

    // first up , get the name of the file
    std::string arg1(argv[1]);
    B_filename = arg1;
    cout << B_filename << endl;
    int startpos = B_filename.find('t');
    int endpos = B_filename.find('.');
    string number = B_filename.substr(B_filename.find('t')+1,B_filename.find('.')-B_filename.find('t')-1);
    double t = atoi(number.c_str());

    griddata griddata;
    griddata.Nx = initialNx;
    griddata.Ny = initialNy;
    griddata.Nz = initialNz;
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    // all major allocations are here
    // the main data storage arrays, contain info associated with the grid
    vector<double>u(Nx*Ny*Nz);
    vector<double>v(Nx*Ny*Nz);
    vector<double>ucvmag(Nx*Ny*Nz);

    uvfile_read_BINARY(u,v,ucvmag,griddata);
    extractinnerlayer(u,v,ucvmag,griddata);
    print_uv(u,v,ucvmag,t,griddata);

    return 0;
}

/*************************Functions for knot initialisation*****************************/

void print_marked( vector<int>&marked,int shelllabel, const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,n;
    stringstream ss;
    ss << shelllabel<< "marked.vtk";
    ofstream uvout (ss.str().c_str());

    uvout << "# vtk DataFile Version 3.0\nUV fields\nASCII\nDATASET STRUCTURED_POINTS\n";
    uvout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    uvout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    uvout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    uvout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    uvout << "SCALARS marked float\nLOOKUP_TABLE default\n";


    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                uvout << marked[n] << '\n';
            }
        }
    }
    uvout.close();
}
void print_uv( vector<double>&u, vector<double>&v,vector<double>&ucvmag, double t, const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,n;
    stringstream ss;
    ss << "uv_plot_stripped_" << t << ".vtk";
    ofstream uvout (ss.str().c_str(),std::ios::binary | std::ios::out);

    uvout << "# vtk DataFile Version 3.0\nUV fields\nBINARY\nDATASET STRUCTURED_POINTS\n";
    uvout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    uvout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    uvout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    uvout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    uvout << "SCALARS u float\nLOOKUP_TABLE default\n";


    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                float val =  FloatSwap(u[n]);
                uvout.write((char*) &val, sizeof(float));
            }
        }
    }

    uvout << "\n" << "SCALARS v float\nLOOKUP_TABLE default\n";


    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                float val =  FloatSwap(v[n]);
                uvout.write( (char*) &val, sizeof(float));
            }
        }
    }

    uvout << "\n" << "SCALARS ucrossv float\nLOOKUP_TABLE default\n";

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                //float val = FloatSwap(sqrt(ucvx[n]*ucvx[n] + ucvy[n]*ucvy[n] + ucvz[n]*ucvz[n]));
                float val = FloatSwap(ucvmag[n]);
                uvout.write( (char*) &val, sizeof(float));
            }
        }
    }

    uvout.close();
}

int uvfile_read_BINARY(vector<double>&u, vector<double>&v,vector<double>&ucvmag,const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    string temp,buff;
    stringstream ss;
    ifstream fin (B_filename.c_str(), std::ios::in  | std::ios::binary);
    int i,j,k,n;

    for(i=0;i<10;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n=pt(i,j,k,griddata);
                char* memblock;
                char* swapped;
                memblock = new char [sizeof(float)];
                swapped = new char [sizeof(float)];
                fin.read(memblock,sizeof(float));
                ByteSwap(memblock, swapped);
                float value = 12;
                memcpy(&value, swapped, 4);
                u[n] = value;
                delete[] memblock;
                delete[] swapped;
            }
        }
    }

    for(i=0;i<3;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n=pt(i,j,k,griddata);
                char* memblock;
                char* swapped;
                memblock = new char [sizeof(float)];
                swapped = new char [sizeof(float)];
                fin.read(memblock,sizeof(float));
                ByteSwap(memblock, swapped);
                float value = 12;
                memcpy(&value, swapped, 4);
                v[n] = value;
                delete[] memblock;
                delete[] swapped;
            }
        }
    }

    for(i=0;i<3;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n=pt(i,j,k,griddata);
                char* memblock;
                char* swapped;
                memblock = new char [sizeof(float)];
                swapped = new char [sizeof(float)];
                fin.read(memblock,sizeof(float));
                ByteSwap(memblock, swapped);
                float value = 12;
                memcpy(&value, swapped, 4);
                ucvmag[n] = value;
                delete[] memblock;
                delete[] swapped;
            }
        }
    }

    fin.close();

    return 0;
}

void extractinnerlayer(vector<double>&u,vector<double>&v,vector<double>&ucvmag,griddata& oldgriddata)
{
    // the plan is as follows: begin with points of u>ucrit, and ucv>0.1. From these points, grow a shell, i.e connect the dots outward. After this shell is grown, simply find all the other u>ucrit points which are NOT marked as being on this shell, and delete them
    int Nx = oldgriddata.Nx;
    int Ny = oldgriddata.Ny;
    int Nz = oldgriddata.Nz;
    double ucrit = -1.2;

    // first of all, go thorugh the array, and mark everything above ucrit, and with ucv>0.1  (ie actually on the knot) with a "to be evaluated" integer
    std::vector<int>marked(u.size(),0);
    for(int n = 0; n<u.size();n++)
    {
        if(u[n]>ucrit && ucvmag[n]>0.1){marked[n]=-1;}
    }

    // okay , grow the shell, starting from this inner layer
    growshell(u,marked,ucrit, oldgriddata);
    // now simply remove the u>ucrit which does not sit on this inner layer
    for(int n = 0; n<u.size();n++)
    {
        if(marked[n]!=-2 && u[n]>ucrit) { u[n] = -1.03; v[n] = -0.66;}
    }
}
void growshell(vector<double>&u,vector<int>& marked,double ucrit, const griddata& griddata)
{
    bool stillboundaryleft = true;
    while(stillboundaryleft)
    {
        grow(u,marked,ucrit,griddata);
        stillboundaryleft = false;
        for(int n = 0; n<u.size();n++)
        {
            if(marked[n]==-1) stillboundaryleft =true;
        }

    }
    // okay we have our marked points - they are marked with a 2 in the marked array. lets set all the uv values we find their to the resting state values
}
void grow(const vector<double>&u,vector<int>&marked,double ucrit,const griddata& griddata)
{
    // the marked array has the following values
    // 0 - not evaluated 
    // -1 - a boundary, to be grown 
    // -2 - the interrior, already grown
    // -3 - a temporary state, marked as a boundary during the update
    // positive numbers - layers of shells already marked
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    for(int i=0;i<Nx;i++)
    {
        for(int j=0; j<Ny; j++)
        {
            for(int k=0; k<Nz; k++)   //Central difference
            {
                int n = pt(i,j,k,griddata);
                if(marked[n] ==-1)
                {

                    for(int iinc=-1;iinc<=1;iinc++)
                    {
                        for(int jinc=-1; jinc<=1; jinc++)
                        {
                            for(int kinc=-1; kinc<=1; kinc++)   //Central difference
                            {
                                int neighboringn = pt(incabsorb(i,iinc,Nx),incabsorb(j,jinc,Ny),incabsorb(k,kinc,Nz),griddata);
                                if(marked[neighboringn] == 0 && u[neighboringn] > ucrit) marked[neighboringn] = -3;
                            }
                        }
                    }
                }
            }
        }
    }
    for(int n = 0; n<u.size();n++)
    {
        if(marked[n]==-1){marked[n] =-2;}
        if(marked[n]==-3){marked[n] =-1;}
    }
}

inline int circularmod(int i, int N)    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
{
    if(i<0) return N - ((-i)%N);
    else return i%N;
}
// inlined functions for incrementing things respecting boundaries
inline int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}
inline int incabsorb(int i, int p, int N)    //increment with reflecting boundary between -1 and 0 and N-1 and N
{
    if(i+p<0) return (0);
    if(i+p>N-1) return (N-1);
    return (i+p);
}
// this function is specifically designed to incremenet, in the direction specified, respecting the boundary conditions, which are global enums
inline int gridinc(int i, int p, int N, int direction )    //increment with reflecting boundary between -1 and 0 and N-1 and N
{
    return incp(i,p,N);
}
inline double x(int i,const griddata& griddata)
{
    return (i+0.5-griddata.Nx/2.0)*h;
}
inline double y(int i,const griddata& griddata)
{
    return (i+0.5-griddata.Ny/2.0)*h;
}
inline double z(int i,const griddata& griddata)
{
    return (i+0.5-griddata.Nz/2.0)*h;
}
inline  int pt( int i,  int j,  int k,const griddata& griddata)       //convert i,j,k to single index
{
    return (i*griddata.Ny*griddata.Nz+j*griddata.Nz+k);
}
inline int sign(int i)
{
    if(i==0) return 0;
    else return i/abs(i);
}
float FloatSwap( float f )
{    
    union
    {
        float f;
        char b[4];
    } dat1, dat2;

    dat1.f = f;
    dat2.b[0] = dat1.b[3];
    dat2.b[1] = dat1.b[2];
    dat2.b[2] = dat1.b[1];
    dat2.b[3] = dat1.b[0];
    return dat2.f;
}

void ByteSwap(const char* TobeSwapped, char* swapped )
{    
    swapped[0] = TobeSwapped[3];
    swapped[1] = TobeSwapped[2];
    swapped[2] = TobeSwapped[1];
    swapped[3] = TobeSwapped[0];
    return; 
}
