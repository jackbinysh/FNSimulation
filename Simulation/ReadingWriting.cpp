#include "ReadingWriting.h"
#include "FN_Constants.h"
#include "FN_Knot.h"
#include <string.h>

int uvfile_read_BINARY(vector<double>&u, vector<double>&v,const Griddata& griddata)
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

    fin.close();

    return 0;
}

int uvfile_read_ASCII(vector<double>&u, vector<double>&v,const Griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    string temp,buff;
    stringstream ss;
    ifstream fin (B_filename.c_str());
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
                ss.clear();
                ss.str("");
                if(fin.good())
                {
                    if(getline(fin,buff))
                    {
                        ss << buff;
                        ss >> u[n];
                    }
                }
                else
                {
                    cout << "Something went wrong!\n";
                    return 1;
                }
            }
        }
    }

    for(i=0;i<2;i++)
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
                ss.clear();
                ss.str("");
                if(fin.good())
                {
                    if(getline(fin,buff)) ss << buff;
                    ss >> v[n];
                }
                else
                {
                    cout << "Something went wrong!\n";
                    return 1;
                }
            }
        }
    }

    fin.close();

    return 0;
}

int uvfile_read(vector<double>&u, vector<double>&v, vector<double>& ku, vector<double>& kv, vector<double>& ucvx, vector<double>& ucvy,vector<double>& ucvz, vector<double>& ucvmag,Griddata& griddata)
{
    string buff,datatype,dimensions,xdim,ydim,zdim;
    ifstream fin (B_filename.c_str());
    for(int i=0;i<4;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff) &&(i==2)) datatype = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    if(datatype.compare("ASCII")==0)
    {
        uvfile_read_ASCII(u,v,griddata);
    }
    else if(datatype.compare("BINARY")==0)
    {
        uvfile_read_BINARY(u,v,griddata);
    }

    // okay we've read in the file - now, did we want to interpolate?
    if(interpolationflag)
    {
        cout << "interpolating grid \n";

        Griddata interpolatedgriddata;
        interpolatedgriddata.Nx = interpolatedNx;
        interpolatedgriddata.Ny = interpolatedNy;
        interpolatedgriddata.Nz = interpolatedNz;
        interpolatedgriddata.h = ((initialNx-1)*initialh)/(interpolatedNx-1);

        vector<double>interpolatedugrid(interpolatedNx*interpolatedNy*interpolatedNz);
        vector<double>interpolatedvgrid(interpolatedNx*interpolatedNy*interpolatedNz);

        // interpolate u and v
        likely::TriCubicInterpolator interpolatedu(u, initialh, initialNx,initialNy,initialNz);
        likely::TriCubicInterpolator interpolatedv(v, initialh, initialNx,initialNy,initialNz);
        for(int i=0;i<interpolatedNx;i++)
        {
            for(int j=0; j<interpolatedNy; j++)
            {
                for(int k=0; k<interpolatedNz; k++)   //Central difference
                {

                    // get the point in space this gridpoint corresponds to
                    double px= x(i,interpolatedgriddata);
                    double py= y(j,interpolatedgriddata);
                    double pz= z(k,interpolatedgriddata);
                    // interpolate
                    interpolatedugrid[pt(i,j,k,interpolatedgriddata)]= interpolatedu(px,py,pz);
                    interpolatedvgrid[pt(i,j,k,interpolatedgriddata)]= interpolatedv(px,py,pz);
                }
            }
        }

        // resize all arrays, set the new u and v arrays, and set the new griddata
        ucvx.resize(interpolatedNx*interpolatedNy*interpolatedNz);
        ucvy.resize(interpolatedNx*interpolatedNy*interpolatedNz);
        ucvz.resize(interpolatedNx*interpolatedNy*interpolatedNz);
        ucvmag.resize(interpolatedNx*interpolatedNy*interpolatedNz);
        ku.resize(4*interpolatedNx*interpolatedNy*interpolatedNz);
        kv.resize(4*interpolatedNx*interpolatedNy*interpolatedNz);

        u = interpolatedugrid;
        v = interpolatedvgrid;

        griddata=interpolatedgriddata;
    }

    return 0;
}

void print_knot( double t, vector<knotcurve>& knotcurves,const Griddata& griddata)
{
    for( int c=0; c < (knotcurves.size()) ; c++)
    {

        /***Write values to file*******/
        stringstream ss;
        ss << "globaldata" << "_" << c <<  ".txt";
        ofstream wrout (ss.str().c_str(), std::ofstream::app);
        wrout << t << '\t' << knotcurves[c].writhe << '\t' << knotcurves[c].twist << '\t' << knotcurves[c].length << '\n';
        wrout.close();

        ss.str("");
        ss.clear();

        ss << "knotplot" << c << "_" << t <<  ".vtk";
        ofstream knotout (ss.str().c_str());

        int i;
        int n = knotcurves[c].knotcurve.size();

        knotout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        knotout << "POINTS " << n << " float\n";

        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].xcoord << ' ' << knotcurves[c].knotcurve[i].ycoord << ' ' << knotcurves[c].knotcurve[i].zcoord << '\n';
        }

        knotout << "\n\nCELLS " << n << ' ' << 3*n << '\n';

        for(i=0; i<n; i++)
        {
            knotout << 2 << ' ' << i << ' ' << incp(i,1,n) << '\n';
        }

        knotout << "\n\nCELL_TYPES " << n << '\n';

        for(i=0; i<n; i++)
        {
            knotout << "3\n";
        }

        knotout << "\n\nPOINT_DATA " << n << "\n\n";

        knotout << "\nSCALARS Curvature float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].curvature << '\n'; }

        knotout << "\nSCALARS Torsion float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].torsion << '\n';
        }

        knotout << "\nVECTORS A float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].ax << ' ' << knotcurves[c].knotcurve[i].ay << ' ' << knotcurves[c].knotcurve[i].az << '\n';
        }

        knotout << "\nVECTORS V float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].vx << ' ' << knotcurves[c].knotcurve[i].vy << ' ' << knotcurves[c].knotcurve[i].vz << '\n';
        }
        knotout << "\nVECTORS t float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].tx << ' ' << knotcurves[c].knotcurve[i].ty << ' ' << knotcurves[c].knotcurve[i].tz << '\n';
        }
        knotout << "\nVECTORS n float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].nx << ' ' << knotcurves[c].knotcurve[i].ny << ' ' << knotcurves[c].knotcurve[i].nz << '\n';
        }
        knotout << "\nVECTORS b float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].bx << ' ' << knotcurves[c].knotcurve[i].by << ' ' << knotcurves[c].knotcurve[i].bz << '\n';
        }
        knotout << "\nVECTORS vdotn float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].vdotnx << ' ' << knotcurves[c].knotcurve[i].vdotny << ' ' << knotcurves[c].knotcurve[i].vdotnz << '\n';
        }
        knotout << "\nVECTORS vdotb float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].vdotbx << ' ' << knotcurves[c].knotcurve[i].vdotby << ' ' << knotcurves[c].knotcurve[i].vdotbz << '\n';
        }
        knotout << "\n\nCELL_DATA " << n << "\n\n";
        knotout << "\nSCALARS Writhe float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].writhe << '\n';
        }

        knotout << "\nSCALARS Twist float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].twist << '\n';
        }

        knotout << "\nSCALARS Length float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].length << '\n';
        }
        knotout.close();
    }
}

void print_B_phi( vector<double>&phi, const Griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    double h = griddata.h;
    int i,j,k,n;
    string fn = "phi.vtk";

    ofstream Bout (fn.c_str());

    Bout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
    Bout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    Bout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    Bout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    Bout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    Bout << "SCALARS Phi float\nLOOKUP_TABLE default\n";
    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                Bout << phi[n] << '\n';
            }
        }
    }
    Bout.close();
}

void print_uv( vector<double>&u, vector<double>&v, vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz,vector<double>&ucvmag, double t, const Griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    double h = griddata.h;
    int i,j,k,n;
    stringstream ss;
    ss << "uv_plot" << t << ".vtk";
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
