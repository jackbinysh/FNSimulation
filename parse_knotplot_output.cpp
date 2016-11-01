//
//  parse_knotplot_output.cpp
//  
//
//  Created by Carl Whitfield on 01/11/2016.
//
//

/*
This file takes text output from KnotPlot of the form
x0 y0 z0
x1 y1 z1
...
and produces a template surface evolver file with vertices and boundary edges. Futher edges and faces will need to be computed by the user. It also generates a mathematica file to visualise the knot with numbered points.
 
To generate the input: Create a curve in KnotPlot then type in the terminal:
 save XXXXX.txt raw
 
 To generate the mathematica output:
 Copy and paste the contents of XXXXX_mm_input.txt into a mathematica file.
*/


#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

int main (int argc, char *argv[])
{
    unsigned int i, npts=0;
    string filename;
    string buff;
    stringstream ss;
    double xpos, ypos, zpos;
    
    if (argc<2)
    {
        cout << "Please enter a filename to read in\n";
        cin  >> filename;
    }
    else
    {
        if (argc > 2) cout << "Ignoring arguments after filename\n";
        filename = argv[1];
    }
    
    ifstream knotin (filename.c_str());
    string filehead; //stores the bit before the point
    string feoutname;
    string mmoutname;
    /*Some text manipulation*/
    ss.clear();
    ss.str("");
    ss << filename;
    getline(ss,filehead,'.');
    ss.clear();
    ss.str("");
    ss << filehead << ".fe";
    ss >> feoutname;
    ss.clear();
    ss.str("");
    ss << filehead << "_mm_input.txt";
    ss >> mmoutname;
    /***********************/
    ofstream feout (feoutname.c_str());   //file for Surface Evolver input
    ofstream mmout (mmoutname.c_str());   //file for mathematica commands
    
    feout << "// Surface Evolver template for knotplot generated file " << filehead << "\n\n";
    feout << "vertices\n\n";
    mmout << "ptsize = 0.2;     (*Size of spheres*)\n\n";
    mmout << "knotpts = {";
    
    while(knotin.good())
    {
        ss.clear();
        ss.str("");
        if(getline(knotin,buff))
        {
            npts++;
            ss << buff;
            ss >> xpos >> ypos >> zpos;
            feout << npts << ' ' << xpos << ' ' << ypos << ' ' << zpos << " fixed\n";
            if(npts>1) mmout << ',';
            mmout << '{' << xpos << ',' << ypos << ',' << zpos << "}";
        }
    }
    
    knotin.close();
    
    feout << "\nedges\n\n";
    
    for(i=1;i<npts;i++)
    {
        feout << i << ' ' << i << ' ' << (i+1) << " fixed\n";
    }
    feout << npts << ' ' << npts << ' ' << 1 << " fixed\n\n\n\n\n";
    
    feout << "read \n\n//Good command to use\n gogo1 := { K 30; u; g 10; K 30; u; g 10; r; u; g 10; u; }\n";
    
    feout.close();
    
    mmout << "};   (*Matrix of point positions*)\n\n";
    
    mmout << "Show[Table[Graphics3D[{Sphere[knotpts[[i]],ptsize],Text[i,knotpts[[i]]]}],{i,1," << npts << "}], Graphics3D[{Line[knotpts],Line[{knotpts[[" << npts << "]], knotpts[[1]]}]}]]      (*Numbered plot of knot*) \n";

    mmout.close();
    
    return 0;
}
