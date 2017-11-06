function [] = CurveWrite(knotplot, filename )
    % write a knotplot file to disk
    n = length(knotplot.POINTS);
    fileID = fopen(filename,'w');
    fprintf(fileID,'# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET UNSTRUCTURED_GRID\n');
    
    fprintf(fileID,'POINTS %d float\n',n);
    fprintf(fileID,'%f %f %f\n',knotplot.POINTS');
    
    fprintf(fileID,'\n\nCELLS %d %d\n',n,3*n);
    for i = 0:n-2
        fprintf(fileID,'%d %d %d\n',2,i,incp(i,1,n));
    end
    fprintf(fileID,'%d %d %d\n',2,n-1,0);
    
    fprintf(fileID,'\n\nCELL_TYPES %d\n',n);
    for i = 0:n-1
        fprintf(fileID,'%d\n',3);
    end
    
    fprintf(fileID,'\n\nPOINT_DATA %d\n\n',n);
    
    fprintf(fileID,'\nSCALARS Curvature float\nLOOKUP_TABLE default\n');
    fprintf(fileID,'%f\n',knotplot.Curvature);
    
    fprintf(fileID,'\nSCALARS Torsion float\nLOOKUP_TABLE default\n');
    fprintf(fileID,'%f\n',knotplot.Torsion);
    
    fprintf(fileID,'\nVECTORS A float\n');
    fprintf(fileID,'%f %f %f\n',knotplot.A');

    %fprintf(fileID,'\nVECTORS V float\n');
    %fprintf(fileID,'%f %f %f\n',knotplot.V);

    fprintf(fileID,'\nVECTORS t float\n');
    fprintf(fileID,'%f %f %f\n',knotplot.t');
    
    fprintf(fileID,'\nVECTORS n float\n');
    fprintf(fileID,'%f %f %f\n',knotplot.n');
    
    fprintf(fileID,'\nVECTORS b float\n');
    fprintf(fileID,'%f %f %f\n',knotplot.b');
    
    fprintf(fileID,'\nVECTORS vdotn float\n');
    fprintf(fileID,'%f %f %f\n',knotplot.vdotn');
    
    fprintf(fileID,'\nVECTORS vdotb float\n');
    fprintf(fileID,'%f %f %f\n',knotplot.vdotb');

    fprintf(fileID,'\n\nCELL_DATA %d\n\n',n');
    
    fprintf(fileID,'\nSCALARS Writhe float\nLOOKUP_TABLE default\n');
    fprintf(fileID,'%f\n',knotplot.Writhe);
   
    fprintf(fileID,'\nSCALARS Twist float\nLOOKUP_TABLE default\n');
    fprintf(fileID,'%f\n',knotplot.Twist);
    
    fprintf(fileID,'\nSCALARS Length float\nLOOKUP_TABLE default\n');
    fprintf(fileID,'%f\n',knotplot.Length);
    
    
    fclose(fileID);

end
