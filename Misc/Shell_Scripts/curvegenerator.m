
% trefoil
N = 100
points = [];
for i = 0:N-1
    x = sin(((2*pi)/N)*i)+ 2*sin((2*(2*pi)/N)*i);
    y = cos(((2*pi)/N)*i)- 2*cos((2*(2*pi)/N)*i);
    z =-sin((3*(2*pi)/N)*i);
    
    points(end+1,:) = [x y z];
end

% pringle
clear all;
N = 1000
points = [];
for i = 0:N-1
    x = sin(((2*pi)/N)*i);
    y = cos(((2*pi)/N)*i);
    z =0.2*sin( ( (4*pi)/N ) *i); 
    points(end+1,:) = [x y z];
end
%scatter3(points(:,1),points(:,2) ,points(:,3) )
dlmwrite('Pringle.txt',points,'delimiter',' ','precision',10)


clear all;
 x0 = 0.024; 
 y0 = 0.036;
 z0 = -0.022;
r= 0.57;
p=1;
 q=2;
 
NP = 800
points = [];
for s = 0:NP-1    
      x = x0 + (1.24+r*cos(2.0*pi*q*s/NP))*cos(2.0*pi*p*s/NP);
      y = y0 + (1.24+r*cos(2.0*pi*q*s/NP))*sin(2.0*pi*p*s/NP);
      z = z0 + r*sin(2.0*pi*q*s/NP);
    points(end+1,:) = [x y z];
end
dlmwrite('Torus_0.txt',points,'delimiter',' ','precision',10)

points = [];
for s = 0:NP-1    
      x = x0 + (1.24-r*cos(2.0*pi*q*s/NP))*cos(2.0*pi*p*s/NP);
      y = y0 - (1.24-r*cos(2.0*pi*q*s/NP))*sin(2.0*pi*p*s/NP);
      z = z0 + r*sin(2.0*pi*q*s/NP);
    points(end+1,:) = [x y z];
end
dlmwrite('Torus_1.txt',points,'delimiter',' ','precision',10)