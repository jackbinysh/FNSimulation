%13/04/2017, Jack Binysh

% an analysis of the 4_1 data.
% the script reads in knotplot files, gets their Hasimoto complex
% curvatures, and plots them.

clear all;

% read in parameters
minvalue = 9003.4;
maxvalue = 9984.4;
increment = 10.9;

%9.1124e+03;

% big matrices to hold the curvatures and torsions.
Curvatures=[];
Torsions=[];
Lengths = [];
vdotns = [];
vdotbs = [];
vs = {};
points = {};
Phis = []; 
%this variable holds the length up to point i on the curve - recalling that
%i is essentially an affine paramter for the curve, and doesnt correspond
%to physical arc length. NOTE: we should imagine this vector as 0 indexed. ie
%the first element corresponds to the length at point 1, which is 1 along
%the curve, and so on.
IntegratedLengths = [];

% Read in the data.
for index = minvalue:increment:maxvalue
    filename = strcat('/home/jackbinysh/Code/Knot_surface_calc/4_1_results/knotplot0_',num2str(index),'.vtk') ; 
    
    knotplot = CurveRead(filename);
    Curvatures=horzcat(Curvatures,knotplot.Curvature);
    Torsions=horzcat(Torsions,knotplot.Torsion);
    Lengths=horzcat(Lengths,knotplot.Length);
    points=horzcat(points,num2cell(knotplot.POINTS,2));
    
    
    % get size of vdotn and vdotb
    n = knotplot.n;
    vn = knotplot.vdotn;
    magvn = arrayfun(@(i)dot(vn(i,:),n(i,:)),1:length(vn))';
    b = knotplot.b;
    vb = knotplot.vdotb;
    magvb = arrayfun(@(i)dot(vb(i,:),b(i,:)),1:length(vb))';
    
    % get the velocity vector itself
    v = num2cell(vn + vb,2);
    
    vs=horzcat(vs,v);
    vdotns=horzcat(vdotns,magvn);
    vdotbs=horzcat(vdotbs,magvb);
    IntegratedLengths=horzcat(IntegratedLengths,cumsum(knotplot.Length));
end

% Two things:
%(1) Hasimoto allows us a phase factor freedom in his complex curvature  -
%this is just shifting the origin of the curve. This is annoying, I want to
%overlay them so we remove this degree of freedom
%(2) Remeber, we arent paramerized by arc lenght with i ! It would be nice
%to standadrize everytinig so we knot we are measuring by arc length.

%adressing (1)
CorrectOffset =0;
offsetK = [];
offsetTau = [];
offsetLength = [];
offsetvn = [];
offsetvb = [];
offsetr = {};
previousstart=zeros(1,3);
for i = 1:size(Curvatures,2)
    
    ArcLength = IntegratedLengths(:,i);
    TotalLength=ArcLength(end);
    Tau = Torsions(:,i);
    K = Curvatures(:,i);
    vn =vdotns(:,i);
    vb = vdotbs(:,i);
    r = cell2mat(points(:,i));
    x = r(:,1);
    y = r(:,2);
    z = r(:,3);
    
    ExtendedLength = vertcat(ArcLength-TotalLength, ArcLength, TotalLength + ArcLength);
    ExtendedK = vertcat(K,K,K);
    ExtendedTau = vertcat(Tau,Tau,Tau);
    Extendedvn = vertcat(vn,vn,vn);
    Extendedvb = vertcat(vb,vb,vb);
    Extendedx = vertcat(x,x,x);
    Extendedy = vertcat(y,y,y);
    Extendedz = vertcat(z,z,z);
    
    
    myfunc = @(offset)CurveOffsetError(offset, Curvatures(:,1),IntegratedLengths(:,1), Curvatures(:,i),IntegratedLengths(:,i));
     trialoffsets = linspace(-0.6*IntegratedLengths(end,1),0.6*IntegratedLengths(end,1),100);
  
    minima = arrayfun(@(x)fminsearch(myfunc,x),trialoffsets);
    % ok so we have these candidate minima. But some will be local and
    % others global, and only one should match closely to the physical
    % start position of the previous curve.
    
    % lets just coarsely strip off the ones that are clearly wrong
    minima = minima(myfunc(minima)<(min(myfunc(minima))+ max(myfunc(minima)))/2);
    
    % of these guys, only one set should match physical space well. lets
    % check!  
    trial = zeros(1,3);
    Error = zeros(length(minima),1);
    for i = 1:length(minima)

        trial(1)=  interp1(ExtendedLength, Extendedx,- minima(i));
        trial(2)=  interp1(ExtendedLength, Extendedy,- minima(i));
        trial(3)=  interp1(ExtendedLength, Extendedz,- minima(i));
        Error(i) = norm(trial - previousstart);

    end
    [~,location] = min(Error);
    offset = minima(location);

    
%     offsetx=  interp1(ExtendedLength, Extendedx, ArcLength - minima);
%     offsety=  interp1(ExtendedLength, Extendedy, ArcLength - offset);
%     offsetz=  interp1(ExtendedLength, Extendedz, ArcLength - offset);
%     offsetr = horzcat(offsetr,num2cell([offsetx,offsety,offsetz],2));
%     
%     [startpoint,startindex] = min(myfunc(trialoffsets));
%     
%     
%     startingoffset = trialoffsets(startindex);
%     
%     CorrectOffset = fminsearch(myfunc,startingoffset);
    

    
    offsetK = horzcat(offsetK, interp1(ExtendedLength, ExtendedK, ArcLength - offset));
    offsetTau = horzcat(offsetTau, interp1(ExtendedLength, ExtendedTau, ArcLength - offset));
    offsetLength = horzcat(offsetLength, interp1(ExtendedLength, ExtendedLength, ArcLength - offset));
    offsetvn = horzcat(offsetvn, interp1(ExtendedLength, Extendedvn, ArcLength - offset));
    offsetvb= horzcat(offsetvb, interp1(ExtendedLength, Extendedvb, ArcLength - offset));
    offsetx=  interp1(ExtendedLength, Extendedx, ArcLength - offset);
    offsety=  interp1(ExtendedLength, Extendedy, ArcLength - offset);
    %Phis = horzcat(Phis,HasimotoTransform(offsetK(:,end),offsetTau(:,end),Lengths(:,i),0));
    offsetz=  interp1(ExtendedLength, Extendedz, ArcLength - offset);
    offsetr = horzcat(offsetr,num2cell([offsetx,offsety,offsetz],2));
    
   % keep track of the first element of offsetr.
    previousstart = [offsetx(1),offsety(1),offsetz(1)]
    
    
    % okay, get the hasimoto transform
    %Phis = horzcat(Phis,HasimotoTransform(offsetK(:,end),offsetTau(:,end),Lengths(:,i),0));
end

hold on
comrs = [];
comvs = [];
for l = 1:size(Curvatures,2)

    ri = points(:,l);
    mi = Lengths(:,l);
    miri = zeros(length(ri),3);
    for i = 1:length(points)
        miri(i,:) = ri{i} *mi(i);
    end
    comr = sum(miri)/sum(mi);
    comrs = vertcat(comrs,comr);

%     vi = vs(:,l);
%     mi = Lengths(:,l);
%     mivi = zeros(length(ri),3);
%     for i = 1:length(points)
%         mivi(i,:) = vi{i} *mi(i);
%     end
%     comv = sum(mivi)/sum(mi);
%     comvs = vertcat(comvs,comv);

end
% 
% x = cell2mat(vi);
% scatter3(x(:,1),x(:,2),x(:,3))
%     figure(1); hold on; scatter3(comrs(1),comrs(2),comrs(3))
%     figure(2); hold on; scatter3(comvs(1),comvs(2),comvs(3))
%         figure(2); hold on; scatter3(comvs(:,1),comvs(:,2),comvs(:,3))


hold on
for i = 1:length(offsetr)
    traj = cell2mat(offsetr(i,:)');
    correctedtraj = traj - comrs
     centredtraj = correctedtraj - repmat(mean(correctedtraj), length(correctedtraj),1);
     scatter3(centredtraj(:,1),centredtraj(:,2),centredtraj(:,3))
end
 
 x = cross(centredtraj(1,:),centredtraj(10,:))
r =vrrotvec( x, [0 0 1])
a =vrrotvec2mat(r) 


q = arrayfun(@(i)a*centredtraj(i,:)',(1:length(correctedtraj)),'UniformOutput',false)
q = cell2mat(q)'

scatter(q(:,1),q(:,2))
 
 
 scatter3(traj(:,1),traj(:,2),traj(:,3))
scatter3(correctedtraj(:,1),correctedtraj(:,2),correctedtraj(:,3))
scatter3(centredtraj(:,1),centredtraj(:,2),centredtraj(:,3))


for i = 1: length(points)
    for j = 1:length(points)
        omega(i,j) = num2cell(cross(vi{i},vi{j})/dot(vi{i},ri{j}),2);
    end
end



scatter3(q(:,1),q(:,2),q(:,3))
scatter(q(:,1),q(:,2))




ri = points(:,l);
mi = Lengths(:,l);
miri = zeros(length(ri),3);
for i = 1:length(points)
    miri(i,:) = ri{i} *mi(i);
end
comr = sum(miri)/sum(mi);

vi = vs(:,l);
mi = Lengths(:,l);
mivi = zeros(length(ri),3);
for i = 1:length(points)
    mivi(i,:) = vi{i} *mi(i);
end
comv = sum(mivi)/sum(mi);

% right now we have the velocities with the com velocity taken off. lets get the axis of rotation:
%omega = cell(length(points),length(points));
%for i = 1: length(points)
%    for j = 1:length(points)
%        omega(i,j) = num2cell(cross(vi{i},vi{j})/dot(vi{i},ri{j}),2);
%    end
%end


