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
    
    % get size of vdotn and vdotb
    vn = knotplot.vdotn;
    magvn = arrayfun(@(i)norm(vn(i,:)),1:length(vn))';
    vb = knotplot.vdotb;
    magvb = arrayfun(@(i)norm(vb(i,:)),1:length(vb))';
    
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
for i = 1:size(Curvatures,2)
    myfunc = @(offset)CurveOffsetError(offset, Curvatures(:,1),IntegratedLengths(:,1), Curvatures(:,i),IntegratedLengths(:,i));
    
    trialoffsets = -0.45*IntegratedLengths(end,1):0.1:0.45*IntegratedLengths(end,1);
    [startpoint,startindex] = min(myfunc(trialoffsets));
    startingoffset = trialoffsets(startindex);
    
    CorrectOffset = fminsearch(myfunc,startingoffset);
    
    ArcLength = IntegratedLengths(:,i);
    TotalLength=ArcLength(end);
    offset = CorrectOffset;
    Tau = Torsions(:,i);
    K = Curvatures(:,i);
    vn =vdotns(:,i);
    vb = vdotbs(:,i);
    
    ExtendedLength = vertcat(ArcLength-TotalLength, ArcLength, TotalLength + ArcLength);
    ExtendedK = vertcat(K,K,K);
    ExtendedTau = vertcat(Tau,Tau,Tau);
    Extendedvn = vertcat(vn,vn,vn);
    Extendedvb = vertcat(vb,vb,vb);
    
    offsetK = horzcat(offsetK, interp1(ExtendedLength, ExtendedK, ArcLength - offset));
    offsetTau = horzcat(offsetTau, interp1(ExtendedLength, ExtendedTau, ArcLength - offset));
    offsetLength = horzcat(offsetLength, interp1(ExtendedLength, ExtendedLength, ArcLength - offset));
    offsetvn = horzcat(offsetvn, interp1(ExtendedLength, Extendedvn, ArcLength - offset));
    offsetvb= horzcat(offsetvb, interp1(ExtendedLength, Extendedvb, ArcLength - offset));
    
    % okay, get the hasimoto transform
    Phis = horzcat(Phis,HasimotoTransform(offsetK(:,end),offsetTau(:,end),Lengths(:,i),0));
end
