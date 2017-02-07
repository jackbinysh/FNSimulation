%a script to perform curve smoothing ,analysis etc of filaments from
% fitzhugh nagumo model. Created 20/12/2016, Jack Binysh


% read in the vtk files. 
M = importdata('knotplot9.28.vtk');
M = M.textdata;

start = find(strcmp(M,'POINTS'));
finish = find(strcmp(M,'CELLS'));

oldpoints = str2double(M(start+1:finish-1,1:3));


M = importdata('knotplot9.44.vtk');
M = M.textdata;

start = find(strcmp(M,'POINTS'));
finish = find(strcmp(M,'CELLS'));

points = str2double(M(start+1:finish-1,1:3));

%%% okay we have the points! now lets do some smoothing. first make the
%%% window

knotcurveold = windowconv(oldpoints);
knotcurve = windowconv(points);



%% now do the geometry on the smoothed curves %%

Length = length(knotcurve);
OldLength = length(knotcurveold);

% get offset
delta = knotcurve - repmat(knotcurveold(1,:), length(knotcurve),1);
result = arrayfun(@(idx) norm(delta(idx,:)), 1:size(delta,1))';
[value,offset]= min(result);


% rather than mess around with indexing, make periodic
periodicknotcurveold = [knotcurveold;knotcurveold];
periodicknotcurve = [knotcurve;knotcurve;knotcurve;knotcurve;knotcurve];

for s = 1:OldLength
    intersection = 0;
    m = s+offset+2*Length;
    stepnum = 0;
    
    while(intersection ~= 1)
        
        n = periodicknotcurveold(s+1,:) - periodicknotcurveold(s,:); 
        V0 =  periodicknotcurveold(s,:);
        P0 =  periodicknotcurve(m,:);
        P1 = periodicknotcurve(m+1,:);
        [intersectionpoint, intersection] = plane_line_intersect(n,V0,P0,P1);
        stepnum = stepnum +1 ;
        if(mod(stepnum,2) == 0)
            m = m -stepnum;
        else
            m = m+stepnum;
        end   
    end
    
    v(s,:) = (intersectionpoint - knotcurveold(s,:))/0.16;
 
end
