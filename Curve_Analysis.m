% a script to perform curve smoothing ,analysis etc of filaments from
% % fitzhugh nagumo model. Created 20/12/2016, Jack Binysh
% 
% % read in the vtk files. 
% M = importdata('knotplot20.vtk');
% M = M.textdata;
% 
% start = find(strcmp(M,'POINTS'));
% finish = find(strcmp(M,'CELLS'));
% 
% points = str2double(M(start+1:finish-1,1:3));
% 
% %%% okay we have the points! now lets do some smoothing. first make the
% %%% window
% 
% %cutoff freq
% fc = 0.002;
% %width of attenuation
% b = 0.001;
% 
% N = ceil(4/b);
% if mod(N,2) == 0
%     N = N+1;
% end
% w = blackman(N);
% h = sinc(2*fc*((0:N-1) - (N-1)/2))'.* w;
% h = h/sum(h);
% 
% % make signal periodic
% numreps = ( 2.*round(((N/length(points))+1)/2)-1  )+ 2;
% Length = length(points);
% reppoints = repmat(points,numreps,1);
% 
% %convolve
% for i = 1:3
%     x(:,i) = conv(reppoints(:,i),h,'same');
%     knotcurve(:,i) = x(floor(numreps/2)*Length:ceil(numreps/2)*Length,i);
% end
% 
% %% now do the geometry on the smoothed curves %%

Length = length(knotcurve);
OldLength = length(knotcurveold);

% get offset
delta = knotcurve - repmat(knotcurveold(1,:), length(knotcurve),1);
result = arrayfun(@(idx) norm(delta(idx,:)), 1:size(delta,1))';
[value,offset]= min(result);

for s = 1:OldLength
    intersection = 0;
    m = s+offset;
    stepnum = 0;
    
    while(intersection ~= 1)
        
         if (m == 0) m = Length; end;
         
        Oldstart = mod(s,OldLength); if (Oldstart == 0) Oldstart = OldLength; end;
        Oldend = mod(s+1,OldLength); if (Oldend == 0) Oldend = OldLength; end;
        Start = mod(m,Length); if (Start == 0) Start = Length; end;
        End = mod(m+1,Length); if (End == 0) End = Length; end;
        
        n = knotcurveold(Oldend,:) - knotcurveold(Oldstart,:); 
        V0 =  knotcurveold(Oldstart,:);
        P0 =  knotcurve(Start,:);
        P1 = knotcurve(End,:);
        [intersectionpoint, intersection] = plane_line_intersect(n,V0,P0,P1);
        stepnum = stepnum +1 ;
        if(mod(stepnum,2) == 0)
            m = mod(m-stepnum,Length);
        else
            m = mod(m+stepnum,Length);
        end   
    end
    
    v(s,:) = intersectionpoint - knotcurveold(s,:);
 
end





