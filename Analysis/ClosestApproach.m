function [ closestapproach ] = ClosestApproach( knotdata, excludedlength)
%CURVEREAD For each point on the knot, outside of an excluded region on the
%filament, find the point of closest approach

    data = knotdata.POINTS;
    lengths = knotdata.Length;
    N = length(data);
    closestapproach = zeros(N,2);
    
    
    nc = num2cell(data,2);
    [x,y] = ndgrid(1:length(data));
    distancematrix = cellfun(@(x,y)sqrt(sum((x-y).^2)),nc(x),nc(y));
    
 
    for i=1:N
        distances = distancematrix(i,:);
        
        % now, working from i outwards, get the guys to be excluded
        % go up
        jhigh=i;
        totallength = 0;
        while (totallength < excludedlength)
            
            totallength = totallength + lengths(jhigh);
            distances(jhigh)=inf;
            jhigh = incp(jhigh,1,N);
        end
        
        jlow=i;
        totallength = 0;
        while totallength<excludedlength
            jlow = incp(jlow,-1,N);
            totallength = totallength + lengths(jlow);
            distances(jlow)=inf;          
        end
        
        [value,index] = min(distances);
        closestapproach(i,1) = value;
        closestapproach(i,2) = index;
    end
end
       
        
        
        
        
    
    

