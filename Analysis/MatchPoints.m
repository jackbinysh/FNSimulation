function [ closestintersection ] = MatchPoints(oldpoints,points)
% output the points of closest intersection between old and new curve, by
% looking in the normal plane of each points of the old curve.

 NPold =length(oldpoints);
 NP =length(points);
% the two curves are basically identical, offset from one another slightly.
% begin by aligningthe two curves

deltas = (points - repmat(oldpoints(1,:),NP,1));
norms = sqrt(sum(deltas.^2,2));
[value,offsetlocation]=min(norms);

%offsetlocation contains the poisition on the new curve closest to the old
%curve. I will just shift them round to bring things to alignment.
points = circshift(points,-(offsetlocation-1));

% okay now point 1 on the old and new curves are closely aligned. Now we
% run round the curve, and work outwards from the corresponding point on
% the new curve ,looking for a segment intersection.
closestintersectionnewmethod = zeros(NPold,3);

for s = 1:length(oldpoints)
    
    distance = 0;
    check = 0;
    while check ~= 1
    
        V0 = oldpoints(s,:);
        n = oldpoints(incp(s,1,NPold),:) - oldpoints(s,:);
        
        t = incp(s,distance,NP);       
        P0 = points(t,:);
        P1 = points(incp(t,1,NP),:);
        [I,check]=plane_line_intersect(n,V0,P0,P1);
        
        if(check ==1) break; end
        
        t = incp(s,-distance,NP);       
        P0 = points(t,:);
        P1 = points(incp(t,-1,NP),:);
        [I,check]=plane_line_intersect(n,V0,P0,P1);
        
        if(check ==1) break; end
        
        distance = distance +1 ;
        
    end
    closestintersectionnewmethod(s,:) = I;
end
        
        
        
% % output the points of closest intersection between old and new curve, by
% % looking in the normal plane of each points of the old curve.
% 
% NPold =length(oldpoints);
% NP =length(points);
% 
% closestintersection = zeros(NPold,3);
% for s = 1:length(oldpoints)
%     
%     closestdistance = Inf;
%     
%     for t = 1:length(points)
%         V0 = oldpoints(s,:);
%         n = oldpoints(incp(s,1,NPold),:) - oldpoints(s,:);
%         P0 = points(t,:);
%         P1 = points(incp(t,1,NP),:);
%         [I,check]=plane_line_intersect(n,V0,P0,P1);
%         
%         if(check ==1)
%             intersectiondistance = norm(I - V0);
%             if(intersectiondistance< closestdistance)
%                 closestdistance = intersectiondistance;
%                 closestintersection(s,:) = I;
%             end
%         end
%     end
% end




