function [ KnotplotWithGeometry ] = ComputeGeometry(knotplot)
    % recomputes the geometry - curvature , Frenet frame, of curve from points
    NP = length(knotplot.POINTS);
    
    %compute lengths
    points = knotplot.POINTS;
    nextpoints = circshift(points,-1);
    dr = nextpoints - points;
    lengths = sqrt(sum(dr.^2,2));
    
    %compute tangents
    
    tangents = zeros(NP,3);
    for i=1:NP
        dsp = lengths(i);
        dsm = lengths(incp(i,-1,NP));
        tangents(i,:) = ( dsm/(dsp*(dsp+dsm)) )* points(incp(i,1,NP),:) + ((dsp-dsm)/(dsp*dsm))* points(i,:)- (dsp/(dsm*(dsp+dsm)))* points(incp(i,-1,NP),:);
    end
    
    %compute Curvatures and Normal
    N = zeros(NP,3);
    curvatures = zeros(NP,1);
    for i=1:NP
        dsp = lengths(i);
        dsm = lengths(incp(i,-1,NP));
        kappaN(i,:) = ( dsm/(dsp*(dsp+dsm)) )* tangents(incp(i,1,NP),:) + ((dsp-dsm)/(dsp*dsm))* tangents(i,:)- (dsp/(dsm*(dsp+dsm)))* tangents(incp(i,-1,NP),:);
        curvatures(i) = norm(kappaN(i,:));
        N(i,:) = kappaN(i,:)/curvatures(i);
    end
    
    B = zeros(NP,3);
    % compute binormals
    for i=1:NP
       B(i,:) = cross(tangents(i,:),N(i,:));
    end  
    
    % compute twists
    A = knotplot.A;
    for i=1:NP
        dAds(i,:) = ( dsm/(dsp*(dsp+dsm)) )* A(incp(i,1,NP),:) + ((dsp-dsm)/(dsp*dsm))* A(i,:)- (dsp/(dsm*(dsp+dsm)))* A(incp(i,-1,NP),:);
        Twist(i) = dot(tangents(i,:), cross( A(i,:),dAds(i,:) ) );
    end  
    
    % okay we have everything
    
    KnotplotWithGeometry = knotplot;
    
    KnotplotWithGeometry.POINTS = points;
    KnotplotWithGeometry.Length = lengths;
    KnotplotWithGeometry.Curvature = curvatures;
    KnotplotWithGeometry.n = N;
    KnotplotWithGeometry.b = B;
    KnotplotWithGeometry.Twist = Twist;

end