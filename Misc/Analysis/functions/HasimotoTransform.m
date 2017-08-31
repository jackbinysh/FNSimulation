function [ phi ] = HasimotoTransform( K,Tau,lengths,offset)
%HASIMOTOTRANSFORM compute the complex curvature of a curve, a la hasimoto

    clear j;
    ArcLength = cumsum(lengths);

    % lets construct an extended version of our periodic curve
    TotalLength = ArcLength(end);

    ExtendedLength = vertcat(ArcLength-TotalLength, ArcLength, TotalLength + ArcLength);
    ExtendedK = vertcat(K,K,K);
    ExtendedTau = vertcat(Tau,Tau,Tau);
    Extendedlengths = vertcat(lengths,lengths,lengths);

    offsetK = interp1(ExtendedLength, ExtendedK, ArcLength - offset);
    offsetTau = interp1(ExtendedLength, ExtendedTau, ArcLength - offset);
    
    IntegratedTau = cumsum(offsetTau.*lengths);
    phi = offsetK.*exp(IntegratedTau*j);   
   

end

