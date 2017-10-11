function [ result ] = incp(i,p,N)
%INCP returns the periodic incremenet of i by p, for a thing of period N.
% remeber, matlab is 1 indexed, so its a little different to C

    if(i+p<=0)
        result = N+i+p;
    else
        result = 1+mod(i+p-1,N);
    end
end
    