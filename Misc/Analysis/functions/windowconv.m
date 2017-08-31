function knotcurve = windowconv(points)

%cutoff freq
fc = 0.002;
%width of attenuation
b = 0.001;

N = ceil(4/b);
if mod(N,2) == 0
    N = N+1;
end
w = blackman(N);
h = sinc(2*fc*((0:N-1) - (N-1)/2))'.* w;
h = h/sum(h);

% make signal periodic
numreps = ( 2.*round(((N/length(points))+1)/2)-1  )+ 2;
Length = length(points);
reppoints = repmat(points,numreps,1);

%convolve
for i = 1:3
    x(:,i) = conv(reppoints(:,i),h,'same');
    knotcurve(:,i) = x(floor(numreps/2)*Length:ceil(numreps/2)*Length,i);
end

end