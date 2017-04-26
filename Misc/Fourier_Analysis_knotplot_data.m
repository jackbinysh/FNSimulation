data = VarName4;
plot(data)

% do it on first 1000 timesteps
data = data(1:800)
Fs = 1;
T = 1/Fs;
L = length(data);
Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:(end-1)) = 2*P1(2:(end-1));
f = Fs*(0:(L/2))/L;
figure
plot(f,P1) 

%do it on the last 300 timesteps
data = VarName4;
data = data(4000:end);
Fs = 1;
T = 1/Fs;
L = length(data);
Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:(end-1)) = 2*P1(2:(end-1));
f = Fs*(0:(L/2))/L;
figure
plot(f,P1) 


% okay lets kill the oscillations in length of these frequencies
%ASSUME DATA IS OF ODD LENGTH
data = VarName4;
plot(data)

Fs = 1;
T = 1/Fs;
L = length(data);
Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:(L+1)/2);
P1(2:(end-1)) = 2*P1(2:(end-1));
f = Fs*(0:((L)/2))/L;
figure
plot(f,P1) 

% Y is the fft, its conjugate symettric
% lets cut the peak with the oscillation frequency in it out
% we located it at f = 0.0875 and beyond
mask = f>0.08;
% we need to double it, due to the conjugate symettry
doubledmask = [mask(1:end-1) fliplr(mask)]
filteredtransform = (1-doubledmask)'.*Y;
filtereddata = ifft(filteredtransform);




