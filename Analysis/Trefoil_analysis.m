%13/04/2017
% A script which reads in the info from a knotplot file, in ascii vtk
% format, takes individual segmets velocities and plots them

%% Read in the data
clear all
curvatures=[];
twistsquared=[];
vdotn=[];
vdotb=[];

times = [];

minvalue = 168;
maxvalue = 168;

for index = minvalue:11.2:maxvalue

    % first up, clear the knotplotdata
    knotplotdata = [];
    % grabb the filename we want to read in
    filename = strcat('/home/jackbinysh/Desktop/Trefoil/knotplot0_',num2str(index),'.vtk') ; 
    % if it exists, read it in
    if (exist (filename, 'file'))
        knotplotdata = importdata(filename);
            % have a vector of the times
            times(end+1) = index;
    end
    if( ~isempty(knotplotdata))
    
        knotplotdata = knotplotdata.textdata;
        % okay, lets get the number of points
        name ='POINTS' ;
        y = cellfun(@(x) strfind(name,x), knotplotdata,'UniformOutput',false);
        location = ~cellfun('isempty', y);
        [row,col] = find(location);
        points = knotplotdata{row,2};
        points = str2num(points);

        % now grab the data, starting from 2 lines below each field, extending
        % 'points' way down

        % scalar data
        names = {'Curvature','Torsion','Twist','Writhe','Length'};
        for i = 1:length(names)
         name = names{i};
         y = cellfun(@(x) strfind(x,name), knotplotdata,'UniformOutput',false);
        location = ~cellfun('isempty', y);
        [row,col] = find(location);
        data = knotplotdata((row+2):(row+2+points-1),1);
        data = str2double(data);
        knotplot.(name) = data;
        end

        % Vector data
        names = {'POINTS','vdotn','vdotb'};
        for i = 1:length(names)
         name = names{i};
         y = cellfun(@(x) strfind(x,name), knotplotdata,'UniformOutput',false);
        location = ~cellfun('isempty', y);
        [row,col] = find(location);
        data = knotplotdata((row+1):(row+1+points-1),1:3);
        data = str2double(data);
        knotplot.(name) = data;
        end

        % great, we have a struct, knotplot, with fields with all the data we want.
        %lets collate it all into one big matrix so we can run stats
        curvatures = vertcat(curvatures,knotplot.Curvature);
        twistsquared = vertcat(twistsquared,knotplot.Twist.^2);

        vdotnmagnitudes = sqrt(dot(knotplot.vdotn,knotplot.vdotn,2));
        vdotn = vertcat(vdotn,vdotnmagnitudes);
        
        vdotbmagnitudes = sqrt(dot(knotplot.vdotb,knotplot.vdotb,2));
        vdotb = vertcat(vdotb,vdotbmagnitudes);
    
    end
end

%% okay, from here we have some proposed linear relationships. Lets 
% perform some exploratory analysis: 

% the simplest proposed relationsship is linear in curvature and square of
% twist
figure
plot3(curvatures,twistsquared,vdotn,'o');
title ('vdotn against curve and twist squared');

figure
plot(vdotn);
title ('vdotn against length along curve');

figure
plot3(curvatures,twistsquared,vdotb,'o');
title ('vdotb against curve and twist squared');

figure
plot(vdotb);
title ('vdotb against length along curve');

% from these plots, we see that vdotb appears to have more stucture in it
% than vdotn. Accordingly, a linear relation may do a better job for vdot n
% than it does for vdotb. For example, lets look at the residuals for any
% proposed linear relation for vdotb

%normal residuals
 [sf,goodness,output] = fit([curvatures,twistsquared],vdotn,'poly11');
 % normal surface fit
 figure
 plot(sf,[curvatures,twistsquared],vdotn);
 title ('vdotn fit');
 figure
 plot(output.residuals)
 title ('vdotn residuals along length');
 
 
% binormal residuals
 [sf,goodness,output] = fit([curvatures,twistsquared],vdotb,'poly11');
 % binormal surface fit
 figure
 plot(sf,[curvatures,twistsquared],vdotb);
 title ('vdotb fit');
  figure
 plot(output.residuals)
 title ('vdotb residuals along length');
 
 %% the next level up is to include dw/ds as a predictor.
 % first of all, we'd better get a decent version of dw/ds
%  Twist = knotplot.Twist;
%  LaggedTwist = Twist;
% LaggedTwist(2:end) = Twist(1:end-1);
% LaggedTwist(1)=Twist(end);
% deltaw = (Twist - LaggedTwist);
% deltawoverdeltas = deltaw./knotplot.Length;
% this is kind of rough though. lets smooth it out in the freq domain:
data = Twist;
Fs = 1;
T = 1/Fs;
L = length(data);
Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% ok, now we want to filter frequencies greater than some cutoff. since our
% signal is even, we want to frequency mask to respect conjugate symettry
% for 2:end
cutoff = 0.1;
filteredtransform = Y;
for i = 1:(L+1)/2
    if(f(i)>cutoff) 
        filteredtransform(i) = 0;
        filteredtransform(end-i+2) = 0; % hit the conjugate symettric part too
    end
end
filtereddata = ifft(filteredtransform);
plot((filtereddata));



% Y is the fft, its conjugate symettric
% lets cut the peak with the oscillation frequency in it out
% we located it at f = 0.0875 and beyond
mask = f>0.7;
% we need to double it, due to the conjugate symettry
doubledmask = [mask(1:end-1) fliplr(mask)];

% now construct the filtered transform
filteredtransform = zeros(L,1);
for i = 1:L
    if(doubledmask(i) ==0) 
        filteredtransform(i) = Y(i);
    end
end
filtereddata = ifft(filteredtransform);
figure
plot(real(filtereddata));
Twist = filtereddata;

LaggedTwist = Twist
LaggedTwist(2:end) = Twist(1:end-1);
LaggedTwist(1)=Twist(end);
deltaw = (Twist - LaggedTwist);
dwds = deltaw./knotplot.Length

% okay, so we've got a smoothed out twist and a smoothed out dwds.
% lets see, with these predictors, whether we get a better fit.
% the proposed model is, say, 
%Vb = alpha k + beta w^2 + gamma dw/ds
Predictors = [curvatures,twistsquared,dw/ds];

model = fitlm(Predictors,vdotn,'linear');
 







% 
%         hold on
%         magnitudes = sqrt(dot(knotplot.vdotn,knotplot.vdotn,2));
%         plot3(knotplot.Curvature,knotplot.Twist.^2,magnitudes,'o');
% 
 %[sf,goodness,output] = fit([curvatures,twistsquared],vdotb,'poly11');
%plot(sf,[curvatures,twistsquared],vdotn);
 %plot(sf,[knotplot.Curvature,knotplot.Twist.^2],vdotb,'Style','Residuals');
% % scatter3(points(:,1),points(:,2),points(:,3),1,output.residuals);
% plot3(curvatures,twistsquared,vdotb,'o');
%  [sf,goodness,output] = fit([curvatures,twistsquared],vdotb,'poly11');
% plot(sf,[curvatures,twistsquared],vdotb);
% 
%  
%  
% plot3(curvatures,twistsquared,vdotn,'o');
% plot(sf,[curvatures,twistsquared],vdotb);



% 
% magnitudes = sqrt(dot(knotplot.vdotn,knotplot.vdotn,2));
% figure
% hold on
% for i = 1: length(magnitudes)
%      plot(knotplot.Curvature(i),magnitudes(i),'o');
%      hold on
% end
% 
% prediction = 1.3465* knotplot.Curvature;
% plot(knotplot.Curvature,(1.3465* knotplot.Curvature)+0.0022);
% 
%         % now lets take averages and so on
% 
%         % average curvature
%         averagecurvature(end+1) = sum(knotplot.Curvature)/length(knotplot.Curvature);
%         % average velocity in the normal directoin
%         magnitudes = sqrt(dot(knotplot.vdotn,knotplot.vdotn,2));
%         averagevdotn(end+1) = sum(magnitudes)/length(magnitudes);
%         % average velocity in the binormal directoin
%         magnitudes = sqrt(dot(knotplot.vdotb,knotplot.vdotb,2));
% %         averagevdotb(end+1) = sum(magnitudes)/length(magnitudes);
%  LaggedTwist = Twist
% LaggedTwist(2:end) = Twist(1:end-1);
% LaggedTwist(1)=Twist(end);
% deltaw = (Twist - LaggedTwist);
% deltawoverdeltas = deltaw./knotplot.Length
% 
% smooth(deltawoverdeltas)
% 
%  scatter3(curvatures,twistsquared,deltawoverdeltas,vdotb);






    