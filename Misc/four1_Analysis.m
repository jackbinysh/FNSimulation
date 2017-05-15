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

minvalue = 9990;
maxvalue = 9990;

for index = minvalue:11.2:maxvalue

    % first up, clear the knotplotdata
    knotplotdata = [];
    % grabb the filename we want to read in
    filename = strcat('/home/jackbinysh/Desktop/4_1/four1/knotplot0_',num2str(index),'.vtk') ; 
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

%% time for a little geometry. lets take the point cloud, get the covariance and diagonalize it.
% this will give the symettry axes.
ri = knotplot.POINTS;
mi = knotplot.Length;
miri = ri;
for i = 1:length(points)
    miri(i,:) = ri(i,:) *knotplot.Length(i);
end
com = sum(miri)/sum(mi);
% centre them
for i = 1:length(points)
    ri(i,:) = ri(i,:) -com;
end
r = cov(miri);
[U,S, V] = svd(r);

rotatedpoints = ri;
for i = 1:length(points)
    rotatedpoints(i,:) = (V' *ri(i,:)')';
end
scatter3(rotatedpoints(:,1),rotatedpoints(:,2),rotatedpoints(:,3))
scatter3(points(:,1),points(:,2),points(:,3))


segmentlengths = knotplot.Length;
arclength = zeros(length(segmentlengths),1)';
for i = 1:length(segmentlengths)-1
arclength(i+1) = arclength(i) + segmentlengths(i);
end

%% run fourier transforms of the curvature, position etc, to pick out frequencies.
data = rotatedpoints(:,1);
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

