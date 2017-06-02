%13/04/2017
% A script which reads in the info from a knotplot file, in ascii vtk
% format, takes averages and plots them. this is useful for data from the
% unknot.


%% Read in the data
clear all
averagecurvature=[];
averagevdotn=[];
averagevdotb=[];
times = [];

minvalue = 67.2;
maxvalue = 1052.8;

for index = minvalue:11.2:maxvalue

    % first up, clear the knotplotdata
    knotplotdata = [];
    % grabb the filename we want to read in
    filename = strcat('/home/jackbinysh/Code/Knot_surface_calc/knotplot0_',num2str(index),'.vtk') ; 
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
        names = {'Curvature','Torsion','Twist','Writhe'};
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
        names = {'vdotn','vdotb'};
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
        % now lets take averages and so on

        % average curvature
        averagecurvature(end+1) = sum(knotplot.Curvature)/length(knotplot.Curvature);
        % average velocity in the normal directoin
        magnitudes = sqrt(dot(knotplot.vdotn,knotplot.vdotn,2));
        averagevdotn(end+1) = sum(magnitudes)/length(magnitudes);
        % average velocity in the binormal directoin
        magnitudes = sqrt(dot(knotplot.vdotb,knotplot.vdotb,2));
        averagevdotb(end+1) = sum(magnitudes)/length(magnitudes);
    end
end

%% fitting to the curves and plotting
% vdotn against curvature
x = linspace(0,max(averagecurvature),10);
vdotnfit = polyfit(averagecurvature,averagevdotn,1);
plot(averagecurvature,averagevdotn);
hold on
plot(x, polyval(vdotnfit,x))
title ('vdotn against curvature');

% vdotb against curvature
figure;
x = linspace(0,max(averagecurvature),10);
vdotbfit = polyfit(averagecurvature,averagevdotb,1);
plot(averagecurvature,averagevdotb);
hold on
plot(x, polyval(vdotbfit,x));
title ('vdotb against curvature');

% lets look at the radius against time
% the propsed relation is r^2(t) = r^2(0) - 2(alpha)t, where alpha is the
% gradient of the vdot n against curvature fit
figure
plot(times,(1./averagecurvature).^2);
curvaturefit = polyfit(times,(1./averagecurvature).^2,1);
y1 = polyval(curvaturefit,times);
hold on
plot(times,y1);
title ('r^2 against time');

% we can also look at the other quantities against time
figure
plot(times,averagevdotb);
title ('vdotb against time');

% we can also look at the other quantities against time
figure
plot(times,averagevdotn);
title ('vdotn against time');

% radius against time
figure
plot(times,(1./averagecurvature));
curvaturefit = polyfit(times,(1./averagecurvature),1);
y1 = polyval(curvaturefit,times);
hold on
plot(times,y1);
title ('r against time');






 

    
    