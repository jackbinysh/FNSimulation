function [ knotplot ] = CurveRead( filename )
%CURVEREAD Read in a knotplot file

    % if it exists, read it in
    if (exist (filename, 'file'))
        knotplotdata = importdata(filename);
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
        names = {'POINTS','n','b','vdotn','vdotb'};
        for i = 1:length(names)
         name = names{i};
        % y = cellfun(@(x) strcmp(x,name), knotplotdata,'UniformOutput',false);
        y = cellfun(@(x) strcmp(x,name), knotplotdata);
        %location = ~cellfun('isempty', y);
        [row,col] = find(y);
        data = knotplotdata((row+1):(row+1+points-1),1:3);
        data = str2double(data);
        knotplot.(name) = data;
        end
end

