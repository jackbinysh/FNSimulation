function [ VirtualFilament ] = ConstructSmoothedCurves(start,finish,windowsize,datadir)
    % for all the knotplots between start and end, construct the smoothed
    % filament, and write to file

    % read in initial chunk
    chunkofknotplots = {};
    chunk = [];
    for index = start:1:(start+windowsize-1)
        filename = strcat(datadir,'knotplot0_',num2str(index),'.vtk') ; 
        knotplot = CurveRead(filename);  
        chunkofknotplots = [chunkofknotplots {knotplot}];
        points = knotplot.POINTS;
        chunk = [chunk {points}];

    end
    VirtualFilament = ConstructVirtualFilament(chunk,windowsize);

    middleindex = (windowsize/2);
    middleTime = middleindex + (start-1)
    knotplotofinterest = chunkofknotplots{middleindex};
    knotplotofinterest.POINTS = VirtualFilament;
    knotplotofinterest = ComputeGeometry(knotplotofinterest);
    CurveWrite(knotplotofinterest,strcat(datadir,'Smoothedknotplot0_',num2str(middleTime),'.vtk'))


    %remove one at the start, add one at the end.
    for index = (start+windowsize):1:finish

        filename = strcat(datadir,'/knotplot0_',num2str(index),'.vtk') ; 
        knotplot = CurveRead(filename);  
        points = knotplot.POINTS;

        chunkofknotplots = [chunkofknotplots {knotplot}];
        chunkofknotplots = chunkofknotplots(2:end);
        chunk = [chunk {points}];
        chunk = chunk(2:end);
        VirtualFilament = ConstructVirtualFilament(chunk,windowsize);

        middleTime = middleTime+1;
        knotplotofinterest = chunkofknotplots{middleindex};
        knotplotofinterest.POINTS = VirtualFilament;
        knotplotofinterest = ComputeGeometry(knotplotofinterest);
        CurveWrite(knotplotofinterest,strcat(datadir,'Smoothedknotplot0_',num2str(middleTime),'.vtk'));
        middleTime
    end

end