function [ VirtualFilament ] = ConstructVirtualFilament(listofcurves,windowsize)
    % a function which takes a curve, and a window a curves taken at
    % neighboring time points, and from the trajectories of the points
    % constructs a smoothed "virtual filament" for the original curve.
    % the list of curves shoud have the structure xxx ... xxC MIDDLE xx...xxx, where 
    % C is the curve of interest - i.e EVEN length and curve one behind the
    % "middle"

    traj = [];

    % first of all, construct the trajectory data.
    curve = listofcurves{(windowsize)/2};

    % do the ones ahead of the current curve
    oldpoints = curve;
    for index = ((windowsize)/2):1:windowsize

        points = listofcurves{index};  
        closestapproach = MatchPoints(oldpoints,points);

        traj = [traj {closestapproach}];
        oldpoints=closestapproach;
    end

    % now do the ones behind
    oldpoints = curve;
    for index = ((windowsize/2)-1):-1:1

        points = listofcurves{index} ; 
        closestapproach = MatchPoints(oldpoints,points);
        traj = [ {closestapproach} traj];
        oldpoints=closestapproach;
    end
    
    %reshape this array a little so each cell is a trajectory, not the
    %curve.
    trajectories = {};
    for i = 1:length(traj{1})
        trajectory = [];
        for j = 1:length(traj)
            points = traj{j}(i,:);
            trajectory = [trajectory; points];
        end
        trajectories{i} = trajectory;
    end
  
    % okay now we have the trajectories as a list of cells. Lets do the
    % smoothing
    
    % grab each trajectory, and smooth it using a filter
    
    filteredtrajectories = trajectories;
    for trajindex = 1:length(trajectories)
        for q = 1:3
            data = trajectories{trajindex}(:,q);
            Fs = 1;
            T = 1/Fs;
            L = length(data);
            Y = fft(data);
            f = Fs*(0:(L/2))/L;

            % a butterworth filter
            cutoff = 0.3/11.2;
            Yprime = zeros(length(Y),1);
            for i = 1:(( length(data)/2 ) + 1) % stupid matlab 1 indexing.
                filter = 1/sqrt( 1+ ( f(i)/cutoff )^6 );
                Yprime(i) = filter * Y(i);
            end
            n= length(data);
            % enforce conjugacy
            for i = 2:( length(data)/2 )
                Yprime(n-(i-2))=conj(Yprime(i));
            end
            filteredtrajectories{trajindex}(:,q) =  ifft(Yprime);
        end
    end
    
% finally, grab the middle filament
    VirtualFilament = filteredtrajectories{1};
    for i = 1:length(filteredtrajectories)
        VirtualFilament(i,:) = filteredtrajectories{i}((windowsize)/2,:);
    end
    
end