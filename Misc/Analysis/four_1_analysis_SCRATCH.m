%13/04/2017, Jack Binysh


clear all;

% read in parameters
minvalue = 9001;
maxvalue = 9100;
increment = 1;

knotplots = {};
for index = minvalue:increment:maxvalue
    filename = strcat('/home/jackbinysh/Code/Knot_surface_calc/4_1_results/knotplot0_',num2str(index),'.vtk') ; 
    knotplot = CurveRead(filename);
    knotplots = horzcat(knotplot,knotplots);
end

CorrectOffset =0;
offsetK = [];
offsetTau = [];
offsetLength = [];
offsetr = {};
offsetA = {};
previousstart=zeros(1,3);
TotalLengths=[];
for i = 1:length(knotplots)
 
    if(i==1)
        knotplot = knotplots(i);
        FirstArcLength = cumsum(knotplot.Length);
        FirstTotalLength=FirstArcLength(end);
        FirstK = knotplot.Curvature/FirstTotalLength;
        FirstArcLength = FirstArcLength/FirstTotalLength; 
        r = knotplot.POINTS/FirstTotalLength;
        x = r(:,1);
        y = r(:,2);
        z = r(:,3);
        previousstart = [x(1),y(1),z(1)];
        TotalLengths(1) = FirstTotalLength;
    end  
    
    knotplot = knotplots(i);
    ArcLength = cumsum(knotplot.Length);
    TotalLength=ArcLength(end);
    TotalLengths(i) = TotalLength;
    % from here, EVERYTHING IS SCALED
    Tau = knotplot.Torsion/TotalLength;
    K = knotplot.Curvature/TotalLength;
    r = knotplot.POINTS/TotalLength;
    x = r(:,1);
    y = r(:,2);
    z = r(:,3);
    A = knotplot.A;
    Ax = A(:,1);
    Ay = A(:,2);
    Az = A(:,3);
    ArcLength = ArcLength/TotalLength; 
    
    ExtendedLength = vertcat(ArcLength-1, ArcLength, ArcLength + 1);
    ExtendedK = vertcat(K,K,K);
    ExtendedTau = vertcat(Tau,Tau,Tau);
    Extendedx = vertcat(x,x,x);
    Extendedy = vertcat(y,y,y);
    Extendedz = vertcat(z,z,z);
    ExtendedAx = vertcat(Ax,Ax,Ax);
    ExtendedAy = vertcat(Ay,Ay,Ay);
    ExtendedAz = vertcat(Az,Az,Az);

    % error function
    myfunc = @(offset)CurveOffsetError(offset,FirstK,FirstArcLength, K,ArcLength);
    trialoffsets = linspace(-0.6,0.6,100);
    minima = arrayfun(@(q)fminsearch(myfunc,q),trialoffsets);    

    % ok so we have these candidate minima. But some will be local and
    % others global, and only one should match closely to the physical
    % start position of the previous curve.
    
    % lets just coarsely strip off the ones that are clearly wrong
    minima = minima(myfunc(minima)<(min(myfunc(minima))+ max(myfunc(minima)))/2);
    
    % of these guys, only one set should match physical space well. lets
    % check!  
    trial = zeros(1,3);
    Error = zeros(length(minima),1);
    for i = 1:length(minima)

        trial(1)=  interp1(ExtendedLength, Extendedx,- minima(i));
        trial(2)=  interp1(ExtendedLength, Extendedy,- minima(i));
        trial(3)=  interp1(ExtendedLength, Extendedz,- minima(i));
        Error(i) = norm(trial - previousstart);

    end
    [~,location] = min(Error);
    offset = minima(location);

    
    offsetK = horzcat(offsetK, interp1(ExtendedLength, ExtendedK, FirstArcLength - offset));
    offsetTau = horzcat(offsetTau, interp1(ExtendedLength, ExtendedTau, FirstArcLength - offset));
    offsetLength = horzcat(offsetLength, interp1(ExtendedLength, ExtendedLength, FirstArcLength - offset));
    offsetx=  interp1(ExtendedLength, Extendedx, FirstArcLength - offset);
    offsety=  interp1(ExtendedLength, Extendedy, FirstArcLength - offset);
    offsetz=  interp1(ExtendedLength, Extendedz, FirstArcLength - offset);
    offsetr = horzcat(offsetr,num2cell([offsetx,offsety,offsetz],2));
    offsetAx=  interp1(ExtendedLength, ExtendedAx, FirstArcLength - offset);
    offsetAy=  interp1(ExtendedLength, ExtendedAy, FirstArcLength - offset);
    offsetAz=  interp1(ExtendedLength, ExtendedAz, FirstArcLength - offset);
    offsetA = horzcat(offsetA,num2cell([offsetAx,offsetAy,offsetAz],2));
    
  % keep track of the first element of offsetr.
    previousstart = [offsetx(1),offsety(1),offsetz(1)];
    
    
end

figure
hold on
for i= 1:size(offsetK,2)
 plot(1:size(offsetK,2),offsetK(i,:)*TotalLengths(i))
end
figure
 plot(1:size(offsetK,2),offsetK(end-3,:)*TotalLengths(end-3))

figure
hold on
for i= 1:size(offsetK,2)
data = offsetK(i,:);
Fs = 1;
T = 1/Fs;
L = length(data);
Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:(L/2+1));
P1(2:(end-1)) = 2*P1(2:(end-1));
f = Fs*(0:(L/2))/L;
% I will normalise by the first peaks size
x = P1(f>0.05 & f < 0.12);
plot(f,P1/max(x) )
xlim([0.05 0.3])
end

figure
hold on
for i= 1
    curve = cell2mat(offsetr(:,i));
    scatter3( curve(:,1), curve(:,2), curve(:,3))
end
scatter3( curve(end-3,1), curve(end-3,2), curve(end-3,3))

hold on
for i= size(offsetr,1)-3
    traj =cell2mat(offsetr(i,:)');
   traj(:,1) =traj(:,1).*TotalLengths';
    traj(:,2)=traj(:,2).*TotalLengths';
   traj(:,3) =traj(:,3).*TotalLengths';
   plot3(traj(:,1),traj(:,2),traj(:,3))
% data =  traj(1:200,1);
% data = data - data(1);
% norm = sum(data.^2);
% plot(data) ;
% % Fs = 1;
% % T = 1/Fs;
% % L = length(data);
% % Y = fft(data);
% % P2 = abs(Y/L);
% % P1 = P2(1:L/2+1);
% % P1(2:(end-1)) = 2*P1(2:(end-1));
% % f = Fs*(0:(L/2))/L;
% % P1 = P1
% % plot(f,P1/max(P1)) ;
end

hold on
for i= 1:12
 plot(1:length(offsetK(:,i)),offsetK(:,i)*TotalLengths(i))
 pause(0.5)
end

hold on
for i= 1:12
 plot(FirstArcLength,offsetK(:,i))
 ylim([1*10^-4 6*10^-4 ])
 pause(0.5)
end



for i= 1:12
 plot(FirstArcLength,offsetK(:,i))
 ylim([1*10^-4 6*10^-4 ])
 pause(0.5)
end




% 
% % 
% data = TotalLengths;
% Fs = 1;
% T = 1/Fs;
% L = length(data);
% Y = fft(data);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:(end-1)) = 2*P1(2:(end-1));
% f = Fs*(0:(L/2))/L;
% figure
% plot(f,P1) 

% dA = {};
% for i= 1:(size(offsetA,2)-1)
%     
%     A = cell2mat(offsetA(:,i));
%     nextA = cell2mat(offsetA(:,i+1)); 
%     deltaA = nextA - A;
%     dA = horzcat(dA,num2cell(deltaA,2));
% end
% 
% 
% 
% summednorms = 0;
% hold on
% for i= 1:140
%     A = cell2mat(dA(i,:)');
%     norms = arrayfun(@(idx) norm(A(idx,:)), 1:size(A,1));
%     plot(1:length(norms),norms)
% end
% 
% 
% 
% summednorms = 0;
% hold on
% for i= 1:50
%     curve = cell2mat(offsetr(:,i));
%     A = cell2mat(dA(:,i));
%     norms = arrayfun(@(idx) norm(A(idx,:)), 1:size(A,1));
%     summednorms = norms+summednorms;  
% end
% plot(FirstArcLength,summednorms/50)
% 
% hold on
% for i= 1:2
%     curve = cell2mat(offsetr(:,i));
%     A = cell2mat(offsetA(:,i));
%     quiver3(curve(:,1),curve(:,2),curve(:,3),A(:,1),A(:,2),A(:,3),0.7)
%     xlim([-0.3 -0.1 ])
%     ylim([-0.12 0.04 ])
%     zlim([-0.28 -0.08 ])
%     view([1 0 0])
%     pause(0.5)
% end
% 
% hold on
% for i= 1:size(offsetr,2)
%     A = cell2mat(offsetA(:,i));
%  plot(offsetK(:,i),A(:,1))
% pause(0.5)
% end
% % 
% % % % 
% % 
% hold on
% for i= 1:size(offsetr,2)
%  plot(FirstArcLength*TotalLengths(i),offsetK(:,i))
%  ylim([1*10^-4 6*10^-4 ])
% pause(0.5)
% end
%  hold on
% for i= 1:size(offsetr,2)
%  plot(FirstArcLength,offsetTau(:,i))
% pause(0.5)
%  end
% % % % 
% % % % hold on
% % % % for i= 1:size(offsetr,1)
% % % %  plot(1:51,offsetK(i,:))
% % % % pause(0.5)
% % % % end
% % % % 
% % 
% hold on
% for i= 1:11
%     curve = cell2mat(offsetr(:,i));
%     plot3(curve(:,1),curve(:,2),curve(:,3))
%     pause(0.5)
% end


% 
% hold on
% for i= 1:size(offsetr,1)
%     traj =cell2mat(offsetr(i,:)');
%    traj(:,1) =traj(:,1).*TotalLengths';
%     traj(:,2)=traj(:,2).*TotalLengths';
%    traj(:,3) =traj(:,3).*TotalLengths';
%    plot3(traj(:,1),traj(:,2),traj(:,3))
% % data =  traj(1:200,1);
% % data = data - data(1);
% % norm = sum(data.^2);
% % plot(data) ;
% % % Fs = 1;
% % % T = 1/Fs;
% % % L = length(data);
% % % Y = fft(data);
% % % P2 = abs(Y/L);
% % % P1 = P2(1:L/2+1);
% % % P1(2:(end-1)) = 2*P1(2:(end-1));
% % % f = Fs*(0:(L/2))/L;
% % % P1 = P1
% % % plot(f,P1/max(P1)) ;
% end
% end
% % 
% % 
% % 
% % 
% % 
% % hold on
% % curve = cell2mat(offsetr(:,1));
% % plot3(curve(:,1),curve(:,2),curve(:,3))

% 
% loc3 = floor(interp1(FirstArcLength,1:length(FirstArcLength), 0.1));
% loc4 = floor(interp1(FirstArcLength,1:length(FirstArcLength), 0.78));

