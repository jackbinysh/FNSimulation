clear all;

directories = dir;
legendinfo= []; 
unknottingtime = [];

%% 
for i = 3:length(directories)
    myName = directories(i).name;
    
    if(isempty(strfind(myName, '.m')) && ~strcmp(myName(end),'u'))
       
        filename = ['/home/jackbinysh/Desktop/knotunknotsimulations_periodic_bc/' ,myName '/globaldata_0.txt' ];
       fileid = fopen(filename);
       T = textscan(fileid,'%f %f %f %f');
       data = cell2mat(T);        
        % get lengths
        lengths = data(:,2);
        %get normalised lenghts
        normalisedlengths =lengths/lengths(1);
        
        % get the unknotting times
        legendinfo{end+1}=[myName]; 
        hold on;
        figure(1)
        plot(lengths,'*');
        figure(2)
        hold on;
        plot(normalisedlengths,'*');    
    end
end
figure(1)
legend(legendinfo)
figure(2)
legend(legendinfo)


%% 

for i = 3:length(directories)
    myName = directories(i).name;
    
    if(isempty(strfind(myName, '.m')) && ~strcmp(myName(end),'u') && ~strcmp(myName,'eight11u')&& ~strcmp(myName,'eight20u'))
       
        filename = ['/home/jackbinysh/Desktop/knotunknotsimulations/' ,myName '/fistcomponentwrithe.txt' ];
        T = table2array(readtable(filename,'Delimiter', ' '));
        T = T(:,1:4);
        [values,order] = sort(T(:,1));
        data = T(order,:);
        
        % get lengths
        lengths = data(:,4);
        %get normalised lenghts
        normalisedlengths =lengths/length(1);
        
        % get the unknotting times
        
        if(~isempty(find(lengths<50, 1,'first')))
            unknottingtime(end+1,1) = find(lengths<50, 1,'first');
            unknottingtime(end,2) = lengths(1);
            legendinfo{end+1}=[myName]; 
            hold on;
            figure(1)
            plot(lengths);
            figure(2)
            hold on;
            plot(lengths/lengths(1));
        end       
    end
end
figure(1)
legend(legendinfo)
figure(2)
legend(legendinfo)

% % 
% % 
% % 
figure; 
hold on;
for i = 1:length(unknottingtime(:,1))
plot((unknottingtime(i,2)),(unknottingtime(i,1))/(unknottingtime(i,2)),'o')
end
legend(legendinfo)

figure; 
hold on;
for i = 1:length(unknottingtime(:,1))
plot((unknottingtime(i,2)),(unknottingtime(i,1)),'o')
end
legend(legendinfo)
% 
x = linspace(0,max(unknottingtime(:,2)),10);
p = polyfit(unknottingtime(:,2),unknottingtime(:,1),1)
plot(unknottingtime(:,2), polyval(p,unknottingtime(:,2) ) )
q = polyfit(unknottingtime(:,2),unknottingtime(:,1),2);
plot(x, polyval(q,x))
% % 
% % %%%% zykov thinks:




