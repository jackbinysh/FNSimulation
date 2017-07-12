
% IMPORT THE ADJ MATRIX

data =  adjacencymatrix;

 % to reaplace nans with zeroes
 data(isnan(data))=0
 
 % to get edge list
 edgelist = adj2edgeL(data);
 
 % write it
 csvwrite('edges.csv',edgelist)
 
 % NOW IMPORT THE LABELS
 
 
 labels = adjacencymatrix;
 
 augmentedlabels = horzcat(num2cell(1:length(labels))', labels)
 
 % to write the csv
fid = fopen('nodes.csv','wt');
 if fid>0
     for k=1:size(augmentedlabels,1)
         fprintf(fid,'%i,%s\n',augmentedlabels{k,:});
     end
     fclose(fid);
 end