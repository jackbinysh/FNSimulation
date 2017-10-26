#takes an obj file for a 1-D curve,made in e.g. Blender
#which is an unordered list of vertices + an edge list
# and produces an ordered an list of vertices
vertexlist = []
correctedvertexlist = []
edgelist = []

name = "component1"

TheFile = open(name+".obj","r")
Lines = TheFile.readlines()

for line in Lines[:]:
    # strip off newline
    filteredline = line.split('\n')[0]
    filteredline = filteredline.split(' ')

    if filteredline[0]=='v':
        vertexlist.append(filteredline[1:4])
    if filteredline[0]=='l':
        edgelist.append(filteredline[1:4])

# convert to edgelist to ints
edgelist = [list(map(int, x)) for x in edgelist]

# we will search this object for where to go next
row0 = [row[0] for row in edgelist] 
row1 = [row[1] for row in edgelist] 

vertex = edgelist[0][0]; 
for i in range(len(edgelist)):



    # the graph is undirected, which makes traversing it not totally trivial.

    numvertexoccurences = row0.count(vertex)

    # numvertexoccurences can be 0 or 2 (as well as boring 1) because the graph is undirected.
    # if its 2 (which will only be possible on the first one) we will just pick a direction. 

    # if its 0 then there are 2 occurnces in the 2nd column of edgelist.
    # we just traversed one of these edges - we should flip the other and then continue

    if numvertexoccurences == 0:
        indices = [i for i, x in enumerate(row1) if x == vertex];
        edge0 = edgelist[indices[0]]
        edge1 = edgelist[indices[1]]
        # these two edges  will look like, eg [60,46],[45,46]. one of these edges
        # will be the one we just traversed. We flip the untraversed edge and traverse it.
        indextoflip = (indices[1] if edge0==currentedge else indices[0])
        # ok flip it
        temp =  row0[indextoflip]
        row0[indextoflip] = row1[indextoflip]
        row1[indextoflip] = temp;
        edgelist[indextoflip][0] = row0[indextoflip]
        edgelist[indextoflip][1] = row1[indextoflip]

    currentindex = row0.index(vertex);
    currentedge = edgelist[currentindex]
    vertex = currentedge[1]
    # add the guy we just found to the corrected vertex list,remebering to 0 index
    correctedvertexlist.append(vertexlist[vertex-1])

# write the vertex list to file
outputfile = open(name+".txt", "w")

for point in correctedvertexlist: 
    string = " ".join(point)
    outputfile.write(string)
    outputfile.write("\n")

outputfile.close();



