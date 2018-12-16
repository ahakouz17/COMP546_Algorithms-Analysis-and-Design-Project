import networkx as nx
import numpy
import MySQLdb as mdb


con = mdb.connect(host='', user='', passwd='', db='')
cur = con.cursor()

#
hotspot_count = {}
try:
    cur.execute("select * from templates_hotspots")
    row = cur.fetchall()
except Exception as e:
    print e

for entry in row:
    hotspot_count[entry[0]] = entry[1]

# read the list of interfaces that was unscored by rosetta to discard if needed
unscored_interfaces = []
with open("../unscored.txt",'r') as unscored_file:
    for line in unscored_file:
        line = line.strip()
        unscored_interfaces.append(line) 


#print len(hotspot_count)
#print hotspot_count['5aviAC']

def pickRepresentative(longestClique):
    maxTemp = ''
    maxcount = 0
    unscored = False
    for t in longestClique:
        if  (hotspot_count[t] > maxcount and not unscored) or (unscored):
            if t in unscored_interfaces and maxcount == 0:
                unscored = True
                maxcount = hotspot_count[t]
                maxTemp = t
            if t not in unscored_interfaces:
                unscored = False
                maxcount = hotspot_count[t]
                maxTemp = t
    print("picked " + maxTemp + " to represent " + str(len(longestClique)) + " templates. Because it has " + str(maxcount) + " hotspot, unscored: "+ str(maxTemp in unscored_interfaces) +"\n")
    return maxTemp

adjfile = open("../adjacencyMatrices.txt", 'r')
current_entry = ""
tempIndex = []
total =0
repTotal = 0
allRepresentatives = []
outFile = open("../cliques.txt", 'w')
clique_lengths = []
while True:
    line = adjfile.readline()
    if len(line) <= 0 or line == "":
        break
    if line[:6] == "entry:":
        current_entry = line[6:].strip()
        print("Clustering " + current_entry + " group ...")
        outFile = open("cliques.txt", 'a')
        outFile.write("Entry:" + current_entry + '\n')
        outFile.close()
        continue
    if line[:10] == "templates:":
        tempIndex = line[10:].strip().split("\t")
        total += len(tempIndex)
        #print(tempIndex)
        continue
    n = len(tempIndex)
    graph = []
    row = line
    for i in range(n):
        if row == '':
            break
        row = row.strip().split("\t")
        row = [float(r.strip()) for r in row]
        graph.append(row)
        if i != n-1:
            row = adjfile.readline()
    if current_entry not in ['']:  # you can add the exceptions here and specify their threshold in the else statement
        graph = [[1 if graph[i][j] >= 0.8 else 0 for j in range(len(graph))] for i in range(len(graph))]
    else:
        graph = [[1 if graph[i][j] >= 0.8 else 0 for j in range(len(graph))] for i in range(len(graph))]

    g = numpy.matrix(graph)
    G = nx.from_numpy_matrix(g)
    cliques = []
    outFile = open("../cliques.txt", 'a')
    representatives = []
    while True:
        cliquesList = list(nx.find_cliques(G))
        if len(cliquesList) <= 0:
            break
        longestClique = max(cliquesList, key=len)
        #print(longestClique)
        cliques.append(longestClique)
        #print("A clique found with length: " + str(len(longestClique)))
        for n in longestClique:
            clique_lengths.append(len(longestClique))
            G.remove_node(n)
        rep = pickRepresentative([tempIndex[t] for t in longestClique])
        representatives.append(rep)
        outFile.write(rep + "\t1:" + str(len(longestClique)) + "\t" + "\t".join(
            [tempIndex[c] for c in longestClique]) + "\n")

    outFile.close()
    #cliques = [[tempIndex[t] for t in c] for c in cliques]
    #representatives = [c[0] for c in cliques]
    allRepresentatives = allRepresentatives + representatives
    repTotal += len(representatives)
    #outFile = open("../cliques.txt", 'a')
    #outFile.write("Entry:" + current_entry + '\n')
    #outFile.write("Templates:\t" + "\t".join(tempIndex))
    #outFile.write("\nRepresentatives(" + str(len(representatives))+ " of "+ str(len(tempIndex)) + "):\t" + "\t".join(map(str, representatives)) + "\n")
    #outFile.close()
    #print("Representatives (" + str(len(representatives))+ " of "+ str(len(tempIndex)) + "):\t" + "\t".join(map(str, representatives)))
    #print("cliques:")
    #print(cliques)
    #print("Clustered " + str(repTotal) + " from " + str(total))
    #print("All Representatives: ")
    #print(allRepresentatives)
    #print("\n")
outFile = open("../cliques.txt", 'a')
outFile.write("Clustered " + str(repTotal) + " from " + str(total)+"\n")
outFile.write("\n\nAll Representatives:"+ "\t".join(map(str, allRepresentatives)) + "\n")
outFile.write("\nAverage Clique Length= "+ str(sum(clique_lengths)/len(clique_lengths)))
outFile.close()
#print(clique_lengths)
#clique_lengthsnp = numpy.array(clique_lengths)
#unique_elements, counts_elements = numpy.unique(clique_lengthsnp, return_counts=True)
#print("\n".join(map(str, unique_elements)))
#print("\n\n\n" + "\n".join(map(str, counts_elements)))
print("Clustered " + str(repTotal) + " from " + str(total))
print("\nAverage Clique Length= "+ str(sum(clique_lengths)/len(clique_lengths)))
