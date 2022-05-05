import networkx as nx
import time
import math

def modularity(Pt, G, m):
    deg_in = 0
    deg_tot = 0
    P = G.subgraph(Pt)
    for x in P.nodes:
        deg_tot = deg_tot + G.degree(x, weight = 'weight')
        deg_in = deg_in + P.degree(x, weight = 'weight')

    return (deg_in/(2 * m)) - ( (deg_tot/(2 * m)) * (deg_tot/(2 * m)) )

def sharedEdges(sg1, sg2, G):

    sub1 = G.subgraph(sg1)
    sub2 = G.subgraph(sg2)
    combo = sg1 + sg2
    combog = G.subgraph(combo)
    hsub = list(sub1.edges()) + list(sub2.edges())
    
    hcombo = list(combog.edges())
    
            
    lastr = list(set(hcombo)^set(hsub))
    last = []
    for l in lastr:
        if l[0] < l[1]:
            last.append((l[1], l[0]))
        else:
            last.append((l[0], l[1]))
    last = set(last)        
    
    retSum = 0
    
    for i in last:
        retSum += G[i[0]][i[1]]['weight']
        
    return retSum

def phase2(groups, G):
    retLis = []
    for i in range(0, len(groups)):
        for j in range(i, len(groups)):
            weight = sharedEdges(groups[i], groups[j], G)
            if weight == 0:
                continue
            else:
                add = [i, j, weight]
                retLis.append(add)

    return retLis



def modGain(C, G, i, m):
    tC = C.copy()
    tC.append(i)

    after = modularity(tC, G, m)
    before = modularity(C, G, m) + modularity([i], G, m)

    return after - before

def modDiff(D, G, i, m):
    tD = D.copy()
    tD.remove(i)

    after = modularity(tD, G, m) + modularity([i], G, m)
    before = modularity(D, G, m)

    return after - before

def phase1(C, D, G, i, m):

    
    gain = modGain(C, G, i, m)
    diff = modDiff(D, G, i, m)
    
    return gain + diff



 

    


def avg_short_path(graphHold):
    pathlens = 0

    for i in list(graphHold.nodes):
        paths_dict = nx.single_source_dijkstra_path_length(graphHold,i)
        pathlens += sum(list(paths_dict.values()))

    return pathlens/(graphHold.number_of_nodes() * (graphHold.number_of_nodes() - 1))

def avg_cluster(graphHold):
    cluster_dict = nx.clustering(graphHold, weight = "weight")
    return sum(list(cluster_dict.values()))/graphHold.number_of_nodes()

#graph 1 
#is an undirected graph because there are 2 pairs of the same edge in the file
#this can be checked by looking at the number of lines in the file (len(contents)) and comparing it to the number of edges in the graph
#graph is also weighted
#num nodes is 638
#num edges is 18625
file_list = ['macaque-visual-cortex.txt', 'cat-cortex.txt', 'macaque-cortex.txt', 'macaque-large-scale.txt','macaque-cerebral-cortex.txt', 'coactivation-matrix.txt']
for name in file_list:
    print(name)
    t = time.process_time()
    with open(name) as f:
        if name == 'macaque-cerebral-cortex.txt':
            file = f.readlines()
            hold = []
            another_hold = []
            for line in file:
                if 'CASE' in line:
                    continue
                line_split = line.split()
                s_t_flne = []
                s_t_flne.append(line_split[2] + line_split[3])
                
                s_t_flne.append(float(line_split[4]))
                s_t_flne.append(line_split[2])
                s_t_flne.append(line_split[3])
                hold.append(s_t_flne)
                another_hold.append(line_split[2] + line_split[3])
                #print(s_t_flne)
            a_hold = set(another_hold)
            new_hold = []
            for i in a_hold:
                counter = 0
                sumf = 0
                source = ''
                target = ''
                temp = []
                for x in hold:
                    
                    if x[0] == i:
                        counter += 1
                        sumf += x[1]
                        source = x[2]
                        target = x[3]
                temp.append(source)
                temp.append(target)
                temp.append(math.log(sumf/counter) * -1)
                new_hold.append(temp)


            G = nx.Graph()

            for z in new_hold:
                G.add_edge(z[0], z[1], weight = z[2])
        else:
            contents = f.readlines()
            
            #if name == 'coactivation-matrix.txt':
            G = nx.Graph()
            #else:
                #G = nx.DiGraph()
            sum_weight = 0
            for line in contents:
                if "#tail" in line:
                    continue
                contents_list = line.split()
                #print(contents_list[0] + " - " + contents_list[1] + " - " + contents_list[2])
                sum_weight += float(contents_list[2])
                G.add_edge(contents_list[0], contents_list[1], weight = float(contents_list[2]))

        


        currMod = -1
        counter = 0
        while currMod != 0:
            m = G.size(weight = 'weight')
            groups = []
            for h in G.nodes:
                groups.append([h])
            
            graphModules = G.nodes
            done = False
            count = 0
            while not done:
                
                flag = True
                for n in graphModules:
                    incident = list(G.edges(n))
                    
                    use = []
                    
                    for d in groups:
                        if n in d:
                            use = d
                            
                    max = 0
                    hold = None
                    for s in incident: 
                        temp = []
                        if s[1] in use:
                            continue
                        for f in groups:
                            if s[1] in f:
                                
                                temp = f
                        
                        p1 = phase1(use, temp, G, s[1], m)
                        #print(p1)
                        if max < p1:
                            max = p1
                            hold = s[1]
                    #print('\n\n\n')    
                    if hold is None or max == 0:
                        count += 1
                        if count == G.number_of_nodes():
                            flag = True
                            break
                        
                        continue
                    else:
                        flag = False
                        count = 0
                        for c in groups: 
                            if hold in c:
                                c.remove(hold)
                                
                            if n in c:
                                c.append(hold)

                        groups = list(filter(None, groups))
                if flag:
                    done = True
                res = list(filter(None, groups))
                summod = 0
                for r in res:
                    summod += modularity(r, G, m)
                #print(summod)
                #print("LOOPED")
                


            res = list(filter(None, groups))
            
            summod = 0
            for r in res:
                summod += modularity(r, G, m)
            print("Number of Nodes in " + name + ": \t" + str(G.number_of_nodes()))
            print("level: " + str(counter))
            print("Modularity: " + str(summod))
            print("Modularity: " + str(nx.algorithms.community.modularity(G, res)))
            print("Communities: " + str(len(res)))
            print()
            p2 = phase2(res, G)
            print(p2)
            currMod = summod
            newG = nx.Graph()
            for s in p2:
                newG.add_edge(s[0], s[1], weight = float(s[2]))
            G = newG

            
            counter += 1

    print("Exectution time: " + str(time.process_time() - t))
    print('\n\n')


        
                
                





        #print("Number of Nodes in " + name + ": \t" + str(G.number_of_nodes()))

        #print("Number of Edges in " + name + ": \t" + str(G.number_of_edges()))

        #print("Number of Lines in " + name + ": \t" + str(len(contents)))

        #print("Weight Ratio -> closer to 1 means unweighted")
        #print(sum_weight/len(contents))
            
        #print("Avg Shortest Path in " + name + ": \t" + str(avg_short_path(G)))
        
        #print("Builtin Avg Shortest Path in " + name + ": \t" + str(nx.average_shortest_path_length(G, weight = "weight")))
        
        #print("Avg Cluster Coeff in " + name + ": \t" + str(avg_cluster(G)))

        #print("Builtin Avg Cluster Coeff in " + name + ": \t" + str(nx.average_clustering(G, weight = "weight")))

