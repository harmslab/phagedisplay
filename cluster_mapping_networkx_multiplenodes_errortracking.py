## Map sequences to clusters and protein data sets -- for generating network files. Takes inputs of sequences that
## belong to a cluster and tests for origin from A5 or A6 reference dataset, then builds a networkx object. 
import numpy as np
import pickle
import glob
import networkx as nx


#bring in the S100A5 and S100A6 reference files and build dictionaries of the sequences, for searching. 
A5_0_file = "hA5-0_reference.ref"
A6_0_file = "hA6-0_reference.ref"
A6_1_file = "hA6-1_reference.ref"
aA5A6_0_file = "aA5A6-0_reference.ref"
aA5A6_1_file = "aA5A6-1_reference.ref"
biotin_file = "biotin.ref"

with open(A5_0_file, "r") as A5_0_input:
    A5list_0 = [line.strip() for line in A5_0_input.readlines()]
#print(A5list_0)
with open(A6_0_file, "r") as A6_0_input:
    A6list_0 = [line.strip() for line in A6_0_input.readlines()]
#print(A6list_0)
with open(A6_1_file, "r") as A6_1_input:
    A6list_1 = [line.strip() for line in A6_1_input.readlines()]
#print(A6list_1)
with open(aA5A6_0_file, "r") as aA5A6_0_input:
    aA5A6list_0 = [line.strip() for line in aA5A6_0_input.readlines()]
#print(aA5A6list_0)
with open(aA5A6_1_file, "r") as aA5A6_1_input:
    aA5A6list_1 = [line.strip() for line in aA5A6_1_input.readlines()]
#print(aA5A6list_1)
with open(biotin_file, "r") as biotin_input:
    biotinlist = [line.strip() for line in biotin_input.readlines()]
#print(aA5A6list_1)

###define a function to check membership of cluster sequences to protein reference files and output to files
def mapping(index, filename):
    
    #lists to hold sequences that contribure to edges
    A5_0_edges = list()
    A6_0_edges = list()
    A6_1_edges = list()
    aA5A6_0_edges = list()
    aA5A6_1_edges = list()
    biotin_edges = list() 
        
    with open(filename, 'r') as f:
        output = [line.strip() for line in f.readlines()]
    #print(output)

    all_seqs = A5list_0+A6list_0+A6list_1+aA5A6list_0+aA5A6list_1+biotinlist
    
    for key in output:
        if key in A5list_0:
            A5_0_edges.append(key)

        if key in A6list_0:
            A6_0_edges.append(key)

        if key in A6list_1:
            A6_1_edges.append(key)

        if key in aA5A6list_0:
            aA5A6_0_edges.append(key)

        if key in aA5A6list_1:
            aA5A6_1_edges.append(key)
            
        if key in biotinlist:
            biotin_edges.append(key)
            
        if key not in all_seqs:
            missed.append(filename + key) 
    
            

    attributes = {
        "members" : output,
        "a5_0_members": A5_0_edges,
        "a6_0_members": A6_0_edges,
        "a6_1_members": A6_1_edges,
        "aA5A6_0_members": aA5A6_0_edges,
        "aA5A6_1_members": aA5A6_1_edges, 
        "biotin_members": biotin_edges
    }
    
    print()
    # Construct edges, between the clusters and the A5 (index 0) and A6 (index 1) dummy nodes. 
    edges = list()
    try:
        edges.append((index, 0, 1.0/len(attributes["a5_0_members"])))
    except:
        pass
    try:
        edges.append((index, 1, 1.0/len(attributes["a6_0_members"])))
    except:
        pass
    try:
        edges.append((index, 2, 1.0/len(attributes["a6_1_members"])))
    except:
        pass
    try:
        edges.append((index, 3, 1.0/len(attributes["aA5A6_0_members"])))
    except:
        pass
    try:
        edges.append((index, 4, 1.0/len(attributes["aA5A6_1_members"])))
    except:
        pass
    try:
        edges.append((index, 5, 1.0/len(attributes["biotin_members"])))
    except:
        pass
    
    return attributes, edges

#open an empty networkx graph object, then append the nodes and edges from inside the mapping function
file_extension = "*.fasta"
G = nx.Graph()
G.add_node(0, name="A5_0")
G.add_node(1, name="A6_0")
G.add_node(2, name="A6_1")
G.add_node(3, name="aA5A6_0")
G.add_node(4, name="aA5A6_1")
G.add_node(5, name="biotin")

missed = list() #open a list to put missed sequences/errors into
counter = 6

for file in glob.glob(file_extension)[0:]:
    
    node_attr, edges = mapping(counter, file) 
    G.add_node(counter, **node_attr)
    G.add_weighted_edges_from(edges)
    counter += 1
if len(missed) > 0:  
	np.savetxt("missed_sequences.missed", missed, fmt="%s")
#saves dict as pickle
nx.write_gpickle(G, open( "network_graph.p", "wb"))

###Plotting

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

pos=nx.spring_layout(G)
nodelist=list(G)
#print(nodelist)
plt.figure(figsize=[10,10])
nx.draw_networkx_nodes(G,pos,
                       nodelist=nodelist[6:],
                       node_color='gray',
                       node_size=100,
                   alpha=0.8)
nx.draw_networkx_nodes(G,pos,
                       nodelist=nodelist[0:1],
                       node_color='m',
                       node_size=1000,
                   alpha=0.8)
nx.draw_networkx_nodes(G,pos,
                       nodelist=nodelist[1:3],
                       node_color='orange',
                       node_size=1000,
                   alpha=0.8)
nx.draw_networkx_nodes(G,pos,
                       nodelist=nodelist[3:5],
                       node_color='lightgreen',
                       node_size=1000,
                   alpha=0.8)
nx.draw_networkx_nodes(G,pos,
                       nodelist=nodelist[5:6],
                       node_color='white',
                       node_size=1000,
                   alpha=0.8)

# edges
nx.draw_networkx_edges(G,pos,width=0.5,alpha=0.5,edge_color='darkgray')

# labels
labels={}
labels[0]=r'A50'
labels[1]=r'A60'
labels[2]=r'A61'
labels[3]=r'A5A60'
labels[4]=r'A5A61'
labels[5]=r'biotin'

nx.draw_networkx_labels(G,pos,labels,font_size=8)

plt.axis('off')
plt.savefig("network_diagram.pdf") # save as pdf
