import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import pygraphviz as pgv

# Read the data
file_path = "nodes.dmp"
df = pd.read_csv(file_path, delimiter = r"\s+\|\s+", usecols=[0, 1, 2], engine = "python", header=None, names = ["child", "parent", "rank"]) 
print('The head of the dataframe: \n', df.head())

# Add ranks
attrs = dict(zip(df['child'], df['rank']))

# Create the graph
G = nx.from_pandas_edgelist(df, source='parent', target='child', create_using=nx.DiGraph())
nx.set_node_attributes(G, attrs, name="rank")

# Check self-loops and remove these edges
print('Nodes with self-loops: ', list(nx.nodes_with_selfloops(G)))
G.remove_edge(1,1)

# Type of the graph
print('Initial graph:')
print('Is it a tree? ', nx.is_tree(G))
print('Is it rooted? ', nx.is_arborescence(G)) # to check if it is rooted

# The number of nodes
print('The number of nodes: ', G.number_of_nodes())

#----------------------------------------------------------------------------------------------------------------------------------------------

# Create the restricted graph that represents the NCBI taxonomy
valid_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

# Nodes to remove
not_nodes = [node for node, attrs in G.nodes(data=True) if not attrs['rank'] in valid_ranks]

restricted_taxonomy = G.copy()
restricted_taxonomy.add_edge(1,1)

def delete_nodes_preserve_neighbors(graph, nodes_to_delete):
    for node in nodes_to_delete:
        list_succes = list(graph.successors(node)) 
        list_pred = list(graph.predecessors(node))
        if node != 1:
            graph.remove_node(node)
        for source in list_pred:
            for target in list_succes:
                if source != target and not graph.has_edge(source, target):
                    graph.add_edge(source, target)

# We keep the node 1 (root) in order to have a structure of the tree 

delete_nodes_preserve_neighbors(restricted_taxonomy, not_nodes)
restricted_taxonomy.remove_edge(1,1)
print('Restricted taxonomy is created')

print('Is it a tree? ', nx.is_tree(restricted_taxonomy))

# The number of nodes
print('The number of nodes in the restricted version: ', restricted_taxonomy.number_of_nodes())
print('The number of edges in the restricted version: ', restricted_taxonomy.number_of_edges())

#----------------------------------------------------------------------------------------------------------------------------------------------

# Read the mappings of sequence reads to NCBI taxonomic identifiers
mapping_file = 'mapping.txt'
mapping = {}

with open(mapping_file, 'r') as file:
    for line in file:
        parts = line.strip().split()
        read_id = parts[0]
        taxonomic_nodes = [int(node) for node in parts[1:]]  # Convert nodes to integers
        mapping[read_id] = taxonomic_nodes

#To check if the mapping file is correct, and all nodes are the leaves of the restricted tree     
"""
for i in range(10):
    for node_id in list(mapping.values())[i]:
        attributes = restricted_taxonomy.nodes[node_id]
        if attributes['rank'] != "species":
            print("Incorrect mapping.txt file!")
            break
    else:
        continue  # only executed if the inner loop did NOT break
    break  # only executed if the inner loop DID break
"""

# Define a function to find the full lineage of a taxonomic node
def find_full_lineage(taxonomy_graph, node):
    lineage = []
    while node in taxonomy_graph:
        lineage.append(node)
        parents = list(taxonomy_graph.predecessors(node))
        if not parents:
            break
        node = parents[0]
    return lineage

# Full lineage of every node of every readID is stored in the dictionnary lineages
lineages = {}
for read_id, taxonomic_nodes in mapping.items():
    #print(f"Read ID: {read_id}")
    seq = []
    for node in taxonomic_nodes:
        full_lineage = find_full_lineage(restricted_taxonomy, node)
        full_lineage.reverse()
        #print(f"    Node {node}: {full_lineage}")
        seq.append(full_lineage)
    lineages[read_id] = seq


# LCAs for each read
LCAs = {}
for read_id, seq in lineages.items():
    min_length = min(len(lst) for lst in seq) 
    lst1 = seq[0]
    all_same = all(lst1[0] == lst[0] for lst in seq[1:]) #this value is always True (the root)
    i = 0
    while all_same and (i+1) < min_length:
        LCAs[read_id] = lst1[i]
        all_same = all(lst1[i+1] == lst[i+1] for lst in seq[1:])

print("{The LCAs for each read: \n", LCAs)

#Matching with taxonomy ranks (superkingdoms)
#2759 -> eukaryota
#2 -> bacteria
#10239 -> viruses
#2157 -> archaea
sk_matching = {2759: 'eukaryota', 2: 'bacteria', 10239: 'viruses', 2157: 'archaea'}

# Superkingdoms for each read
Superkingdoms = {}
for read_id, seq in lineages.items():
    skingdoms_for_read ={}
    skingdoms_for_read['others'] = 0
    for l in seq:
        if l[1] in sk_matching:
            if sk_matching[l[1]] in skingdoms_for_read:
                skingdoms_for_read[sk_matching[l[1]]] += 1
            else:
                skingdoms_for_read[sk_matching[l[1]]] = 1
        else:
            skingdoms_for_read['others'] += 1
    Superkingdoms[read_id] = skingdoms_for_read
print("The Superkingdoms for each read: \n")
for read_id, sk_number in Superkingdoms.items():
    print(f"The Superkingdoms for {read_id} read: {Superkingdoms[read_id]}")




# LCA skeletons for each sequence read
Skeletons = {}
for read_id, seq in lineages.items():
    # Create a set of all unique nodes
    all_nodes = {item for sublist in seq for item in sublist}

    # Construct a directed graph (DAG)
    G = nx.DiGraph()
    G.add_nodes_from(all_nodes)

    # Connect nodes in the lineages
    for lineage in seq:
        for i in range(len(lineage) - 1):
            G.add_edge(lineage[i], lineage[i+1])
    Skeletons[read_id] = G
    # The number of nodes
    print(f"For the {read_id} sequence read: ")
    print('The number of nodes: ', Skeletons[read_id].number_of_nodes())
    #print('The number of edges: ', Skeletons[read_id].number_of_edges())
    #print('Is it a tree? ', nx.is_tree(Skeletons[read_id]))

#Drawing for one of the reads
"""
G = Skeletons['R00010']
pos = nx.spring_layout(G)

plt.figure(figsize=(10, 8))
nx.draw_networkx(G, pos, with_labels=True, node_size=100, node_color='lightblue', font_size=6, arrows=True)
plt.title("LCA for the R00010 read")
plt.show
#plt.savefig("LCA_R00010.png")
#plt.close()
"""

#Drawing all of them
from networkx.drawing.nx_pydot import graphviz_layout
import matplotlib.pyplot as plt
G = nx.DiGraph()
for read_id, graph in Skeletons.items():
    pos = graphviz_layout(graph, prog="dot")
    plt.figure(figsize=(10, 8))
    nx.draw(graph, pos)
    plt.title(read_id)
    plt.savefig(f"LCA_{read_id}.png")
    plt.close()

