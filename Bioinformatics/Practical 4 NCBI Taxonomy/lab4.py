import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

# Test the code on a smaller tree

# Read the data
file_path = "test.dmp"
df = pd.read_csv(file_path, delimiter = r"\s+\|\s+", usecols=[0, 1, 2], engine = "python", header=None, names = ["child", "parent", "rank"]) 
print('Test small tree: \n', df.head())

# Add ranks
attrs = dict(zip(df['child'], df['rank']))

# Create the graph
G = nx.from_pandas_edgelist(df, source='parent', target='child', create_using=nx.DiGraph())
nx.set_node_attributes(G, attrs, name="rank")

# Check self-loops and remove these edges
print('Nodes with self-loops: ', list(nx.nodes_with_selfloops(G)))
G.remove_edge(1,1)

# Type of the graph
print('Is it a tree? ', nx.is_tree(G))
print('Is it directed? ', nx.is_directed(G))
print('Is it rooted? ', nx.is_arborescence(G))
print('Is it a DAG? ', nx.is_directed_acyclic_graph(G))

# The number of nodes
print('The number of nodes: ', G.number_of_nodes())

# Draw the tree
pos = nx.spring_layout(G)
nx.draw_networkx(G, pos, with_labels=True, node_size=500, node_color='lightblue', font_size=10, font_weight='bold', arrows=True)
plt.title("Rooted Tree with Ranks")
plt.show()


# Nodes to remove
not_nodes = [node for node, attrs in G.nodes(data=True) if attrs['rank'] != 'rank']

# Nodes to remove
restricted_nodes = [node for node, attrs in G.nodes(data=True) if attrs['rank'] == 'rank']

restricted_taxonomy = G.copy()

#Lists of successors and predecessors for the nodes wich will be removed
list_succes = [list(G.successors(i)) for i in not_nodes]
list_pred = [list(G.predecessors(i)) for i in not_nodes]

# Add edges from the destination node to the children
restricted_taxonomy.add_edges_from([(list_pred[i][0], list_succes[i][j]) for i in range(1,len(list_succes)) for j in range(len(list_succes[i]))])

# Restrict
restricted_taxonomy = nx.subgraph(restricted_taxonomy, restricted_nodes)


# Draw the restricted graph
pos = nx.spring_layout(restricted_taxonomy)
nx.draw_networkx(restricted_taxonomy, pos, with_labels=True, node_size=500, node_color='lightblue', font_size=10, font_weight='bold', arrows=True)
plt.title("Restricted graph")
plt.show()


#--------------------------------------------------------------------------------------

# Work with the whole dataset

# Read the data
file_path = "nodes.dmp"
df = pd.read_csv(file_path, delimiter = r"\s+\|\s+", usecols=[0, 1, 2], engine = "python", header=None, names = ["child", "parent", "rank"]) 
print('Test small tree: \n', df.head())

# Add ranks
attrs = dict(zip(df['child'], df['rank']))

# Create the graph
G = nx.from_pandas_edgelist(df, source='parent', target='child', create_using=nx.DiGraph())
nx.set_node_attributes(G, attrs, name="rank")

"""
list(G.successors(10239))
G.nodes[2840056]["rank"]
"""

# Check self-loops and remove these edges
print('Nodes with self-loops: ', list(nx.nodes_with_selfloops(G)))
G.remove_edge(1,1)

# Type of the graph
print('Is it a tree? ', nx.is_tree(G))
print('Is it directed? ', nx.is_directed(G))
print('Is it rooted? ', nx.is_arborescence(G)) # to check if it is rooted
print('Is it a DAG? ', nx.is_directed_acyclic_graph(G))

# The number of nodes
print('The number of nodes: ', G.number_of_nodes())

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

print('Is it a tree? ', nx.is_tree(restricted_taxonomy))

# The number of nodes
print('The number of nodes in the restricted version: ', restricted_taxonomy.number_of_nodes())
print('The number of edges in the restricted version: ', restricted_taxonomy.number_of_edges())


# Search for the top taxonomic rank in the NCBI taxonomy
print('Is it a DAG? ', nx.is_directed_acyclic_graph(restricted_taxonomy))

long_path = nx.dag_longest_path(restricted_taxonomy)
print('The rank of the top nodes is: ', restricted_taxonomy.nodes[long_path[1]]["rank"])

