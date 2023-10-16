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
print('Is it rooted? ', nx.is_arborescence(G))
print('Is it a DAG? ', nx.is_directed_acyclic_graph(G))

# The number of nodes
print('The number of nodes: ', G.number_of_nodes())

valid_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# Nodes to remove
not_nodes = [node for node, attrs in G.nodes(data=True) if not attrs['rank'] in valid_ranks]

# Nodes to remove
restricted_nodes = [node for node, attrs in G.nodes(data=True) if attrs['rank'] in valid_ranks]

restricted_taxonomy = G.copy()

# Remove the info about the node 1 from these lists because it's a root
not_nodes[0]
not_nodes = not_nodes[1:]

#Lists of successors and predecessors for the nodes wich will be removed
list_succes = [list(G.successors(i)) for i in not_nodes]
list_pred = [list(G.predecessors(i)) for i in not_nodes]

# Add edges from the destination node to the children
restricted_taxonomy.add_edges_from([(list_pred[i][0], list_succes[i][j]) for i in range(1,len(list_succes)) for j in range(len(list_succes[i]))])

# Restrict
restricted_taxonomy = nx.subgraph(restricted_taxonomy, restricted_nodes)

# The number of nodes
print('The number of nodes in the restricted version: ', restricted_taxonomy.number_of_nodes())

# The nodes with in-degree 0 and their ranks
search_top = [(node, restricted_taxonomy.nodes[node]["rank"])  for node in restricted_taxonomy.nodes() if restricted_taxonomy.in_degree(node) == 0]
search_top_list = [x[1] for x in search_top]
unique_elements = [x for i, x in enumerate(search_top_list) if x not in search_top_list[:i]]