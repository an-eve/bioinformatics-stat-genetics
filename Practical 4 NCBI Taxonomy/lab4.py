import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

file_path = "nodes.dmp"

df = pd.read_csv(file_path, delimiter = r"\s+\|\s+", usecols=[0, 1, 2], engine = "python", header=None, names = ["child", "parent", "rank"]) 
df.head()

attrs = dict(zip(df['child'], df['rank']))


# Find rows where values in the first and second columns are the same
duplicate_rows = df[df['parent'] == df['child']]

# Delete the duplicate rows by index
df = df.drop(duplicate_rows.index)

# Reset the index if needed
df = df.reset_index(drop=True)

G = nx.from_pandas_edgelist(df, source='parent', target='child', create_using=nx.DiGraph())

nx.set_node_attributes(G, attrs, name="rank")

list(nx.nodes_with_selfloops(G))
G.remove_edge(1,1)

nx.is_tree(G)
nx.is_directed(G)
nx.is_arborescence(G)
nx.is_directed_acyclic_graph(G)

G.number_of_nodes()

valid_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
restricted_nodes = [node for node in G.nodes() if not G.nodes[node]['rank'] in valid_ranks]

G.edges()
G.nodes(data=True)

G.nodes[1]["rank"] = "rank"

not_tree_nodes = [node for node, attrs in G.nodes(data=True) if attrs['rank'] != 'rank']
restricted_tree_nodes = [node for node, attrs in G.nodes(data=True) if attrs['rank'] == 'rank']
restricted_tree = nx.subgraph(restricted_tree, restricted_tree_nodes)

not_tree_nodes = [node for node, attrs in G.nodes(data=True) if not attrs['rank'] in valid_ranks]
restricted_tree_nodes = [node for node, attrs in G.nodes(data=True) if attrs['rank'] in valid_ranks]
restricted_tree = nx.subgraph(restricted_tree, restricted_tree_nodes)

restricted_tree.nodes(data=True)

[attrs for node, attrs in G.nodes(data=True)]

"""
# Draw the directed graph
pos = nx.spring_layout(G)
nx.draw_networkx(G, pos, with_labels=True, node_size=500, node_color='lightblue', font_size=10, font_weight='bold', arrows=True)
plt.title("Rooted Tree with Ranks")
plt.show()
"""

"""
# Draw the directed graph
pos = nx.spring_layout(restricted_tree)
nx.draw_networkx(restricted_tree, pos, with_labels=True, node_size=500, node_color='lightblue', font_size=10, font_weight='bold', arrows=True)
plt.title("Rooted Tree with Ranks")
plt.show()
"""
ls = [list(G.successors(i)) for i in not_tree_nodes]
lp = [list(G.predecessors(i)) for i in not_tree_nodes]
list(restricted_tree.predecessors(2))

restricted_tree = G.copy()
# Add edges from the destination node to the children
restricted_tree.add_edges_from([(lp[i][0], ls[i][j]) for i in range(1,len(ls)) for j in range(len(ls[i]))])

# Remove edges from the source node to the children
restricted_tree.remove_edges_from([(source_node, child) for child in children])

restricted_tree.number_of_nodes()

nx.is_tree(restricted_tree)
nx.is_arborescence(restricted_tree)

[node for node in restricted_tree.nodes() if restricted_tree.in_degree(node) == 0]

restricted_tree.nodes[1090]["rank"]

ls