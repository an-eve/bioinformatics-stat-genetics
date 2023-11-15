**Assignments:**

1. Given a file **nodes.dmp** for the NCBI taxonomy and a file **mapping.txt** of mappings of sequence reads
to NCBI taxonomic identifiers, write a Python script to find the lineage of the mapped nodes for each sequence
read. Give the code of your Python script as your answer to this question, using the $\LaTeX$ package **listings**.
2. What is the taxonomic rank of the LCA of the mapped nodes for each of the sequence reads?
3. Do the sequence reads come from archaea, bacteria, eukaryota, or viruses?
4. Given a file **nodes.dmp** for the NCBI taxonomy and a file **mapping.txt** of mappings of sequence
reads to NCBI taxonomic identifiers, write a Python script to build the LCA skeleton tree for each sequence read.
Give the code of your Python script as your answer to this question, using the $\LaTeX$ package **listings**.
5. How many nodes are there in the LCA skeleton tree for each of the sequence reads?
6. What is the taxonomic rank of the root of the LCA skeleton tree for each of the sequence reads?

**Definitions:**

* **The lineage** of a node refers to the sequence of nodes that are ancestors of the given node, starting from the root of the tree and going up to the immediate parent of the node itself. These ancestors are the nodes that form the path leading from the root to the target node.

* **The Lowest Common Ancestor (LCA)** of multiple nodes in a tree is the deepest node that is a common ancestor to all of those nodes.

* **The LCA skeleton tree** is a simplified tree structure that highlights the lowest common ancestor of a specified group of taxa or sequences, providing a more focused view of their evolutionary relationships while reducing complexity compared to a full phylogenetic tree.












