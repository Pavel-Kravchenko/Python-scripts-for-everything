---
name: algo-graph-representations
description: "1. Understand what graphs are and why they matter in bioinformatics 2. Implement three core graph representations: adjacency matrix, adjacency list, edge list 3. Analyze space/time trade-offs for each"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/01_graph_representations.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# 🔗 Graph Representations

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/01_graph_representations.ipynb`*


## Learning Objectives

By the end of this notebook, you will be able to:

1. Understand what graphs are and why they matter in bioinformatics
2. Implement three core graph representations: adjacency matrix, adjacency list, edge list
3. Analyze space/time trade-offs for each representation
4. Choose the right representation for biological network analysis

---

## Why Graphs in Bioinformatics?

Graphs are **everywhere** in biology:

| Application | Vertices | Edges |
|------------|----------|-------|
| Protein-protein interaction | Proteins | Physical interactions |
| Metabolic networks | Metabolites | Reactions |
| Gene regulatory networks | Genes/TFs | Regulatory relationships |
| Sequence assembly (de Bruijn) | k-mers | Overlaps |
| Phylogenetic trees | Species/sequences | Evolutionary relationships |

Understanding graph representations is the foundation for analyzing these networks.

## 1. Graph Fundamentals

A **graph** G = (V, E) consists of:
- **V**: Set of vertices (nodes)
- **E**: Set of edges connecting vertices

### Types of Graphs

```python
Undirected          Directed           Weighted
A --- B             A --> B            A --5-- B
|     |             |     |            |       |
C --- D             C <-- D            C --3-- D
```python

```python
# Let's define a sample biological network:
# A protein-protein interaction network

proteins = ['TP53', 'MDM2', 'BRCA1', 'ATM', 'CHEK2']
interactions = [
    ('TP53', 'MDM2'),    # p53-MDM2 interaction
    ('TP53', 'BRCA1'),   # p53-BRCA1 interaction
    ('ATM', 'TP53'),     # ATM phosphorylates p53
    ('ATM', 'CHEK2'),    # ATM activates CHEK2
    ('CHEK2', 'TP53'),   # CHEK2 phosphorylates p53
    ('BRCA1', 'ATM'),    # BRCA1-ATM complex
]

print(f"Network has {len(proteins)} proteins and {len(interactions)} interactions")
```python

## 2. Adjacency Matrix

A 2D matrix where `matrix[i][j] = 1` if there's an edge from vertex i to vertex j.

**Properties:**
- Space: O(V²)
- Edge lookup: O(1)
- Finding neighbors: O(V)
- Best for: Dense graphs, when V is small

```python
import numpy as np

class AdjacencyMatrix:
    """Graph representation using adjacency matrix."""
    
    def __init__(self, vertices):
        self.vertices = vertices
        self.vertex_to_idx = {v: i for i, v in enumerate(vertices)}
        self.n = len(vertices)
        self.matrix = np.zeros((self.n, self.n), dtype=int)
    
    def add_edge(self, u, v, directed=False):
        """Add edge between vertices u and v."""
        i, j = self.vertex_to_idx[u], self.vertex_to_idx[v]
        self.matrix[i][j] = 1
        if not directed:
            self.matrix[j][i] = 1
    
    def has_edge(self, u, v):
        """Check if edge exists - O(1)."""
        i, j = self.vertex_to_idx[u], self.vertex_to_idx[v]
        return self.matrix[i][j] == 1
    
    def neighbors(self, u):
        """Get all neighbors of vertex u - O(V)."""
        i = self.vertex_to_idx[u]
        return [self.vertices[j] for j in range(self.n) if self.matrix[i][j] == 1]
    
    def degree(self, u):
        """Get degree of vertex u."""
        i = self.vertex_to_idx[u]
        return int(np.sum(self.matrix[i]))
    
    def display(self):
        """Pretty print the matrix."""
        print("     " + "  ".join(f"{v[:4]:>4}" for v in self.vertices))
        for i, v in enumerate(self.vertices):
            row = "  ".join(f"{self.matrix[i][j]:>4}" for j in range(self.n))
            print(f"{v[:4]:>4} {row}")
```python

```python
# Build our PPI network as adjacency matrix
ppi_matrix = AdjacencyMatrix(proteins)

for p1, p2 in interactions:
    ppi_matrix.add_edge(p1, p2)  # undirected by default

ppi_matrix.display()
```python

```python
# Query the network
print(f"Does TP53 interact with MDM2? {ppi_matrix.has_edge('TP53', 'MDM2')}")
print(f"Does MDM2 interact with ATM? {ppi_matrix.has_edge('MDM2', 'ATM')}")
print(f"\nTP53 interacts with: {ppi_matrix.neighbors('TP53')}")
print(f"TP53 degree (# interactions): {ppi_matrix.degree('TP53')}")
```python

## 3. Adjacency List

Each vertex stores a list of its neighbors.

**Properties:**
- Space: O(V + E)
- Edge lookup: O(degree) — can be O(1) with sets
- Finding neighbors: O(1) to access, O(degree) to iterate
- Best for: Sparse graphs, most biological networks

```python
from collections import defaultdict

class AdjacencyList:
    """Graph representation using adjacency list."""
    
    def __init__(self):
        self.graph = defaultdict(set)  # Using sets for O(1) edge lookup
    
    def add_vertex(self, v):
        """Add isolated vertex."""
        if v not in self.graph:
            self.graph[v] = set()
    
    def add_edge(self, u, v, directed=False):
        """Add edge between u and v."""
        self.graph[u].add(v)
        if not directed:
            self.graph[v].add(u)
    
    def has_edge(self, u, v):
        """Check if edge exists - O(1) with sets."""
        return v in self.graph.get(u, set())
    
    def neighbors(self, u):
        """Get neighbors of u - O(1) access."""
        return self.graph.get(u, set())
    
    def degree(self, u):
        """Get degree of vertex."""
        return len(self.graph.get(u, []))
    
    def vertices(self):
        """Get all vertices."""
        return list(self.graph.keys())
    
    def edges(self):
        """Get all edges (for undirected, each edge appears once)."""
        seen = set()
        for u in self.graph:
            for v in self.graph[u]:
                if (v, u) not in seen:
                    seen.add((u, v))
        return list(seen)
    
    def display(self):
        """Pretty print the adjacency list."""
        for v in sorted(self.graph.keys()):
            neighbors = sorted(self.graph[v])
            print(f"{v}: {neighbors}")
```python

```python
# Build PPI network as adjacency list
ppi_list = AdjacencyList()

for protein in proteins:
    ppi_list.add_vertex(protein)

for p1, p2 in interactions:
    ppi_list.add_edge(p1, p2)

ppi_list.display()
```python

```python
# Same queries, more efficient for sparse graphs
print(f"TP53 neighbors: {ppi_list.neighbors('TP53')}")
print(f"All edges: {ppi_list.edges()}")
```python

## 4. Edge List

Simply a list of all edges.

**Properties:**
- Space: O(E)
- Edge lookup: O(E)
- Finding neighbors: O(E)
- Best for: Edge-centric algorithms (Kruskal's MST), simple storage

```python
class EdgeList:
    """Graph representation using edge list."""
    
    def __init__(self):
        self.edges = []
        self._vertices = set()
    
    def add_edge(self, u, v, weight=1):
        """Add edge (u, v) with optional weight."""
        self.edges.append((u, v, weight))
        self._vertices.add(u)
        self._vertices.add(v)
    
    def has_edge(self, u, v):
        """Check if edge exists - O(E)."""
        return any((e[0] == u and e[1] == v) or (e[0] == v and e[1] == u) 
                   for e in self.edges)
    
    def neighbors(self, u):
        """Get neighbors - O(E)."""
        result = set()
        for e in self.edges:
            if e[0] == u:
                result.add(e[1])
            elif e[1] == u:
                result.add(e[0])
        return result
    
    def vertices(self):
        return self._vertices
    
    def sorted_edges(self):
        """Return edges sorted by weight (for Kruskal's)."""
        return sorted(self.edges, key=lambda e: e[2])
    
    def display(self):
        for u, v, w in self.edges:
            print(f"{u} --({w})-- {v}")
```python

```python
# Example: Weighted PPI network (weights = confidence scores)
ppi_edges = EdgeList()

weighted_interactions = [
    ('TP53', 'MDM2', 0.99),
    ('TP53', 'BRCA1', 0.85),
    ('ATM', 'TP53', 0.92),
    ('ATM', 'CHEK2', 0.88),
    ('CHEK2', 'TP53', 0.76),
    ('BRCA1', 'ATM', 0.81),
]

for p1, p2, conf in weighted_interactions:
    ppi_edges.add_edge(p1, p2, conf)

ppi_edges.display()
```python

```python
# Edges sorted by confidence (useful for filtering)
print("\nHighest confidence interactions:")
for u, v, w in ppi_edges.sorted_edges()[::-1]:
    print(f"  {u} <-> {v}: {w:.2f}")
```python

## 5. Comparison Summary

| Operation | Adjacency Matrix | Adjacency List | Edge List |
|-----------|-----------------|----------------|------------|
| Space | O(V²) | O(V + E) | O(E) |
| Add edge | O(1) | O(1) | O(1) |
| Remove edge | O(1) | O(degree) | O(E) |
| Check edge | O(1) | O(degree) | O(E) |
| Get neighbors | O(V) | O(degree) | O(E) |
| Iterate all edges | O(V²) | O(E) | O(E) |

### When to use each:

- **Adjacency Matrix**: Dense graphs, frequent edge queries, matrix operations
- **Adjacency List**: Most biological networks (sparse), BFS/DFS traversals
- **Edge List**: MST algorithms, simple file I/O, edge-centric operations

## Loading Real Graph Data

The examples above use hard-coded protein networks. Real-world graphs are typically stored as **edge lists** — one edge per line with optional weights. Here is how to parse the format used in `Course/Assets/data/graph1.txt`:

```python
G E 7
G F 8
E D 5
...
```python

Each line is `node_a  node_b  weight`. Loading this into an adjacency list gives you a real weighted graph to experiment with.

```python
from pathlib import Path
from collections import defaultdict

def load_edge_list(filepath: str) -> dict[str, list[tuple[str, float]]]:
    """Parse a weighted edge-list file into an undirected adjacency list.

    File format: one edge per line as  'node_a  node_b  weight'
    Returns: {node: [(neighbor, weight), ...]}
    """
    graph: dict[str, list[tuple[str, float]]] = defaultdict(list)
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            u, v, w = parts[0], parts[1], float(parts[2])
            graph[u].append((v, w))
            graph[v].append((u, w))   # undirected
    return dict(graph)


# Path to the sample graph file (relative to this notebook)
DATA_FILE = Path("../../assets/data/graph1.txt")

g = load_edge_list(str(DATA_FILE))

print("Nodes:", sorted(g.keys()))
print()
for node in sorted(g.keys()):
    neighbors = ", ".join(f"{v}({w:.0f})" for v, w in sorted(g[node]))
    print(f"  {node}: {neighbors}")
```python

## 🧬 Exercise 1: Gene Regulatory Network

Build a **directed** gene regulatory network where edges represent "gene A regulates gene B".

```python
Regulatory relationships:
- MYC → CCND1 (activates)
- MYC → CDKN1A (represses)
- TP53 → CDKN1A (activates)
- TP53 → BAX (activates)
- E2F1 → CCND1 (activates)
- RB1 → E2F1 (represses)
```python

1. Implement using an adjacency list
2. Find all genes regulated by TP53
3. Find all regulators of CDKN1A

```python
# Your solution here
regulations = [
    ('MYC', 'CCND1'),
    ('MYC', 'CDKN1A'),
    ('TP53', 'CDKN1A'),
    ('TP53', 'BAX'),
    ('E2F1', 'CCND1'),
    ('RB1', 'E2F1'),
]

# Implement directed graph using AdjacencyList
# Hint: Use directed=True when adding edges
```python

## 🧬 Exercise 2: Memory Analysis

A typical human PPI network has ~20,000 proteins and ~300,000 interactions.

Calculate:
1. Memory usage for adjacency matrix (assuming 1 byte per cell)
2. Memory usage for adjacency list (assuming 8 bytes per edge entry)
3. Which representation is more efficient?

```python
V = 20000   # proteins
E = 300000  # interactions

# Memory usage analysis
# Adjacency matrix: V x V entries x 8 bytes (int64 per cell)
matrix_bytes = V * V * 8

# Adjacency list: roughly 2*E pointers (undirected) x 8 bytes each
# plus V list-header objects (~56 bytes each in Python)
list_bytes = 2 * E * 8 + V * 56

print("Adjacency Matrix:", round(matrix_bytes / 1e6, 1), "MB")
print("Adjacency List:  ", round(list_bytes / 1e6, 1), "MB")
print("Ratio:", round(matrix_bytes / list_bytes, 1), "x")
print()
print("For sparse graphs (E << V*V), adjacency lists use far less memory.")
print("Adjacency matrices are only practical when E is close to V*V.")
```python

---

## Summary

In this notebook, you learned:

✅ Three fundamental graph representations  
✅ Space and time complexity trade-offs  
✅ How to choose the right representation for biological networks  
✅ Implementation patterns in Python

**Next:** [02_bfs_dfs.ipynb](02_bfs_dfs.ipynb) - Graph traversal algorithms

---

## Key Takeaway

> Most biological networks (PPI, GRN, metabolic) are **sparse** — adjacency lists are almost always the right choice. Reserve adjacency matrices for small, dense networks or when you need fast edge lookups.

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
