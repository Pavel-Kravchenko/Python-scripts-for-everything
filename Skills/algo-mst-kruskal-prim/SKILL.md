---
name: algo-mst-kruskal-prim
description: "1. Understand MST and its applications in biological networks 2. Implement Kruskal's algorithm with Union-Find 3. Implement Prim's algorithm with priority queue 4. Apply MST to phylogenetic tree const"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/04_mst_kruskal_prim.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# 🌳 Minimum Spanning Trees: Kruskal's and Prim's Algorithms

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/04_mst_kruskal_prim.ipynb`*


## Learning Objectives

1. Understand MST and its applications in biological networks
2. Implement Kruskal's algorithm with Union-Find
3. Implement Prim's algorithm with priority queue
4. Apply MST to phylogenetic tree construction

---

## What is a Minimum Spanning Tree?

A **spanning tree** connects all vertices with V-1 edges and no cycles. A **minimum spanning tree** has minimum total edge weight.

**Biological applications:**
- Phylogenetic tree construction (minimum evolution)
- Gene co-expression network simplification
- Clustering (single-linkage hierarchical)

## MST Properties

An MST has exactly **V-1 edges** for a graph with V vertices. If all edge weights are distinct the MST is **unique**. When weights are equal, multiple MSTs may exist.

The **cut property** explains why greedy algorithms work: for any partition of the vertices into two sets (a "cut"), the minimum-weight edge crossing that cut must belong to the MST. This is the theoretical guarantee behind both Kruskal's and Prim's approaches.

```python
import heapq

class UnionFind:
    """Disjoint Set Union with path compression and union by rank."""
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n
    
    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # Path compression
        return self.parent[x]
    
    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px == py:
            return False
        if self.rank[px] < self.rank[py]:
            px, py = py, px
        self.parent[py] = px
        if self.rank[px] == self.rank[py]:
            self.rank[px] += 1
        return True
```python

```python
def kruskal(vertices, edges):
    """
    Kruskal's MST algorithm.
    edges: list of (u, v, weight)
    Time: O(E log E)
    """
    vertex_idx = {v: i for i, v in enumerate(vertices)}
    uf = UnionFind(len(vertices))
    
    # Sort edges by weight
    sorted_edges = sorted(edges, key=lambda e: e[2])
    
    mst = []
    total_weight = 0
    
    for u, v, w in sorted_edges:
        if uf.union(vertex_idx[u], vertex_idx[v]):
            mst.append((u, v, w))
            total_weight += w
            if len(mst) == len(vertices) - 1:
                break
    
    return mst, total_weight

# Example: Distance matrix between species
species = ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Mouse']
distances = [
    ('Human', 'Chimp', 1.2),
    ('Human', 'Gorilla', 1.6),
    ('Human', 'Orangutan', 3.1),
    ('Human', 'Mouse', 8.5),
    ('Chimp', 'Gorilla', 1.5),
    ('Chimp', 'Orangutan', 3.0),
    ('Chimp', 'Mouse', 8.4),
    ('Gorilla', 'Orangutan', 2.8),
    ('Gorilla', 'Mouse', 8.3),
    ('Orangutan', 'Mouse', 8.1),
]

mst, total = kruskal(species, distances)
print("Minimum Spanning Tree (simplified phylogeny):")
for u, v, w in mst:
    print(f"  {u} -- {v}: {w}")
print(f"Total evolutionary distance: {total}")
```python

## Step-by-Step Kruskal

The verbose version below prints each edge decision — accepted when the two endpoints are in different components, rejected when they are already connected (adding it would create a cycle).

```python
def kruskal_verbose(vertices, edges):
    """Kruskal's with step-by-step output."""
    vertex_idx = {v: i for i, v in enumerate(vertices)}
    uf = UnionFind(len(vertices))
    sorted_edges = sorted(edges, key=lambda e: e[2])
    mst = []

    for u, v, w in sorted_edges:
        if uf.union(vertex_idx[u], vertex_idx[v]):
            mst.append((u, v, w))
            print(f"  ACCEPT {u}--{v} (weight {w})")
        else:
            print(f"  reject {u}--{v} (weight {w}) -- would create cycle")
        if len(mst) == len(vertices) - 1:
            break
    return mst

print("Kruskal step-by-step on species distances:")
kruskal_verbose(species, distances)
```python

```python
def prim(vertices, adj):
    """
    Prim's MST algorithm.
    adj: dict {vertex: [(neighbor, weight), ...]}
    Time: O((V + E) log V)
    """
    start = vertices[0]
    visited = {start}
    mst = []
    total_weight = 0
    
    # Min-heap: (weight, u, v)
    heap = [(w, start, v) for v, w in adj[start]]
    heapq.heapify(heap)
    
    while heap and len(visited) < len(vertices):
        w, u, v = heapq.heappop(heap)
        if v in visited:
            continue
        
        visited.add(v)
        mst.append((u, v, w))
        total_weight += w
        
        for neighbor, weight in adj[v]:
            if neighbor not in visited:
                heapq.heappush(heap, (weight, v, neighbor))
    
    return mst, total_weight

# Build adjacency list
from collections import defaultdict
adj = defaultdict(list)
for u, v, w in distances:
    adj[u].append((v, w))
    adj[v].append((u, w))

mst_prim, total_prim = prim(species, adj)
print("\nPrim's MST (same result):")
for u, v, w in mst_prim:
    print(f"  {u} -- {v}: {w}")
```python

## Connection to Phylogenetics

Minimum spanning trees are related to minimum evolution phylogenetic methods. While real phylogenetics uses more sophisticated tree-building (neighbor-joining, maximum likelihood), MST gives a rough approximation: connect species by shortest evolutionary distances without cycles.

The MST backbone reveals which species are closest neighbors and which groups cluster together, matching well-known primate evolutionary relationships.

## Connection to Single-Linkage Clustering

Single-linkage hierarchical clustering produces the same dendrogram as the MST. The merge order in the dendrogram corresponds exactly to the order in which Kruskal's adds edges.

Cutting the **k-1 longest edges** from the MST splits the data into **k clusters** — a fast alternative to running full hierarchical clustering.

```python
def mst_clustering(vertices, edges, k=2):
    """Cluster by removing k-1 longest edges from MST."""
    mst, _ = kruskal(vertices, edges)
    # Sort MST edges by weight descending, remove k-1 longest
    mst_sorted = sorted(mst, key=lambda e: e[2], reverse=True)
    removed = mst_sorted[:k-1]
    kept = mst_sorted[k-1:]

    print(f"Removed edges for {k} clusters:")
    for u, v, w in removed:
        print(f"  {u}--{v} (weight {w})")

    # Find components using remaining edges
    uf = UnionFind(len(vertices))
    vertex_idx = {v: i for i, v in enumerate(vertices)}
    for u, v, w in kept:
        uf.union(vertex_idx[u], vertex_idx[v])

    from collections import defaultdict
    clusters = defaultdict(list)
    for v in vertices:
        clusters[uf.find(vertex_idx[v])].append(v)
    return list(clusters.values())

clusters = mst_clustering(species, distances, k=2)
for i, c in enumerate(clusters, 1):
    print(f"Cluster {i}: {c}")
```python

## Exercises

### Exercise 1 (1 star): Gene Co-expression MST

Given a gene co-expression network with correlation weights, find the MST backbone. Use at least 8 genes. Print the resulting backbone connections and total weight.

```python
# Gene co-expression network: (gene_a, gene_b, correlation_distance)
# Use 1 - |correlation| as the edge weight so highly correlated genes are close
genes = ['BRCA1', 'BRCA2', 'TP53', 'ATM', 'CHEK2', 'RAD51', 'PALB2', 'MDM2']
coexpression = [
    # YOUR CODE HERE: add edges as (gene_a, gene_b, distance)
]

# YOUR CODE HERE: call kruskal and print the MST backbone
```python

### Exercise 2 (2 stars): MST Edge Replacement

Implement a function that, given an MST and a new candidate edge `(u, v, w)`, determines whether adding this edge would improve the tree. When a new edge is added to an MST it creates exactly one cycle. If the new edge is lighter than the heaviest edge in that cycle, replace the heavy edge with the new one.

```python
def mst_update(vertices, mst_edges, new_edge):
    """
    Determine if new_edge improves the MST.
    Returns the updated MST if improvement is possible, otherwise returns the original.

    Parameters
    ----------
    vertices   : list of vertex labels
    mst_edges  : list of (u, v, weight) forming the current MST
    new_edge   : (u, v, weight) candidate edge

    Returns
    -------
    (mst, improved) where improved is True if the tree was updated
    """
    # YOUR CODE HERE
    # Hint: build an adjacency structure from mst_edges, find the path
    # between new_edge endpoints using BFS/DFS, then find the max-weight
    # edge on that path. If new_edge weight < max path weight, swap them.
    pass
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
