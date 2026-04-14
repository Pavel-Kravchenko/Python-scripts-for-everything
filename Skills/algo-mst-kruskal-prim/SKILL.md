---
name: algo-mst-kruskal-prim
description: "Minimum spanning trees — Kruskal's (Union-Find, O(E log E)) and Prim's (heap, O((V+E) log V)), MST clustering"
tool_type: python
primary_tool: Python
---

# Minimum Spanning Trees

An MST connects all V vertices with V-1 edges at minimum total weight. Unique when all edge weights are distinct.

**Cut property:** The minimum-weight edge crossing any partition of vertices must be in the MST. This is why greedy works.

## Union-Find (Disjoint Set)

```python
class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # path compression
        return self.parent[x]

    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px == py: return False
        if self.rank[px] < self.rank[py]: px, py = py, px
        self.parent[py] = px
        if self.rank[px] == self.rank[py]: self.rank[px] += 1
        return True
```

## Kruskal's Algorithm — O(E log E)

Sort edges by weight, greedily add if no cycle (Union-Find check).

```python
def kruskal(vertices, edges):
    """edges: list of (u, v, weight). Returns (mst_edges, total_weight)."""
    vertex_idx = {v: i for i, v in enumerate(vertices)}
    uf = UnionFind(len(vertices))
    mst, total = [], 0
    for u, v, w in sorted(edges, key=lambda e: e[2]):
        if uf.union(vertex_idx[u], vertex_idx[v]):
            mst.append((u, v, w))
            total += w
            if len(mst) == len(vertices) - 1: break
    return mst, total
```

## Prim's Algorithm — O((V + E) log V)

Grow MST from a start vertex using a min-heap of candidate edges.

```python
import heapq
from collections import defaultdict

def prim(vertices, adj):
    """adj: {vertex: [(neighbor, weight), ...]}. Returns (mst_edges, total_weight)."""
    start = vertices[0]
    visited = {start}
    mst, total = [], 0
    heap = [(w, start, v) for v, w in adj[start]]
    heapq.heapify(heap)
    while heap and len(visited) < len(vertices):
        w, u, v = heapq.heappop(heap)
        if v in visited: continue
        visited.add(v)
        mst.append((u, v, w))
        total += w
        for neighbor, weight in adj[v]:
            if neighbor not in visited:
                heapq.heappush(heap, (weight, v, neighbor))
    return mst, total
```

## MST-Based Clustering

Cut the k-1 longest MST edges to get k clusters. Equivalent to single-linkage hierarchical clustering.

```python
def mst_clustering(vertices, edges, k=2):
    mst, _ = kruskal(vertices, edges)
    mst_sorted = sorted(mst, key=lambda e: e[2], reverse=True)
    kept = mst_sorted[k-1:]
    uf = UnionFind(len(vertices))
    vertex_idx = {v: i for i, v in enumerate(vertices)}
    for u, v, w in kept:
        uf.union(vertex_idx[u], vertex_idx[v])
    clusters = defaultdict(list)
    for v in vertices:
        clusters[uf.find(vertex_idx[v])].append(v)
    return list(clusters.values())
```

## Kruskal vs Prim

| | Kruskal | Prim |
|---|---|---|
| Time | O(E log E) | O((V+E) log V) |
| Best for | Sparse graphs | Dense graphs |
| Data structure | Union-Find | Min-heap |

## Pitfalls

- Kruskal requires Union-Find; a naive cycle check makes it O(VE)
- Prim's heap may contain stale entries for already-visited vertices — always check `if v in visited: continue`
- Equal-weight edges can produce different (but equally valid) MSTs
