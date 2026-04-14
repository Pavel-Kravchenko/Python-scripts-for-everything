---
name: graphs-dynamic-programming
description: Graph algorithms (BFS, DFS, Dijkstra, MST, topological sort) and dynamic programming (memoization, tabulation, knapsack, sequence alignment)
primary_tool: NumPy
---

# Graph Algorithms & Dynamic Programming

## When to Use

**Graphs:** BFS/DFS for traversal/reachability. Dijkstra for weighted shortest path. Kruskal/Prim for MST. Topological sort for DAG ordering. Adjacency list for sparse bio networks.

**DP:** Overlapping subproblems + optimal substructure. NW/SW for sequence alignment. Knapsack for resource selection under constraint.

## Complexity Reference

| Algorithm | Time | Space | Notes |
|---|---|---|---|
| BFS / DFS | O(V+E) | O(V) | |
| Dijkstra | O((V+E) log V) | O(V) | Non-negative weights only |
| Kruskal | O(E log E) | O(V) | Better for sparse |
| Prim | O((V+E) log V) | O(V) | Better for dense |
| Topo sort (Kahn) | O(V+E) | O(V) | Detects cycles |
| Knapsack 0/1 | O(nW) | O(W) | Rolling array |
| NW / Edit distance | O(mn) | O(mn) | O(min(m,n)) with rolling row |

**Graph representations:** Adjacency matrix O(V^2) space for dense + frequent lookups. Adjacency list O(V+E) for sparse (bio networks: PPI ~20k nodes, ~300k edges -- list saves ~170x memory). Edge list for MST/I/O.

## DP Framework

1. **State:** What does `dp[i][j]` represent?
2. **Base case:** Smallest subproblem with known answer
3. **Recurrence:** How does `dp[i]` depend on smaller states?
4. **Answer extraction:** Corner, max cell, or traceback?

**Memoization vs tabulation:** Memoization (top-down, `@lru_cache`) computes only needed subproblems but risks stack overflow. Tabulation (bottom-up) enables rolling-array space optimization.

## Sequence Alignment Notes

- **NW (global):** Init borders with gap penalties; traceback from bottom-right
- **SW (local):** Floor cells at 0; traceback from max cell (stops at 0)
- **Affine gaps:** Three matrices M/X/Y; open cost once, extend per extra base
- Log-space for max-product paths: `min(-log(conf))` via Dijkstra

## Code Templates

### Adjacency List + BFS + DFS
```python
from collections import defaultdict, deque

graph = defaultdict(list)  # {node: [(neighbor, weight), ...]}
graph[u].append((v, w)); graph[v].append((u, w))

def bfs(graph, start):
    dist = {start: 0}
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in graph[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist

def dfs(graph, start):
    visited, stack, order = set(), [start], []
    while stack:
        u = stack.pop()
        if u not in visited:
            visited.add(u); order.append(u)
            stack.extend(v for v in graph[u] if v not in visited)
    return order
```

### Dijkstra
```python
import heapq
def dijkstra(graph, start):
    dist, heap, visited, pred = {start: 0}, [(0, start)], set(), {start: None}
    while heap:
        d, u = heapq.heappop(heap)
        if u in visited: continue
        visited.add(u)
        for v, w in graph[u]:
            if d + w < dist.get(v, float('inf')):
                dist[v] = d + w; pred[v] = u
                heapq.heappush(heap, (dist[v], v))
    return dist, pred
```

### Union-Find + Kruskal MST
```python
class UF:
    def __init__(self, n): self.p = list(range(n)); self.rank = [0]*n
    def find(self, x):
        if self.p[x] != x: self.p[x] = self.find(self.p[x])
        return self.p[x]
    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px == py: return False
        if self.rank[px] < self.rank[py]: px, py = py, px
        self.p[py] = px
        if self.rank[px] == self.rank[py]: self.rank[px] += 1
        return True

def kruskal(vertices, edges):
    idx = {v: i for i, v in enumerate(vertices)}
    uf = UF(len(vertices)); mst, total = [], 0
    for u, v, w in sorted(edges, key=lambda e: e[2]):
        if uf.union(idx[u], idx[v]):
            mst.append((u, v, w)); total += w
            if len(mst) == len(vertices) - 1: break
    return mst, total
```

### Topological Sort (Kahn)
```python
def topo_kahn(vertices, adj):
    in_deg = {v: 0 for v in vertices}
    for u in adj:
        for v in adj[u]: in_deg[v] += 1
    queue = deque(v for v in vertices if in_deg[v] == 0)
    result = []
    while queue:
        u = queue.popleft(); result.append(u)
        for v in adj[u]:
            in_deg[v] -= 1
            if in_deg[v] == 0: queue.append(v)
    return result if len(result) == len(vertices) else None  # None = cycle
```

### 0/1 Knapsack (space-optimized + traceback variant)
```python
def knapsack(weights, values, capacity):
    dp = [0] * (capacity + 1)
    for w, v in zip(weights, values):
        for c in range(capacity, w - 1, -1):  # backwards prevents reuse
            dp[c] = max(dp[c], dp[c - w] + v)
    return dp[capacity]
```

### Needleman-Wunsch (global alignment)
```python
import numpy as np
def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-2):
    m, n = len(s1), len(s2)
    dp = np.zeros((m+1, n+1), dtype=int)
    dp[:,0] = [i*gap for i in range(m+1)]
    dp[0,:] = [j*gap for j in range(n+1)]
    for i in range(1, m+1):
        for j in range(1, n+1):
            s = match if s1[i-1]==s2[j-1] else mismatch
            dp[i][j] = max(dp[i-1][j-1]+s, dp[i-1][j]+gap, dp[i][j-1]+gap)
    # Traceback
    a1, a2 = [], []; i, j = m, n
    while i>0 or j>0:
        if i>0 and j>0:
            s = match if s1[i-1]==s2[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1]+s:
                a1.append(s1[i-1]); a2.append(s2[j-1]); i-=1; j-=1; continue
        if i>0 and dp[i][j] == dp[i-1][j]+gap:
            a1.append(s1[i-1]); a2.append('-'); i-=1
        else:
            a1.append('-'); a2.append(s2[j-1]); j-=1
    return dp[m][n], ''.join(reversed(a1)), ''.join(reversed(a2))
```

### Smith-Waterman (local alignment)
```python
def smith_waterman(s1, s2, match=2, mismatch=-1, gap=-1):
    m, n = len(s1), len(s2)
    dp = np.zeros((m+1, n+1), dtype=int)
    best, pos = 0, (0,0)
    for i in range(1, m+1):
        for j in range(1, n+1):
            s = match if s1[i-1]==s2[j-1] else mismatch
            dp[i][j] = max(0, dp[i-1][j-1]+s, dp[i-1][j]+gap, dp[i][j-1]+gap)
            if dp[i][j] > best: best=dp[i][j]; pos=(i,j)
    a1, a2 = [], []; i, j = pos
    while dp[i][j] > 0:
        s = match if s1[i-1]==s2[j-1] else mismatch
        if dp[i][j]==dp[i-1][j-1]+s: a1.append(s1[i-1]); a2.append(s2[j-1]); i-=1; j-=1
        elif dp[i][j]==dp[i-1][j]+gap: a1.append(s1[i-1]); a2.append('-'); i-=1
        else: a1.append('-'); a2.append(s2[j-1]); j-=1
    return best, ''.join(reversed(a1)), ''.join(reversed(a2))
```

## Pitfalls

- **Dijkstra with negative weights** fails silently -- use Bellman-Ford
- **Topo sort on cyclic graph:** Kahn returns partial list; check `len(result) == len(vertices)`
- **DFS cycle detection (undirected):** track parent to avoid false positives
- **Knapsack iterate order:** 0/1 backwards; unbounded forwards
- **Knapsack traceback** requires full 2D table, not space-optimized 1D
- **SW traceback** stops at 0, not at corner -- do not reuse NW traceback logic
- **Affine gap alignment** needs 3 matrices (M/X/Y); single matrix gives wrong gap-run scores
- **MST on disconnected graph** produces a spanning forest, not a single tree

## Bioinformatics Applications

| Algorithm | Application |
|---|---|
| BFS k-hop | Proteins within 2 interactions of a hub (PPI) |
| Dijkstra log-space | Most reliable signaling path: `min(-log(conf))` |
| Kruskal MST | Simplified phylogenetic tree from pairwise distances |
| MST clustering | Single-linkage: cut k-1 longest edges for k clusters |
| Topo sort | Gene regulatory cascade order; pipeline scheduling |
| NW | Global alignment of same-length orthologs |
| SW | Domain search, short-read mapping |
| 0/1 Knapsack | Gene panel selection within sequencing budget |
