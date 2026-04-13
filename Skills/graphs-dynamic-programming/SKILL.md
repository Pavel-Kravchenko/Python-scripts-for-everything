---
name: graphs-dynamic-programming
description: Graph algorithms (BFS, DFS, Dijkstra, MST, topological sort) and dynamic programming (memoization, tabulation, knapsack, sequence alignment)
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Graph Algorithms & Dynamic Programming

## When to Use

**Graphs:**
- Network traversal, reachability, connected components → BFS/DFS
- Shortest path in unweighted graph → BFS
- Shortest path with non-negative weights → Dijkstra
- Minimum spanning tree (clustering, phylogeny) → Kruskal or Prim
- Ordering tasks with dependencies → Topological sort
- PPI, metabolic networks, gene regulatory networks → adjacency list

**Dynamic Programming:**
- Problem has overlapping subproblems + optimal substructure → DP
- Counting paths/alignments, minimum edits → tabulation
- Sequence alignment (Needleman-Wunsch, Smith-Waterman) → 2D DP table
- Resource selection under constraint → knapsack

---

## Quick Reference

### Graph Representation Trade-offs

| | Adjacency Matrix | Adjacency List | Edge List |
|---|---|---|---|
| Space | O(V²) | O(V+E) | O(E) |
| Edge lookup | O(1) | O(degree) | O(E) |
| Neighbors | O(V) | O(degree) | O(E) |
| Best for | Dense, frequent lookups | Sparse (bio networks) | MST, I/O |

Biological networks (PPI ~20k proteins, ~300k edges) are sparse → adjacency list saves ~170x memory.

### Algorithm Complexities

| Algorithm | Time | Space | Notes |
|---|---|---|---|
| BFS / DFS | O(V+E) | O(V) | |
| Dijkstra | O((V+E) log V) | O(V) | Non-negative weights only |
| Kruskal | O(E log E) | O(V) | Better for sparse |
| Prim | O((V+E) log V) | O(V) | Better for dense |
| Topo sort (Kahn) | O(V+E) | O(V) | Detects cycles |
| Knapsack 0/1 | O(nW) | O(W) | Rolling array |
| Edit distance / NW | O(mn) | O(mn) | O(min(m,n)) with rolling row |

---

## Key Patterns

### DP Framework (4 questions)
1. **State:** What does `dp[i]` or `dp[i][j]` represent?
2. **Base case:** Smallest subproblem with known answer
3. **Recurrence:** How does `dp[i]` depend on smaller states?
4. **Answer extraction:** Where is the final answer? (corner, max cell, traceback)

### Memoization vs Tabulation

| | Memoization (top-down) | Tabulation (bottom-up) |
|---|---|---|
| Style | Recursion + cache | Iteration + table |
| Subproblems | Only needed | All |
| Stack overflow | Possible | No |
| Space opt | Harder | Easier (rolling array) |

Use `@functools.lru_cache(maxsize=None)` to add memoization with one line.

### Sequence Alignment (DP)
- **NW (global):** Init borders with gap penalties; traceback from bottom-right
- **SW (local):** Floor cells at 0; traceback from max cell
- **Affine gaps:** Three matrices M/X/Y; open cost once, extend per extra base
- Use log-space for max-product paths: `min(-log(conf))` via Dijkstra

---

## Code Templates

### Adjacency List (weighted, undirected)
```python
from collections import defaultdict
graph = defaultdict(list)  # {node: [(neighbor, weight), ...]}
graph[u].append((v, w)); graph[v].append((u, w))
```python

### BFS — level order, shortest hops
```python
from collections import deque
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
```python

### DFS — iterative
```python
def dfs(graph, start):
    visited, stack, order = set(), [start], []
    while stack:
        u = stack.pop()
        if u not in visited:
            visited.add(u); order.append(u)
            stack.extend(v for v in graph[u] if v not in visited)
    return order
```python

### Dijkstra
```python
import heapq
def dijkstra(graph, start):
    dist = {start: 0}
    heap = [(0, start)]
    visited = set()
    pred = {start: None}
    while heap:
        d, u = heapq.heappop(heap)
        if u in visited: continue
        visited.add(u)
        for v, w in graph[u]:
            if d + w < dist.get(v, float('inf')):
                dist[v] = d + w; pred[v] = u
                heapq.heappush(heap, (dist[v], v))
    return dist, pred

def path(pred, end):
    p = []; v = end
    while v is not None: p.append(v); v = pred.get(v)
    return p[::-1]
```python

### Union-Find (for Kruskal)
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
```python

### Kruskal's MST
```python
def kruskal(vertices, edges):  # edges: [(u, v, w)]
    idx = {v: i for i, v in enumerate(vertices)}
    uf = UF(len(vertices))
    mst, total = [], 0
    for u, v, w in sorted(edges, key=lambda e: e[2]):
        if uf.union(idx[u], idx[v]):
            mst.append((u, v, w)); total += w
            if len(mst) == len(vertices) - 1: break
    return mst, total
```python

### Prim's MST
```python
def prim(vertices, adj):  # adj: {v: [(neighbor, w)]}
    start = vertices[0]; visited = {start}; heap = []; mst = []; total = 0
    for v, w in adj[start]: heapq.heappush(heap, (w, start, v))
    while heap and len(visited) < len(vertices):
        w, u, v = heapq.heappop(heap)
        if v in visited: continue
        visited.add(v); mst.append((u, v, w)); total += w
        for nb, nw in adj[v]:
            if nb not in visited: heapq.heappush(heap, (nw, v, nb))
    return mst, total
```python

### Topological Sort — Kahn's Algorithm
```python
from collections import defaultdict, deque
def topo_kahn(vertices, adj):  # adj: {u: [v, ...]}
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
```python

### 0/1 Knapsack (space-optimized)
```python
def knapsack(weights, values, capacity):
    dp = [0] * (capacity + 1)
    for w, v in zip(weights, values):
        for c in range(capacity, w - 1, -1):  # backwards prevents reuse
            dp[c] = max(dp[c], dp[c - w] + v)
    return dp[capacity]

def knapsack_items(weights, values, names, capacity):
    n = len(weights)
    dp = [[0]*(capacity+1) for _ in range(n+1)]
    for i in range(1, n+1):
        for c in range(capacity+1):
            dp[i][c] = dp[i-1][c]
            if weights[i-1] <= c:
                dp[i][c] = max(dp[i][c], dp[i-1][c-weights[i-1]] + values[i-1])
    selected = []; c = capacity
    for i in range(n, 0, -1):
        if dp[i][c] != dp[i-1][c]:
            selected.append(names[i-1]); c -= weights[i-1]
    return dp[n][capacity], selected[::-1]
```python

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
```python

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
```python

---

## Common Pitfalls

- **Dijkstra with negative weights:** fails silently — use Bellman-Ford instead
- **Topological sort on cyclic graph:** Kahn returns partial list; check `len(result) == len(vertices)`
- **DFS cycle detection (undirected):** track parent to avoid false positives from the edge you came from
- **Knapsack iterate order:** 0/1 → backwards; unbounded → forwards
- **Knapsack traceback:** requires full 2D table, not the space-optimized 1D version
- **SW traceback:** stops at 0, not at corner — don't use NW traceback logic
- **Affine gap alignment:** needs 3 matrices (M/X/Y); single matrix gives wrong scores for gap runs
- **MST on disconnected graph:** Kruskal will produce a spanning forest, not a single tree

---

## Bioinformatics Connections

| Algorithm | Application |
|---|---|
| BFS k-hop neighborhood | Find all proteins within 2 interactions of a hub gene (PPI network) |
| Dijkstra (log-space) | Most reliable signaling path: `min(-log(confidence))` |
| Kruskal MST | Simplified phylogenetic tree from pairwise evolutionary distances |
| MST clustering | Single-linkage hierarchical clustering; cut k-1 longest edges → k clusters |
| Topological sort | Gene regulatory cascade order; RNA-seq pipeline scheduling |
| Critical path (DAG) | Minimum bioinformatics pipeline runtime with unlimited parallelism |
| Needleman-Wunsch | Global alignment of same-length orthologs |
| Smith-Waterman | Local alignment — domain search, short-read mapping |
| Edit distance | Sequence divergence metric; traceback gives exact mutation operations |
| LCS | Conservation score: `LCS_len / max(len(s1), len(s2))` — fast pre-filter before full alignment |
| 0/1 Knapsack | Gene panel selection within sequencing budget |
| Subset sum | Check if exact sequencing cost target is achievable |

---

## Related Skills

- `biopython-databases` — BioPython pairwise aligner (wraps NW/SW), BLAST
- `numpy-pandas-wrangling` — NumPy arrays for DP matrices
- `python-collections-regex` — `heapq`, `defaultdict`, `functools.lru_cache`
- `complexity-sorting-searching` — complexity classes, Big-O analysis for graph/DP algorithms
