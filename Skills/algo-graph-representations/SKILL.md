---
name: algo-graph-representations
description: Graph representations — adjacency matrix, adjacency list, edge list — with complexity tables, implementation patterns, and biological network use cases.
tool_type: python
primary_tool: NumPy
---

# Graph Representations

## Complexity Comparison

| Operation | Adjacency Matrix | Adjacency List | Edge List |
|-----------|-----------------|----------------|------------|
| Space | O(V²) | O(V + E) | O(E) |
| Add edge | O(1) | O(1) | O(1) |
| Remove edge | O(1) | O(degree) | O(E) |
| Check edge | O(1) | O(degree) w/ list; O(1) w/ set | O(E) |
| Get neighbors | O(V) | O(1) access + O(degree) iterate | O(E) |
| Iterate all edges | O(V²) | O(E) | O(E) |

## Decision Table

| Scenario | Use |
|----------|-----|
| Dense graph (E ≈ V²) | Adjacency Matrix |
| Sparse graph (most bio networks) | Adjacency List |
| Frequent edge existence queries | Adjacency Matrix |
| BFS / DFS traversal | Adjacency List |
| MST (Kruskal's), file I/O | Edge List |
| Small V, matrix operations needed | Adjacency Matrix |

**Most biological networks (PPI, GRN, metabolic) are sparse — adjacency lists are almost always the right choice.**

## Bio Network Vertex/Edge Mapping

| Application | Vertices | Edges |
|------------|----------|-------|
| PPI network | Proteins | Physical interactions |
| Metabolic network | Metabolites | Reactions |
| Gene regulatory | Genes / TFs | Regulatory relationships |
| de Bruijn graph | k-mers | Overlaps |
| Phylogenetic tree | Species / sequences | Evolutionary relationships |

## Implementations

### Adjacency Matrix
```python
import numpy as np

class AdjacencyMatrix:
    def __init__(self, vertices: list):
        self.vertices = vertices
        self.idx = {v: i for i, v in enumerate(vertices)}
        self.n = len(vertices)
        self.matrix = np.zeros((self.n, self.n), dtype=int)

    def add_edge(self, u, v, directed: bool = False):
        i, j = self.idx[u], self.idx[v]
        self.matrix[i][j] = 1
        if not directed:
            self.matrix[j][i] = 1

    def has_edge(self, u, v) -> bool:          # O(1)
        return bool(self.matrix[self.idx[u]][self.idx[v]])

    def neighbors(self, u) -> list:            # O(V)
        i = self.idx[u]
        return [self.vertices[j] for j in range(self.n) if self.matrix[i][j]]

    def degree(self, u) -> int:
        return int(np.sum(self.matrix[self.idx[u]]))
```

### Adjacency List
```python
from collections import defaultdict

class AdjacencyList:
    def __init__(self):
        self.graph: dict = defaultdict(set)  # sets give O(1) edge lookup

    def add_edge(self, u, v, directed: bool = False):
        self.graph[u].add(v)
        if not directed:
            self.graph[v].add(u)

    def has_edge(self, u, v) -> bool:         # O(1) with sets
        return v in self.graph.get(u, set())

    def neighbors(self, u) -> set:            # O(1) access
        return self.graph.get(u, set())

    def degree(self, u) -> int:
        return len(self.graph.get(u, []))

    def edges(self) -> list:                  # undirected: each edge once
        seen = set()
        result = []
        for u in self.graph:
            for v in self.graph[u]:
                key = tuple(sorted([str(u), str(v)]))
                if key not in seen:
                    seen.add(key)
                    result.append((u, v))
        return result
```

### Edge List
```python
class EdgeList:
    def __init__(self):
        self.edges: list[tuple] = []
        self._vertices: set = set()

    def add_edge(self, u, v, weight: float = 1.0):
        self.edges.append((u, v, weight))
        self._vertices.update([u, v])

    def sorted_by_weight(self) -> list:       # for Kruskal's MST
        return sorted(self.edges, key=lambda e: e[2])
```

### Loading from file
```python
from pathlib import Path
from collections import defaultdict

def load_edge_list(filepath: str) -> dict[str, list[tuple]]:
    """Parse 'node_a  node_b  weight' file into adjacency list."""
    graph: dict = defaultdict(list)
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            u, v, w = line.split()
            graph[u].append((v, float(w)))
            graph[v].append((u, float(w)))  # undirected
    return dict(graph)
```

## Memory Analysis (Human PPI: 20,000 proteins, 300,000 interactions)

```python
V, E = 20_000, 300_000
matrix_mb = V * V * 8 / 1e6              # ~3,200 MB
list_mb   = (2 * E * 8 + V * 56) / 1e6  # ~5.9 MB
# Ratio: ~542x — adjacency list wins for sparse graphs
```

## Pitfalls

- **Adjacency matrix for large sparse graphs**: 20k proteins → 3.2 GB matrix for a network that fits in 6 MB as an adjacency list.
- **List vs set for adjacency list**: using a `list` makes `has_edge` O(degree); use `set` for O(1) lookup if duplicate edges are not needed.
- **Undirected graphs and `edges()`**: naively iterating gives each edge twice. Track seen pairs with a set or only iterate `u < v`.
- **Vertex not in graph**: `graph[v]` on a `defaultdict` silently creates an empty entry. Use `graph.get(v, set())` to avoid spurious vertices.
- **Edge list for BFS/DFS**: finding neighbors requires scanning all edges — O(E) per node. Use an adjacency list instead.
