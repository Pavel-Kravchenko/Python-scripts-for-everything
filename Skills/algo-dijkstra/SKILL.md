---
name: algo-dijkstra
description: "1. Understand weighted graphs and their biological applications 2. Implement Dijkstra's algorithm with a priority queue 3. Apply shortest path finding to confidence-weighted PPI networks"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/03_dijkstra.ipynb"
---

# ⚡ Dijkstra's Algorithm: Shortest Paths in Weighted Graphs

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/03_dijkstra.ipynb`*

# ⚡ Dijkstra's Algorithm: Shortest Paths in Weighted Graphs

## Learning Objectives

1. Understand weighted graphs and their biological applications
2. Implement Dijkstra's algorithm with a priority queue
3. Apply shortest path finding to confidence-weighted PPI networks

---

```python
import heapq
from collections import defaultdict

class WeightedGraph:
    """Weighted graph using adjacency list."""
    def __init__(self):
        self.adj = defaultdict(list)  # vertex -> [(neighbor, weight), ...]
    
    def add_edge(self, u, v, weight, directed=False):
        self.adj[u].append((v, weight))
        if not directed:
            self.adj[v].append((u, weight))
    
    def vertices(self):
        return list(self.adj.keys())
```

## 1. Dijkstra's Algorithm

Finds shortest paths from a source to all vertices in a graph with **non-negative** edge weights.

**Key idea:** Greedily select the vertex with minimum distance, update neighbors.

**Time complexity:** O((V + E) log V) with binary heap

```python
def dijkstra(graph, start):
    """
    Dijkstra's algorithm using min-heap.

    Works correctly for both undirected and directed graphs, including
    directed graphs where some vertices only appear as edge targets
    (no outgoing edges). Those vertices are discovered lazily.

    Time complexity: O((V + E) log V) with binary heap
    Space complexity: O(V)

    Parameters
    ----------
    graph : WeightedGraph
        Weighted graph with non-negative edge weights.
    start : hashable
        Source vertex.

    Returns
    -------
    distances : dict
        Shortest distance from start to every reachable vertex.
    predecessors : dict
        predecessors[v] = the vertex immediately before v on the
        shortest path from start.  Used by reconstruct_path().
    """
    # Initialise known vertices; vertices with no outgoing edges that are
    # only reachable as targets are added lazily via distances.get().
    distances = {v: float('inf') for v in graph.vertices()}
    distances[start] = 0
    predecessors = {start: None}

    # Min-heap: (tentative_distance, vertex)
    heap = [(0, start)]
    visited = set()

    while heap:
        dist, u = heapq.heappop(heap)

        # Skip stale entries (we may push the same vertex multiple times)
        if u in visited:
            continue
        visited.add(u)

        for v, weight in graph.adj[u]:
            if v in visited:
                continue
            new_dist = dist + weight
            # Use .get() with inf default so vertices that only appear as
            # targets in a directed graph (absent from graph.vertices()) are
            # handled correctly without a KeyError.
            if new_dist < distances.get(v, float('inf')):
                distances[v] = new_dist
                predecessors[v] = u
                heapq.heappush(heap, (new_dist, v))

    return distances, predecessors


def reconstruct_path(predecessors, end):
    """Reconstruct shortest path from predecessors dictionary.

    Parameters
    ----------
    predecessors : dict
        Output of dijkstra().
    end : hashable
        Destination vertex.

    Returns
    -------
    list
        Ordered list of vertices from source to end, or empty list if
        end is not reachable.
    """
    if end not in predecessors:
        return []
    path = []
    current = end
    while current is not None:
        path.append(current)
        current = predecessors.get(current)
    return path[::-1]
```

## Step-by-Step Walkthrough

The verbose version prints the priority queue decisions at each step, making it easy to trace why each distance estimate gets updated.

```python
def dijkstra_verbose(graph, start):
    """Dijkstra's algorithm with step-by-step output."""
    distances = {v: float('inf') for v in graph.vertices()}
    distances[start] = 0
    predecessors = {start: None}
    heap = [(0, start)]
    visited = set()
    step = 0

    while heap:
        dist, u = heapq.heappop(heap)
        if u in visited:
            continue
        visited.add(u)
        step += 1
        print(f"Step {step}: Visit {u} (distance {dist:.2f})")

        for v, weight in graph.adj[u]:
            if v in visited:
                continue
            new_dist = dist + weight
            # Use .get() to safely handle directed-graph target vertices
            if new_dist < distances.get(v, float('inf')):
                old = distances.get(v, float('inf'))
                distances[v] = new_dist
                predecessors[v] = u
                heapq.heappush(heap, (new_dist, v))
                old_str = f"{old:.2f}" if old != float('inf') else "inf"
                print(f"  Update {v}: {old_str} -> {new_dist:.2f} via {u}")

    return distances, predecessors

# Small demo graph
demo = WeightedGraph()
for u, v, w in [('A', 'B', 1.0), ('A', 'C', 4.0), ('B', 'C', 2.0), ('B', 'D', 5.0), ('C', 'D', 1.0)]:
    demo.add_edge(u, v, w)

print("Dijkstra walkthrough from A:")
dijkstra_verbose(demo, 'A')
```

```python
# Build weighted PPI network (weights = interaction confidence)
# Lower weight = higher confidence (we minimize distance)
ppi = WeightedGraph()
interactions = [
    ('TP53', 'MDM2', 0.01),   # Very high confidence
    ('TP53', 'BRCA1', 0.15),
    ('TP53', 'ATM', 0.08),
    ('MDM2', 'RB1', 0.20),
    ('BRCA1', 'CHEK2', 0.12),
    ('ATM', 'CHEK2', 0.05),
    ('RB1', 'E2F1', 0.10),
    ('CHEK2', 'CDC25A', 0.18),
]

for u, v, conf in interactions:
    # Convert confidence to distance: higher confidence = lower weight
    ppi.add_edge(u, v, 1 - conf)

# Find shortest paths from TP53
distances, predecessors = dijkstra(ppi, 'TP53')

print("Shortest path distances from TP53:")
for gene, dist in sorted(distances.items(), key=lambda x: x[1]):
    path = reconstruct_path(predecessors, gene)
    print(f"  {gene}: {dist:.2f} via {' -> '.join(path)}")
```

## Why Non-Negative Weights?

Dijkstra's greedy approach relies on a key invariant: once a vertex is popped from the priority queue and marked visited, its distance is final. This holds only when all edge weights are non-negative — a negative edge could create a shorter path to an already-finalized vertex, violating the invariant.

For graphs with negative weights, use **Bellman-Ford** (O(VE)) which relaxes all edges V-1 times and can also detect negative cycles.

## 🧬 Exercise: Most Reliable Path

Instead of minimizing distance, find the path that **maximizes** product of confidence scores (most reliable signaling pathway).

```python
import math

def most_reliable_path(graph, start, end, confidences):
    """
    Find path maximizing product of edge confidences.
    Hint: max(∏ conf) = max(∑ log(conf)) = min(∑ -log(conf))
    """
    # Convert to log space and use Dijkstra
    log_graph = WeightedGraph()
    for u, neighbors in graph.adj.items():
        for v, _ in neighbors:
            conf = confidences.get((u, v), confidences.get((v, u), 0.5))
            log_graph.add_edge(u, v, -math.log(conf), directed=True)
    
    distances, preds = dijkstra(log_graph, start)
    path = reconstruct_path(preds, end)
    reliability = math.exp(-distances[end])
    
    return path, reliability

# Test with original confidence scores
conf_scores = {('TP53', 'MDM2'): 0.99, ('TP53', 'BRCA1'): 0.85, ('TP53', 'ATM'): 0.92,
               ('MDM2', 'RB1'): 0.80, ('BRCA1', 'CHEK2'): 0.88, ('ATM', 'CHEK2'): 0.95,
               ('RB1', 'E2F1'): 0.90, ('CHEK2', 'CDC25A'): 0.82}

path, rel = most_reliable_path(ppi, 'TP53', 'CDC25A', conf_scores)
print(f"Most reliable path to CDC25A: {' -> '.join(path)}")
print(f"Path reliability: {rel:.4f}")
```

## Exercise 2 (2 stars): Network Diameter

The **diameter** of a weighted graph is the longest shortest path between any two vertices. It measures how "spread out" the network is — in a PPI network it tells you the maximum signaling distance.

Implement the function below and apply it to the PPI network.

```python
def graph_diameter(graph):
    """
    Find the diameter of a weighted graph.
    Returns (diameter, source, target) for the furthest-apart vertex pair.

    Parameters
    ----------
    graph : WeightedGraph

    Returns
    -------
    (diameter, source, target, path)
    """
    # YOUR CODE HERE
    # Hint: run Dijkstra from every vertex, track the maximum finite distance.
    pass

# Apply to the PPI network
# result = graph_diameter(ppi)
# print(f"PPI network diameter: {result[0]:.2f}")
# print(f"Furthest pair: {result[1]} -> {result[2]}")
# print(f"Path: {' -> '.join(result[3])}")
```
