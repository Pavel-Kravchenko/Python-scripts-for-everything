---
name: algo-bfs-dfs
description: "1. Implement Breadth-First Search (BFS) and Depth-First Search (DFS) 2. Understand when to use each traversal algorithm 3. Apply graph traversals to biological network analysis 4. Find connected compo"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/02_bfs_dfs.ipynb"
---

# 🔍 Graph Traversals: BFS and DFS

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/02_bfs_dfs.ipynb`*

# 🔍 Graph Traversals: BFS and DFS

## Learning Objectives

By the end of this notebook, you will be able to:

1. Implement Breadth-First Search (BFS) and Depth-First Search (DFS)
2. Understand when to use each traversal algorithm
3. Apply graph traversals to biological network analysis
4. Find connected components and shortest paths

---

```python
from collections import defaultdict, deque
import matplotlib.pyplot as plt
import numpy as np

class Graph:
    """Simple undirected graph using adjacency list."""
    def __init__(self):
        self.adj = defaultdict(set)
    
    def add_edge(self, u, v):
        self.adj[u].add(v)
        self.adj[v].add(u)
    
    def neighbors(self, v):
        return self.adj[v]
    
    def vertices(self):
        return list(self.adj.keys())
```

## 1. Breadth-First Search (BFS)

BFS explores vertices level by level, using a **queue** (FIFO).

```
       A           Level 0: A
      /|\ 
     B C D         Level 1: B, C, D
    /|   |
   E F   G         Level 2: E, F, G
```

**Use cases:**
- Shortest path in unweighted graphs
- Finding all nodes within k hops (network neighborhood)
- Level-order traversal

```python
def bfs(graph, start):
    """
    Breadth-First Search from start vertex.
    Returns: (visited order, distances from start)
    Time: O(V + E), Space: O(V)
    """
    visited = set([start])
    queue = deque([start])
    order = []
    distances = {start: 0}
    
    while queue:
        vertex = queue.popleft()
        order.append(vertex)
        
        for neighbor in graph.neighbors(vertex):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
                distances[neighbor] = distances[vertex] + 1
    
    return order, distances

# Build a sample PPI network
ppi = Graph()
edges = [
    ('TP53', 'MDM2'), ('TP53', 'BRCA1'), ('TP53', 'ATM'),
    ('MDM2', 'RB1'), ('BRCA1', 'CHEK2'), ('ATM', 'CHEK2'),
    ('RB1', 'E2F1'), ('CHEK2', 'CDC25A')
]
for u, v in edges:
    ppi.add_edge(u, v)

order, distances = bfs(ppi, 'TP53')
print("BFS traversal from TP53:", order)
print("\nDistances from TP53:")
for gene, dist in sorted(distances.items(), key=lambda x: x[1]):
    print(f"  {gene}: {dist} hops")
```

```python
def find_neighbors_within_k(graph, start, k):
    """Find all vertices within k edges of start (network neighborhood)."""
    _, distances = bfs(graph, start)
    return [v for v, d in distances.items() if d <= k and v != start]

# Find genes within 2 hops of TP53
neighborhood = find_neighbors_within_k(ppi, 'TP53', 2)
print(f"Genes within 2 hops of TP53: {neighborhood}")
```

## 2. Depth-First Search (DFS)

DFS explores as far as possible along each branch before backtracking, using a **stack** (LIFO) or recursion.

**Use cases:**
- Cycle detection
- Topological sorting
- Finding connected components
- Path finding (any path, not necessarily shortest)

```python
def dfs_iterative(graph, start):
    """
    Depth-First Search using explicit stack.
    Time: O(V + E), Space: O(V)
    """
    visited = set()
    stack = [start]
    order = []
    
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            order.append(vertex)
            # Add neighbors in reverse for consistent ordering
            for neighbor in sorted(graph.neighbors(vertex), reverse=True):
                if neighbor not in visited:
                    stack.append(neighbor)
    
    return order

def dfs_recursive(graph, start, visited=None):
    """DFS using recursion."""
    if visited is None:
        visited = set()
    
    visited.add(start)
    order = [start]
    
    for neighbor in sorted(graph.neighbors(start)):
        if neighbor not in visited:
            order.extend(dfs_recursive(graph, neighbor, visited))
    
    return order

print("DFS iterative from TP53:", dfs_iterative(ppi, 'TP53'))
print("DFS recursive from TP53:", dfs_recursive(ppi, 'TP53'))
```

## 3. Finding Connected Components

In biological networks, connected components represent functional modules or pathways.

```python
def find_connected_components(graph):
    """Find all connected components using DFS."""
    visited = set()
    components = []
    
    for vertex in graph.vertices():
        if vertex not in visited:
            # Start new component
            component = []
            stack = [vertex]
            
            while stack:
                v = stack.pop()
                if v not in visited:
                    visited.add(v)
                    component.append(v)
                    for neighbor in graph.neighbors(v):
                        if neighbor not in visited:
                            stack.append(neighbor)
            
            components.append(component)
    
    return components

# Add a disconnected component
ppi.add_edge('KRAS', 'BRAF')
ppi.add_edge('BRAF', 'MEK1')

components = find_connected_components(ppi)
print(f"Found {len(components)} connected components:")
for i, comp in enumerate(components, 1):
    print(f"  Component {i}: {comp}")
```

## 4. BFS vs DFS Comparison

| Feature | BFS | DFS |
|---------|-----|-----|
| Data structure | Queue | Stack/Recursion |
| Explores | Level by level | Branch by branch |
| Shortest path | ✅ Yes (unweighted) | ❌ No |
| Memory | O(width) | O(depth) |
| Cycle detection | Possible | Natural fit |
| Topological sort | No | Yes |

## 🧬 Exercise 1: Shortest Path Between Genes

Implement a function to find the shortest path between two genes using BFS.

```python
def shortest_path(graph, start, end):
    """Find shortest path from start to end using BFS."""
    if start == end:
        return [start]
    
    visited = {start}
    queue = deque([(start, [start])])
    
    while queue:
        vertex, path = queue.popleft()
        
        for neighbor in graph.neighbors(vertex):
            if neighbor == end:
                return path + [neighbor]
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return None  # No path exists

# Test: Find path from E2F1 to ATM
path = shortest_path(ppi, 'E2F1', 'ATM')
print(f"Shortest path E2F1 → ATM: {' → '.join(path)}")
```

## 🧬 Exercise 2: Cycle Detection

Implement cycle detection using DFS. In biological networks, cycles can indicate feedback loops.

```python
def has_cycle(graph):
    """Detect if undirected graph has a cycle using DFS."""
    visited = set()
    
    def dfs(vertex, parent):
        visited.add(vertex)
        for neighbor in graph.neighbors(vertex):
            if neighbor not in visited:
                if dfs(neighbor, vertex):
                    return True
            elif neighbor != parent:
                return True  # Found cycle
        return False
    
    for v in graph.vertices():
        if v not in visited:
            if dfs(v, None):
                return True
    return False

# Add a cycle
ppi.add_edge('ATM', 'BRCA1')  # Creates TP53-ATM-BRCA1-TP53 cycle
print(f"Network has cycle: {has_cycle(ppi)}")
```

---

## Summary

✅ BFS uses a queue, explores level-by-level, finds shortest paths  
✅ DFS uses a stack/recursion, explores depth-first, good for cycles  
✅ Both have O(V + E) time complexity  
✅ Connected components reveal network modules

**Next:** [03_dijkstra.ipynb](03_dijkstra.ipynb) - Shortest paths with weighted edges
