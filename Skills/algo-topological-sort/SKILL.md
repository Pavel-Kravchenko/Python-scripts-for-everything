---
name: algo-topological-sort
description: "Topological sort — DFS and Kahn's BFS algorithms for DAG ordering, cycle detection, critical path analysis"
tool_type: python
primary_tool: Python
---

# Topological Sort

Linear ordering of vertices in a DAG such that for every edge u->v, u appears before v. Returns None if graph has a cycle.

## DFS-Based — O(V + E)

```python
from collections import defaultdict, deque

class DirectedGraph:
    def __init__(self):
        self.adj = defaultdict(list)
        self.vertices = set()
    def add_edge(self, u, v):
        self.adj[u].append(v)
        self.vertices.add(u)
        self.vertices.add(v)

def topological_sort_dfs(graph):
    """Returns ordering or None if cycle exists."""
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {v: WHITE for v in graph.vertices}
    result = []
    def dfs(v):
        color[v] = GRAY
        for neighbor in graph.adj[v]:
            if color[neighbor] == GRAY: return False  # cycle
            if color[neighbor] == WHITE:
                if not dfs(neighbor): return False
        color[v] = BLACK
        result.append(v)
        return True
    for v in graph.vertices:
        if color[v] == WHITE:
            if not dfs(v): return None
    return result[::-1]
```

## Kahn's Algorithm (BFS) — O(V + E)

Process vertices with in-degree 0 first. If result has fewer vertices than graph, a cycle exists.

```python
def topological_sort_kahn(graph):
    in_degree = {v: 0 for v in graph.vertices}
    for u in graph.adj:
        for v in graph.adj[u]:
            in_degree[v] += 1
    queue = deque(v for v in graph.vertices if in_degree[v] == 0)
    result = []
    while queue:
        v = queue.popleft()
        result.append(v)
        for neighbor in graph.adj[v]:
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)
    return result if len(result) == len(graph.vertices) else None
```

Use a min-heap instead of deque to get the lexicographically smallest valid ordering.

## Critical Path Analysis

Longest path in a weighted DAG = minimum execution time with unlimited parallelism.

```python
def critical_path(graph, durations):
    """Returns (makespan, critical_nodes)."""
    order = topological_sort_kahn(graph)
    earliest = {v: 0 for v in graph.vertices}
    for v in order:
        for neighbor in graph.adj[v]:
            earliest[neighbor] = max(earliest[neighbor],
                                     earliest[v] + durations.get(v, 0))
    makespan = max(earliest[v] + durations.get(v, 0) for v in graph.vertices)
    # Backward pass for slack
    latest = {v: makespan - durations.get(v, 0) for v in graph.vertices}
    for v in reversed(order):
        for neighbor in graph.adj[v]:
            latest[v] = min(latest[v], latest[neighbor] - durations.get(v, 0))
    critical_nodes = [v for v in order if earliest[v] == latest[v]]
    return makespan, critical_nodes
```

## Pitfalls

- DFS topo sort detects cycles via GRAY nodes (back edges); Kahn's detects cycles via incomplete result
- A DAG can have many valid topological orderings — neither algorithm guarantees a unique result
- DFS approach uses O(V) stack space for deep graphs
