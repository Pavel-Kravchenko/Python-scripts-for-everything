---
name: algo-topological-sort
description: "1. Understand directed acyclic graphs (DAGs) 2. Implement topological sort using DFS and Kahn's algorithm 3. Apply to biological pathway analysis"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/05_topological_sort.ipynb"
---

# 📐 Topological Sort

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/05_topological_sort.ipynb`*

# 📐 Topological Sort

## Learning Objectives

1. Understand directed acyclic graphs (DAGs)
2. Implement topological sort using DFS and Kahn's algorithm
3. Apply to biological pathway analysis

---

## Applications in Biology

- **Gene regulatory networks:** Order of gene activation
- **Metabolic pathways:** Order of enzymatic reactions
- **Workflow scheduling:** Bioinformatics pipeline execution

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
```

```python
def topological_sort_dfs(graph):
    """
    Topological sort using DFS.
    Returns ordering or None if cycle exists.
    """
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {v: WHITE for v in graph.vertices}
    result = []
    
    def dfs(v):
        color[v] = GRAY
        for neighbor in graph.adj[v]:
            if color[neighbor] == GRAY:
                return False  # Cycle detected
            if color[neighbor] == WHITE:
                if not dfs(neighbor):
                    return False
        color[v] = BLACK
        result.append(v)
        return True
    
    for v in graph.vertices:
        if color[v] == WHITE:
            if not dfs(v):
                return None
    
    return result[::-1]

# Gene regulatory cascade
grn = DirectedGraph()
regulations = [
    ('Growth_Factor', 'RAS'),
    ('RAS', 'RAF'),
    ('RAF', 'MEK'),
    ('MEK', 'ERK'),
    ('ERK', 'Transcription_Factors'),
    ('Transcription_Factors', 'Cell_Cycle_Genes'),
]

for u, v in regulations:
    grn.add_edge(u, v)

order = topological_sort_dfs(grn)
print("Signal transduction order:")
for i, gene in enumerate(order, 1):
    print(f"  {i}. {gene}")
```

## Multiple Valid Orderings

A DAG can have many valid topological orderings — any permutation of vertices that is consistent with all edge directions is valid. For the linear cascade above there is only one ordering, but branching DAGs have many.

Kahn's algorithm with a **priority queue** (min-heap) instead of a plain queue yields the lexicographically smallest valid ordering. Enumerating all valid orderings requires backtracking.

```python
def all_topological_sorts(graph):
    """Find all valid topological orderings (small graphs only)."""
    in_degree = {v: 0 for v in graph.vertices}
    for u in graph.adj:
        for v in graph.adj[u]:
            in_degree[v] += 1

    results = []
    visited = set()
    current = []

    def backtrack():
        if len(current) == len(graph.vertices):
            results.append(current[:])
            return
        for v in sorted(graph.vertices):
            if v not in visited and in_degree[v] == 0:
                visited.add(v)
                current.append(v)
                for neighbor in graph.adj[v]:
                    in_degree[neighbor] -= 1
                backtrack()
                current.pop()
                visited.remove(v)
                for neighbor in graph.adj[v]:
                    in_degree[neighbor] += 1

    backtrack()
    return results

# Small branching DAG: A->C, B->C, C->D, C->E
small = DirectedGraph()
for u, v in [('A', 'C'), ('B', 'C'), ('C', 'D'), ('C', 'E')]:
    small.add_edge(u, v)

orderings = all_topological_sorts(small)
print(f"Found {len(orderings)} valid topological orderings:")
for o in orderings:
    print(" -> ".join(o))
```

```python
def topological_sort_kahn(graph):
    """
    Kahn's algorithm: BFS-based topological sort.
    Also detects cycles.
    """
    in_degree = {v: 0 for v in graph.vertices}
    for u in graph.adj:
        for v in graph.adj[u]:
            in_degree[v] += 1
    
    # Start with nodes having no incoming edges
    queue = deque([v for v in graph.vertices if in_degree[v] == 0])
    result = []
    
    while queue:
        v = queue.popleft()
        result.append(v)
        
        for neighbor in graph.adj[v]:
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)
    
    if len(result) != len(graph.vertices):
        return None  # Cycle exists
    
    return result

print("\nKahn's algorithm order:")
for i, gene in enumerate(topological_sort_kahn(grn), 1):
    print(f"  {i}. {gene}")
```

## Critical Path Analysis

In a pipeline DAG where each node has an estimated runtime, the **critical path** is the longest path from source to sink. It determines the minimum total execution time even with unlimited parallelism — no matter how many machines you have, the steps on the critical path must execute sequentially.

Algorithm: process vertices in topological order, propagating the earliest possible start time forward.

```python
def critical_path(graph, durations):
    """
    Find the critical path (longest path) in a weighted DAG.

    Parameters
    ----------
    graph     : DirectedGraph
    durations : dict {vertex: runtime}

    Returns
    -------
    (makespan, critical_nodes)
    """
    order = topological_sort_kahn(graph)
    # earliest[v] = earliest time v can START
    earliest = {v: 0 for v in graph.vertices}

    for v in order:
        finish_v = earliest[v] + durations.get(v, 0)
        for neighbor in graph.adj[v]:
            if finish_v > earliest[neighbor]:
                earliest[neighbor] = finish_v

    # Makespan = max finish time across all nodes
    makespan = max(earliest[v] + durations.get(v, 0) for v in graph.vertices)

    # Trace back: a node is on the critical path if removing it would reduce makespan
    # Simple approach: find all nodes whose slack == 0
    # latest[v] = latest time v can START without delaying makespan
    latest = {v: makespan - durations.get(v, 0) for v in graph.vertices}
    for v in reversed(order):
        for neighbor in graph.adj[v]:
            candidate = latest[neighbor] - durations.get(v, 0)
            if candidate < latest[v]:
                latest[v] = candidate

    critical_nodes = [v for v in order if earliest[v] == latest[v]]
    return makespan, critical_nodes


# Pipeline with estimated runtimes (minutes)
pipeline_dag = DirectedGraph()
pipeline_steps = [
    ('FastQC', 'Trimmomatic'),
    ('Trimmomatic', 'STAR_Alignment'),
    ('STAR_Alignment', 'FeatureCounts'),
    ('STAR_Alignment', 'BAM_QC'),
    ('FeatureCounts', 'DESeq2'),
    ('DESeq2', 'Pathway_Analysis'),
    ('BAM_QC', 'MultiQC'),
    ('FastQC', 'MultiQC'),
]
for u, v in pipeline_steps:
    pipeline_dag.add_edge(u, v)

runtimes = {
    'FastQC': 5, 'Trimmomatic': 20, 'STAR_Alignment': 45,
    'FeatureCounts': 10, 'BAM_QC': 8, 'DESeq2': 15,
    'Pathway_Analysis': 12, 'MultiQC': 3,
}

makespan, critical = critical_path(pipeline_dag, runtimes)
print(f"Minimum pipeline runtime: {makespan} minutes")
print(f"Critical path: {' -> '.join(critical)}")
```

## 🧬 Exercise: Pipeline Scheduling

Given a bioinformatics pipeline with dependencies, find valid execution order.

```python
pipeline = DirectedGraph()
steps = [
    ('FastQC', 'Trimmomatic'),
    ('Trimmomatic', 'STAR_Alignment'),
    ('STAR_Alignment', 'FeatureCounts'),
    ('STAR_Alignment', 'BAM_QC'),
    ('FeatureCounts', 'DESeq2'),
    ('DESeq2', 'Pathway_Analysis'),
    ('BAM_QC', 'MultiQC'),
    ('FastQC', 'MultiQC'),
]

for u, v in steps:
    pipeline.add_edge(u, v)

execution_order = topological_sort_kahn(pipeline)
print("Pipeline execution order:")
for i, step in enumerate(execution_order, 1):
    print(f"  {i}. {step}")
```

## Exercise 2 (2 stars): Parallel Scheduling

Given the pipeline DAG with runtimes, determine which steps can run in parallel at each time point. Print a Gantt-chart-style schedule showing when each step starts and finishes, and which steps overlap.

```python
def parallel_schedule(graph, durations):
    """
    Compute the earliest start and finish time for each step.
    Print a Gantt-chart style schedule.

    Parameters
    ----------
    graph     : DirectedGraph
    durations : dict {vertex: runtime}
    """
    # YOUR CODE HERE
    # Hint: use topological order to compute earliest start times,
    # then group steps by start time to show parallel slots.
    pass

# Test with the pipeline above
parallel_schedule(pipeline_dag, runtimes)
```
