# Graph Algorithms

This module covers fundamental graph algorithms and data structures essential for solving connectivity, pathfinding, and optimization problems.

## Topics Covered

### 1. Graph Representations
- **Adjacency Matrix** - O(1) edge lookup, O(V²) space
- **Adjacency List** - O(V+E) space, efficient for sparse graphs
- **Edge List** - Simple representation, O(E) space

### 2. Graph Traversals
- **Breadth-First Search (BFS)** - Level-order traversal, shortest path in unweighted graphs
- **Depth-First Search (DFS)** - Explore as far as possible, backtrack

### 3. Shortest Path Algorithms
- **Dijkstra's Algorithm** - Single-source shortest path, non-negative weights
- **Bellman-Ford** - Handles negative weights, detects negative cycles
- **Floyd-Warshall** - All-pairs shortest paths

### 4. Minimum Spanning Trees
- **Kruskal's Algorithm** - Greedy edge selection with Union-Find
- **Prim's Algorithm** - Greedy vertex selection with Priority Queue

### 5. Advanced Topics
- **Topological Sort** - Ordering of DAG vertices
- **Cycle Detection** - Finding cycles in directed/undirected graphs
- **Strongly Connected Components** - Kosaraju's/Tarjan's algorithm

## Complexity Reference

| Algorithm | Time Complexity | Space Complexity | Use Case |
|-----------|-----------------|------------------|----------|
| BFS | O(V + E) | O(V) | Shortest path (unweighted) |
| DFS | O(V + E) | O(V) | Connectivity, cycles |
| Dijkstra (heap) | O((V+E) log V) | O(V) | Shortest path (non-negative) |
| Bellman-Ford | O(V × E) | O(V) | Shortest path (negative weights) |
| Floyd-Warshall | O(V³) | O(V²) | All-pairs shortest paths |
| Kruskal | O(E log E) | O(V) | MST (sparse graphs) |
| Prim | O((V+E) log V) | O(V) | MST (dense graphs) |
| Topological Sort | O(V + E) | O(V) | Task scheduling |

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_graph_representations.ipynb](01_graph_representations.ipynb) | Adjacency matrix, list, and edge list |
| [02_bfs_dfs.ipynb](02_bfs_dfs.ipynb) | Breadth-first and depth-first traversals |
| [03_dijkstra.ipynb](03_dijkstra.ipynb) | Shortest path with priority queue |
| [04_mst_kruskal_prim.ipynb](04_mst_kruskal_prim.ipynb) | Minimum spanning tree algorithms |
| [05_topological_sort.ipynb](05_topological_sort.ipynb) | DAG ordering and cycle detection |

## Interactive Visualizations

- [BFS & DFS Visualization](../interactive/graphs/bfs-dfs.html) - Watch graph traversals step by step
- [Dijkstra's Algorithm](../interactive/graphs/dijkstra.html) - See shortest paths computed live
- [MST Visualization](../interactive/graphs/mst.html) - Kruskal's and Prim's algorithms

## Key Concepts

### When to Use BFS vs DFS

| Scenario | Use |
|----------|-----|
| Shortest path in unweighted graph | BFS |
| Detecting cycles | DFS |
| Topological sorting | DFS |
| Connected components | Either |
| Minimum spanning tree | Neither (use Kruskal/Prim) |
| Finding any path | DFS (often simpler) |
| Level-order traversal | BFS |

### Graph Types

- **Directed vs Undirected** - Edges with/without direction
- **Weighted vs Unweighted** - Edges with/without costs
- **Cyclic vs Acyclic** - Contains/doesn't contain cycles
- **Connected vs Disconnected** - All vertices reachable from any vertex
- **Dense vs Sparse** - E ≈ V² vs E ≈ V
