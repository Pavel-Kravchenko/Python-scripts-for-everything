# Data Structure Operation Complexity Reference

Quick-reference card for time complexity of common operations.

## Core Data Structures

| Data Structure | Search | Insert | Delete | Notes |
|---|---|---|---|---|
| **Unsorted Array** | O(n) | O(1) amortized | O(n) | Insert at end is O(1); delete requires shifting |
| **Sorted Array** | O(log n) | O(n) | O(n) | Binary search for lookup; insert/delete require shifting |
| **Singly Linked List** | O(n) | O(1) head | O(n) | Insert at head is O(1); delete requires traversal to find predecessor |
| **Doubly Linked List** | O(n) | O(1) | O(1) given node | Delete is O(1) if you have a pointer to the node |
| **Stack** | O(n) | O(1) push | O(1) pop | LIFO access only |
| **Queue** | O(n) | O(1) enqueue | O(1) dequeue | FIFO access only |

## Tree Structures

| Data Structure | Search | Insert | Delete | Notes |
|---|---|---|---|---|
| **Binary Search Tree** | O(h) | O(h) | O(h) | h = height; worst case O(n) for skewed tree |
| **AVL Tree** | O(log n) | O(log n) | O(log n) | Guaranteed balanced; at most 2 rotations per insert |
| **Red-Black Tree** | O(log n) | O(log n) | O(log n) | Looser balance than AVL; used in most standard libraries |
| **Binary Heap** | O(n) | O(log n) | O(log n) | Extract-min/max is O(log n); build heap is O(n) |

## Hash-Based Structures

| Data Structure | Search | Insert | Delete | Notes |
|---|---|---|---|---|
| **Hash Table (chaining)** | O(1 + alpha) | O(1 + alpha) | O(1 + alpha) | alpha = load factor = n/m |
| **Hash Table (open addressing)** | O(1/(1-alpha)) | O(1/(1-alpha)) | O(1/(1-alpha)) | Degrades as load factor approaches 1 |
| **Bloom Filter** | O(k) | O(k) | N/A | k = number of hash functions; probabilistic; no false negatives |

## Sorting Algorithms

| Algorithm | Best | Average | Worst | Space | Stable |
|---|---|---|---|---|---|
| **Bubble Sort** | O(n) | O(n^2) | O(n^2) | O(1) | Yes |
| **Selection Sort** | O(n^2) | O(n^2) | O(n^2) | O(1) | No |
| **Insertion Sort** | O(n) | O(n^2) | O(n^2) | O(1) | Yes |
| **Shell Sort** | O(n log n) | depends on gap | O(n^2) | O(1) | No |
| **Merge Sort** | O(n log n) | O(n log n) | O(n log n) | O(n) | Yes |
| **QuickSort** | O(n log n) | O(n log n) | O(n^2) | O(log n) | No |
| **Counting Sort** | O(n + k) | O(n + k) | O(n + k) | O(k) | Yes |
| **Radix Sort** | O(d(n + k)) | O(d(n + k)) | O(d(n + k)) | O(n + k) | Yes |

## Graph Algorithms

| Algorithm | Time | Space | Notes |
|---|---|---|---|
| **BFS** | O(V + E) | O(V) | Shortest path in unweighted graphs |
| **DFS** | O(V + E) | O(V) | Cycle detection, topological sort |
| **Dijkstra (binary heap)** | O((V + E) log V) | O(V) | Non-negative weights only |
| **Bellman-Ford** | O(VE) | O(V) | Handles negative weights; detects negative cycles |
| **Kruskal (MST)** | O(E log E) | O(V) | Edge-centric; uses Union-Find |
| **Prim (MST)** | O((V + E) log V) | O(V) | Vertex-centric; uses priority queue |
| **Topological Sort** | O(V + E) | O(V) | DAGs only |

## String Algorithms

| Algorithm | Preprocessing | Search | Notes |
|---|---|---|---|
| **Naive** | None | O(nm) | Simple but slow |
| **KMP** | O(m) | O(n) | Prefix function; no backtracking on text |
| **Rabin-Karp** | O(m) | O(n) avg, O(nm) worst | Rolling hash; good for multi-pattern |
| **Aho-Corasick** | O(sum of patterns) | O(n + matches) | Multi-pattern automaton |
| **Suffix Array** | O(n log n) | O(m log n) | Space-efficient alternative to suffix trees |
| **Suffix Tree** | O(n) | O(m) | Powerful but memory-intensive |

---

*Based on 5th semester algorithms course materials, Lomonosov Moscow State University*
