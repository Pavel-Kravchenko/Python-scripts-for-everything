# 04. Linear Data Structures

**Tier 4: Algorithms & Data Structures**

The building blocks of algorithm design: linked lists, stacks, queues, and dynamic arrays. Understanding how these structures allocate and move data is essential for reasoning about the performance of higher-level tools. Includes amortized analysis of dynamic array growth -- the same principle that explains Python list performance.

## Topics Covered

- Singly linked lists: insertion, deletion, reversal, cycle detection
- Doubly linked lists: bidirectional traversal, O(1) removal
- Circular linked lists and their applications
- Stacks (LIFO): array-based and linked-list-based implementations
- Stack applications: balanced parentheses, function call simulation, expression evaluation
- Queues (FIFO): circular array, linked-list implementations
- Deques and priority queues
- Dynamic arrays: growth factor, amortized O(1) append analysis

## Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [01_linked_lists.ipynb](01_linked_lists.ipynb) | Singly/doubly/circular linked lists, operations, cycle detection |
| 2 | [02_stacks_queues.ipynb](02_stacks_queues.ipynb) | Stack and queue implementations, applications |
| 3 | [03_dynamic_arrays.ipynb](03_dynamic_arrays.ipynb) | Amortized analysis, resizing strategy |

## Complexity Reference

| Structure | Access | Search | Insert (head/tail) | Delete (head/tail) | Space |
|-----------|--------|--------|--------------------|--------------------|-------|
| Singly Linked List | O(n) | O(n) | O(1) / O(n) | O(1) / O(n) | O(n) |
| Doubly Linked List | O(n) | O(n) | O(1) / O(1) | O(1) / O(1) | O(n) |
| Stack (array) | O(n) | O(n) | O(1) amortized | O(1) | O(n) |
| Queue (circular) | O(n) | O(n) | O(1) | O(1) | O(n) |
| Dynamic Array | O(1) | O(n) | O(1) amortized | O(n) | O(n) |

## Prerequisites

- [01. Complexity Analysis](../01_Complexity_Analysis/README.md) -- amortized analysis uses Big-O

---

[<< 03. Searching Algorithms](../03_Searching_Algorithms/README.md) | [Tier 4 Overview](../README.md) | [05. Tree Structures >>](../05_Tree_Structures/README.md)
