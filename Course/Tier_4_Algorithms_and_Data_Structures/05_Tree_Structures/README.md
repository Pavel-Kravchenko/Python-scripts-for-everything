# 05. Tree Structures

**Tier 4: Algorithms & Data Structures**

Hierarchical data structures that enable O(log n) search, insertion, and deletion. This module progresses from the basic Binary Search Tree through AVL self-balancing trees to Red-Black trees -- the structure underlying Python's `sortedcontainers` and many database indices. Phylogenetic trees in bioinformatics share the same traversal algorithms studied here.

## Topics Covered

- Binary Search Tree (BST): insert, delete, search, traversals
- Inorder, preorder, postorder, level-order traversals
- BST degeneration and worst-case O(n) behavior
- AVL trees: balance factor, LL/RR/LR/RL rotations, height guarantee
- Red-Black tree properties (coloring invariants)
- Red-Black insertion: recoloring and rotation cases
- Interval trees for overlapping range queries (extension)

## Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [01_binary_search_trees.ipynb](01_binary_search_trees.ipynb) | BST operations, all four traversals, deletion cases |
| 2 | [02_avl_trees.ipynb](02_avl_trees.ipynb) | Self-balancing, rotation types, height analysis |
| 3 | [03_red_black_trees.ipynb](03_red_black_trees.ipynb) | RB properties, insertion with recoloring and rotation |

## Complexity Reference

| Operation | BST (avg) | BST (worst) | AVL | Red-Black |
|-----------|-----------|-------------|-----|-----------|
| Search | O(log n) | O(n) | O(log n) | O(log n) |
| Insert | O(log n) | O(n) | O(log n) | O(log n) |
| Delete | O(log n) | O(n) | O(log n) | O(log n) |
| Space | O(n) | O(n) | O(n) | O(n) |

## Bioinformatics Connections

- Phylogenetic tree traversal uses the same inorder/preorder/postorder algorithms from `01_binary_search_trees.ipynb`. See [Phylogenetics](../../Tier_2_Core_Bioinformatics/06_Phylogenetics/01_phylogenetics.ipynb).
- Interval trees (an extension of BSTs) are used to query genomic annotations -- finding all genes that overlap a given chromosomal region.

## Prerequisites

- [04. Linear Data Structures](../04_Linear_Data_Structures/README.md) -- stacks and queues are used in traversal algorithms

---

[<< 04. Linear Data Structures](../04_Linear_Data_Structures/README.md) | [Tier 4 Overview](../README.md) | [06. Hash-Based Structures >>](../06_Hash_Based_Structures/README.md)
