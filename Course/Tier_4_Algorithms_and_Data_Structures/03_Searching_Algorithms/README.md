# 03. Searching Algorithms

**Tier 4: Algorithms & Data Structures**

Efficient search is at the heart of bioinformatics: every sequence alignment, interval query, and quality threshold check is a search problem. This module covers linear scan, binary search and its variants, and the two-pointer technique that turns many quadratic problems into linear ones.

## Topics Covered

- Linear search and its O(n) bound
- Binary search on sorted arrays: iterative and recursive forms
- Binary search variants: lower bound, upper bound, search for condition boundary
- Two-pointer technique: opposite ends, sliding window, fast/slow pointers
- Applications: searching sorted genomic intervals, bisecting quality thresholds

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_linear_binary_search.ipynb](01_linear_binary_search.ipynb) | Linear search, binary search variants, two-pointer techniques (50 cells) |

## Bioinformatics Connections

- Searching sorted genomic intervals (BED, GTF) uses binary search; `bisect` from the Python standard library implements the same logic.
- Bisecting quality thresholds (e.g., finding the minimal PHRED score cutoff that meets a coverage target) is a direct application of binary search on a monotone condition.

## Prerequisites

- [01. Complexity Analysis](../01_Complexity_Analysis/README.md)
- [02. Sorting Algorithms](../02_Sorting_Algorithms/README.md) -- binary search requires sorted input

---

[<< 02. Sorting Algorithms](../02_Sorting_Algorithms/README.md) | [Tier 4 Overview](../README.md) | [04. Linear Data Structures >>](../04_Linear_Data_Structures/README.md)
