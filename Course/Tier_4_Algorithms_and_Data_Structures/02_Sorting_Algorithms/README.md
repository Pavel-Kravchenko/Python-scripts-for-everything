# 02. Sorting Algorithms

**Tier 4: Algorithms & Data Structures**

A systematic study of sorting from O(n²) comparison sorts to O(n) linear-time sorts. Each algorithm is implemented step-by-step with visualizations, making their mechanics and complexity trade-offs concrete. Sorting is one of the most-studied problems in CS and a foundation for understanding algorithm design.

## Topics Covered

- Bubble sort, Selection sort, Insertion sort
- Shell sort and gap sequences
- Merge sort (divide-and-conquer)
- QuickSort (partitioning, pivot selection, worst-case avoidance)
- Counting sort, Radix sort, Bucket sort
- Stability, in-place vs. auxiliary space, adaptive behavior
- Lower bound proof: why comparison sorts cannot beat O(n log n)

## Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [01_comparison_sorts.ipynb](01_comparison_sorts.ipynb) | Bubble, Selection, Insertion, Shell, Merge, QuickSort with step-by-step visualization |
| 2 | [02_linear_sorts.ipynb](02_linear_sorts.ipynb) | Counting sort, Radix sort, Bucket sort |

## Complexity Reference

| Algorithm | Best | Average | Worst | Space | Stable |
|-----------|------|---------|-------|-------|--------|
| Bubble Sort | O(n) | O(n²) | O(n²) | O(1) | Yes |
| Selection Sort | O(n²) | O(n²) | O(n²) | O(1) | No |
| Insertion Sort | O(n) | O(n²) | O(n²) | O(1) | Yes |
| Shell Sort | O(n log n) | depends | O(n²) | O(1) | No |
| Merge Sort | O(n log n) | O(n log n) | O(n log n) | O(n) | Yes |
| QuickSort | O(n log n) | O(n log n) | O(n²) | O(log n) | No |
| Counting Sort | O(n+k) | O(n+k) | O(n+k) | O(k) | Yes |
| Radix Sort | O(nk) | O(nk) | O(nk) | O(n+k) | Yes |

## Bioinformatics Connections

- BAM coordinate sorting uses a merge-sort-based approach; samtools sort relies on the same principles covered in `01_comparison_sorts.ipynb`. See [Variant Calling](../../Tier_3_Applied_Bioinformatics/02_Variant_Calling_and_SNP_Analysis/01_variant_calling_and_snp_analysis.ipynb).
- Variant prioritization (ranking variants by quality score) uses counting and radix sort ideas for integer-keyed data.

## Prerequisites

- [01. Complexity Analysis](../01_Complexity_Analysis/README.md) -- Big-O notation required to read the complexity table

---

[<< 01. Complexity Analysis](../01_Complexity_Analysis/README.md) | [Tier 4 Overview](../README.md) | [03. Searching Algorithms >>](../03_Searching_Algorithms/README.md)
