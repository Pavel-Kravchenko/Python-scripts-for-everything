# 01. Complexity Analysis

**Tier 4: Algorithms & Data Structures**

The mathematical foundation for reasoning about algorithm efficiency. Covers asymptotic notation, time and space complexity analysis, recurrence relations, and the basic algorithmic paradigms of recursion and iteration. Understanding complexity lets you choose the right tool when processing large genomic datasets.

## Topics Covered

- Big-O, Big-Omega, Big-Theta notation and their definitions
- Best, average, and worst-case analysis
- Common growth rate families: O(1), O(log n), O(n), O(n log n), O(n²), O(2ⁿ)
- Analyzing loops, nested loops, and recursive calls
- Recurrence relations and the Master Theorem
- Recursion fundamentals: call stack, base cases, recursive leap of faith
- Basic algorithmic paradigms: divide-and-conquer, decrease-and-conquer

## Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [01_complexity_analysis.ipynb](01_complexity_analysis.ipynb) | Asymptotic notation, growth rates, loop analysis, recurrences |
| 2 | [02_basic_algorithms.ipynb](02_basic_algorithms.ipynb) | Recursion, iteration, algorithm analysis practice |

## Complexity Reference

| Growth Rate | Name | Example |
|-------------|------|---------|
| O(1) | Constant | Hash table lookup |
| O(log n) | Logarithmic | Binary search |
| O(n) | Linear | Linear scan |
| O(n log n) | Linearithmic | Merge sort |
| O(n²) | Quadratic | Bubble sort |
| O(2ⁿ) | Exponential | Recursive subset enumeration |

## Bioinformatics Connections

- Choosing efficient tools when processing large genomes depends directly on understanding complexity trade-offs. All Tier 2--3 notebooks assume this fluency.

## Prerequisites

- Tier 1 Python (functions, loops, recursion basics)

---

[<< 00. Skills Check](../00_Skills_Check/README.md) | [Tier 4 Overview](../README.md) | [02. Sorting Algorithms >>](../02_Sorting_Algorithms/README.md)
