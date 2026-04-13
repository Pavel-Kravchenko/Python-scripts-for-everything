---
name: algo-comparison-sorts
description: "This notebook provides comprehensive coverage of fundamental comparison-based sorting algorithms."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/02_Sorting_Algorithms/01_comparison_sorts.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Comparison-Based Sorting Algorithms

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/02_Sorting_Algorithms/01_comparison_sorts.ipynb`*

# Comparison-Based Sorting Algorithms

This notebook provides comprehensive coverage of fundamental comparison-based sorting algorithms.

## Table of Contents
1. [Bubble Sort](#1-bubble-sort) - Simple but inefficient
2. [Merge Sort](#2-merge-sort) - Divide and conquer, stable
3. [Shell Sort](#3-shell-sort) - Improved insertion sort
4. [QuickSort](#4-quicksort) - Recursive and iterative versions
5. [Performance Comparison](#5-performance-comparison)
6. [Summary](#6-summary)

```python
# Common imports used throughout the notebook
import time
import random
import copy
from typing import List, Optional, Tuple, Callable
import matplotlib.pyplot as plt
import numpy as np

%matplotlib inline
plt.style.use('seaborn-v0_8-whitegrid')
```

---

# 1. Bubble Sort

## Theory

**Bubble Sort** is one of the simplest sorting algorithms. It works by repeatedly stepping through the list, comparing adjacent elements and swapping them if they are in the wrong order. This process continues until the list is sorted.

### How It Works
1. Start at the beginning of the array
2. Compare each pair of adjacent elements
3. If they're in the wrong order, swap them
4. After each pass, the largest unsorted element "bubbles up" to its correct position
5. Repeat until no swaps are needed

### When to Use
- **Educational purposes** - excellent for understanding sorting concepts
- **Small datasets** - acceptable for arrays with < 50 elements
- **Nearly sorted data** - can be efficient with optimized early termination
- **Memory constrained** - requires O(1) extra space

### Why It's Useful
- Simple to understand and implement
- Stable sorting algorithm (preserves relative order of equal elements)
- In-place algorithm (no extra memory needed)
- Can detect already-sorted arrays in O(n)

## ASCII Art Visualization

```
Initial: [5, 3, 8, 1, 2]

Pass 1:
[5, 3, 8, 1, 2]    Compare 5 and 3
 ↑  ↑              5 > 3? Yes → Swap!
[3, 5, 8, 1, 2]    Compare 5 and 8
    ↑  ↑           5 > 8? No → Keep
[3, 5, 8, 1, 2]    Compare 8 and 1
       ↑  ↑        8 > 1? Yes → Swap!
[3, 5, 1, 8, 2]    Compare 8 and 2
          ↑  ↑     8 > 2? Yes → Swap!
[3, 5, 1, 2, 8]    ← 8 bubbled to end! ✓
             ^

Pass 2:
[3, 5, 1, 2, 8]    Compare 3 and 5
 ↑  ↑              3 > 5? No → Keep
[3, 5, 1, 2, 8]    Compare 5 and 1
    ↑  ↑           5 > 1? Yes → Swap!
[3, 1, 5, 2, 8]    Compare 5 and 2
       ↑  ↑        5 > 2? Yes → Swap!
[3, 1, 2, 5, 8]    ← 5 in place! ✓
          ^  ^

Pass 3:
[3, 1, 2, 5, 8]    Compare 3 and 1
 ↑  ↑              3 > 1? Yes → Swap!
[1, 3, 2, 5, 8]    Compare 3 and 2
    ↑  ↑           3 > 2? Yes → Swap!
[1, 2, 3, 5, 8]    ← 3 in place! ✓
       ^  ^  ^

Pass 4:
[1, 2, 3, 5, 8]    Compare 1 and 2
 ↑  ↑              1 > 2? No → Keep
[1, 2, 3, 5, 8]    ← No swaps! Array is sorted! ✓
 ^  ^  ^  ^  ^

Final: [1, 2, 3, 5, 8]
```

## Complexity Analysis

| Case | Time Complexity | Space Complexity | Stable? |
|------|-----------------|------------------|--------|
| Best | O(n) | O(1) | Yes |
| Average | O(n²) | O(1) | Yes |
| Worst | O(n²) | O(1) | Yes |

- **Best case O(n)**: Array is already sorted (with optimization)
- **Worst case O(n²)**: Array is sorted in reverse order
- **Space O(1)**: Only uses a constant amount of extra memory for swapping

```python
def bubble_sort(arr: List[int], verbose: bool = False) -> List[int]:
    """
    Sort an array using the Bubble Sort algorithm.
    
    Bubble Sort repeatedly steps through the list, compares adjacent elements,
    and swaps them if they are in the wrong order. The pass through the list
    is repeated until the list is sorted.
    
    Args:
        arr: List of integers to be sorted
        verbose: If True, print each step of the sorting process
    
    Returns:
        The sorted list (sorted in-place)
    
    Time Complexity:
        - Best: O(n) when array is already sorted
        - Average: O(n²)
        - Worst: O(n²) when array is reverse sorted
    
    Space Complexity: O(1) - in-place sorting
    
    Example:
        >>> bubble_sort([64, 34, 25, 12, 22, 11, 90])
        [11, 12, 22, 25, 34, 64, 90]
    """
    n = len(arr)
    
    if verbose:
        print(f"Initial array: {arr}")
        print("-" * 50)
    
    # Traverse through all array elements
    for i in range(n - 1):
        swapped = False  # Optimization: track if any swaps occurred
        
        if verbose:
            print(f"\nPass {i + 1}:")
        
        # Last i elements are already in place, so we reduce the range
        for j in range(n - 1 - i):
            # Compare adjacent elements
            if arr[j] > arr[j + 1]:
                # Swap if current element is greater than next
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
                swapped = True
                
                if verbose:
                    print(f"  Swapped {arr[j + 1]} and {arr[j]}: {arr}")
        
        # If no swapping occurred in this pass, array is sorted
        if not swapped:
            if verbose:
                print(f"  No swaps needed - array is sorted!")
            break
    
    if verbose:
        print("-" * 50)
        print(f"Final sorted array: {arr}")
    
    return arr
```

```python
# Example 1: Basic usage with step-by-step output
test_array = [5, 3, 8, 1, 2]
print("=== Bubble Sort Example ===")
bubble_sort(test_array.copy(), verbose=True)
```

```python
# Example 2: Already sorted array (best case)
sorted_array = [1, 2, 3, 4, 5]
print("\n=== Best Case: Already Sorted ===")
bubble_sort(sorted_array.copy(), verbose=True)
```

```python
# Example 3: Reverse sorted array (worst case)
reverse_array = [5, 4, 3, 2, 1]
print("\n=== Worst Case: Reverse Sorted ===")
bubble_sort(reverse_array.copy(), verbose=True)
```

---

# 2. Merge Sort

## Theory

**Merge Sort** is an efficient, stable, divide-and-conquer sorting algorithm. It divides the input array into two halves, recursively sorts them, and then merges the sorted halves.

### How It Works
1. **Divide**: Split the array into two halves
2. **Conquer**: Recursively sort each half
3. **Combine**: Merge the two sorted halves into one sorted array

### When to Use
- **Large datasets** - consistent O(n log n) performance
- **Linked lists** - very efficient for linked data structures
- **External sorting** - excellent for sorting data that doesn't fit in memory
- **Stability required** - preserves relative order of equal elements
- **Parallel processing** - easily parallelizable

### Why It's Useful
- Guaranteed O(n log n) time complexity in all cases
- Stable sorting algorithm
- Well-suited for external sorting (files, databases)
- Highly parallelizable

## ASCII Art Visualization

```
                    [38, 27, 43, 3, 9, 82, 10]
                              |
              ┌───────────────┴───────────────┐
              ↓                               ↓
        [38, 27, 43, 3]                 [9, 82, 10]
              |                               |
       ┌──────┴──────┐                 ┌──────┴──────┐
       ↓             ↓                 ↓             ↓
   [38, 27]      [43, 3]           [9, 82]        [10]
       |             |                 |             |
    ┌──┴──┐       ┌──┴──┐          ┌──┴──┐          |
    ↓     ↓       ↓     ↓          ↓     ↓          ↓
  [38]  [27]    [43]   [3]       [9]   [82]       [10]
    ↓     ↓       ↓     ↓          ↓     ↓          ↓
    └──┬──┘       └──┬──┘          └──┬──┘          |
       ↓             ↓                 ↓             |
   [27, 38]      [3, 43]           [9, 82]        [10]
       ↓             ↓                 ↓             ↓
       └──────┬──────┘                 └──────┬──────┘
              ↓                               ↓
       [3, 27, 38, 43]                  [9, 10, 82]
              ↓                               ↓
              └───────────────┬───────────────┘
                              ↓
                [3, 9, 10, 27, 38, 43, 82]


MERGE PROCESS DETAIL:
Merging [27, 38] and [3, 43]:

Step 1: Compare 27 and 3  → Take 3    Result: [3]
        [27, 38]  [43]
         ↑        ↑

Step 2: Compare 27 and 43 → Take 27   Result: [3, 27]
        [38]  [43]
         ↑     ↑

Step 3: Compare 38 and 43 → Take 38   Result: [3, 27, 38]
        []  [43]
             ↑

Step 4: Left exhausted    → Take 43   Result: [3, 27, 38, 43]
```

## Complexity Analysis

| Case | Time Complexity | Space Complexity | Stable? |
|------|-----------------|------------------|--------|
| Best | O(n log n) | O(n) | Yes |
| Average | O(n log n) | O(n) | Yes |
| Worst | O(n log n) | O(n) | Yes |

- **Time O(n log n)**: Dividing takes O(log n), merging takes O(n) at each level
- **Space O(n)**: Requires auxiliary array for merging
- **Stable**: Maintains relative order of equal elements

```python
def merge_sort(arr: List[int], verbose: bool = False, depth: int = 0) -> List[int]:
    """
    Sort an array using the Merge Sort algorithm.
    
    Merge Sort is a divide-and-conquer algorithm that:
    1. Divides the array into two halves
    2. Recursively sorts each half
    3. Merges the sorted halves back together
    
    Args:
        arr: List of integers to be sorted
        verbose: If True, print each step of the sorting process
        depth: Current recursion depth (for visualization)
    
    Returns:
        A new sorted list
    
    Time Complexity: O(n log n) in all cases
    Space Complexity: O(n) for the auxiliary arrays
    
    Example:
        >>> merge_sort([38, 27, 43, 3, 9, 82, 10])
        [3, 9, 10, 27, 38, 43, 82]
    """
    indent = "  " * depth  # For visualization
    
    if verbose:
        print(f"{indent}Splitting: {arr}")
    
    # Base case: arrays with 0 or 1 element are already sorted
    if len(arr) <= 1:
        return arr
    
    # Divide: find the middle point
    mid = len(arr) // 2
    
    # Conquer: recursively sort both halves
    left_half = merge_sort(arr[:mid], verbose, depth + 1)
    right_half = merge_sort(arr[mid:], verbose, depth + 1)
    
    # Combine: merge the sorted halves
    merged = merge(left_half, right_half)
    
    if verbose:
        print(f"{indent}Merged: {left_half} + {right_half} → {merged}")
    
    return merged


def merge(left: List[int], right: List[int]) -> List[int]:
    """
    Merge two sorted arrays into one sorted array.
    
    Args:
        left: First sorted list
        right: Second sorted list
    
    Returns:
        A new sorted list containing all elements from both inputs
    """
    result = []
    i = j = 0
    
    # Compare elements from both arrays and add smaller one
    while i < len(left) and j < len(right):
        if left[i] <= right[j]:  # <= ensures stability
            result.append(left[i])
            i += 1
        else:
            result.append(right[j])
            j += 1
    
    # Add remaining elements from left array (if any)
    while i < len(left):
        result.append(left[i])
        i += 1
    
    # Add remaining elements from right array (if any)
    while j < len(right):
        result.append(right[j])
        j += 1
    
    return result
```

```python
# Example 1: Basic usage with step-by-step output
test_array = [38, 27, 43, 3, 9, 82, 10]
print("=== Merge Sort Example ===")
print(f"Original: {test_array}")
print("\nSorting process:")
result = merge_sort(test_array.copy(), verbose=True)
print(f"\nSorted: {result}")
```

```python
# Example 2: Demonstrating stability
print("\n=== Stability Demonstration ===")
# Using tuples (value, original_index) to show stability
data = [(3, 'a'), (1, 'b'), (3, 'c'), (2, 'd'), (1, 'e')]
print(f"Original: {data}")
# Sort by first element only
sorted_data = sorted(data, key=lambda x: x[0])  # Python's sort is stable like merge sort
print(f"Stable sorted: {sorted_data}")
print("Note: (3, 'a') comes before (3, 'c'), and (1, 'b') comes before (1, 'e')")
```

---

# 3. Shell Sort

## Theory

**Shell Sort** is an optimization of Insertion Sort that allows exchanging items that are far apart. It starts by sorting elements far apart from each other and progressively reduces the gap between elements to be compared.

### How It Works
1. Choose a gap sequence (e.g., n/2, n/4, n/8, ..., 1)
2. For each gap, perform a gapped insertion sort
3. Reduce the gap and repeat
4. Final pass with gap=1 is a standard insertion sort

### When to Use
- **Medium-sized arrays** - good for n < 10,000
- **Memory constrained** - in-place with O(1) extra space
- **Embedded systems** - simple implementation, good cache performance
- **Nearly sorted data** - performs well on partially sorted arrays

### Why It's Useful
- Better than O(n²) in practice for most gap sequences
- In-place sorting (O(1) extra memory)
- Simple to implement
- Good for medium-sized arrays
