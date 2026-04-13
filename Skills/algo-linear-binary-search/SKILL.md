---
name: algo-linear-binary-search
description: "This notebook covers fundamental search algorithms, from the simple linear search to various binary search variants."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/03_Searching_Algorithms/01_linear_binary_search.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Linear and Binary Search Algorithms

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/03_Searching_Algorithms/01_linear_binary_search.ipynb`*


This notebook covers fundamental search algorithms, from the simple linear search to various binary search variants.

## Table of Contents
1. [Linear Search](#1-linear-search)
2. [Binary Search (Classic)](#2-binary-search-classic)
3. [Binary Search Variants](#3-binary-search-variants)
   - Find First Occurrence
   - Find Last Occurrence
   - Lower Bound / Upper Bound
   - Search in Rotated Sorted Array
4. [Performance Comparison](#4-performance-comparison)
5. [Common Pitfalls](#5-common-pitfalls)

```python
# Required imports for all sections
import numpy as np
import time
import random
import matplotlib.pyplot as plt
from typing import List, Optional, Tuple

%matplotlib inline
plt.style.use('seaborn-v0_8-whitegrid')
```python

---
## 1. Linear Search

### Theory

**When to use?**
- When the array is **unsorted** or cannot be sorted
- When the array is **small** (< 10-20 elements)
- When you need to search **only once** (sorting + binary search is overkill)
- When elements are stored in a **linked list** (no random access)

**Prerequisites:**
- None! Works on any collection

**Key insight:**
- Check each element one by one until found or end of array
- Simple but inefficient for large datasets

### ASCII Art Visualization

```python
Find 7 in [3, 1, 4, 1, 5, 9, 2, 6, 7]

Step 1: [3, 1, 4, 1, 5, 9, 2, 6, 7]
         ↑
         3 ≠ 7 → continue

Step 2: [3, 1, 4, 1, 5, 9, 2, 6, 7]
            ↑
            1 ≠ 7 → continue

Step 3: [3, 1, 4, 1, 5, 9, 2, 6, 7]
               ↑
               4 ≠ 7 → continue

...continues checking each element...

Step 9: [3, 1, 4, 1, 5, 9, 2, 6, 7]
                                 ↑
                                 7 = 7 → FOUND at index 8!

Summary:
         [3, 1, 4, 1, 5, 9, 2, 6, 7]
          ↑  ↑  ↑  ↑  ↑  ↑  ↑  ↑  ✓
          1  2  3  4  5  6  7  8  9 comparisons
```python

### Complexity Table

| Metric | Complexity | Notes |
|--------|------------|-------|
| Time (Best) | O(1) | Element at first position |
| Time (Average) | O(n/2) = O(n) | Element in the middle |
| Time (Worst) | O(n) | Element at end or not found |
| Space | O(1) | Only a few variables needed |

```python
def linear_search(arr: List[int], target: int) -> int:
    """
    Search for target in array using linear search.
    
    Args:
        arr: List of elements to search through
        target: Element to find
    
    Returns:
        Index of target if found, -1 otherwise
    
    Time Complexity: O(n)
    Space Complexity: O(1)
    """
    for i in range(len(arr)):
        if arr[i] == target:
            return i
    return -1


def linear_search_with_sentinel(arr: List[int], target: int) -> int:
    """
    Linear search optimized with sentinel value.
    Eliminates bound checking in the loop.
    
    Args:
        arr: List of elements (will be temporarily modified)
        target: Element to find
    
    Returns:
        Index of target if found, -1 otherwise
    """
    if not arr:
        return -1
    
    # Save the last element and place sentinel
    last = arr[-1]
    arr[-1] = target
    
    i = 0
    while arr[i] != target:
        i += 1
    
    # Restore the last element
    arr[-1] = last
    
    # Check if we found the target or just the sentinel
    if i < len(arr) - 1 or arr[-1] == target:
        return i
    return -1
```python

```python
def linear_search_verbose(arr: List[int], target: int) -> int:
    """Linear search with step-by-step output."""
    print(f"Searching for {target} in {arr}")
    print("-" * 50)
    
    for i in range(len(arr)):
        print(f"Step {i+1}: Checking index {i}, value = {arr[i]}", end="")
        if arr[i] == target:
            print(f" ✓ FOUND!")
            return i
        print(f" ✗ Continue")
    
    print(f"Element {target} not found after {len(arr)} comparisons")
    return -1

# Test cases
print("=" * 60)
print("TEST 1: Normal case")
print("=" * 60)
result = linear_search_verbose([3, 1, 4, 1, 5, 9, 2, 6, 7], 7)
print(f"\nResult: index {result}\n")
```python

```python
print("=" * 60)
print("TEST 2: Element at beginning (best case)")
print("=" * 60)
result = linear_search_verbose([7, 1, 4, 1, 5, 9, 2, 6, 3], 7)
print(f"\nResult: index {result}\n")
```python

```python
print("=" * 60)
print("TEST 3: Edge cases")
print("=" * 60)

# Empty array
print("\nEmpty array:")
print(f"linear_search([], 5) = {linear_search([], 5)}")

# Single element - found
print("\nSingle element (found):")
print(f"linear_search([5], 5) = {linear_search([5], 5)}")

# Single element - not found
print("\nSingle element (not found):")
print(f"linear_search([3], 5) = {linear_search([3], 5)}")

# Element not in array
print("\nElement not found:")
print(f"linear_search([1, 2, 3, 4], 10) = {linear_search([1, 2, 3, 4], 10)}")
```python

---
## 2. Binary Search (Classic)

### Theory

**When to use?**
- When the array is **sorted**
- When you need to search **multiple times** (sorting cost is amortized)
- When the array is **large** and O(n) is too slow
- When you need to find insertion points or bounds

**Prerequisites:**
- Array must be **sorted** in ascending (or descending) order
- Random access to elements (arrays, not linked lists)

**Key insight:**
- Divide search space in half with each comparison
- If middle element < target → search right half
- If middle element > target → search left half
- Each step eliminates half of remaining elements

### ASCII Art Visualization

```python
Find 7 in [1, 2, 3, 4, 5, 6, 7, 8, 9]
          indices: 0  1  2  3  4  5  6  7  8

Step 1: [1, 2, 3, 4, 5, 6, 7, 8, 9]
         L           M           R
         ↑           ↑           ↑
       left=0     mid=4      right=8
       
       arr[4] = 5
       5 < 7 → target is in RIGHT half
       Update: left = mid + 1 = 5

Step 2: [1, 2, 3, 4, 5, 6, 7, 8, 9]
                     L  M     R
                     ↑  ↑     ↑
                  left=5  mid=6  right=8
       
       arr[6] = 7
       7 = 7 → FOUND at index 6!

Summary: Only 2 comparisons vs 7 for linear search!

Search space reduction:
Step 1: 9 elements → 4 elements (eliminated 5)
Step 2: 4 elements → FOUND
```python

### Complexity Table

| Metric | Complexity | Notes |
|--------|------------|-------|
| Time (Best) | O(1) | Target is at middle |
| Time (Average) | O(log n) | Halving search space |
| Time (Worst) | O(log n) | Target at extreme or not found |
| Space (Iterative) | O(1) | Only a few variables |
| Space (Recursive) | O(log n) | Call stack depth |

```python
def binary_search_iterative(arr: List[int], target: int) -> int:
    """
    Search for target in sorted array using iterative binary search.
    
    Args:
        arr: Sorted list of elements
        target: Element to find
    
    Returns:
        Index of target if found, -1 otherwise
    
    Time Complexity: O(log n)
    Space Complexity: O(1)
    """
    left, right = 0, len(arr) - 1
    
    while left <= right:
        # Prevent integer overflow: mid = left + (right - left) // 2
        mid = (left + right) // 2
        
        if arr[mid] == target:
            return mid
        elif arr[mid] < target:
            left = mid + 1  # Search right half
        else:
            right = mid - 1  # Search left half
    
    return -1  # Not found


def binary_search_recursive(arr: List[int], target: int, 
                            left: int = None, right: int = None) -> int:
    """
    Search for target in sorted array using recursive binary search.
    
    Args:
        arr: Sorted list of elements
        target: Element to find
        left: Left boundary (inclusive)
        right: Right boundary (inclusive)
    
    Returns:
        Index of target if found, -1 otherwise
    
    Time Complexity: O(log n)
    Space Complexity: O(log n) - call stack
    """
    # Initialize boundaries on first call
    if left is None:
        left = 0
    if right is None:
        right = len(arr) - 1
    
    # Base case: search space exhausted
    if left > right:
        return -1
    
    mid = (left + right) // 2
    
    if arr[mid] == target:
        return mid
    elif arr[mid] < target:
        return binary_search_recursive(arr, target, mid + 1, right)
    else:
        return binary_search_recursive(arr, target, left, mid - 1)
```python

```python
def binary_search_verbose(arr: List[int], target: int) -> int:
    """Binary search with step-by-step visualization."""
    print(f"Searching for {target} in {arr}")
    print("-" * 60)
    
    left, right = 0, len(arr) - 1
    step = 1
    
    while left <= right:
        mid = (left + right) // 2
        
        # Visualization
        print(f"\nStep {step}:")
        print(f"  Array:  {arr}")
        print(f"  Indices: ", end="")
        for i in range(len(arr)):
            print(f"{i:>3}", end="")
        print()
        
        # Show pointers
        pointer_line = "           "
        for i in range(len(arr)):
            if i == left and i == mid and i == right:
                pointer_line += "LMR"
            elif i == left and i == mid:
                pointer_line += "LM "
            elif i == mid and i == right:
                pointer_line += "MR "
            elif i == left:
                pointer_line += "L  "
            elif i == mid:
                pointer_line += "M  "
            elif i == right:
                pointer_line += "R  "
            else:
                pointer_line += "   "
        print(pointer_line)
        
        print(f"  left={left}, mid={mid}, right={right}")
        print(f"  arr[{mid}] = {arr[mid]}")
        
        if arr[mid] == target:
            print(f"  {arr[mid]} == {target} → FOUND at index {mid}!")
            return mid
        elif arr[mid] < target:
            print(f"  {arr[mid]} < {target} → search RIGHT half")
            left = mid + 1
        else:
            print(f"  {arr[mid]} > {target} → search LEFT half")
            right = mid - 1
        
        step += 1
    
    print(f"\nElement {target} not found after {step-1} steps")
    return -1

# Test
print("=" * 70)
print("TEST: Finding 7 in sorted array")
print("=" * 70)
result = binary_search_verbose([1, 2, 3, 4, 5, 6, 7, 8, 9], 7)
```python

```python
print("=" * 70)
print("TEST: Element not in array")
print("=" * 70)
result = binary_search_verbose([1, 2, 3, 4, 5, 6, 8, 9, 10], 7)
```python

```python
print("=" * 60)
print("Edge Cases")
print("=" * 60)

# Empty array
print("\nEmpty array:")
print(f"binary_search_iterative([], 5) = {binary_search_iterative([], 5)}")

# Single element
print("\nSingle element (found):")
print(f"binary_search_iterative([5], 5) = {binary_search_iterative([5], 5)}")

# Single element (not found)
print("\nSingle element (not found):")
print(f"binary_search_iterative([3], 5) = {binary_search_iterative([3], 5)}")

# Compare iterative vs recursive
test_arr = [1, 3, 5, 7, 9, 11, 13, 15]
print(f"\nComparing iterative vs recursive on {test_arr}:")
for target in [1, 7, 15, 6]:
    iter_result = binary_search_iterative(test_arr, target)
    rec_result = binary_search_recursive(test_arr, target)
    print(f"  target={target:2d}: iterative={iter_result:2d}, recursive={rec_result:2d}")
```python

---
## 3. Binary Search Variants

### 3.1 Find First Occurrence

**When to use?**
- When array has **duplicate elements**
- Need the **leftmost** position of target

**Key insight:**
- When we find target, don't stop—continue searching left
- Keep track of the last found position

#### ASCII Art Visualization

```python
Find FIRST 5 in [1, 3, 5, 5, 5, 7, 9]
indices:         0  1  2  3  4  5  6

Step 1: [1, 3, 5, 5, 5, 7, 9]
         L        M        R
       
       arr[3] = 5 = target
       Found! But is it the FIRST?
       Save result=3, search LEFT: right = mid - 1 = 2

Step 2: [1, 3, 5, 5, 5, 7, 9]
         L  M  R
       
       arr[1] = 3 < 5 → search RIGHT
       left = mid + 1 = 2

Step 3: [1, 3, 5, 5, 5, 7, 9]
               LMR
       
       arr[2] = 5 = target
       Found! Update result=2, search LEFT: right = 1

Step 4: left (2) > right (1) → STOP

Answer: First occurrence at index 2
                   ↓
        [1, 3, 5, 5, 5, 7, 9]
               ↑  ↑  ↑
           first  │  last
```python

```python
def find_first_occurrence(arr: List[int], target: int) -> int:
    """
    Find the first (leftmost) occurrence of target in sorted array.
    
    Args:
        arr: Sorted list of elements (may contain duplicates)
        target: Element to find
    
    Returns:
        Index of first occurrence, -1 if not found
    
    Time Complexity: O(log n)
    Space Complexity: O(1)
    """
    left, right = 0, len(arr) - 1
    result = -1
    
    while left <= right:
        mid = (left + right) // 2
        
        if arr[mid] == target:
            result = mid  # Save this position
            right = mid - 1  # Continue searching left
        elif arr[mid] < target:
            left = mid + 1
        else:
            right = mid - 1
    
    return result


def find_first_occurrence_recursive(arr: List[int], target: int,
                                     left: int = None, right: int = None,
                                     result: int = -1) -> int:
    """
    Recursive version of find_first_occurrence.
    """
    if left is None:
        left = 0
    if right is None:
        right = len(arr) - 1
    
    if left > right:
        return result
    
    mid = (left + right) // 2
    
    if arr[mid] == target:
        # Found target, but continue searching left for earlier occurrence
        return find_first_occurrence_recursive(arr, target, left, mid - 1, mid)
    elif arr[mid] < target:
        return find_first_occurrence_recursive(arr, target, mid + 1, right, result)
    else:
        return find_first_occurrence_recursive(arr, target, left, mid - 1, result)


# Test
test_arr = [1, 3, 5, 5, 5, 7, 9]
print(f"Array: {test_arr}")
print(f"Find first 5: index {find_first_occurrence(test_arr, 5)}")
print(f"Find first 7: index {find_first_occurrence(test_arr, 7)}")
print(f"Find first 4: index {find_first_occurrence(test_arr, 4)}")
```python
