---
name: algo-linear-sorts
description: "This notebook covers **non-comparison-based** sorting algorithms that can achieve **O(n)** time complexity under certain conditions."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/02_Sorting_Algorithms/02_linear_sorts.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Linear Time Sorting Algorithms

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/02_Sorting_Algorithms/02_linear_sorts.ipynb`*

# Linear Time Sorting Algorithms

This notebook covers **non-comparison-based** sorting algorithms that can achieve **O(n)** time complexity under certain conditions.

## Table of Contents
1. [Why Can Linear Sorts Beat O(n log n)?](#why-linear)
2. [Counting Sort](#counting-sort)
3. [Radix Sort](#radix-sort)
4. [Bucket Sort](#bucket-sort)
5. [Complexity Comparison](#complexity-comparison)
6. [Key Takeaways](#key-takeaways)

---

<a id="why-linear"></a>
## Why Can These Algorithms Beat O(n log n)?

### The Comparison Sort Lower Bound

**Theorem:** Any comparison-based sorting algorithm requires **Ω(n log n)** comparisons in the worst case.

**Proof Intuition:**
- With `n` elements, there are `n!` possible permutations
- Each comparison gives 1 bit of information (≤ or >)
- To distinguish `n!` orderings, we need at least `log₂(n!)` comparisons
- By Stirling's approximation: `log₂(n!) ≈ n log₂(n)`

### How Linear Sorts Escape This Bound

**They don't compare elements!** Instead, they use:
- **Direct addressing** (Counting Sort): Elements map directly to positions
- **Digit decomposition** (Radix Sort): Sort by individual digits
- **Distribution** (Bucket Sort): Elements distributed to buckets by value

### The Trade-off

| Advantage | Cost |
|-----------|------|
| O(n) time | Requires assumptions about input |
| Stable by design | Extra space for auxiliary structures |
| Predictable performance | Limited to specific data types |

---

<a id="counting-sort"></a>
## 1. Counting Sort

### Theory

**Counting Sort** sorts integers by counting occurrences of each value, then using arithmetic to place elements in their correct positions.

#### Assumptions & Limitations
- Input must be **integers** (or keys that map to integers)
- Range of values `k` must be known and **reasonable** (O(n) or close)
- If `k >> n`, space and time become O(k), making it inefficient

#### When to Use
- Sorting integers with a **small, known range** (e.g., ages 0-120, grades 0-100)
- When **stability** is required
- As a subroutine in **Radix Sort**

### ASCII Art Visualization

```
Input: [4, 2, 2, 8, 3, 3, 1]
Range: 1-8

Step 1: Count occurrences
┌─────────────────────────────────────────┐
│ Index: 0   1   2   3   4   5   6   7   8│
│ Count: 0   1   2   2   1   0   0   0   1│
│            ↑   ↑   ↑   ↑               ↑│
│           (1) (2) (3) (4)             (8)│
└─────────────────────────────────────────┘

Step 2: Cumulative sum (prefix sum)
┌─────────────────────────────────────────┐
│ Index: 0   1   2   3   4   5   6   7   8│
│ Count: 0   1   3   5   6   6   6   6   7│
│            │   │   │   │               │ │
│            ▼   ▼   ▼   ▼               ▼ │
│     "1 goes" "2s go" "3s go" "4 goes" "8 goes"│
│     to pos 0 to 1-2 to 3-4  to pos 5  to pos 6│
└─────────────────────────────────────────┘

Step 3: Build output (right to left for stability)
┌─────────────────────────────────────────┐
│ Processing: 1 ← 3 ← 3 ← 8 ← 2 ← 2 ← 4  │
│                                         │
│ Output: [1, 2, 2, 3, 3, 4, 8]           │
│          ↑  ↑  ↑  ↑  ↑  ↑  ↑            │
│         pos 0  1  2  3  4  5  6         │
└─────────────────────────────────────────┘
```

```python
from typing import List, Optional


def counting_sort(arr: List[int], max_val: Optional[int] = None) -> List[int]:
    """
    Sort an array of non-negative integers using Counting Sort.
    
    This is a stable, non-comparison-based sorting algorithm that works
    by counting the occurrences of each unique value.
    
    Args:
        arr: List of non-negative integers to sort
        max_val: Maximum value in array (computed if not provided)
    
    Returns:
        New sorted list
    
    Time Complexity: O(n + k) where k is the range of values
    Space Complexity: O(k) for the count array + O(n) for output
    
    Example:
        >>> counting_sort([4, 2, 2, 8, 3, 3, 1])
        [1, 2, 2, 3, 3, 4, 8]
    """
    if not arr:
        return []
    
    # Determine the range
    k = max_val if max_val is not None else max(arr)
    
    # Step 1: Count occurrences of each value
    count = [0] * (k + 1)
    for num in arr:
        count[num] += 1
    
    # Step 2: Compute cumulative counts (prefix sum)
    # count[i] now contains the number of elements <= i
    for i in range(1, k + 1):
        count[i] += count[i - 1]
    
    # Step 3: Build output array (iterate right-to-left for stability)
    output = [0] * len(arr)
    for num in reversed(arr):
        # count[num] - 1 gives the correct position for num
        count[num] -= 1
        output[count[num]] = num
    
    return output
```

```python
def counting_sort_verbose(arr: List[int]) -> List[int]:
    """
    Counting Sort with step-by-step visualization.
    
    Args:
        arr: List of non-negative integers to sort
    
    Returns:
        Sorted list
    """
    if not arr:
        return []
    
    print(f"Input array: {arr}")
    print(f"Length: {len(arr)}")
    print("=" * 60)
    
    k = max(arr)
    print(f"\nStep 1: Create count array (size = max + 1 = {k + 1})")
    
    # Count occurrences
    count = [0] * (k + 1)
    for i, num in enumerate(arr):
        count[num] += 1
        print(f"  arr[{i}] = {num} → count[{num}]++")
    
    print(f"\n  Count array: {count}")
    print(f"  Indices:     {list(range(len(count)))}")
    
    # Cumulative sum
    print("\nStep 2: Compute cumulative sum (prefix sum)")
    for i in range(1, k + 1):
        old_val = count[i]
        count[i] += count[i - 1]
        print(f"  count[{i}] = {old_val} + count[{i-1}] = {count[i]}")
    
    print(f"\n  Cumulative:  {count}")
    print(f"  Meaning: count[v] = number of elements ≤ v")
    
    # Build output
    print("\nStep 3: Build output (right to left for stability)")
    output = [0] * len(arr)
    
    for num in reversed(arr):
        pos = count[num] - 1
        count[num] -= 1
        output[pos] = num
        print(f"  Place {num} at position {pos} → output = {output}")
    
    print("\n" + "=" * 60)
    print(f"Sorted array: {output}")
    
    return output
```

```python
# Example: Step-by-step execution
print("COUNTING SORT DEMONSTRATION")
print("=" * 60)

test_array = [4, 2, 2, 8, 3, 3, 1]
result = counting_sort_verbose(test_array)
```

```python
# Test stability: elements with same key maintain relative order
print("\nSTABILITY TEST")
print("=" * 60)

# Using tuples (value, original_index) to track stability
def counting_sort_with_objects(arr: List[tuple]) -> List[tuple]:
    """Sort tuples by first element (key), demonstrating stability."""
    if not arr:
        return []
    
    k = max(item[0] for item in arr)
    count = [0] * (k + 1)
    
    for item in arr:
        count[item[0]] += 1
    
    for i in range(1, k + 1):
        count[i] += count[i - 1]
    
    output = [None] * len(arr)
    for item in reversed(arr):  # Right-to-left for stability!
        count[item[0]] -= 1
        output[count[item[0]]] = item
    
    return output

# Items: (key, identifier)
items = [(3, 'A'), (1, 'B'), (3, 'C'), (2, 'D'), (3, 'E')]
print(f"Input:  {items}")
print("        Note: Three items with key=3: A, C, E")

sorted_items = counting_sort_with_objects(items)
print(f"Output: {sorted_items}")
print("        Items with key=3 maintain order: A, C, E ✓")
```

---

<a id="radix-sort"></a>
## 2. Radix Sort

### Theory

**Radix Sort** sorts numbers by processing individual digits, from least significant to most significant (LSD) or vice versa (MSD).

#### Key Insight
If we use a **stable** sort for each digit position, the overall sort is correct.

#### Assumptions & Limitations
- Works on integers or fixed-length strings
- All numbers must have the same number of digits (pad with leading zeros)
- Requires a stable sorting algorithm for each digit (typically Counting Sort)

#### When to Use
- Sorting **large numbers** with many digits (where comparison would be expensive)
- Sorting **strings** of fixed length
- When `d × (n + k) < n log n` (d = digits, k = base/radix)

### ASCII Art Visualization

```
Input: [170, 45, 75, 90, 802, 24, 2, 66]

Pad to 3 digits: [170, 045, 075, 090, 802, 024, 002, 066]

═══════════════════════════════════════════════════════════════
Pass 1: Sort by 1s digit (rightmost)
═══════════════════════════════════════════════════════════════

  170  045  075  090  802  024  002  066
    ↓    ↓    ↓    ↓    ↓    ↓    ↓    ↓
    0    5    5    0    2    4    2    6

  Buckets:  0: [170, 090]
            2: [802, 002]
            4: [024]
            5: [045, 075]
            6: [066]

  Result: [170, 090, 802, 002, 024, 045, 075, 066]

═══════════════════════════════════════════════════════════════
Pass 2: Sort by 10s digit (middle)
═══════════════════════════════════════════════════════════════

  170  090  802  002  024  045  075  066
   ↓    ↓    ↓    ↓    ↓    ↓    ↓    ↓
   7    9    0    0    2    4    7    6

  Buckets:  0: [802, 002]
            2: [024]
            4: [045]
            6: [066]
            7: [170, 075]
            9: [090]

  Result: [802, 002, 024, 045, 066, 170, 075, 090]

═══════════════════════════════════════════════════════════════
Pass 3: Sort by 100s digit (leftmost)
═══════════════════════════════════════════════════════════════

  802  002  024  045  066  170  075  090
  ↓    ↓    ↓    ↓    ↓    ↓    ↓    ↓
  8    0    0    0    0    1    0    0

  Buckets:  0: [002, 024, 045, 066, 075, 090]
            1: [170]
            8: [802]

  Result: [002, 024, 045, 066, 075, 090, 170, 802]

═══════════════════════════════════════════════════════════════
Final: [2, 24, 45, 66, 75, 90, 170, 802] ✓
═══════════════════════════════════════════════════════════════
```

```python
def counting_sort_by_digit(arr: List[int], exp: int) -> List[int]:
    """
    Helper function: Sort array by a specific digit position using Counting Sort.
    
    Args:
        arr: List of non-negative integers
        exp: The digit position (1 for units, 10 for tens, 100 for hundreds, etc.)
    
    Returns:
        Array sorted by the specified digit
    
    Note:
        Uses base 10, so digit values are 0-9.
    """
    n = len(arr)
    output = [0] * n
    count = [0] * 10  # Digits 0-9
    
    # Count occurrences of each digit at position 'exp'
    for num in arr:
        digit = (num // exp) % 10
        count[digit] += 1
    
    # Cumulative count
    for i in range(1, 10):
        count[i] += count[i - 1]
    
    # Build output (right-to-left for stability)
    for num in reversed(arr):
        digit = (num // exp) % 10
        count[digit] -= 1
        output[count[digit]] = num
    
    return output


def radix_sort(arr: List[int]) -> List[int]:
    """
    Sort an array of non-negative integers using Radix Sort (LSD).
    
    Uses Counting Sort as the stable subroutine for each digit.
    Processes digits from least significant to most significant.
    
    Args:
        arr: List of non-negative integers to sort
    
    Returns:
        New sorted list
    
    Time Complexity: O(d × (n + k)) where:
        - d = number of digits in the maximum number
        - n = number of elements
        - k = base (10 for decimal)
    
    Space Complexity: O(n + k)
    
    Example:
        >>> radix_sort([170, 45, 75, 90, 802, 24, 2, 66])
        [2, 24, 45, 66, 75, 90, 170, 802]
    """
    if not arr:
        return []
    
    # Find the maximum number to determine the number of digits
    max_num = max(arr)
    
    # Make a copy to avoid modifying the original
    result = arr.copy()
    
    # Process each digit position (1s, 10s, 100s, ...)
    exp = 1
    while max_num // exp > 0:
        result = counting_sort_by_digit(result, exp)
        exp *= 10
    
    return result
```

```python
def radix_sort_verbose(arr: List[int]) -> List[int]:
    """
    Radix Sort with step-by-step visualization.
    
    Args:
        arr: List of non-negative integers to sort
    
    Returns:
        Sorted list
    """
    if not arr:
        return []
    
    print(f"Input array: {arr}")
    print("=" * 70)
    
    max_num = max(arr)
    num_digits = len(str(max_num))
    print(f"Maximum value: {max_num} ({num_digits} digits)")
    print(f"Number of passes required: {num_digits}")
    
    result = arr.copy()
    exp = 1
    pass_num = 1
    
    while max_num // exp > 0:
        print(f"\n{'='*70}")
        print(f"Pass {pass_num}: Sort by {'1s' if exp == 1 else f'{exp}s'} digit")
        print(f"{'='*70}")
        
        # Show digits being examined
        digits = [(num // exp) % 10 for num in result]
        print(f"\nCurrent array: {result}")
        print(f"Digits at position {exp}: {digits}")
        
        # Show bucket distribution
        buckets = {i: [] for i in range(10)}
        for num in result:
            digit = (num // exp) % 10
            buckets[digit].append(num)
        
        print(f"\nBucket distribution:")
        for digit, nums in buckets.items():
            if nums:
                print(f"  Bucket {digit}: {nums}")
        
        # Perform the sort
        result = counting_sort_by_digit(result, exp)
        print(f"\nAfter pass {pass_num}: {result}")
        
        exp *= 10
        pass_num += 1
    
    print(f"\n{'='*70}")
    print(f"FINAL SORTED ARRAY: {result}")
    print(f"{'='*70}")
    
    return result
```
