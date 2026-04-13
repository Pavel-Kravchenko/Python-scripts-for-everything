---
name: algo-complexity-analysis
description: "This notebook provides a comprehensive introduction to analyzing algorithm efficiency using Big O notation and related concepts."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/01_Complexity_Analysis/01_complexity_analysis.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# 📊 Algorithmic Complexity Analysis & Big O Notation

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/01_Complexity_Analysis/01_complexity_analysis.ipynb`*

# 📊 Algorithmic Complexity Analysis & Big O Notation

This notebook provides a comprehensive introduction to analyzing algorithm efficiency using Big O notation and related concepts.

---

## Table of Contents
1. [Theory: What is Algorithmic Complexity?](#1-theory-section)
2. [Growth Rate Visualizations](#2-growth-rate-visualizations)
3. [Code Examples by Complexity Class](#3-code-examples-by-complexity-class)
4. [Complexity Analysis Practice](#4-complexity-analysis-practice)
5. [Space Complexity](#5-space-complexity)
6. [Key Takeaways](#6-key-takeaways)

---

# 1. Theory Section

## What is Algorithmic Complexity?

**Algorithmic complexity** measures the resources (time and space) an algorithm requires as a function of the input size. It helps us:

- **Compare algorithms** independently of hardware
- **Predict performance** at scale
- **Make informed decisions** about which algorithm to use

### Why Do We Care About Efficiency?

Consider searching through a dataset:

```
┌─────────────────────────────────────────────────────────────────┐
│  Input Size (n)  │  O(n) Linear  │  O(n²) Quadratic           │
├─────────────────────────────────────────────────────────────────┤
│       10         │    10 ops     │      100 ops               │
│      100         │   100 ops     │    10,000 ops              │
│    1,000         │ 1,000 ops     │ 1,000,000 ops              │
│  1,000,000       │   1M ops      │ 1,000,000,000,000 ops  ⚠️  │
└─────────────────────────────────────────────────────────────────┘
```

**At scale, algorithm choice matters more than hardware speed!**

## Big O, Big Ω, and Big Θ Notation

These mathematical notations describe algorithm behavior:

```
╔═══════════════════════════════════════════════════════════════════════╗
║  Notation   │  Name         │  Meaning                               ║
╠═══════════════════════════════════════════════════════════════════════╣
║  O(f(n))    │  Big O        │  Upper bound (worst case)              ║
║             │               │  "grows no faster than f(n)"           ║
╠═══════════════════════════════════════════════════════════════════════╣
║  Ω(f(n))    │  Big Omega    │  Lower bound (best case)               ║
║             │               │  "grows at least as fast as f(n)"      ║
╠═══════════════════════════════════════════════════════════════════════╣
║  Θ(f(n))    │  Big Theta    │  Tight bound (average case)            ║
║             │               │  "grows exactly as fast as f(n)"       ║
╚═══════════════════════════════════════════════════════════════════════╝
```

### Visual Representation

```
    Time │
         │         ╱ O(f(n)) - Upper bound
         │        ╱
         │       ╱
         │    ══╱══ Θ(f(n)) - Tight bound (actual growth)
         │     ╱
         │    ╱
         │   ╱_____ Ω(f(n)) - Lower bound
         │  ╱
         └────────────────────── n
```

### In Practice

We typically use **Big O** because:
1. We want to know the **worst-case** scenario
2. It gives us an **upper bound guarantee**
3. It's the most commonly discussed in interviews and documentation

## Common Complexity Classes

From fastest to slowest:

```
╔════════════════════════════════════════════════════════════════════════════╗
║  Complexity   │  Name           │  Example Operations                     ║
╠════════════════════════════════════════════════════════════════════════════╣
║  O(1)         │  Constant       │  Array access, hash lookup              ║
║  O(log n)     │  Logarithmic    │  Binary search                          ║
║  O(n)         │  Linear         │  Simple loop, linear search             ║
║  O(n log n)   │  Linearithmic   │  Efficient sorting (merge, heap)        ║
║  O(n²)        │  Quadratic      │  Nested loops, bubble sort              ║
║  O(n³)        │  Cubic          │  Triple nested loops                    ║
║  O(2ⁿ)        │  Exponential    │  Recursive Fibonacci, subsets           ║
║  O(n!)        │  Factorial      │  Permutations, traveling salesman       ║
╚════════════════════════════════════════════════════════════════════════════╝
```

### Rule: Drop Constants and Lower-Order Terms

```
O(2n + 5)       →  O(n)
O(n² + n)       →  O(n²)
O(500)          →  O(1)
O(n/2)          →  O(n)
O(n² + n³)      →  O(n³)
```

**Why?** At large n, only the dominant term matters.

---

# 2. Growth Rate Visualizations

## ASCII Chart: Comparing Growth Rates

```
Time │
     │                                              ⁂ O(2ⁿ)
     │                                         ⁂
     │                                    ⁂
     │                               ⁂
     │                          ⁂
     │                                          ★ O(n²)
     │                     ⁂               ★
     │                               ★
     │                ⁂        ★
     │                    ★           ──────── O(n log n)
     │           ⁂   ★      ──────
     │       ⁂  ★  ─────
     │    ⁂ ★ ────              ═══════════ O(n)
     │  ⁂★────      ═══════════
     │ ★──   ═══════
     │──═════               ················ O(log n)
     │═·············
     │·                     ________________ O(1)
     └───────────────────────────────────────────── n
            10    20    30    40    50    60
```

## Relative Growth Comparison

```
For n = 64:

O(1)        │█                                        = 1
O(log n)    │██████                                   = 6
O(n)        │████████████████████████████████████████ = 64
O(n log n)  │████████████████████████████████████████████████████████████ ≈ 384
O(n²)       │████████████████████████████████████████ ×100 MORE... = 4,096
O(2ⁿ)       │ WOULD NOT FIT ON ANY SCREEN! ≈ 18 quintillion
```

```python
import math

def visualize_growth_rates(n_values=None):
    """
    Print a table showing how different complexity classes scale.
    
    Time Complexity: O(k) where k is len(n_values)
    Space Complexity: O(1)
    """
    if n_values is None:
        n_values = [1, 10, 100, 1000, 10000]
    
    print("\n" + "="*80)
    print(f"{'n':>8} | {'O(1)':>10} | {'O(log n)':>10} | {'O(n)':>10} | {'O(n log n)':>12} | {'O(n²)':>15}")
    print("="*80)
    
    for n in n_values:
        o_1 = 1
        o_log_n = round(math.log2(n), 2) if n > 0 else 0
        o_n = n
        o_n_log_n = round(n * math.log2(n), 2) if n > 0 else 0
        o_n2 = n * n
        
        print(f"{n:>8} | {o_1:>10} | {o_log_n:>10} | {o_n:>10,} | {o_n_log_n:>12,.0f} | {o_n2:>15,}")
    
    print("="*80)

visualize_growth_rates()
```

---

# 3. Code Examples by Complexity Class

Each example includes:
- Function implementation
- Detailed comments explaining the complexity
- Docstring with time and space complexity

## O(1) - Constant Time

The operation takes the same time regardless of input size.

```
Time │
     │ ________________________________
     │
     └────────────────────────────────── n
```

```python
def get_first_element(arr):
    """
    Return the first element of an array.
    
    Time Complexity: O(1) - Direct index access is constant time
    Space Complexity: O(1) - No additional space used
    
    Args:
        arr: List of elements
    Returns:
        First element or None if empty
    """
    # Array index access is O(1) because:
    # - Arrays store elements in contiguous memory
    # - Address calculation: base_address + (index * element_size)
    # - This calculation is constant regardless of array size
    if len(arr) == 0:
        return None
    return arr[0]  # O(1) - same time for array of 10 or 10 million elements


def is_even(number):
    """
    Check if a number is even.
    
    Time Complexity: O(1) - Single arithmetic operation
    Space Complexity: O(1) - No additional space used
    """
    # Modulo operation is O(1) - it's a single CPU instruction
    # The size of the number doesn't affect the operation time
    # (for fixed-width integers)
    return number % 2 == 0


def hash_table_lookup(hash_table, key):
    """
    Look up a value in a hash table.
    
    Time Complexity: O(1) average - Hash computation + direct access
    Space Complexity: O(1) - No additional space used
    
    Note: Worst case is O(n) if all keys hash to the same bucket,
    but with good hash functions this is extremely rare.
    """
    # 1. Compute hash of key - O(1)
    # 2. Access bucket at computed index - O(1)
    # Total: O(1) average case
    return hash_table.get(key)  # Python dict uses hash table internally


# Demonstration
test_array = list(range(1000000))  # 1 million elements
print(f"First element of million-item array: {get_first_element(test_array)}")
print(f"Is 42 even? {is_even(42)}")
```

## O(log n) - Logarithmic Time

Each step eliminates half (or a constant fraction) of the remaining elements.

```
Time │
     │        ·········
     │   ·····
     │  ·
     │ ·
     │·
     └────────────────────────────────── n
```

```python
def binary_search(arr, target):
    """
    Find target in a sorted array using binary search.
    
    Time Complexity: O(log n)
        - Each iteration cuts the search space in half
        - For n elements: log₂(n) iterations maximum
        - Example: 1 million elements → ~20 comparisons max
    
    Space Complexity: O(1) - Only uses a few variables
    
    Args:
        arr: Sorted list of comparable elements
        target: Value to find
    Returns:
        Index of target, or -1 if not found
    """
    left, right = 0, len(arr) - 1
    
    # Loop runs O(log n) times
    # Why? Each iteration: search_space = search_space / 2
    # Starting with n, after k iterations: n / 2^k
    # Stops when n / 2^k = 1, solving: k = log₂(n)
    iterations = 0
    while left <= right:
        iterations += 1
        mid = (left + right) // 2  # O(1) - arithmetic
        
        if arr[mid] == target:      # O(1) - comparison
            print(f"Found in {iterations} iterations (log₂({len(arr)}) ≈ {math.log2(len(arr)):.1f})")
            return mid
        elif arr[mid] < target:
            left = mid + 1  # Eliminate left half
        else:
            right = mid - 1  # Eliminate right half
    
    print(f"Not found after {iterations} iterations")
    return -1


def count_digits(n):
    """
    Count the number of digits in a positive integer.
    
    Time Complexity: O(log n)
        - We divide by 10 each iteration
        - Number of digits ≈ log₁₀(n)
    
    Space Complexity: O(1)
    """
    if n == 0:
        return 1
    
    count = 0
    n = abs(n)
    
    # Each iteration removes one digit (divides by 10)
    # For n with d digits: d = ⌊log₁₀(n)⌋ + 1 iterations
    while n > 0:
        count += 1
        n //= 10  # Remove last digit
    
    return count


# Demonstration
import math
sorted_array = list(range(0, 1000000, 2))  # 500,000 even numbers
print("Binary search for 424242:")
binary_search(sorted_array, 424242)

print(f"\nDigits in 12345678: {count_digits(12345678)}")
```

## O(n) - Linear Time

Time grows directly proportional to input size.

```
Time │
     │                    ╱
     │                  ╱
     │                ╱
     │              ╱
     │            ╱
     │          ╱
     │        ╱
     │      ╱
     │    ╱
     │  ╱
     │╱
     └────────────────────────────────── n
```

```python
def find_maximum(arr):
    """
    Find the maximum element in an unsorted array.
    
    Time Complexity: O(n)
        - Must check every element (could be anywhere)
        - Exactly n comparisons for n elements
    
    Space Complexity: O(1) - Single variable for tracking max
    
    Args:
        arr: Non-empty list of comparable elements
    Returns:
        Maximum element
    """
    if not arr:
        return None
    
    maximum = arr[0]  # O(1)
    
    # Loop executes n-1 times (once per element after first)
    # Each iteration does O(1) work (comparison, possible assignment)
    # Total: O(n)
    for element in arr[1:]:
        if element > maximum:  # O(1) comparison
            maximum = element  # O(1) assignment
    
    return maximum


def linear_search(arr, target):
    """
    Search for target by checking each element.
    
    Time Complexity: O(n)
        - Best case: O(1) - target is first element
        - Worst case: O(n) - target is last or not present
        - Average case: O(n/2) = O(n) - target is in middle
    
    Space Complexity: O(1)
    """
    # In worst case, we examine all n elements
    for i, element in enumerate(arr):
        if element == target:  # O(1) comparison
            return i
    return -1


def sum_array(arr):
    """
    Calculate sum of all elements.
    
    Time Complexity: O(n) - Must visit each element exactly once
    Space Complexity: O(1) - Single accumulator variable
    """
    total = 0  # O(1) space
    
    # n iterations, O(1) work each
    for num in arr:
        total += num  # O(1)
    
    return total


# Demonstration
test_data = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5]
print(f"Maximum: {find_maximum(test_data)}")
print(f"Sum: {sum_array(test_data)}")
print(f"Index of 9: {linear_search(test_data, 9)}")
```
