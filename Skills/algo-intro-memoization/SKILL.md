---
name: algo-intro-memoization
description: "1. Recognize problems that can be solved with dynamic programming 2. Understand the difference between recursion and memoization 3. Implement memoization using Python decorators and dictionaries 4. Ap"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/01_intro_memoization.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# 🔄 Introduction to Dynamic Programming: Memoization

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/01_intro_memoization.ipynb`*


## Learning Objectives

By the end of this notebook, you will be able to:

1. Recognize problems that can be solved with dynamic programming
2. Understand the difference between recursion and memoization
3. Implement memoization using Python decorators and dictionaries
4. Apply DP to sequence problems relevant to bioinformatics

---

## Why Dynamic Programming in Bioinformatics?

Dynamic programming is the **foundation** of many core algorithms:

| Algorithm | Application | DP Type |
|-----------|-------------|----------|
| Needleman-Wunsch | Global sequence alignment | 2D table |
| Smith-Waterman | Local sequence alignment | 2D table |
| Nussinov | RNA secondary structure | 2D table |
| Viterbi | HMM state decoding | 1D with states |
| Gene finding | ORF prediction | Multi-state |

Understanding DP fundamentals will help you understand these algorithms deeply.

## 1. The Problem with Pure Recursion

Let's start with the classic Fibonacci sequence:

```python
F(0) = 0
F(1) = 1
F(n) = F(n-1) + F(n-2) for n > 1
```python

```python
import time

def fib_recursive(n):
    """Pure recursive Fibonacci - O(2^n) time complexity!"""
    if n <= 1:
        return n
    return fib_recursive(n - 1) + fib_recursive(n - 2)

# Test small values
for i in range(10):
    print(f"F({i}) = {fib_recursive(i)}")
```python

```python
# Let's see how long it takes for larger values
def time_function(func, n):
    start = time.time()
    result = func(n)
    elapsed = time.time() - start
    return result, elapsed

for n in [20, 25, 30, 35]:
    result, elapsed = time_function(fib_recursive, n)
    print(f"F({n}) = {result}, took {elapsed:.3f} seconds")
```python

### Why is it so slow?

The recursion tree has **overlapping subproblems**:

```python
                    F(5)
                   /    \
                F(4)     F(3)       ← F(3) computed twice!
               /   \     /   \
            F(3)  F(2)  F(2)  F(1)  ← F(2) computed 3 times!
           /  \   
        F(2)  F(1)
```python

For F(n), we make ~2^n function calls, even though there are only n unique subproblems.

## 2. Memoization: Remember What We've Computed

**Memoization** stores results of expensive function calls and returns the cached result when the same inputs occur again.

```python
def fib_memoized(n, memo=None):
    """Fibonacci with memoization - O(n) time, O(n) space."""
    if memo is None:
        memo = {}
    
    # Base cases
    if n <= 1:
        return n
    
    # Check if already computed
    if n in memo:
        return memo[n]
    
    # Compute and store
    memo[n] = fib_memoized(n - 1, memo) + fib_memoized(n - 2, memo)
    return memo[n]

# Now it's fast!
for n in [20, 25, 30, 35, 100]:
    result, elapsed = time_function(fib_memoized, n)
    print(f"F({n}) = {result}, took {elapsed:.6f} seconds")
```python

### Using Python's `@functools.lru_cache`

Python provides a built-in decorator for memoization:

```python
from functools import lru_cache

@lru_cache(maxsize=None)  # Unlimited cache
def fib_cached(n):
    """Fibonacci with automatic memoization."""
    if n <= 1:
        return n
    return fib_cached(n - 1) + fib_cached(n - 2)

# Clear cache for fair timing
fib_cached.cache_clear()

result, elapsed = time_function(fib_cached, 100)
print(f"F(100) = {result}")
print(f"Time: {elapsed:.6f} seconds")
print(f"Cache info: {fib_cached.cache_info()}")
```python

## 3. The DP Mindset: Identify Subproblems

To apply dynamic programming, ask:

1. **Can I define the problem recursively?** (recurrence relation)
2. **Are there overlapping subproblems?** (same subproblems solved multiple times)
3. **What are the base cases?** (smallest solvable problems)

If yes to all three → DP is likely applicable.

## 4. Example: Climbing Stairs (Counting Paths)

You're climbing a staircase with `n` steps. Each time you can climb 1 or 2 steps. How many distinct ways can you reach the top?

This models many biological scenarios:
- Counting possible RNA folding paths
- Enumerating assembly paths

```python
@lru_cache(maxsize=None)
def climb_stairs(n):
    """
    Count ways to climb n stairs taking 1 or 2 steps at a time.
    
    Recurrence: ways(n) = ways(n-1) + ways(n-2)
    - ways(n-1): last step was 1 stair
    - ways(n-2): last step was 2 stairs
    """
    if n <= 2:
        return n
    return climb_stairs(n - 1) + climb_stairs(n - 2)

for n in range(1, 11):
    print(f"Stairs={n}: {climb_stairs(n)} ways")
```python

## 5. 🧬 Bioinformatics Example: Counting Alignments

Given two sequences of lengths `m` and `n`, how many possible alignments exist?

At each position, we can:
- Match/mismatch (consume both sequences)
- Gap in sequence 1 (consume only sequence 2)
- Gap in sequence 2 (consume only sequence 1)

```python
@lru_cache(maxsize=None)
def count_alignments(m, n):
    """
    Count possible alignments between sequences of length m and n.
    
    This is the number of paths through an (m+1) x (n+1) alignment matrix.
    """
    # Base cases: one sequence exhausted
    if m == 0 or n == 0:
        return 1  # Only gap insertions possible
    
    # Recurrence: sum of three moves
    return (count_alignments(m - 1, n - 1) +  # Match/mismatch
            count_alignments(m - 1, n) +       # Gap in seq2
            count_alignments(m, n - 1))        # Gap in seq1

# Count alignments for small sequences
for m, n in [(3, 3), (5, 5), (10, 10)]:
    count_alignments.cache_clear()
    count = count_alignments(m, n)
    print(f"Sequences of length {m} and {n}: {count:,} possible alignments")
```python

```python
# The growth is exponential - this is why we need scoring!
import matplotlib.pyplot as plt

lengths = list(range(1, 15))
counts = []
for l in lengths:
    count_alignments.cache_clear()
    counts.append(count_alignments(l, l))

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(lengths, counts, 'bo-')
plt.xlabel('Sequence Length')
plt.ylabel('Number of Alignments')
plt.title('Alignment Count (Linear Scale)')

plt.subplot(1, 2, 2)
plt.semilogy(lengths, counts, 'ro-')
plt.xlabel('Sequence Length')
plt.ylabel('Number of Alignments (log)')
plt.title('Alignment Count (Log Scale)')

plt.tight_layout()
plt.show()
```python

## 6. Example: Coin Change (Minimum Operations)

Given coin denominations and a target amount, find the minimum number of coins needed.

This is analogous to:
- Finding minimum edit distance
- Optimal resource allocation

```python
def min_coins(coins, amount, memo=None):
    """
    Find minimum coins needed to make the amount.
    Returns -1 if impossible.
    """
    if memo is None:
        memo = {}
    
    # Base cases
    if amount == 0:
        return 0
    if amount < 0:
        return float('inf')
    
    if amount in memo:
        return memo[amount]
    
    # Try each coin and take minimum
    min_count = float('inf')
    for coin in coins:
        result = min_coins(coins, amount - coin, memo)
        if result != float('inf'):
            min_count = min(min_count, result + 1)
    
    memo[amount] = min_count
    return min_count

# Example: US coins
coins = [1, 5, 10, 25]
for amount in [11, 30, 41, 63]:
    result = min_coins(coins, amount)
    print(f"Amount ${amount/100:.2f}: {result} coins minimum")
```python

## 🧬 Exercise 1: RNA Folding Count

In RNA secondary structure, bases can pair (A-U, G-C). Given an RNA sequence length `n`, count the number of possible non-crossing base pairing configurations.

Simplified rules:
- Position i can pair with position j if j > i + 3 (minimum loop size)
- Pairings cannot cross

This is related to the Catalan numbers.

Implement `count_rna_structures(n)` that returns the count of valid structures.

```python
# Hint: The recurrence involves considering whether position n pairs
# with some earlier position k, and recursively counting structures
# for the resulting subproblems.

@lru_cache(maxsize=None)
def count_rna_structures(n):
    """Count non-crossing RNA secondary structures of length n."""
    # Your solution here
    pass

# Test (these are Catalan-like numbers)
# for n in range(1, 15):
#     print(f"n={n}: {count_rna_structures(n)} structures")
```python

## 🧬 Exercise 2: Longest Increasing Subsequence

Given a sequence of numbers (e.g., gene expression values over time), find the length of the longest increasing subsequence.

Example: `[10, 22, 9, 33, 21, 50, 41, 60, 80]`
LIS: `[10, 22, 33, 50, 60, 80]` → length 6

Implement using memoization.

```python
def longest_increasing_subsequence(arr):
    """Find length of longest increasing subsequence."""
    # Your solution here
    pass

# Test
# test_arr = [10, 22, 9, 33, 21, 50, 41, 60, 80]
# print(f"LIS length: {longest_increasing_subsequence(test_arr)}")
```python

---

## Summary

In this notebook, you learned:

✅ Why pure recursion can be exponentially slow  
✅ How memoization trades space for time  
✅ Python's `@lru_cache` decorator for automatic memoization  
✅ The DP mindset: identify subproblems, recurrence, base cases  
✅ Applications to sequence alignment counting  

**Next:** [02_tabulation.ipynb](02_tabulation.ipynb) - Bottom-up DP with tabulation

---

## Key Takeaway

> **Memoization** is top-down DP: start with the big problem, recursively break it down, cache results. It's the easiest way to add DP to a recursive solution — just add `@lru_cache`!

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
