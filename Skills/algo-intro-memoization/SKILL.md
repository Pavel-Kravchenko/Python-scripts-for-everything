---
name: algo-intro-memoization
description: "Top-down DP with memoization — cache recursive subproblems for exponential-to-linear speedup, @lru_cache shortcut"
tool_type: python
primary_tool: Python
---

# Memoization (Top-Down DP)

## Core Idea

Memoization caches results of overlapping subproblems. Pure recursion on Fibonacci is O(2^n); with memoization it becomes O(n) time, O(n) space.

## DP Applicability Checklist

1. Recursive substructure (recurrence relation exists)
2. Overlapping subproblems (same subproblem solved multiple times)
3. Base cases defined

## Manual Memoization

```python
def fib_memoized(n, memo=None):
    if memo is None:
        memo = {}
    if n <= 1:
        return n
    if n in memo:
        return memo[n]
    memo[n] = fib_memoized(n - 1, memo) + fib_memoized(n - 2, memo)
    return memo[n]
```

## @lru_cache (Preferred)

```python
from functools import lru_cache

@lru_cache(maxsize=None)
def fib_cached(n):
    if n <= 1:
        return n
    return fib_cached(n - 1) + fib_cached(n - 2)

fib_cached.cache_clear()  # reset before timing
print(fib_cached.cache_info())  # hits, misses, maxsize, currsize
```

## Counting Sequence Alignments

```python
@lru_cache(maxsize=None)
def count_alignments(m, n):
    """Count possible alignments between sequences of length m and n."""
    if m == 0 or n == 0:
        return 1
    return (count_alignments(m - 1, n - 1) +  # match/mismatch
            count_alignments(m - 1, n) +       # gap in seq2
            count_alignments(m, n - 1))        # gap in seq1
```

## Coin Change (Minimum Operations)

```python
def min_coins(coins, amount, memo=None):
    if memo is None:
        memo = {}
    if amount == 0:
        return 0
    if amount < 0:
        return float('inf')
    if amount in memo:
        return memo[amount]
    min_count = float('inf')
    for coin in coins:
        result = min_coins(coins, amount - coin, memo)
        if result != float('inf'):
            min_count = min(min_count, result + 1)
    memo[amount] = min_count
    return min_count
```

## Key Takeaway

Memoization = top-down DP. Start with the big problem, recurse down, cache results. Easiest upgrade: add `@lru_cache` to any pure recursive function.

## Pitfalls

- Python's default recursion limit is 1000 — use `sys.setrecursionlimit()` or switch to tabulation for deep recursion
- Mutable default arguments (`memo={}`) persist across calls — use `memo=None` pattern
- `@lru_cache` requires hashable arguments — won't work with list/dict params
