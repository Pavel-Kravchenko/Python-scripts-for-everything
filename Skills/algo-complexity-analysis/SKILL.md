---
name: algo-complexity-analysis
description: Big O notation, complexity classes, space vs time trade-offs, and complexity analysis rules.
tool_type: python
primary_tool: Python
---

## Complexity Classes (fastest to slowest)

| Class | Name | Example |
|-------|------|---------|
| O(1) | Constant | Dict lookup, array index |
| O(log n) | Logarithmic | Binary search |
| O(n) | Linear | Linear scan, single pass |
| O(n log n) | Linearithmic | Merge sort, heap sort |
| O(n²) | Quadratic | Nested loops, bubble sort |
| O(2ⁿ) | Exponential | Brute-force subset enumeration |
| O(n!) | Factorial | Permutation enumeration |

## Simplification Rules

```
O(2n + 5)    →  O(n)      # drop constants
O(n² + n)    →  O(n²)     # drop lower-order terms
O(500)       →  O(1)
O(n² + n³)   →  O(n³)
```

## At Scale (n = 1,000,000)

| Complexity | Operations | Feasible? |
|------------|-----------|-----------|
| O(1) | 1 | Yes |
| O(log n) | ~20 | Yes |
| O(n) | 1,000,000 | Yes |
| O(n log n) | ~20,000,000 | Yes |
| O(n²) | 10¹² | No |
| O(2ⁿ) | ∞ | Never |

## Analysis Patterns

```python
# O(1): direct access
return arr[0]
return hash_table.get(key)

# O(log n): halving the search space each step
while left <= right:
    mid = (left + right) // 2
    ...

# O(n): single traversal
for item in arr:
    ...

# O(n log n): divide-and-conquer
# merge sort, heapq.nlargest(), sorted()

# O(n²): nested loops over same collection
for i in range(n):
    for j in range(n):   # or range(i, n) — still O(n²)
        ...

# O(n²) disguised: string concatenation in a loop
result = ""
for s in items:
    result += s    # creates new string each time → O(n²) total
# Fix: ''.join(items)  → O(n)
```

## Space Complexity

| Pattern | Space |
|---------|-------|
| Fixed variables | O(1) |
| Single copy of input | O(n) |
| Recursion depth d | O(d) stack frames |
| 2D DP table | O(n²) or O(n) with rolling array |

```python
# O(1) space: iterative with fixed variables
def find_max(arr):
    m = arr[0]
    for x in arr[1:]:
        m = max(m, x)
    return m

# O(n) space: storing results
def cumsum(arr):
    result = []          # grows with n
    total = 0
    for x in arr:
        total += x
        result.append(total)
    return result
```

## Amortized Complexity

- Python `list.append()`: O(1) amortized (occasional O(n) resize, but rare)
- Python `dict` lookup: O(1) average; O(n) worst case (all keys collide — rare with good hash)

## Best / Average / Worst

| Algorithm | Best | Average | Worst |
|-----------|------|---------|-------|
| Binary search | O(1) | O(log n) | O(log n) |
| Quicksort | O(n log n) | O(n log n) | O(n²) |
| Merge sort | O(n log n) | O(n log n) | O(n log n) |
| Hash table lookup | O(1) | O(1) | O(n) |
| BFS/DFS | O(V+E) | O(V+E) | O(V+E) |

## Pitfalls

- **Hidden O(n) inside a loop**: `in` on a list is O(n); inside an O(n) loop = O(n²). Use a `set` for O(1) membership.
- **String concatenation**: `s += x` in a loop is O(n²) total. Use `''.join(parts)` or `io.StringIO`.
- **`sorted()` is O(n log n)**: calling it inside a loop makes the loop O(n² log n).
- **Recursion depth**: unbounded recursion on large n hits Python's default 1000-frame limit. Use iterative approach or `sys.setrecursionlimit`.
- **Space vs time trade-off**: memoization trades O(n) space for O(n) → O(1) repeated lookups.
- **Worst-case vs average-case**: quicksort is O(n²) on sorted input; always use random pivot or `timsort` for general data.
