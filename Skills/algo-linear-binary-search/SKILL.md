---
name: algo-linear-binary-search
description: "Linear and binary search — implementations, variants (first/last occurrence), and when to use each"
tool_type: python
primary_tool: Python
---

# Linear and Binary Search

## When to Use Which

| Use Linear Search | Use Binary Search |
|---|---|
| Unsorted data | Sorted data |
| Small arrays (< ~20 elements) | Large arrays |
| Single search (sorting + binary is overkill) | Multiple searches (amortize sort cost) |
| Linked lists (no random access) | Random-access structures (arrays) |

## Complexity

| | Linear | Binary |
|---|---|---|
| Best | O(1) | O(1) |
| Average/Worst | O(n) | O(log n) |
| Space (iterative) | O(1) | O(1) |
| Prerequisite | None | Sorted array |

## Linear Search

```python
def linear_search(arr, target):
    for i in range(len(arr)):
        if arr[i] == target:
            return i
    return -1
```

### Sentinel Optimization

Eliminates bound checking in the inner loop:

```python
def linear_search_with_sentinel(arr, target):
    if not arr:
        return -1
    last = arr[-1]
    arr[-1] = target
    i = 0
    while arr[i] != target:
        i += 1
    arr[-1] = last
    if i < len(arr) - 1 or arr[-1] == target:
        return i
    return -1
```

## Binary Search (Iterative)

```python
def binary_search(arr, target):
    left, right = 0, len(arr) - 1
    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == target:
            return mid
        elif arr[mid] < target:
            left = mid + 1
        else:
            right = mid - 1
    return -1
```

## Find First Occurrence

When array has duplicates, find the leftmost position:

```python
def find_first_occurrence(arr, target):
    left, right = 0, len(arr) - 1
    result = -1
    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == target:
            result = mid
            right = mid - 1  # keep searching left
        elif arr[mid] < target:
            left = mid + 1
        else:
            right = mid - 1
    return result
```

## Pitfalls

- Integer overflow in `(left + right) // 2` — in Python this is safe (arbitrary precision), but in C/Java use `left + (right - left) // 2`
- Off-by-one errors: `left <= right` (inclusive bounds) vs `left < right` (exclusive right) require different update rules
- Binary search on unsorted data silently gives wrong results — no error, just incorrect
- `bisect` module in stdlib: `bisect_left` = lower bound, `bisect_right` = upper bound
