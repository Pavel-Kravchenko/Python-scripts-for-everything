---
name: algo-comparison-sorts
description: Comparison-based sorting algorithms — bubble, merge, shell, quicksort — with complexity tables and implementation patterns.
tool_type: python
primary_tool: Python
---

# Comparison-Based Sorting Algorithms

## Complexity Table

| Algorithm | Best | Average | Worst | Space | Stable? |
|-----------|------|---------|-------|-------|---------|
| Bubble Sort | O(n) | O(n²) | O(n²) | O(1) | Yes |
| Merge Sort | O(n log n) | O(n log n) | O(n log n) | O(n) | Yes |
| Shell Sort | O(n log n) | O(n^1.3) | O(n²) | O(1) | No |
| QuickSort | O(n log n) | O(n log n) | O(n²) | O(log n) | No |

## Decision Table

| Scenario | Use |
|----------|-----|
| Small array (< 50 elements) | Bubble or Insertion |
| Need guaranteed O(n log n) | Merge Sort |
| General purpose, in-place | QuickSort |
| Memory constrained, medium n | Shell Sort |
| Need stability | Merge Sort |
| External / disk sorting | Merge Sort |
| Nearly sorted data | Insertion / Bubble (with early termination) |

## Implementations

### Bubble Sort — O(n²), O(1) space, stable
```python
def bubble_sort(arr: list[int]) -> list[int]:
    n = len(arr)
    for i in range(n - 1):
        swapped = False
        for j in range(n - 1 - i):
            if arr[j] > arr[j + 1]:
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
                swapped = True
        if not swapped:
            break  # already sorted — O(n) best case
    return arr
```

### Merge Sort — O(n log n), O(n) space, stable
```python
def merge_sort(arr: list[int]) -> list[int]:
    if len(arr) <= 1:
        return arr
    mid = len(arr) // 2
    left = merge_sort(arr[:mid])
    right = merge_sort(arr[mid:])
    return _merge(left, right)

def _merge(left: list[int], right: list[int]) -> list[int]:
    result, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        if left[i] <= right[j]:   # <= preserves stability
            result.append(left[i]); i += 1
        else:
            result.append(right[j]); j += 1
    result.extend(left[i:])
    result.extend(right[j:])
    return result
```

### Shell Sort — O(n^1.3) average, O(1) space
```python
def shell_sort(arr: list[int]) -> list[int]:
    n = len(arr)
    gap = n // 2
    while gap > 0:
        for i in range(gap, n):
            temp = arr[i]
            j = i
            while j >= gap and arr[j - gap] > temp:
                arr[j] = arr[j - gap]
                j -= gap
            arr[j] = temp
        gap //= 2
    return arr
```

### QuickSort — O(n log n) average, O(n²) worst
```python
def quicksort(arr: list[int], lo: int = 0, hi: int | None = None) -> list[int]:
    if hi is None:
        hi = len(arr) - 1
    if lo < hi:
        p = _partition(arr, lo, hi)
        quicksort(arr, lo, p - 1)
        quicksort(arr, p + 1, hi)
    return arr

def _partition(arr: list[int], lo: int, hi: int) -> int:
    pivot = arr[hi]
    i = lo - 1
    for j in range(lo, hi):
        if arr[j] <= pivot:
            i += 1
            arr[i], arr[j] = arr[j], arr[i]
    arr[i + 1], arr[hi] = arr[hi], arr[i + 1]
    return i + 1
```

## Pitfalls

- **QuickSort worst case on sorted input**: always O(n²) with last-element pivot; use random pivot or median-of-three for real data.
- **Merge Sort space cost**: O(n) auxiliary is unavoidable for the standard top-down version; use TimSort (Python's built-in) for production.
- **Stability matters for multi-key sorts**: if you sort records by key A then key B, an unstable sort can destroy the A ordering. Bubble and Merge are safe; Quick and Shell are not.
- **Python's `list.sort()` is TimSort**: O(n log n) worst case, stable, O(n) best case. Prefer it over rolling your own unless learning or the comparison is custom.
- **Shell sort gap sequence affects performance**: the naive `n//2` sequence gives O(n²) worst case; Knuth's `(3^k-1)/2` (1, 4, 13, 40…) gives O(n^1.3).
