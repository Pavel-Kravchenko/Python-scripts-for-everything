---
name: complexity-sorting-searching
description: Big-O analysis, comparison and linear sorting algorithms, binary search, and complexity reference tables
---

# Complexity Analysis, Sorting & Searching

## When to Use

- Need to choose between algorithms → check complexity class first
- Input size > 10k and O(n²) in inner loop → switch to O(n log n)
- Values in bounded integer range → counting sort beats comparison sort
- Sorted array + need fast lookup → binary search (O(log n) vs O(n))
- BAM/BED coordinate queries → binary search on sorted position array
- BLAST hits ranking → sort by e-value (Python `sorted` = Timsort)
- k-mer counting → counting/radix sort on fixed-length strings
- Genomic interval overlap → sort by start, binary search for query position

---

## Quick Reference

### Big-O, Big-Ω, Big-Θ

| Notation | Name | Meaning |
|---|---|---|
| O(f(n)) | Big-O | Upper bound (worst case) |
| Ω(f(n)) | Big-Omega | Lower bound (best case) |
| Θ(f(n)) | Big-Theta | Tight bound (both upper and lower) |

**Growth rate hierarchy (slowest → fastest):**
O(1) < O(log n) < O(n) < O(n log n) < O(n²) < O(n³) < O(2ⁿ) < O(n!)

**Rules:** Drop constants and lower-order terms. O(3n² + 5n + 2) = O(n²).

### Complexity Classes at a Glance

| Complexity | Name | n=10 | n=100 | n=1000 | Example |
|---|---|---|---|---|---|
| O(1) | Constant | 1 | 1 | 1 | Array access, hash lookup |
| O(log n) | Logarithmic | 3 | 7 | 10 | Binary search |
| O(n) | Linear | 10 | 100 | 1,000 | Linear search, single loop |
| O(n log n) | Linearithmic | 33 | 664 | 9,966 | Merge sort, heap sort |
| O(n²) | Quadratic | 100 | 10,000 | 1,000,000 | Bubble/insertion sort (worst) |
| O(2ⁿ) | Exponential | 1,024 | (huge) | (∞) | Naive subset enumeration |

### Master Theorem (simplified)

For T(n) = aT(n/b) + f(n):
- f(n) = O(n^(log_b a - ε)) → T(n) = Θ(n^(log_b a))  [recursion-heavy]
- f(n) = Θ(n^(log_b a)) → T(n) = Θ(n^(log_b a) log n)  [balanced]
- f(n) = Ω(n^(log_b a + ε)) → T(n) = Θ(f(n))  [work-heavy]

### Sorting Algorithm Comparison

| Algorithm | Best | Average | Worst | Space | Stable | Notes |
|---|---|---|---|---|---|---|
| Bubble Sort | O(n) | O(n²) | O(n²) | O(1) | Yes | Early exit if no swaps |
| Selection Sort | O(n²) | O(n²) | O(n²) | O(1) | No | Min swaps (n-1) |
| Insertion Sort | O(n) | O(n²) | O(n²) | O(1) | Yes | Fast on small/nearly-sorted |
| Shell Sort | O(n log n) | gap-dependent | O(n²) | O(1) | No | Gap sequence matters |
| Merge Sort | O(n log n) | O(n log n) | O(n log n) | O(n) | Yes | Guaranteed; BAM coordinate sort |
| QuickSort | O(n log n) | O(n log n) | O(n²) | O(log n) | No | Worst on sorted input |
| Heap Sort | O(n log n) | O(n log n) | O(n log n) | O(1) | No | In-place; cache-unfriendly |
| Counting Sort | O(n+k) | O(n+k) | O(n+k) | O(k) | Yes | Integers only; k = value range |
| Radix Sort | O(d(n+k)) | O(d(n+k)) | O(d(n+k)) | O(n+k) | Yes | d = digits; k = base |
| Bucket Sort | O(n+k) | O(n+k) | O(n²) | O(n) | Yes | Uniform distribution assumed |

### Binary Search Complexity

| Variant | Time | Space |
|---|---|---|
| Standard (iterative) | O(log n) | O(1) |
| Standard (recursive) | O(log n) | O(log n) |
| Lower/upper bound | O(log n) | O(1) |
| Rotated array | O(log n) | O(1) |

---

## Key Patterns

### Analyzing Loops
```python
# Single loop → O(n)
# Nested loop → O(n²)
# Loop halving → O(log n)
# Loop + recursive halve → O(n log n)  [merge sort pattern]
for i in range(n):          # O(n)
    j = n
    while j > 1:            # O(log n) inner
        j //= 2
# Total: O(n log n)
```

### When Linear Sorts Beat O(n log n)
- Values are non-negative integers with bounded range k → **counting sort**: O(n+k)
- Fixed-length strings/integers → **radix sort**: O(d·n) where d = key length
- Floats uniformly distributed in [0,1) → **bucket sort**: O(n) average
- Comparison-based sorting has a lower bound of Ω(n log n); linear sorts bypass this by not comparing

---

## Code Templates

### Bubble Sort — O(n²) / O(n) best
```python
def bubble_sort(arr):
    n = len(arr)
    for i in range(n - 1):
        swapped = False
        for j in range(n - 1 - i):
            if arr[j] > arr[j + 1]:
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
                swapped = True
        if not swapped:
            break  # already sorted
    return arr
```

### Insertion Sort — O(n²) / O(n) best
```python
def insertion_sort(arr):
    for i in range(1, len(arr)):
        key = arr[i]
        j = i - 1
        while j >= 0 and arr[j] > key:
            arr[j + 1] = arr[j]
            j -= 1
        arr[j + 1] = key
    return arr
```

### Merge Sort — O(n log n) all cases
```python
def merge_sort(arr):
    if len(arr) <= 1:
        return arr
    mid = len(arr) // 2
    left = merge_sort(arr[:mid])
    right = merge_sort(arr[mid:])
    return _merge(left, right)

def _merge(left, right):
    result, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        if left[i] <= right[j]:  # <= preserves stability
            result.append(left[i]); i += 1
        else:
            result.append(right[j]); j += 1
    return result + left[i:] + right[j:]
```

### QuickSort — O(n log n) avg, O(n²) worst
```python
def quicksort(arr, low=0, high=None):
    if high is None: high = len(arr) - 1
    if low < high:
        pivot_idx = _partition(arr, low, high)
        quicksort(arr, low, pivot_idx - 1)
        quicksort(arr, pivot_idx + 1, high)
    return arr

def _partition(arr, low, high):
    pivot, i = arr[high], low
    for j in range(low, high):
        if arr[j] < pivot:
            arr[i], arr[j] = arr[j], arr[i]; i += 1
    arr[i], arr[high] = arr[high], arr[i]
    return i
```

### Heap Sort — O(n log n) in-place
```python
import heapq
def heap_sort(arr):
    heapq.heapify(arr)                      # O(n) build
    return [heapq.heappop(arr) for _ in range(len(arr))]  # O(n log n)
```

### Counting Sort — O(n+k), stable
```python
def counting_sort(arr, max_val=None):
    if not arr: return []
    k = max_val or max(arr)
    count = [0] * (k + 1)
    for num in arr: count[num] += 1
    for i in range(1, k + 1): count[i] += count[i - 1]  # prefix sum
    output = [0] * len(arr)
    for num in reversed(arr):      # reversed → stable
        count[num] -= 1
        output[count[num]] = num
    return output
```

### Radix Sort — O(d·(n+k)), stable
```python
def radix_sort(arr):
    if not arr: return []
    max_num = max(arr)
    result = arr.copy()
    exp = 1
    while max_num // exp > 0:
        result = _counting_by_digit(result, exp)
        exp *= 10
    return result

def _counting_by_digit(arr, exp):
    count, output = [0]*10, [0]*len(arr)
    for num in arr: count[(num // exp) % 10] += 1
    for i in range(1, 10): count[i] += count[i-1]
    for num in reversed(arr):
        d = (num // exp) % 10; count[d] -= 1; output[count[d]] = num
    return output
```

### Binary Search — standard template
```python
def binary_search(arr, target):
    left, right = 0, len(arr) - 1
    while left <= right:
        mid = left + (right - left) // 2   # overflow-safe
        if arr[mid] == target:   return mid
        elif arr[mid] < target:  left = mid + 1
        else:                    right = mid - 1
    return -1
```

### Lower Bound (first index >= target) — bisect_left equivalent
```python
def lower_bound(arr, target):
    left, right = 0, len(arr)
    while left < right:
        mid = (left + right) // 2
        if arr[mid] < target: left = mid + 1
        else: right = mid
    return left
```

### Upper Bound (first index > target) — bisect_right equivalent
```python
def upper_bound(arr, target):
    left, right = 0, len(arr)
    while left < right:
        mid = (left + right) // 2
        if arr[mid] <= target: left = mid + 1
        else: right = mid
    return left
```

### First / Last Occurrence
```python
def find_first(arr, target):
    left, right, result = 0, len(arr)-1, -1
    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == target: result = mid; right = mid - 1  # keep going left
        elif arr[mid] < target: left = mid + 1
        else: right = mid - 1
    return result

def find_last(arr, target):
    left, right, result = 0, len(arr)-1, -1
    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == target: result = mid; left = mid + 1   # keep going right
        elif arr[mid] < target: left = mid + 1
        else: right = mid - 1
    return result

def count_occurrences(arr, target):
    first = find_first(arr, target)
    return 0 if first == -1 else find_last(arr, target) - first + 1
```

### Search in Rotated Sorted Array
```python
def search_rotated(arr, target):
    left, right = 0, len(arr) - 1
    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == target: return mid
        if arr[left] <= arr[mid]:  # left half sorted
            if arr[left] <= target < arr[mid]: right = mid - 1
            else: left = mid + 1
        else:                      # right half sorted
            if arr[mid] < target <= arr[right]: left = mid + 1
            else: right = mid - 1
    return -1
```

### Python stdlib shortcuts
```python
import bisect
bisect.bisect_left(arr, x)   # lower_bound
bisect.bisect_right(arr, x)  # upper_bound
sorted(iterable, key=fn)     # Timsort O(n log n), stable
arr.sort(key=fn, reverse=True)
# Sort BLAST hits by e-value:
hits.sort(key=lambda h: h['evalue'])
# Sort genomic records by (chrom, position):
records.sort(key=lambda r: (r['chrom'], r['start']))
```

---

## Common Pitfalls

- **`mid = (left + right) // 2`** overflows in C/Java; use `left + (right - left) // 2`
- **`while left < right` vs `<=`**: use `<=` for exact match; use `<` for boundary (lower/upper bound)
- **`left = mid` instead of `mid + 1`** causes infinite loop when `left == mid`
- **QuickSort on sorted input** hits O(n²); use random pivot or 3-way partition
- **Counting sort requires non-negative integers**; can offset for negatives
- **Radix sort with large k base** uses extra memory; base 10 is safe default
- **Bucket sort degrades to O(n²)** if distribution is not uniform (all fall in one bucket)
- **Merge sort on linked lists** is O(n log n) with O(1) extra space; on arrays needs O(n)
- **Stability matters** when sorting by secondary key after primary (e.g., chrom then position)

---

## Bioinformatics Connections

| Problem | Algorithm | Why |
|---|---|---|
| BAM coordinate sorting | Merge sort | Stable, guaranteed O(n log n), used by samtools |
| BLAST hit ranking | `sorted()` / Timsort | Sort by e-value ascending; stable |
| k-mer lexicographic order | Radix sort | Fixed-length strings; O(k·n) beats O(n log n) for short k |
| DNA base frequency sort | Counting sort | Alphabet size = 4 (A,C,G,T); O(n) |
| Gene record stable sort | Timsort / merge sort | Preserve original order for ties |
| VCF/BED interval query | Binary search on sorted positions | Find first overlapping interval: `lower_bound(starts, query)` |
| Prefix sum on coverage | Cumulative count array | Counting sort step; O(n) range queries |

---

## Related Skills

- `linear-tree-hash-structures` — heap/priority queue internals, BST O(h) operations
- `graphs-dynamic-programming` — graph traversal complexity (BFS/DFS O(V+E))
- `string-algorithms` — KMP O(n), suffix array O(n log n) for pattern search
- `sequence-alignment` — DP O(nm) Smith-Waterman, complexity trade-offs
