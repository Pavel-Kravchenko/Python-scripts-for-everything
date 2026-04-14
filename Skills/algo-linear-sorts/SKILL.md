---
name: algo-linear-sorts
description: "Counting sort, radix sort, bucket sort — O(n) non-comparison sorts that bypass the Omega(n log n) lower bound"
tool_type: python
primary_tool: Python
---

# Linear Time Sorting Algorithms

## Why O(n) is Possible

Comparison sorts need Omega(n log n) (decision tree with n! leaves). Linear sorts bypass this by not comparing elements — they use direct addressing, digit decomposition, or distribution.

**Trade-off:** O(n) time but requires assumptions about input (integers, bounded range, etc.) and extra space.

## Counting Sort — O(n + k)

k = range of values. Stable. Used as subroutine in radix sort.

```python
def counting_sort(arr, max_val=None):
    if not arr:
        return []
    k = max_val if max_val is not None else max(arr)
    count = [0] * (k + 1)
    for num in arr:
        count[num] += 1
    for i in range(1, k + 1):
        count[i] += count[i - 1]  # prefix sum
    output = [0] * len(arr)
    for num in reversed(arr):  # right-to-left for stability
        count[num] -= 1
        output[count[num]] = num
    return output
```

**When to use:** Integers with small known range (ages 0-120, grades 0-100). If k >> n, space becomes O(k) and it's inefficient.

## Radix Sort — O(d(n + k))

Sort by individual digits from least significant to most significant, using a stable sort (counting sort) per digit. d = number of digits, k = base (10 for decimal).

```python
def counting_sort_by_digit(arr, exp):
    output = [0] * len(arr)
    count = [0] * 10
    for num in arr:
        digit = (num // exp) % 10
        count[digit] += 1
    for i in range(1, 10):
        count[i] += count[i - 1]
    for num in reversed(arr):
        digit = (num // exp) % 10
        count[digit] -= 1
        output[count[digit]] = num
    return output

def radix_sort(arr):
    if not arr:
        return []
    max_num = max(arr)
    result = arr.copy()
    exp = 1
    while max_num // exp > 0:
        result = counting_sort_by_digit(result, exp)
        exp *= 10
    return result
```

**When to use:** Large numbers with many digits, fixed-length strings, when d x (n+k) < n log n.

## When to Choose Which

| Algorithm | Time | Space | Stable | Input Requirement |
|---|---|---|---|---|
| Counting Sort | O(n + k) | O(n + k) | Yes | Non-negative integers, small range |
| Radix Sort | O(d(n + k)) | O(n + k) | Yes | Integers or fixed-length strings |
| Bucket Sort | O(n) avg | O(n + k) | Yes | Uniformly distributed values |

## Pitfalls

- Counting sort with negative integers requires offset: shift all values by `abs(min)`
- Radix sort on negative numbers: sort absolute values, then reverse negatives and prepend
- Right-to-left iteration in counting sort is essential for stability — left-to-right breaks stable ordering
