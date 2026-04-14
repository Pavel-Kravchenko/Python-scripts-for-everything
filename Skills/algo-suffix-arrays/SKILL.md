---
name: algo-suffix-arrays
description: "Suffix arrays — space-efficient sorted suffix index, O(n log n) construction, O(m log n) pattern search"
tool_type: python
primary_tool: Python
---

# Suffix Arrays

SA[i] = starting position of the i-th lexicographically smallest suffix. Space-efficient alternative to suffix trees (just an int array).

## Example: "banana$"

```
SA = [6, 5, 3, 1, 0, 4, 2]
  SA[0]=6 -> "$"
  SA[1]=5 -> "a$"
  SA[2]=3 -> "ana$"
  SA[3]=1 -> "anana$"
  SA[4]=0 -> "banana$"
  SA[5]=4 -> "na$"
  SA[6]=2 -> "nana$"
```

## Naive Construction — O(n^2 log n)

```python
def build_suffix_array_naive(text):
    """Sort all suffixes. O(n^2 log n) due to O(n) string comparisons."""
    n = len(text)
    suffixes = [(text[i:], i) for i in range(n)]
    suffixes.sort(key=lambda x: x[0])
    return [pos for _, pos in suffixes]
```

## Manber-Myers — O(n log n)

Sort by doubling prefix length. Each round uses ranks from previous round as O(1) comparison keys.

```python
def build_suffix_array_manber_myers(text):
    """O(n log n) via iterative doubling of comparison prefix length."""
    n = len(text)
    if n <= 1: return list(range(n))
    sa = list(range(n))
    rank = [ord(c) for c in text]
    tmp = [0] * n
    k = 1
    while k < n:
        def key(i): return (rank[i], rank[i + k] if i + k < n else -1)
        sa.sort(key=key)
        tmp[sa[0]] = 0
        for i in range(1, n):
            tmp[sa[i]] = tmp[sa[i-1]] + (1 if key(sa[i]) > key(sa[i-1]) else 0)
        rank, tmp = tmp, rank
        if rank[sa[-1]] == n - 1: break  # all ranks unique
        k *= 2
    return sa
```

## Pattern Search — O(m log n)

Binary search for the pattern as a prefix of sorted suffixes.

```
Find "ana" in SA of "banana$":
  Binary search finds range [2, 4) in SA
  -> SA[2]=3, SA[3]=1 -> positions 1 and 3
```

## Construction Complexity

| Method | Time | Space |
|---|---|---|
| Naive | O(n^2 log n) | O(n^2) |
| Manber-Myers | O(n log n) | O(n) |
| SA-IS / DC3 | O(n) | O(n) |

## Pitfalls

- Always append a terminal character ('$') smaller than any character in the alphabet — ensures no suffix is a prefix of another
- Naive construction stores all suffix strings: O(n^2) space. Manber-Myers only uses rank arrays: O(n) space
- For LCP (longest common prefix) array construction, use Kasai's algorithm in O(n) after building SA
