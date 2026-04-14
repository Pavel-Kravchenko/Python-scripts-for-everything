---
name: algo-kmp-algorithm
description: "KMP string matching — O(n+m) pattern search using failure function prefix table"
tool_type: python
primary_tool: Python
---

# Knuth-Morris-Pratt (KMP) Algorithm

## Key Insight

After a partial match fails, the prefix function tells us how far to shift the pattern without re-comparing characters we already know match. The text pointer never moves backward.

## Prefix Function (Failure Function)

pi[i] = length of the longest proper prefix of P[0..i] that is also a suffix.

Example: Pattern "ABABACA" -> pi = [0, 0, 1, 2, 3, 0, 1]

## Build Prefix Table — O(n)

```python
def build_prefix(pattern: str) -> list[int]:
    """O(n) prefix function using fallback on previously computed values."""
    n = len(pattern)
    if n == 0:
        return []
    pi = [0] * n
    j = 0
    for i in range(1, n):
        while j > 0 and pattern[j] != pattern[i]:
            j = pi[j - 1]
        if pattern[j] == pattern[i]:
            j += 1
        pi[i] = j
    return pi
```

Amortized O(n): j increases at most n times total, so total decreases <= n.

## KMP Search — O(n + m)

```python
def kmp_search(text: str, pattern: str) -> list[int]:
    """Find all occurrences of pattern in text. O(n + m) time, O(m) space."""
    if not pattern or len(pattern) > len(text):
        return []
    pi = build_prefix(pattern)
    matches = []
    j = 0
    for i in range(len(text)):
        while j > 0 and text[i] != pattern[j]:
            j = pi[j - 1]
        if text[i] == pattern[j]:
            j += 1
        if j == len(pattern):
            matches.append(i - j + 1)
            j = pi[j - 1]  # look for next (possibly overlapping) match
    return matches
```

## Complexity

| | Time | Space |
|---|---|---|
| Preprocessing | O(m) | O(m) |
| Search | O(n) | O(1) |
| **Total** | **O(n + m)** | **O(m)** |

## Pitfalls

- The prefix table build looks O(n^2) due to nested while-in-for, but the amortized argument (j can only decrease as many times as it increased) guarantees O(n)
- KMP finds overlapping matches by default — `j = pi[j-1]` after a full match continues searching
- For multiple patterns of different lengths, Aho-Corasick is more appropriate than running KMP per pattern
