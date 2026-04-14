---
name: algo-naive-pattern-matching
description: "Naive O(nm) brute-force string matching — sliding window baseline, when it's acceptable, motivation for KMP/Rabin-Karp"
tool_type: python
primary_tool: Python
---

# Naive Pattern Matching

## Algorithm

Slide pattern over text, comparing character by character at each position.

```python
def naive_search(text: str, pattern: str) -> list[int]:
    n, m = len(text), len(pattern)
    positions = []
    if m == 0 or m > n:
        return positions
    for i in range(n - m + 1):
        match = True
        for j in range(m):
            if text[i + j] != pattern[j]:
                match = False
                break
        if match:
            positions.append(i)
    return positions
```

## Complexity

| Case | Time | When |
|---|---|---|
| Best | O(n) | First char mismatches quickly |
| Worst | O(n x m) | Many partial matches (e.g., T="AAA...B", P="AA...B") |
| Average | O(n) for large alphabets | Random text, mismatches found early |
| Space | O(1) auxiliary | |

## When Naive is Fine

- **Small patterns:** m is small (e.g., word in a document), so O(nm) ~ O(n)
- **Large alphabet:** More characters means mismatches found quickly (average ~1-2 comparisons per position)
- **One-off search:** No preprocessing amortization needed

## Why It's Inefficient

After a partial match fails, it discards everything learned and restarts from scratch at the next position. KMP/Rabin-Karp/Boyer-Moore fix this.

| Algorithm | Time | Key Idea |
|---|---|---|
| KMP | O(n + m) | Failure function skips redundant comparisons |
| Rabin-Karp | O(n + m) avg | Rolling hash, O(1) window comparison |
| Boyer-Moore | O(n/m) best | Bad character rule skips large portions |

## Pitfalls

- Finds overlapping matches by default — this is correct behavior
- Python `str.find()` / `in` uses a mix of Boyer-Moore and Horspool internally and is faster than hand-written naive search
- Worst case is triggered by low-entropy text (e.g., DNA with small alphabet "ACGT")
