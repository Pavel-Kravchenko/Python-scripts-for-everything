---
name: algo-rabin-karp
description: "Rabin-Karp hash-based string matching — rolling hash for O(n+m) average, excels at multi-pattern search"
tool_type: python
primary_tool: Python
---

# Rabin-Karp Algorithm

## Key Insight

Compare hash values instead of strings. If `hash(window) != hash(pattern)`, definitely no match. Only verify character-by-character when hashes match (may be a collision).

## Rolling Hash Formula

Treat string as a number in base d (typically 256 for ASCII), mod a prime q:

```
h(s[i..i+m-1]) = s[i]*d^(m-1) + s[i+1]*d^(m-2) + ... + s[i+m-1]
```

To slide window from position i to i+1 in O(1):
```
new_hash = (d * (old_hash - s[i] * d^(m-1)) + s[i+m]) mod q
```

## Implementation

```python
def rabin_karp_search(text: str, pattern: str, d: int = 256, q: int = 101) -> list[int]:
    """O(n+m) average, O(nm) worst case. Always verify on hash match."""
    n, m = len(text), len(pattern)
    if m == 0 or m > n:
        return []

    # Precompute d^(m-1) mod q
    h = 1
    for _ in range(m - 1):
        h = (h * d) % q

    # Initial hashes
    pattern_hash = window_hash = 0
    for i in range(m):
        pattern_hash = (d * pattern_hash + ord(pattern[i])) % q
        window_hash = (d * window_hash + ord(text[i])) % q

    matches = []
    for i in range(n - m + 1):
        if pattern_hash == window_hash:
            if text[i:i+m] == pattern:  # verify to handle collisions
                matches.append(i)
        if i < n - m:
            window_hash = (d * (window_hash - ord(text[i]) * h) + ord(text[i + m])) % q
            if window_hash < 0:
                window_hash += q
    return matches
```

## Complexity

| Case | Time | Notes |
|---|---|---|
| Average | O(n + m) | Few hash collisions |
| Worst (many collisions) | O(nm) | Every position needs verification |
| Space | O(1) | |

## Multi-Pattern Search (Where Rabin-Karp Shines)

For k patterns of the same length: compute one rolling hash over text, check against a set of k pattern hashes. Average O(n + km) vs naive O(knm).

## Why Prime q?

- Distributes hash values uniformly, reducing collision probability
- Without modulo, hash values overflow (256^10 is astronomically large)
- Common choices: 101 (small demo), 1_000_000_007 (production)

## Pitfalls

- Always verify on hash match — collisions are inevitable, not bugs
- Smaller q = more collisions = more O(m) verifications = slower
- Python handles negative modulo correctly, but other languages may not — add `if hash < 0: hash += q`
- For different-length patterns, Aho-Corasick is more appropriate
