---
name: string-algorithms
description: Pattern matching algorithms — naive, KMP (failure function), Rabin-Karp (rolling hash), and DFA-based matching for sequence search
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# String Matching Algorithms

## When to Use
- **Naive**: short patterns (m small), large alphabet, or one-off search; no setup cost
- **KMP**: guaranteed O(n+m), single pattern, repetitive pattern structure (e.g., ATATATG), streaming input
- **Rabin-Karp**: multiple patterns of same length (hash all, scan once), plagiarism/duplicate detection
- **DFA**: small fixed alphabet (DNA: |Σ|=4), many texts with the same pattern, need O(1) per character with no backtracking

---

## Quick Reference

### Complexity Comparison

| Algorithm | Preprocessing | Search | Space | Notes |
|-----------|--------------|--------|-------|-------|
| Naive     | O(1)         | O(n×m) | O(1)  | Best: O(n) with large alphabet |
| KMP       | O(m)         | O(n)   | O(m)  | Never backtracks in text |
| Rabin-Karp| O(m)         | O(n+m) avg, O(n×m) worst | O(1) | k patterns: O(n + k×m) |
| DFA       | O(m×\|Σ\|)   | O(n)   | O(m×\|Σ\|) | O(1) per char, no fallback logic |

### When to Choose

| Scenario | Use |
|----------|-----|
| Single pattern, worst-case guarantee | KMP |
| Multiple same-length patterns | Rabin-Karp (hash all at once) |
| DNA/binary alphabet, many texts | DFA |
| Short pattern, large alphabet | Naive (first-char mismatch kills cost) |

---

## Key Patterns

### KMP Failure Function
`sp[i]` = length of longest proper prefix of `p[0:i+1]` that is also a suffix.

```python
Pattern:  A  T  G  C  A  T  G
Index:    0  1  2  3  4  5  6
sp:       0  0  0  0  1  2  3
```python

After a mismatch at position `j`, jump to `sp[j-1]` — skip re-examining what we already know matched.

### Rabin-Karp Rolling Hash (polynomial, mod prime)
```python
hash(s[i:i+m]) = (s[i]·b^(m-1) + s[i+1]·b^(m-2) + ... + s[i+m-1]) mod q
roll:  hash(s[i+1:i+m+1]) = (hash(s[i:i+m]) - s[i]·b^(m-1)) · b + s[i+m]  mod q
```python
- `b` = base (e.g., 256 for ASCII, 4 for DNA)
- `q` = large prime (e.g., 101, 1_000_000_007)
- Precompute `b^(m-1) mod q` once; each roll is O(1)

### DFA Construction via KMP
For each state `i` (representing matched `pattern[:i]`) and each character `c`:
```python
transition[i][c] = len(longest prefix of pattern that is a suffix of pattern[:i] + c)
```python
Uses `prefix_length(pattern, pattern[:i] + c)` — same KMP prefix trick.

---

## Code Templates

### Naive Search
```python
def naive_search(text: str, pattern: str) -> list[int]:
    n, m = len(text), len(pattern)
    return [i for i in range(n - m + 1) if text[i:i+m] == pattern]
```python

### KMP — Failure Function + Search
```python
def prefix_function(p: str) -> list[int]:
    sp = [0] * len(p)
    j = 0
    for i in range(1, len(p)):
        while j >= 0 and p[j] != p[i]:
            j = sp[j - 1] if j - 1 >= 0 else -1
        j += 1
        sp[i] = j
    return sp

def kmp_search(text: str, pattern: str) -> list[int]:
    matches = []
    f = prefix_function(pattern)
    n, m = len(text), len(pattern)
    j = 0
    for i in range(n):
        while j >= 0 and text[i] != pattern[j]:
            j = f[j - 1] if j - 1 >= 0 else -1
        j += 1
        if j == m:
            matches.append(i - m + 1)
            j = f[m - 1]   # allow overlapping matches
    return matches
```python

**Alternate (prefix_search on concatenated string):**
```python
def prefix_search(text: str, pattern: str) -> list[int]:
    """Run prefix function on pattern + '#' + text; positions where sp[i] == m."""
    combined = pattern + "#" + text
    m = len(pattern)
    sp = [0] * len(combined)
    j = 0
    for i in range(1, len(combined)):
        while j > 0 and combined[i] != combined[j]:
            j = sp[j - 1]
        if combined[i] == combined[j]:
            j += 1
        sp[i] = j
    # match at combined index i means text index i - (m + 1) - m + 1 = i - 2m
    # but combined[m+1+k] corresponds to text[k], so text index = i - m - 1 - m + 1 = i - 2m
    return [i - 2 * m for i in range(len(combined)) if sp[i] == m]
```python

### Rabin-Karp (polynomial rolling hash)
```python
def rabin_karp_all(text: str, pattern: str, base: int = 256, mod: int = 101) -> list[int]:
    n, m = len(text), len(pattern)
    if m > n:
        return []
    h = pow(base, m - 1, mod)          # base^(m-1) mod q
    ph = th = 0
    for i in range(m):
        ph = (base * ph + ord(pattern[i])) % mod
        th = (base * th + ord(text[i])) % mod
    matches = []
    for i in range(n - m + 1):
        if th == ph and text[i:i+m] == pattern:   # hash match → verify (avoids spurious hits)
            matches.append(i)
        if i < n - m:
            th = (base * (th - ord(text[i]) * h) + ord(text[i + m])) % mod
    return matches

def rabin_karp_multi(text: str, patterns: list[str], base: int = 256, mod: int = 101) -> dict[str, list[int]]:
    """Search k same-length patterns in one pass. All patterns must have equal length."""
    if not patterns:
        return {}
    m = len(patterns[0])
    pattern_hashes = {}
    h = pow(base, m - 1, mod)
    for p in patterns:
        ph = 0
        for c in p:
            ph = (base * ph + ord(c)) % mod
        pattern_hashes.setdefault(ph, []).append(p)
    results = {p: [] for p in patterns}
    th = 0
    for c in text[:m]:
        th = (base * th + ord(c)) % mod
    for i in range(len(text) - m + 1):
        if th in pattern_hashes:
            window = text[i:i+m]
            for p in pattern_hashes[th]:
                if window == p:
                    results[p].append(i)
        if i < len(text) - m:
            th = (base * (th - ord(text[i]) * h) + ord(text[i + m])) % mod
    return results
```python

### DFA Matching
```python
def prefix_length(pattern: str, probe: str) -> int:
    """Longest prefix of pattern that is also a suffix of probe (KMP trick)."""
    combined = pattern + "#" + probe + "$"
    sp = [0] * len(combined)
    j = 0
    for i in range(1, len(combined) - 1):
        while j > 0 and combined[i] != combined[j]:
            j = sp[j - 1]
        if combined[i] == combined[j]:
            j += 1
        sp[i] = j
    return sp[-2]  # last char before '$'

def build_automaton(pattern: str, alphabet: str) -> list[dict[str, int]]:
    """Returns list of transition dicts; state len(pattern) is accept."""
    return [
        {c: prefix_length(pattern, pattern[:i] + c) for c in alphabet}
        for i in range(len(pattern) + 1)
    ]

def dfa_search(text: str, automaton: list[dict[str, int]]) -> list[int]:
    accept = len(automaton) - 1
    state, matches = 0, []
    for i, c in enumerate(text):
        state = automaton[state][c]
        if state == accept:
            matches.append(i - accept + 1)
    return matches

# Usage
alphabet = "ATGC"
pattern  = "ATTCTGATTT"
dfa = build_automaton(pattern, alphabet)
hits = dfa_search("AATGCCGTATTCTATTCTGATTTCTGAATTCTGATTTTTAGT", dfa)
```python

---

## Common Pitfalls

- **KMP overlapping matches**: after a full match, set `j = f[m-1]`, not `j = 0` — otherwise overlapping occurrences like `AA` in `AAAA` are missed.
- **Rabin-Karp spurious hits**: a hash match is not a guarantee; always verify with `text[i:i+m] == pattern`. Skip verification only if collision probability is provably negligible.
- **Rolling hash negative values**: `(th - ord(text[i]) * h) % mod` can be negative in Python — it wraps correctly, but in C/Java you must add `mod` before taking `%`.
- **DFA alphabet completeness**: every character in the text must have a transition defined; unrecognized characters raise `KeyError`. Explicitly handle or restrict to known alphabet.
- **DFA preprocessing cost**: O(m×|Σ|) build time — not worth it for |Σ|=256 (ASCII) with small `m`; prefer KMP there.
- **Naive with `text[i:i+m]`**: creates a new string object per position (O(m) space each); use explicit char comparison for truly O(1) space.
- **KMP `j = -1` sentinel**: the reference implementation uses `j = -1` as a signal to advance without matching, equivalent to the "no prefix" state. Don't confuse with array indexing.

---

## Bioinformatics Connections

| Application | Algorithm | Notes |
|-------------|-----------|-------|
| Restriction site finding (EcoRI: `GAATTC`) | KMP or DFA | Single fixed pattern; DFA fast for streaming FASTQ |
| Motif scanning (TFBS, k-mer search) | Rabin-Karp | Hash all motif variants, single text pass |
| BLAST seed-and-extend | Naive | Short seed (11-mer default); large alphabet → fast mismatch |
| Tandem repeat detection | KMP prefix function | sp[i] reveals internal repetition period |
| Multiple restriction enzymes | Rabin-Karp multi | All enzymes same length → one scan |
| Long-read mapping seeds | DFA | Fixed seed pattern, millions of reads |

**KMP prefix function reveals period:**
```python
# If sp[m-1] > 0 and m % (m - sp[m-1]) == 0, pattern has period (m - sp[m-1])
# "ATGATGATG" -> period 3
```python

---

## Related Skills
- `graphs-dynamic-programming` — edit distance, Smith-Waterman use DP on character grids
- `advanced-string-structures` — O(n log n) suffix array construction enables all-against-all substring queries
- `numpy-pandas-wrangling` — vectorized k-mer counting with `np.frombuffer` for large genomes
