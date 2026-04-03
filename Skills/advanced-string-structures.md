---
name: advanced-string-structures
description: Tries, Aho-Corasick multi-pattern matching, suffix arrays with LCP, and suffix trees for genome indexing
---

# Advanced String Structures

## When to Use
- **Trie**: prefix search, autocomplete, k-mer membership testing, IP longest-prefix match
- **Aho-Corasick**: find all occurrences of k patterns simultaneously in one O(n) pass — restriction enzyme sites, primer binding, motif scanning
- **Suffix Array + LCP**: genome-wide pattern search, longest repeated substring/motif, k-mer frequency counting, de Bruijn graph construction
- **Suffix Tree**: O(m) pattern match, longest common substring of two sequences, tandem repeat detection, LZ compression

---

## Quick Reference

### Complexity Table

| Structure | Build | Pattern Search | Space | Notes |
|-----------|-------|---------------|-------|-------|
| Trie | O(n·m) | O(m) | O(n·m) | n words, m avg length |
| Aho-Corasick | O(M) | O(n + M + z) | O(M) | M = total pattern len, z = matches |
| Suffix Array (naive) | O(n² log n) | O(m log n) | O(n) | |
| Suffix Array (Manber-Myers) | O(n log n) | O(m log n) | O(n) | doubling trick |
| Suffix Array (DC3/SA-IS) | O(n) | O(m log n) | O(n) | optimal |
| LCP Array (Kasai) | O(n) | — | O(n) | after SA |
| Suffix Tree (naive) | O(n²) | O(m) | O(n) | |
| Suffix Tree (Ukkonen) | O(n) | O(m + k) | O(n) | k = occurrences |

---

## Key Patterns

### Trie
- Each node = one character; path root→node = prefix
- `is_end` flag marks word boundaries; prefix queries need only reach a node (not `is_end`)
- Delete: unset `is_end`; remove node only if it has no children

### Aho-Corasick
- Build trie over all patterns; then BFS to attach **failure links** (suffix links)
- `fail[v]` = node for the longest proper suffix of `str(v)` that is a trie prefix
- Depth-1 nodes always fail to root
- **Output (dictionary) links**: `out[v] = direct_patterns + out[fail[v]]` — propagated during BFS so search only checks `node.patterns`
- Search: walk text character by character; on mismatch follow `fail` links until match or root

### Suffix Array
- `SA[i]` = start index of the i-th lexicographically smallest suffix
- Append sentinel `$` (lexicographically smallest) to force all suffixes to be distinct and leaf-terminated
- Pattern search: binary search for leftmost and rightmost position where suffix starts with pattern → O(m log n)
- **LCP[i]** = length of longest common prefix between `SA[i-1]` and `SA[i]`; Kasai's algorithm builds it in O(n)
- Longest repeated substring = `S[SA[k] : SA[k] + max(LCP)]`
- Distinct substrings = n(n+1)/2 − sum(LCP)

### Suffix Tree
- Compressed trie of all suffixes: n leaves, ≤ n−1 internal nodes, O(n) total
- Edge labels are substrings (store as index pairs, not copies)
- Ukkonen's algorithm builds it online in O(n) using suffix links + active point + remainder
- LCS of two strings: concatenate as `S1$S2#`, build generalized suffix tree, find deepest internal node whose subtree has leaves from both strings

---

## Code Templates

### Trie

```python
class TrieNode:
    def __init__(self):
        self.children = {}   # char -> TrieNode
        self.is_end = False

class Trie:
    def __init__(self):
        self.root = TrieNode()

    def insert(self, word: str) -> None:
        node = self.root
        for ch in word:
            if ch not in node.children:
                node.children[ch] = TrieNode()
            node = node.children[ch]
        node.is_end = True

    def search(self, word: str) -> bool:
        node = self._find(word)
        return node is not None and node.is_end

    def starts_with(self, prefix: str) -> bool:
        return self._find(prefix) is not None

    def _find(self, prefix: str):
        node = self.root
        for ch in prefix:
            if ch not in node.children:
                return None
            node = node.children[ch]
        return node

    def get_all_with_prefix(self, prefix: str) -> list[str]:
        results = []
        node = self._find(prefix)
        if node:
            self._dfs(node, prefix, results)
        return results

    def _dfs(self, node, word, results):
        if node.is_end:
            results.append(word)
        for ch, child in sorted(node.children.items()):
            self._dfs(child, word + ch, results)
```

### Aho-Corasick (reference implementation)

```python
from collections import defaultdict

class AhoNode:
    def __init__(self):
        self.out: dict[str, "AhoNode"] = {}
        self.patterns: list[str] = []
        self.fail: "AhoNode | None" = None

def build_automaton(patterns: list[str]) -> AhoNode:
    root = AhoNode()
    for pat in patterns:                      # build trie
        node = root
        for ch in pat:
            if ch not in node.out:
                node.out[ch] = AhoNode()
            node = node.out[ch]
        node.patterns.append(pat)

    from collections import deque
    queue = deque()
    for child in root.out.values():           # depth-1 → root
        child.fail = root
        queue.append(child)

    while queue:                              # BFS failure links
        cur = queue.popleft()
        for ch, child in cur.out.items():
            queue.append(child)
            fb = cur.fail
            while fb is not None and ch not in fb.out:
                fb = fb.fail
            child.fail = fb.out[ch] if fb else root
            child.patterns += child.fail.patterns   # propagate output links
    return root

def find(text: str, patterns: list[str]) -> dict[str, list[int]]:
    root = build_automaton(patterns)
    results: dict[str, list[int]] = defaultdict(list)
    node = root
    for i, ch in enumerate(text):
        while node is not None and ch not in node.out:
            node = node.fail
        if node is None:
            node = root
            continue
        node = node.out[ch]
        for pat in node.patterns:
            results[pat].append(i - len(pat) + 1)
    return dict(results)
```

### Suffix Array (naive O(n² log n) + Kasai LCP)

```python
def build_suffix_array(text: str) -> list[int]:
    """Naive O(n² log n). Append '$' before calling."""
    n = len(text)
    return sorted(range(n), key=lambda i: text[i:])

def build_lcp_array(text: str, sa: list[int]) -> list[int]:
    """Kasai's algorithm, O(n)."""
    n = len(text)
    rank = [0] * n
    for i, s in enumerate(sa):
        rank[s] = i
    lcp = [0] * n
    k = 0
    for i in range(n):
        if rank[i] == 0:
            k = 0
            continue
        j = sa[rank[i] - 1]
        while i + k < n and j + k < n and text[i + k] == text[j + k]:
            k += 1
        lcp[rank[i]] = k
        if k:
            k -= 1
    return lcp

def search_suffix_array(text: str, sa: list[int], pattern: str) -> list[int]:
    """Binary search, O(m log n). Returns sorted list of start positions."""
    n, m = len(text), len(pattern)
    left, right = 0, n
    while left < right:
        mid = (left + right) // 2
        if text[sa[mid]:sa[mid]+m] < pattern:
            left = mid + 1
        else:
            right = mid
    lo = left
    left, right = lo, n
    while left < right:
        mid = (left + right) // 2
        if text[sa[mid]:sa[mid]+m] <= pattern:   # note <=
            left = mid + 1
        else:
            right = mid
    return sorted(sa[lo:left])
```

### Suffix Array Manber-Myers (O(n log n))

```python
def build_suffix_array_mm(text: str) -> list[int]:
    n = len(text)
    sa = sorted(range(n), key=lambda i: text[i])
    rank = [0] * n
    for i in range(1, n):
        rank[sa[i]] = rank[sa[i-1]] + (text[sa[i]] != text[sa[i-1]])
    k = 1
    while k < n:
        key = lambda i: (rank[i], rank[i + k] if i + k < n else -1)
        sa = sorted(range(n), key=key)
        new_rank = [0] * n
        for i in range(1, n):
            new_rank[sa[i]] = new_rank[sa[i-1]] + (key(sa[i]) != key(sa[i-1]))
        rank = new_rank
        if rank[sa[-1]] == n - 1:
            break
        k *= 2
    return sa
```

---

## Common Pitfalls

- **Sentinel character**: suffix arrays/trees require a unique terminator (`$`) smaller than all alphabet characters; omitting it causes ties and incorrect sort order
- **Aho-Corasick output links**: forgetting to propagate `out[fail[v]]` patterns to `out[v]` during BFS causes missed matches for shorter patterns embedded in longer ones (e.g., "he" inside "she")
- **Trie delete**: unmarking `is_end` is correct only if the node still has children; remove the node only when it becomes childless and is not `is_end` anywhere up the chain
- **LCP array indexing**: `LCP[0]` is undefined (set to 0); LCP values are between consecutive entries in sorted order, not consecutive positions in text
- **Suffix tree edge labels**: store `(start, end)` index pairs, not substrings — copying creates O(n²) space

---

## Bioinformatics Connections

| Application | Structure | Detail |
|-------------|-----------|--------|
| Short-read alignment (BWA, Bowtie) | Suffix Array + BWT/FM-index | BWT = `text[SA[i]-1]`; FM-index enables O(m) backward search with compressed SA |
| Multi-motif scanning (restriction sites, TF binding) | Aho-Corasick | One pass over genome finds all k patterns simultaneously |
| K-mer counting / genome assembly | Trie or Suffix Array | Binary search on SA counts k-mer occurrences; de Bruijn graph nodes are (k−1)-mers stored in trie |
| Longest repeated motif / tandem repeats | Suffix Array + LCP | `max(LCP)` gives LRS length; LCP runs of length ≥ k find tandem candidates |
| Longest common substring of two sequences | Suffix Tree (generalized) | Concatenate `seq1 $ seq2 #`; deepest internal node with leaves from both strings |
| Distinct k-mer enumeration | Suffix Array | Walk SA with LCP; adjacent entries sharing k-prefix belong to same k-mer |
| Sequence compression | Suffix Tree | Backbone of LZ77/LZ78; edge labels encode repeated factors |

---

## Related Skills
- `string-algorithms` — KMP/Rabin-Karp/DFA prerequisite patterns used in trie and Aho-Corasick construction
- `linux-git-bash` — indexing tools (bwa index, samtools faidx) that use these structures internally
- `python-advanced-sql` — implementing trie/SA-based k-mer tables backed by SQLite for large datasets
