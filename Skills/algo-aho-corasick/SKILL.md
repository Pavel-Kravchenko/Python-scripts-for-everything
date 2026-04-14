---
name: algo-aho-corasick
description: "Multi-pattern string matching in O(n + m + z) via a trie augmented with KMP-style failure links."
tool_type: python
primary_tool: Python
---

# Aho-Corasick Algorithm

## When to Use

| Approach | Time | Use when |
|----------|------|----------|
| Brute force | O(n × m) | one pattern, tiny text |
| KMP × k | O(n × k + m) | few patterns |
| **Aho-Corasick** | **O(n + m + z)** | many patterns, single pass |

n = text length, m = total pattern length, k = pattern count, z = match count.

## Key Concepts

- **Trie**: insert all patterns; each node = prefix
- **Failure link** `fail[v]`: longest proper suffix of the string at `v` that is also a trie prefix. Computed BFS level-by-level.
- **Output list** `out[v]`: all patterns ending at `v`, including those reachable via failure chain (`out[v] += out[fail[v]]`).

## Build (BFS for failure links)

```
fail[root] = root
depth-1 nodes: fail[c] = root
for each node u (BFS order):
    for char, child in u.children:
        f = fail[u]
        while f != root and char not in f.children:
            f = fail[f]
        child.fail = f.children[char] if char in f.children else root
        child.out = child.out + child.fail.out   # merge output lists
```

## Search

```
state = root
for i, char in enumerate(text):
    while state != root and char not in state.children:
        state = state.fail
    if char in state.children:
        state = state.children[char]
    for pattern in state.out:
        yield pattern, i - len(pattern) + 1
```

## Implementation

```python
from collections import deque

class AhoCorasickNode:
    def __init__(self):
        self.children = {}
        self.fail = None
        self.output = []


class AhoCorasick:
    def __init__(self):
        self.root = AhoCorasickNode()
        self.root.fail = self.root
        self._built = False

    def add_pattern(self, pattern: str) -> None:
        if self._built:
            raise RuntimeError("Cannot add patterns after build()")
        node = self.root
        for ch in pattern:
            node = node.children.setdefault(ch, AhoCorasickNode())
        node.output.append(pattern)

    def build(self) -> None:
        queue = deque()
        for child in self.root.children.values():
            child.fail = self.root
            queue.append(child)
        while queue:
            cur = queue.popleft()
            for ch, child in cur.children.items():
                queue.append(child)
                f = cur.fail
                while f is not self.root and ch not in f.children:
                    f = f.fail
                child.fail = f.children[ch] if ch in f.children else self.root
                if child.fail is child:          # avoid self-loop at root's children
                    child.fail = self.root
                child.output = child.output + child.fail.output
        self._built = True

    def search(self, text: str) -> list[tuple[str, int]]:
        if not self._built:
            raise RuntimeError("Call build() before search()")
        results = []
        node = self.root
        for i, ch in enumerate(text):
            while node is not self.root and ch not in node.children:
                node = node.fail
            if ch in node.children:
                node = node.children[ch]
            for pattern in node.output:
                results.append((pattern, i - len(pattern) + 1))
        return results
```

## Complexity

| Phase | Time | Space |
|-------|------|-------|
| Build trie | O(m) | O(m × σ) |
| Failure links (BFS) | O(m) | O(m) |
| Search | O(n + z) | O(1) extra |
| **Total** | **O(n + m + z)** | **O(m)** |

σ = alphabet size.

## Pitfalls

- **Must call `build()` before `search()`**: failure links are not set during `add_pattern`.
- **Output list merging is O(z) not O(m)**: merging via `out[v] += out[fail[v]]` during BFS is safe because each match is reported once; skipping it causes missed overlapping patterns (e.g., "he" inside "she").
- **Root's self-loop**: `root.fail = root`; without this, the while loop in search has no termination condition.
- **Adding patterns after build raises**: rebuild the automaton from scratch if new patterns are needed.
- **Case sensitivity**: Aho-Corasick is case-sensitive by default; normalize text and patterns to the same case before inserting.
