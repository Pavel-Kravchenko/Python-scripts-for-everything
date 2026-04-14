---
name: algo-suffix-trees
description: "Suffix trees — compressed trie of all suffixes, O(m) pattern search, O(n) construction via Ukkonen's algorithm"
tool_type: python
primary_tool: Python
---

# Suffix Trees

Compressed trie containing all suffixes of a string. O(n) nodes, O(m) pattern search, O(m+k) to find all k occurrences.

## Key Properties

- Exactly n leaves (one per suffix), at most n-1 internal nodes
- Each internal node (except root) has >= 2 children
- Edge labels are substrings, not single characters
- Path from root to leaf i spells suffix starting at position i

## Suffix Tree vs Suffix Array

| | Suffix Tree | Suffix Array |
|---|---|---|
| Space | O(n) but large constant | O(n) with small constant |
| Pattern search | O(m) | O(m log n) |
| Construction | O(n) via Ukkonen | O(n) via SA-IS |
| Implementation | Complex | Simpler |

## Naive Construction — O(n^2)

Build trie of all suffixes, then compress single-child paths.

```python
class SuffixTreeNode:
    def __init__(self):
        self.children = {}      # edge_label -> SuffixTreeNode
        self.suffix_index = -1  # leaf: starting position; internal: -1

    def is_leaf(self):
        return len(self.children) == 0

class SuffixTree:
    def __init__(self, text, terminal='$'):
        self.text = text + terminal
        self.root = SuffixTreeNode()
        self._build_trie()
        self._compress(self.root)

    def _build_trie(self):
        n = len(self.text)
        for i in range(n):
            current = self.root
            for j in range(i, n):
                char = self.text[j]
                if char not in current.children:
                    current.children[char] = SuffixTreeNode()
                current = current.children[char]
            current.suffix_index = i

    def _compress(self, node):
        changed = True
        while changed:
            changed = False
            for key in list(node.children.keys()):
                child = node.children[key]
                if len(child.children) == 1:
                    gk = list(child.children.keys())[0]
                    node.children[key + gk] = child.children[gk]
                    del node.children[key]
                    changed = True
        for child in node.children.values():
            self._compress(child)

    def search(self, pattern):
        """Return all starting positions where pattern occurs."""
        node, remaining = self.root, pattern
        while remaining:
            found = False
            for edge, child in node.children.items():
                if edge[0] == remaining[0]:
                    match_len = min(len(edge), len(remaining))
                    if edge[:match_len] != remaining[:match_len]:
                        return []
                    remaining = remaining[match_len:]
                    node = child
                    found = True
                    break
            if not found:
                return []
        positions = []
        self._collect_leaves(node, positions)
        return sorted(positions)

    def _collect_leaves(self, node, positions):
        if node.is_leaf():
            positions.append(node.suffix_index)
        else:
            for child in node.children.values():
                self._collect_leaves(child, positions)
```

## Operations Complexity

| Operation | Naive Build | Ukkonen Build |
|---|---|---|
| Construction | O(n^2) | O(n) |
| Pattern search | O(m) | O(m) |
| All k occurrences | O(m + k) | O(m + k) |
| Longest repeated substring | O(n) | O(n) |
| LCS of two strings | O(n + m) | O(n + m) |

## Longest Common Substring

Concatenate strings with unique separators ("ABAB$BABA#"), build suffix tree, find deepest internal node with leaves from both strings.

## Pitfalls

- Always append a unique terminal character ('$') — without it, some suffixes that are prefixes of others won't get their own leaf
- Naive construction is O(n^2); for production use, implement Ukkonen's O(n) algorithm or use a suffix array instead
- Edge labels should store (start, end) indices into the original string, not actual substrings, to keep O(n) space
