---
name: algo-tries
description: "Prefix tree for O(m) string insert/search and O(p+k) prefix queries; ideal for autocomplete and dictionary membership."
tool_type: python
primary_tool: Python
---

# Trie (Prefix Tree)

Each node = one character. Path root→node = prefix. `is_end` marks word boundaries.

## Complexity

| Operation | Time | Notes |
|-----------|------|-------|
| Insert | O(m) | m = word length |
| Search (exact) | O(m) | |
| Prefix check | O(p) | p = prefix length |
| Prefix collect | O(p + k) | k = results |
| Delete | O(m) | |

**vs Hash Table**: trie has O(p + k) prefix search vs O(n × m) scan; hash is usually lower memory for exact lookup.

## Implementation

```python
class TrieNode:
    def __init__(self):
        self.children: dict[str, "TrieNode"] = {}
        self.is_end = False


class Trie:
    def __init__(self):
        self.root = TrieNode()

    def _find_node(self, prefix: str) -> TrieNode | None:
        node = self.root
        for ch in prefix:
            if ch not in node.children:
                return None
            node = node.children[ch]
        return node

    def insert(self, word: str) -> None:
        node = self.root
        for ch in word:
            node = node.children.setdefault(ch, TrieNode())
        node.is_end = True

    def search(self, word: str) -> bool:
        node = self._find_node(word)
        return node is not None and node.is_end

    def starts_with(self, prefix: str) -> bool:
        return self._find_node(prefix) is not None

    def get_all_with_prefix(self, prefix: str) -> list[str]:
        results: list[str] = []
        node = self._find_node(prefix)
        if node is None:
            return results
        self._collect(node, prefix, results)
        return results

    def _collect(self, node: TrieNode, word: str, results: list[str]) -> None:
        if node.is_end:
            results.append(word)
        for ch, child in sorted(node.children.items()):
            self._collect(child, word + ch, results)

    def delete(self, word: str) -> bool:
        def _del(node: TrieNode, depth: int) -> bool:
            if depth == len(word):
                if not node.is_end:
                    return False
                node.is_end = False
                return not node.children  # safe to prune if leaf
            ch = word[depth]
            if ch not in node.children:
                return False
            if _del(node.children[ch], depth + 1):
                del node.children[ch]
                return not node.children and not node.is_end
            return False
        return _del(self.root, 0)

    def get_all_words(self) -> list[str]:
        return self.get_all_with_prefix("")
```

## Common Pattern: Autocomplete

```python
def autocomplete(trie: Trie, prefix: str, max_results: int = 5) -> list[str]:
    return trie.get_all_with_prefix(prefix)[:max_results]
```

## Pitfalls

- **`search` vs `starts_with`**: "ca" returns `True` for `starts_with` even if only "cat" was inserted; `search` checks `is_end`.
- **Deleting a prefix of another word**: never remove nodes that have children; only clear `is_end`. The delete implementation above handles this — prune only when the node becomes a childless non-endpoint.
- **Memory vs hash map**: a trie node with a `dict` per node uses more memory than a flat hash map for large, sparse alphabets. Use an array of size |Σ| only when the alphabet is small and dense (e.g., DNA: 4 chars).
- **Case sensitivity**: trie is case-sensitive by default; normalize input before inserting/searching.
