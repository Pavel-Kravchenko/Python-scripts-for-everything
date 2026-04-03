"""
Aho-Corasick multi-pattern string matching algorithm.

Builds a trie over all patterns, then adds suffix (failure) links via BFS
so that when a character fails to match in the trie, search falls back to
the longest proper suffix that is still a prefix of some pattern.

Main entry point: find(text, patterns) -> dict mapping each pattern to the
list of 0-based start positions where it occurs in text.
"""
from __future__ import annotations

from collections import defaultdict


# ---------------------------------------------------------------------------
# Queue (used for BFS during suffix-link construction)
# ---------------------------------------------------------------------------

class _Queue:
    def __init__(self):
        self._items: list = []

    def is_empty(self) -> bool:
        return len(self._items) == 0

    def put(self, x) -> None:
        self._items.append(x)

    def get(self):
        return self._items.pop(0)

    @property
    def size(self) -> int:
        return len(self._items)


# ---------------------------------------------------------------------------
# Trie node
# ---------------------------------------------------------------------------

class AhoNode:
    """Node in the Aho-Corasick trie."""

    def __init__(self):
        self.out: dict[str, "AhoNode"] = {}  # transitions: char -> child node
        self.patterns: list[str] = []        # patterns that end at this node
        self.fail: "AhoNode | None" = None   # suffix (failure) link


# ---------------------------------------------------------------------------
# Build trie and suffix links
# ---------------------------------------------------------------------------

def _build_trie(patterns: list[str]) -> AhoNode:
    """Insert all patterns into a trie and return the root."""
    root = AhoNode()
    for pattern in patterns:
        node = root
        for char in pattern:
            if char not in node.out:
                node.out[char] = AhoNode()
            node = node.out[char]
        node.patterns.append(pattern)
    return root


def build_automaton(patterns: list[str]) -> AhoNode:
    """Build the Aho-Corasick automaton: trie + BFS suffix links.

    Returns the root node of the automaton.
    """
    root = _build_trie(patterns)
    queue = _Queue()

    # All depth-1 nodes fail back to root
    for node in root.out.values():
        queue.put(node)
        node.fail = root

    # BFS to compute suffix links for deeper nodes
    while not queue.is_empty():
        current = queue.get()
        for char, child in current.out.items():
            queue.put(child)
            # Walk up suffix links until we find a node with a char transition
            fallback = current.fail
            while fallback is not None and char not in fallback.out:
                fallback = fallback.fail
            if fallback:
                child.fail = fallback.out[char]
            else:
                child.fail = root
            # Propagate output (dictionary) links
            child.patterns += child.fail.patterns

    return root


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def find(text: str, patterns: list[str]) -> dict[str, list[int]]:
    """Find all occurrences of every pattern in text.

    Returns a dict: pattern -> [list of 0-based start positions].
    """
    root = build_automaton(patterns)
    results: dict[str, list[int]] = defaultdict(list)
    node = root

    for i, char in enumerate(text):
        # Follow suffix links until a matching transition is found or we reach root
        while node is not None and char not in node.out:
            node = node.fail
        if node is None:
            node = root
            continue
        node = node.out[char]
        for pattern in node.patterns:
            results[pattern].append(i - len(pattern) + 1)

    return dict(results)


if __name__ == "__main__":
    # Example with biological-style strings
    patterns = ["ATG", "TGC", "ATGC", "GC"]
    text = "AATGCATGC"
    print(f"Text:     {text}")
    print(f"Patterns: {patterns}")
    matches = find(text, patterns)
    for pat, positions in sorted(matches.items()):
        print(f"  '{pat}': {positions}")

    print()

    # Original example from source
    patterns2 = ["m", "mama", "mkdir"]
    text2 = "mamammmkdirmmmamamdir"
    print(f"Text:     {text2}")
    print(f"Patterns: {patterns2}")
    matches2 = find(text2, patterns2)
    for pat, positions in sorted(matches2.items()):
        print(f"  '{pat}': {positions}")
