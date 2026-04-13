---
name: algo-aho-corasick
description: "The Aho-Corasick algorithm solves the **multi-pattern matching problem**: given a text of length `n` and `k` patterns with total length `m`, find all occurrences of all patterns in the text."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/02_aho_corasick.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Aho-Corasick Algorithm

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/02_aho_corasick.ipynb`*

# Aho-Corasick Algorithm

## Multi-Pattern String Matching

The Aho-Corasick algorithm solves the **multi-pattern matching problem**: given a text of length `n` and `k` patterns with total length `m`, find all occurrences of all patterns in the text.

---

## 1. The Multi-Pattern Matching Problem

### Problem Statement

**Input:**
- Text `T` of length `n`
- Set of `k` patterns `P = {p₁, p₂, ..., pₖ}` with total length `m = |p₁| + |p₂| + ... + |pₖ|`

**Output:**
- All pairs `(i, j)` where pattern `pⱼ` occurs at position `i` in text `T`

### Naive Approaches

| Approach | Time Complexity | Description |
|----------|-----------------|-------------|
| Brute force | O(n × m) | For each pattern, check each position |
| KMP × k times | O(n × k + m) | Run KMP for each pattern separately |
| **Aho-Corasick** | **O(n + m + z)** | Single pass through text |

Where `z` = total number of matches (output size).

### Why Aho-Corasick?

- **Single pass**: Process the text only once, regardless of the number of patterns
- **Optimal**: Linear in input size plus output size
- **Real-world applications**: 
  - Intrusion detection systems (matching malware signatures)
  - Spam filters (detecting forbidden phrases)
  - DNA sequence analysis (finding motifs)
  - Text editors (search/replace multiple terms)

---

## 2. Trie + Failure Links: The Key Insight

### The Idea

Aho-Corasick combines two concepts:

1. **Trie (prefix tree)**: Store all patterns in a trie for efficient prefix matching
2. **Failure links (like KMP)**: When a mismatch occurs, jump to the longest proper suffix that is also a prefix of some pattern

### From KMP to Aho-Corasick

```
KMP:           Single pattern → Prefix function on that pattern
Aho-Corasick:  Multiple patterns → Trie + failure links between nodes
```

### Failure Link Definition

For a node `v` representing string `S`:

**`fail[v]`** = node representing the **longest proper suffix** of `S` that is also a **prefix of some pattern** (exists in the trie)

If no such suffix exists, `fail[v]` = root.

---

## 3. ASCII Art: Trie with Failure Links

**Patterns: "he", "she", "his", "hers"**

```
              (0) root
             /    \
           h        s
           |        |
          (1)      (6)
         /   \      |
        e     i     h
        |     |     |
       (2)*  (4)   (7)
        |     |     |
        r     s     e
        |     |     |
       (3)   (5)*  (8)*
        |
        s
        |
       (9)*

* = pattern end node (output node)

States:
  (0) = ""      (root)
  (1) = "h"
  (2) = "he"    [outputs: "he"]
  (3) = "her"
  (4) = "hi"
  (5) = "his"   [outputs: "his"]
  (6) = "s"
  (7) = "sh"
  (8) = "she"   [outputs: "she", "he"]  ← also outputs "he" via dictionary link!
  (9) = "hers"  [outputs: "hers"]
```

### Failure Links (dashed arrows)

```
fail[1] = 0    "h"   → no proper suffix is prefix of pattern → root
fail[2] = 0    "he"  → "e" not in trie → root
fail[3] = 0    "her" → "er", "r" not in trie → root
fail[4] = 0    "hi"  → "i" not in trie → root
fail[5] = 0    "his" → "is", "s" not in trie → root  (wait, "s" IS in trie!)
fail[5] = 6    "his" → "s" is in trie → state 6
fail[6] = 0    "s"   → no proper suffix → root
fail[7] = 1    "sh"  → "h" is in trie → state 1
fail[8] = 2    "she" → "he" is in trie → state 2  ← This is key!
fail[9] = 6    "hers"→ "s" is in trie → state 6
```

**Key insight**: `fail[8] = 2` means when we're at "she", we also know "he" matches!

---

## 4. Building the Automaton: BFS for Failure Links

### Algorithm

We compute failure links level by level using BFS:

```
1. fail[root] = root
2. For all children c of root: fail[c] = root  (depth-1 nodes)
3. BFS from depth-1 nodes:
   For each node u with child v via edge 'a':
     - Start at f = fail[u]
     - While f ≠ root AND f has no child 'a':
         f = fail[f]
     - If f has child 'a':
         fail[v] = f.child['a']
     - Else:
         fail[v] = root
```

### Example: Computing Failure Links

```
Patterns: "he", "she", "his", "hers"

═══════════════════════════════════════════════════════════════
Level 0: root
═══════════════════════════════════════════════════════════════
  fail[root] = root

═══════════════════════════════════════════════════════════════
Level 1: Process children of root (states 1, 6)
═══════════════════════════════════════════════════════════════
  State 1 ('h'): fail[1] = root (by definition, depth-1)
  State 6 ('s'): fail[6] = root (by definition, depth-1)

═══════════════════════════════════════════════════════════════
Level 2: Process children of level 1 (states 2, 4, 7)
═══════════════════════════════════════════════════════════════
  State 2 ('e' from state 1, representing "he"):
    parent = 1, parent_fail = fail[1] = root
    Does root have 'e' child? NO
    → fail[2] = root

  State 4 ('i' from state 1, representing "hi"):
    parent = 1, parent_fail = fail[1] = root
    Does root have 'i' child? NO
    → fail[4] = root

  State 7 ('h' from state 6, representing "sh"):
    parent = 6, parent_fail = fail[6] = root
    Does root have 'h' child? YES → state 1
    → fail[7] = 1

═══════════════════════════════════════════════════════════════
Level 3: Process children of level 2 (states 3, 5, 8)
═══════════════════════════════════════════════════════════════
  State 3 ('r' from state 2, representing "her"):
    parent = 2, parent_fail = fail[2] = root
    Does root have 'r' child? NO
    → fail[3] = root

  State 5 ('s' from state 4, representing "his"):
    parent = 4, parent_fail = fail[4] = root
    Does root have 's' child? YES → state 6
    → fail[5] = 6

  State 8 ('e' from state 7, representing "she"):
    parent = 7, parent_fail = fail[7] = 1
    Does state 1 have 'e' child? YES → state 2
    → fail[8] = 2  ★ Key: "she" links to "he"!

═══════════════════════════════════════════════════════════════
Level 4: Process children of level 3 (state 9)
═══════════════════════════════════════════════════════════════
  State 9 ('s' from state 3, representing "hers"):
    parent = 3, parent_fail = fail[3] = root
    Does root have 's' child? YES → state 6
    → fail[9] = 6
```

---

## 5. Dictionary Links (Output Links)

### The Problem with Failure Links Alone

When we reach state 8 ("she"), we detect pattern "she". But what about "he"?

Following `fail[8] = 2`, we find state 2 which outputs "he". But this requires extra traversal.

### Dictionary Links

**`dict_link[v]`** = the nearest node reachable via failure links that is a pattern endpoint.

Alternative approach: **Output list** - precompute all patterns that end at each state (including via failure chain).

```
out[8] = ["she", "he"]  ← "she" is direct, "he" comes from fail[8]=2
out[2] = ["he"]
out[5] = ["his"]
out[9] = ["hers"]
```

### Computing Output Lists (during BFS)

```
For each node v (in BFS order):
    out[v] = patterns ending at v
    out[v] += out[fail[v]]  ← inherit from failure link
```

This ensures all shorter patterns "contained" in longer ones are reported.

---

## 6. Search Phase

### Algorithm

```
state = root
for i = 0 to n-1:
    c = text[i]
    while state ≠ root AND state has no edge c:
        state = fail[state]  # follow failure links
    if state has edge c:
        state = state.child[c]
    # Report all patterns ending here
    for pattern in out[state]:
        report match at position (i - len(pattern) + 1)
```

### Search Example

```
Text: "ushers"
Patterns: "he", "she", "hers"

Position 0: 'u'
  state = root, no 'u' edge
  → state = root

Position 1: 's'
  state = root, has 's' edge → state 6
  out[6] = [] → no match

Position 2: 'h'
  state = 6, has 'h' edge → state 7
  out[7] = [] → no match

Position 3: 'e'
  state = 7, has 'e' edge → state 8
  out[8] = ["she", "he"]
  ★ MATCH: "she" at position 3-3+1 = 1
  ★ MATCH: "he" at position 3-2+1 = 2

Position 4: 'r'
  state = 8, no 'r' edge
  → follow fail[8] = 2
  state = 2, has 'r' edge → state 3
  out[3] = [] → no match

Position 5: 's'
  state = 3, has 's' edge → state 9
  out[9] = ["hers"]
  ★ MATCH: "hers" at position 5-4+1 = 2

═══════════════════════════════════════════════════════════════
Results:
  "she"  found at position 1
  "he"   found at position 2
  "hers" found at position 2
═══════════════════════════════════════════════════════════════
```

---

## 7. Complexity Analysis

| Phase | Time Complexity | Space Complexity |
|-------|-----------------|------------------|
| Build trie | O(m) | O(m × σ) or O(m) with dict |
| Compute failure links | O(m) | O(m) |
| Search | O(n + z) | O(1) extra |
| **Total** | **O(n + m + z)** | **O(m)** |

Where:
- `n` = text length
- `m` = total pattern length (sum of all pattern lengths)
- `z` = number of pattern occurrences (output size)
- `σ` = alphabet size

### Why is Search O(n + z)?

**Amortized analysis**: Each character either:
1. Advances the state down the trie (at most `n` times total)
2. Follows a failure link (decreases "depth" in the trie)

Since we can only go up as many times as we went down, total transitions = O(n).
Output reporting adds O(z) for z matches.

---

## 8. Implementation

```python
from collections import deque
from typing import List, Dict, Tuple, Optional


class AhoCorasickNode:
    """A node in the Aho-Corasick automaton.
    
    Attributes:
        children: Dictionary mapping characters to child nodes.
        fail: Failure link - longest proper suffix that is a prefix of some pattern.
        output: List of patterns that end at this node (including via failure chain).
        depth: Depth of this node in the trie (length of string it represents).
    """
    
    def __init__(self, depth: int = 0):
        self.children: Dict[str, 'AhoCorasickNode'] = {}
        self.fail: Optional['AhoCorasickNode'] = None
        self.output: List[str] = []  # Patterns ending at this node
        self.depth: int = depth
    
    def __repr__(self) -> str:
        return f"Node(depth={self.depth}, children={list(self.children.keys())}, output={self.output})"
```

```python
class AhoCorasick:
    """Aho-Corasick automaton for multi-pattern string matching.
    
    The algorithm finds all occurrences of multiple patterns in a text
    in O(n + m + z) time, where:
        n = text length
        m = total pattern length
        z = number of matches
    
    Example:
        >>> ac = AhoCorasick()
        >>> ac.add_pattern("he")
        >>> ac.add_pattern("she")
        >>> ac.add_pattern("hers")
        >>> ac.build()
        >>> ac.search("ushers")
        [('she', 1), ('he', 2), ('hers', 2)]
    """
    
    def __init__(self):
        """Initialize the automaton with an empty root node."""
        self.root = AhoCorasickNode(depth=0)
        self.root.fail = self.root  # Root's failure link points to itself
        self._built = False
        self._patterns: List[str] = []
    
    def add_pattern(self, pattern: str) -> None:
        """Add a pattern to the automaton.
        
        Args:
            pattern: The pattern string to add.
        
        Raises:
            ValueError: If pattern is empty.
            RuntimeError: If automaton has already been built.
        """
        if not pattern:
            raise ValueError("Pattern cannot be empty")
        if self._built:
            raise RuntimeError("Cannot add patterns after building the automaton")
        
        self._patterns.append(pattern)
        node = self.root
        
        # Traverse/create path for this pattern
        for char in pattern:
            if char not in node.children:
                node.children[char] = AhoCorasickNode(depth=node.depth + 1)
            node = node.children[char]
        
        # Mark pattern endpoint
        node.output.append(pattern)
    
    def add_patterns(self, patterns: List[str]) -> None:
        """Add multiple patterns to the automaton.
        
        Args:
            patterns: List of pattern strings to add.
        """
        for pattern in patterns:
            self.add_pattern(pattern)
    
    def build(self) -> None:
        """Build the automaton by computing failure links using BFS.
        
        This method must be called after adding all patterns and before searching.
        Time complexity: O(m) where m is the total length of all patterns.
        """
        if self._built:
            return
        
        queue = deque()
        
        # Initialize: children of root have fail link to root
        for child in self.root.children.values():
            child.fail = self.root
            queue.append(child)
        
        # BFS to compute failure links
        while queue:
            current = queue.popleft()
            
            for char, child in current.children.items():
                queue.append(child)
                
                # Find failure link for child
                fail_node = current.fail
                while fail_node != self.root and char not in fail_node.children:
                    fail_node = fail_node.fail
                
                if char in fail_node.children:
                    child.fail = fail_node.children[char]
                else:
                    child.fail = self.root
                
                # Merge output lists (dictionary/suffix links)
                # This ensures we report all patterns including shorter ones
                child.output = child.output + child.fail.output
        
        self._built = True
    
    def search(self, text: str) -> List[Tuple[str, int]]:
        """Search for all pattern occurrences in the text.
        
        Args:
            text: The text to search in.
        
        Returns:
            List of tuples (pattern, position) for each match found.
            Position is the starting index (0-based) of the match.
        
        Raises:
            RuntimeError: If automaton has not been built yet.
        
        Time complexity: O(n + z) where n is text length, z is number of matches.
        """
        if not self._built:
            raise RuntimeError("Automaton must be built before searching")
        
        results = []
        node = self.root
        
        for i, char in enumerate(text):
            # Follow failure links until we find a match or reach root
            while node != self.root and char not in node.children:
                node = node.fail
            
            # Transition on the current character
            if char in node.children:
                node = node.children[char]
            # else: stay at root
            
            # Report all matches at this position
            for pattern in node.output:
                start_pos = i - len(pattern) + 1
                results.append((pattern, start_pos))
        
        return results
    
    def search_positions_only(self, text: str) -> List[int]:
        """Search and return only the starting positions of matches.
        
        Args:
            text: The text to search in.
        
        Returns:
            List of starting positions (0-based) where any pattern was found.
        """
        if not self._built:
            raise RuntimeError("Automaton must be built before searching")
        
        positions = []
        node = self.root
        
        for i, char in enumerate(text):
            while node != self.root and char not in node.children:
                node = node.fail
            
            if char in node.children:
                node = node.children[char]
            
            for pattern in node.output:
                positions.append(i - len(pattern) + 1)
        
        return positions
    
    def visualize(self) -> None:
        """Print a simple visualization of the trie structure."""
        def _print_node(node, prefix="", label="root"):
            output_str = f" → outputs: {node.output}" if node.output else ""
            fail_info = ""
            if node.fail and node.fail != self.root and node != self.root:
                fail_info = f" (fail→depth {node.fail.depth})"
            print(f"{prefix}{label}{output_str}{fail_info}")
            
            children = list(node.children.items())
            for i, (char, child) in enumerate(children):
                is_last = (i == len(children) - 1)
                child_prefix = prefix + ("    " if is_last else "│   ")
                connector = "└── " if is_last else "├── "
                _print_node(child, prefix + connector[:-4], connector + f"'{char}'")
        
        print("\nTrie Structure:")
        _print_node(self.root)
    
    @property
    def patterns(self) -> List[str]:
        """Return the list of patterns added to the automaton."""
        return self._patterns.copy()
```
