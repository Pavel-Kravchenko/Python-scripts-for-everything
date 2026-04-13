---
name: algo-suffix-trees
description: "A **Suffix Tree** is one of the most powerful data structures for string processing. It's essentially a compressed trie of all suffixes of a given string, enabling extremely efficient pattern matching"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/04_suffix_trees.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Suffix Trees

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/04_suffix_trees.ipynb`*


A **Suffix Tree** is one of the most powerful data structures for string processing. It's essentially a compressed trie of all suffixes of a given string, enabling extremely efficient pattern matching and many other string operations.

## 1. What is a Suffix Tree?

A **Suffix Tree** for a string `S` of length `n` is a compressed trie containing all `n` suffixes of `S`.

### Key Properties

1. **Exactly `n` leaves** - one for each suffix (positions 0 to n-1)
2. **At most `n-1` internal nodes** (excluding root) - total O(n) nodes
3. **Each internal node (except root) has ≥2 children**
4. **Edge labels are substrings** (not single characters like in a trie)
5. **Path from root to each leaf spells out one suffix**
6. **No two edges from the same node start with the same character**

### Suffix Tree vs Suffix Trie

| Feature | Suffix Trie | Suffix Tree |
|---------|-------------|-------------|
| Edge labels | Single characters | Substrings |
| Number of nodes | O(n²) | O(n) |
| Space complexity | O(n²) | O(n) |
| Internal nodes | May have 1 child | Always ≥2 children |

## 2. ASCII Art Visualization

### Suffix Tree for "banana$"

First, let's list all suffixes:
```python
Position 0: banana$
Position 1: anana$
Position 2: nana$
Position 3: ana$
Position 4: na$
Position 5: a$
Position 6: $
```python

The suffix tree structure:
```python
                        (root)
                    /    |     \
                   $    a       banana$
                  [6]   |          [0]
                        |
                       (a)
                      /   \
                    $     na
                   [5]     |
                          (na)
                         /    \
                        $     na$
                       [3]    [1]
                        
                        |
                    (also from root)
                        |
                       na
                        |
                       (na)
                      /    \
                     $     na$
                    [4]    [2]

Numbers in brackets [i] = starting position of that suffix
```python

**Full tree structure:**
```python
                            (root)
                    /     /    \      \
                   $    a      banana$  na
                  [6]   |        [0]     |
                       (a)              (na)
                      /   \            /    \
                    $     na          $     na$
                   [5]     |         [4]    [2]
                          (na)
                         /    \
                        $     na$
                       [3]    [1]
```python

### Pattern Search Example: Finding "ana"

```python
Step 1: Start at root, find edge starting with 'a'
        
            (root)
               |
               a    ← follow 'a'
               |
              (a)
             /   \
           $     na    ← follow 'n', then 'a' 
          [5]     |
                 (na)
                /    \
               $     na$
              [3]    [1]

Step 2: After following 'a' from root, we're at node (a)
Step 3: Find edge starting with 'n' → edge "na"
Step 4: Match 'n' and 'a' against edge "na" → complete match!

Pattern "ana" found as: root → 'a' → "na" (first 2 chars)

Step 5: Collect all leaves below current position:
        - Leaf [3] → suffix "ana$" starts at position 3
        - Leaf [1] → suffix "anana$" starts at position 1

Result: Pattern "ana" occurs at positions 1 and 3

Verification: "banana"
              0123456
              - Position 1: "anana" contains "ana" ✓
              - Position 3: "ana" matches exactly ✓
```python

### Longest Common Substring (LCS) Visualization

```python
S1 = "ABAB"
S2 = "BABA"

Combined string: "ABAB$BABA#"
($ terminates S1, # terminates S2)

Build suffix tree and mark leaves:
- Leaves with positions 0-4 belong to S1 (green)
- Leaves with positions 5-9 belong to S2 (red)
- Internal nodes with both green and red descendants are "blue"

Conceptual tree:

                    (root)
                   /  |  |  \
                  $  AB  BA  #
                 [4]  |   |  [9]
                     ...  ...

Find deepest internal node that has descendants from BOTH strings:
- Node with path "ABA" has leaves from both S1 and S2 → "blue"
- Node with path "BAB" has leaves from both S1 and S2 → "blue"

Deepest blue node depth = 3

LCS = "ABA" or "BAB" (both length 3)
```python

## 3. Construction Methods

### Naive Construction - O(n²)

**Algorithm:**
1. Build a suffix trie by inserting all suffixes
2. Compress paths with single children into single edges

**Time Complexity:** O(n²) - inserting all suffixes takes quadratic time

### Ukkonen's Algorithm - O(n)

**Key Ideas:**
1. Build tree incrementally, left to right
2. Use suffix links to avoid redundant traversals
3. Implicit suffix tree becomes explicit with terminal character

**Optimizations:**
- **Suffix links:** Jump between nodes representing related suffixes
- **Edge compression:** Store (start, end) indices instead of substrings
- **Skip/count trick:** Jump over edges when length is known

### Complexity Comparison

| Operation | Naive | Ukkonen |
|-----------|-------|--------|
| Construction | O(n²) | O(n) |
| Space | O(n) | O(n) |
| Pattern search | O(m) | O(m) |
| All k occurrences | O(m + k) | O(m + k) |
| Longest repeated substring | O(n) | O(n) |
| LCS of two strings | O(n + m) | O(n + m) |

## 4. Implementation

### Node and Tree Classes

```python
class SuffixTreeNode:
    """
    A node in the suffix tree.
    
    Attributes:
        children: Dictionary mapping edge labels (strings) to child nodes
        suffix_index: For leaves, the starting position of the suffix (-1 for internal nodes)
        string_id: Identifier for which string this leaf belongs to (for LCS)
    """
    
    def __init__(self):
        self.children = {}      # edge_label -> SuffixTreeNode
        self.suffix_index = -1  # Starting position of suffix (-1 for internal nodes)
        self.string_id = None   # Which string this suffix belongs to (for LCS)
    
    def is_leaf(self):
        """Returns True if this is a leaf node (represents a complete suffix)."""
        return len(self.children) == 0
    
    def __repr__(self):
        if self.is_leaf():
            return f"Leaf[{self.suffix_index}]"
        return f"Internal(children={list(self.children.keys())})"
```python

```python
class SuffixTree:
    """
    Suffix Tree implementation with naive O(n²) construction.
    
    The tree is built by first constructing a suffix trie (all suffixes
    inserted character by character), then compressing single-child paths
    into single edges with substring labels.
    
    Attributes:
        text: The original text (with terminal character)
        root: The root node of the suffix tree
    """
    
    def __init__(self, text, terminal='$'):
        """
        Build a suffix tree for the given text.
        
        Args:
            text: The input string
            terminal: Terminal character to append (default '$')
        """
        self.text = text + terminal
        self.root = SuffixTreeNode()
        self._build_trie()
        self._compress(self.root)
    
    def _build_trie(self):
        """
        Build a suffix trie by inserting all suffixes.
        This is O(n²) as we insert n suffixes of average length n/2.
        """
        n = len(self.text)
        for i in range(n):
            # Insert suffix starting at position i
            current = self.root
            for j in range(i, n):
                char = self.text[j]
                if char not in current.children:
                    current.children[char] = SuffixTreeNode()
                current = current.children[char]
            # Mark the leaf with its starting position
            current.suffix_index = i
    
    def _compress(self, node):
        """
        Compress single-child paths into single edges.
        This converts the trie into a proper suffix tree.
        """
        # Keep compressing until no more single-child internal nodes
        changed = True
        while changed:
            changed = False
            keys = list(node.children.keys())
            for key in keys:
                child = node.children[key]
                # If child has exactly one child, merge the edges
                if len(child.children) == 1:
                    grandchild_key = list(child.children.keys())[0]
                    grandchild = child.children[grandchild_key]
                    # Create new merged edge label
                    new_key = key + grandchild_key
                    node.children[new_key] = grandchild
                    del node.children[key]
                    changed = True
        
        # Recursively compress children
        for child in node.children.values():
            self._compress(child)
    
    def search(self, pattern):
        """
        Search for all occurrences of a pattern in the text.
        
        Args:
            pattern: The pattern to search for
            
        Returns:
            List of starting positions where pattern occurs
        """
        positions = []
        node = self.root
        remaining = pattern
        
        while remaining:
            # Find edge starting with first character of remaining pattern
            found = False
            for edge_label, child in node.children.items():
                if edge_label[0] == remaining[0]:
                    found = True
                    # Check how much of edge matches remaining pattern
                    match_len = min(len(edge_label), len(remaining))
                    if edge_label[:match_len] == remaining[:match_len]:
                        remaining = remaining[match_len:]
                        node = child
                    else:
                        # Mismatch within edge
                        return []
                    break
            
            if not found:
                return []
        
        # Collect all leaf positions below current node
        self._collect_leaves(node, positions)
        return sorted(positions)
    
    def _collect_leaves(self, node, positions):
        """Collect all suffix indices from leaves under this node."""
        if node.is_leaf():
            positions.append(node.suffix_index)
        else:
            for child in node.children.values():
                self._collect_leaves(child, positions)
    
    def visualize(self, node=None, prefix="", is_last=True, edge_label=""):
        """
        Print a visual representation of the suffix tree.
        """
        if node is None:
            node = self.root
            print("(root)")
        else:
            connector = "└── " if is_last else "├── "
            node_info = f"[{node.suffix_index}]" if node.is_leaf() else ""
            print(f"{prefix}{connector}\"{edge_label}\" {node_info}")
        
        children = list(node.children.items())
        for i, (edge, child) in enumerate(children):
            extension = "    " if is_last else "│   "
            self.visualize(child, prefix + extension, i == len(children) - 1, edge)
    
    def get_suffix_array(self):
        """
        Extract suffix array from the tree via lexicographic DFS traversal.
        
        Returns:
            List of suffix starting positions in lexicographic order
        """
        suffix_array = []
        
        def dfs(node):
            if node.is_leaf():
                suffix_array.append(node.suffix_index)
            else:
                # Visit children in lexicographic order of edge labels
                for edge_label in sorted(node.children.keys()):
                    dfs(node.children[edge_label])
        
        dfs(self.root)
        return suffix_array
```python

### Basic Usage Examples

```python
# Build suffix tree for "banana"
text = "banana"
st = SuffixTree(text)

print(f"Suffix Tree for '{text}$':")
print()
st.visualize()
```python

```python
# Pattern search examples
patterns = ["ana", "ban", "nan", "xyz", "a", "na"]

print(f"Pattern searches in '{text}':")
print("=" * 40)
for p in patterns:
    positions = st.search(p)
    if positions:
        print(f"  '{p}' found at positions: {positions}")
    else:
        print(f"  '{p}' not found")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
