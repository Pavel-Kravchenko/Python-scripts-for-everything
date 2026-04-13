---
name: algo-tries
description: "A **Trie** (pronounced "try") is a tree-like data structure used for efficient storage and retrieval of strings. The name comes from "re**TRIE**val" - it was designed for fast information retrieval."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/01_tries.ipynb"
---

# Trie (Prefix Tree)

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/01_tries.ipynb`*

# Trie (Prefix Tree)

## What is a Trie?

A **Trie** (pronounced "try") is a tree-like data structure used for efficient storage and retrieval of strings. The name comes from "re**TRIE**val" - it was designed for fast information retrieval.

### Key Properties:
- Each node represents a **single character**
- The **path** from root to any node represents a **prefix**
- Words are marked with a special **end-of-word flag** (`is_end`)
- All descendants of a node share a common prefix

### Why Use a Trie?
- **Fast prefix operations**: O(m) where m is the word/prefix length
- **Efficient autocomplete**: Find all words with a given prefix
- **Memory efficient for shared prefixes**: "cat", "car", "card" share "ca"

---

## Trie Structure Visualization

**Trie containing words: "cat", "car", "card", "dog", "do"**

```
           (root)
          /      \
         c        d
         |        |
         a        o
        / \       |*
       t*  r*     g*
           |
           d*

* = end of word (is_end = True)

Path "c-a-t" = "cat"
Path "c-a-r" = "car"
Path "c-a-r-d" = "card"
Path "d-o" = "do"
Path "d-o-g" = "dog"
```

Notice how:
- "cat" and "car" share the prefix "ca"
- "car" and "card" share the prefix "car"
- "do" and "dog" share the prefix "do"

---

## Trie Operations

### 1. Insert Operation

**Insert "cats" into the trie:**

```
Before:              After:
    (root)              (root)
       |                   |
       c                   c
       |                   |
       a                   a
      / \                 / \
     t*  r*              t*  r*
                         |
                         s*
```

**Algorithm:**
1. Start at root
2. For each character in word:
   - If child exists → move to it
   - Else → create new node, then move to it
3. Mark final node as end-of-word

---

### 2. Search Operation

**Search for "car":**

```
    (root)
       ↓ 'c'
       c
       ↓ 'a'
       a
      / ↓ 'r'
     t   r* ← is_end=True → FOUND!
```

**Algorithm:**
1. Start at root
2. For each character:
   - If child doesn't exist → return False
   - Else → move to child
3. Return True only if `is_end` is True

---

### 3. Prefix Search

**Find all words with prefix "ca":**

```
    (root)
       ↓
       c
       ↓
       a  ← reached prefix
      / \
     t*  r*  ← collect all words from here
         |
         d*

Result: ["cat", "car", "card"]
```

**Algorithm:**
1. Navigate to prefix end node
2. DFS/BFS to collect all words with `is_end=True`

---

### 4. Delete Operation

**Delete "car" (but keep "card"):**

```
Before:              After:
       a                   a
      / \                 / \
     t*  r*              t*  r   ← is_end=False now
         |                   |
         d*                  d*
```

Just unmark `is_end` - don't remove nodes that have children!

---

## Complexity Analysis

| Operation | Time Complexity | Space Complexity |
|-----------|-----------------|------------------|
| Insert | O(m) | O(m) worst case |
| Search | O(m) | O(1) |
| Prefix Search | O(p + k) | O(k) for results |
| Delete | O(m) | O(1) |
| Build from n words | O(n × m) | O(n × m) worst |

Where:
- **m** = length of the word
- **p** = length of the prefix
- **k** = number of matching results

### Comparison with Hash Table

| Aspect | Trie | Hash Table |
|--------|------|------------|
| Search | O(m) | O(m) average |
| Prefix search | O(p + k) | O(n × m) |
| Sorted iteration | Natural | Requires sort |
| Memory | Can be higher | Usually lower |

---

## Implementation

```python
class TrieNode:
    """
    A node in the Trie structure.
    
    Attributes:
        children (dict): Maps characters to child TrieNode objects
        is_end (bool): True if this node marks the end of a valid word
    
    Based on Node class pattern from Algorithms_HW9.ipynb
    """
    
    def __init__(self):
        self.children = {}  # char -> TrieNode
        self.is_end = False  # marks end of word
    
    def __repr__(self):
        return f"TrieNode(children={list(self.children.keys())}, is_end={self.is_end})"
```

```python
class Trie:
    """
    A Trie (prefix tree) data structure for efficient string operations.
    
    Supports:
        - Insert words
        - Search for exact words
        - Check if prefix exists
        - Get all words with a given prefix
        - Delete words
        - Count words with prefix
    
    Example:
        >>> trie = Trie()
        >>> trie.insert("hello")
        >>> trie.search("hello")
        True
        >>> trie.starts_with("hel")
        True
    """
    
    def __init__(self):
        """Initialize an empty Trie with a root node."""
        self.root = TrieNode()
    
    def insert(self, word: str) -> None:
        """
        Insert a word into the trie.
        
        Args:
            word: The string to insert
            
        Time Complexity: O(m) where m is word length
        Space Complexity: O(m) worst case (no shared prefix)
        """
        node = self.root
        for char in word:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
        node.is_end = True
    
    def search(self, word: str) -> bool:
        """
        Check if word exists in the trie.
        
        Args:
            word: The string to search for
            
        Returns:
            True if the exact word exists, False otherwise
            
        Time Complexity: O(m)
        """
        node = self._find_node(word)
        return node is not None and node.is_end
    
    def starts_with(self, prefix: str) -> bool:
        """
        Check if any word in the trie starts with the given prefix.
        
        Args:
            prefix: The prefix string to check
            
        Returns:
            True if any word starts with prefix, False otherwise
            
        Time Complexity: O(p) where p is prefix length
        """
        return self._find_node(prefix) is not None
    
    def _find_node(self, prefix: str) -> TrieNode:
        """
        Navigate to the node representing the end of prefix.
        
        Args:
            prefix: The prefix to navigate to
            
        Returns:
            The TrieNode at end of prefix, or None if prefix doesn't exist
        """
        node = self.root
        for char in prefix:
            if char not in node.children:
                return None
            node = node.children[char]
        return node
    
    def get_all_with_prefix(self, prefix: str) -> list:
        """
        Get all words in the trie that start with the given prefix.
        
        Args:
            prefix: The prefix to search for
            
        Returns:
            List of all words starting with the prefix
            
        Time Complexity: O(p + k) where p is prefix length, k is number of results
        """
        results = []
        node = self._find_node(prefix)
        
        if node is None:
            return results
        
        # DFS to collect all words from this node
        self._collect_words(node, prefix, results)
        return results
    
    def _collect_words(self, node: TrieNode, current_word: str, results: list) -> None:
        """
        Recursively collect all words from a given node using DFS.
        
        Args:
            node: Current TrieNode
            current_word: The word built so far
            results: List to append found words to
        """
        if node.is_end:
            results.append(current_word)
        
        for char, child in sorted(node.children.items()):
            self._collect_words(child, current_word + char, results)
    
    def delete(self, word: str) -> bool:
        """
        Delete a word from the trie.
        
        Args:
            word: The word to delete
            
        Returns:
            True if word was deleted, False if word didn't exist
            
        Time Complexity: O(m)
        """
        def _delete_helper(node: TrieNode, word: str, depth: int) -> bool:
            if depth == len(word):
                if not node.is_end:
                    return False  # Word doesn't exist
                node.is_end = False
                return len(node.children) == 0  # Can delete if no children
            
            char = word[depth]
            if char not in node.children:
                return False  # Word doesn't exist
            
            should_delete_child = _delete_helper(node.children[char], word, depth + 1)
            
            if should_delete_child:
                del node.children[char]
                return len(node.children) == 0 and not node.is_end
            
            return False
        
        _delete_helper(self.root, word, 0)
        return True
    
    def count_words_with_prefix(self, prefix: str) -> int:
        """
        Count how many words start with the given prefix.
        
        Args:
            prefix: The prefix to count
            
        Returns:
            Number of words starting with prefix
        """
        node = self._find_node(prefix)
        if node is None:
            return 0
        
        count = [0]
        
        def _count_words(node: TrieNode):
            if node.is_end:
                count[0] += 1
            for child in node.children.values():
                _count_words(child)
        
        _count_words(node)
        return count[0]
    
    def get_all_words(self) -> list:
        """
        Get all words stored in the trie.
        
        Returns:
            List of all words in alphabetical order
        """
        return self.get_all_with_prefix("")
```

---

## Trie Visualization Helper

```python
def visualize_trie(trie: Trie) -> None:
    """
    Print a visual representation of the trie structure.
    
    Args:
        trie: The Trie object to visualize
    """
    def _visualize(node: TrieNode, prefix: str, is_last: bool, indent: str):
        # Determine the connector
        connector = "└── " if is_last else "├── "
        
        # Print current node
        if prefix:  # Don't print for root
            end_marker = "*" if node.is_end else ""
            print(f"{indent}{connector}{prefix[-1]}{end_marker}")
        else:
            print("(root)")
        
        # Update indent for children
        if prefix:
            indent += "    " if is_last else "│   "
        
        # Process children
        children = sorted(node.children.items())
        for i, (char, child) in enumerate(children):
            is_last_child = (i == len(children) - 1)
            _visualize(child, prefix + char, is_last_child, indent)
    
    _visualize(trie.root, "", True, "")
    print("\n* = end of word")
```

### Example 1: Basic Operations

```python
# Create a trie and insert words
trie = Trie()
words = ["cat", "car", "card", "dog", "do"]

for word in words:
    trie.insert(word)
    print(f"Inserted: '{word}'")

print("\nTrie structure:")
visualize_trie(trie)
```

```python
# Search operations
print("Search operations:")
print(f"  search('cat')  = {trie.search('cat')}   # exact word exists")
print(f"  search('ca')   = {trie.search('ca')}  # prefix only, not a word")
print(f"  search('dogs') = {trie.search('dogs')}  # word doesn't exist")

print("\nPrefix operations:")
print(f"  starts_with('ca')  = {trie.starts_with('ca')}   # prefix exists")
print(f"  starts_with('xyz') = {trie.starts_with('xyz')}  # prefix doesn't exist")
```

### Example 2: Autocomplete

```python
# Build autocomplete trie
autocomplete_trie = Trie()
dictionary = [
    "apple", "application", "apply", "apt",
    "banana", "band", "bandana",
    "cat", "car", "card", "care", "careful", "careless",
    "dog", "door", "doom"
]

for word in dictionary:
    autocomplete_trie.insert(word)

print("Dictionary loaded with", len(dictionary), "words")
print("\nTrie structure:")
visualize_trie(autocomplete_trie)
```

```python
def autocomplete(trie: Trie, prefix: str, max_results: int = 5) -> list:
    """
    Get autocomplete suggestions for a prefix.
    
    Args:
        trie: The Trie containing words
        prefix: User input to autocomplete
        max_results: Maximum suggestions to return
        
    Returns:
        List of suggested completions
    """
    suggestions = trie.get_all_with_prefix(prefix)
    return suggestions[:max_results]

# Test autocomplete
test_prefixes = ["app", "car", "do", "b", "xyz"]

print("Autocomplete Suggestions:")
print("=" * 40)
for prefix in test_prefixes:
    suggestions = autocomplete(autocomplete_trie, prefix)
    print(f"  '{prefix}' → {suggestions}")
```
