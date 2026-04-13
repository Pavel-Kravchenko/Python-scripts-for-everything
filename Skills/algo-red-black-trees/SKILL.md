---
name: algo-red-black-trees
description: "A **Red-Black Tree** is a self-balancing binary search tree where each node has an extra bit for color (RED or BLACK). The tree uses these colors to ensure the tree remains approximately balanced duri"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/05_Tree_Structures/03_red_black_trees.ipynb"
---

# Red-Black Trees

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/05_Tree_Structures/03_red_black_trees.ipynb`*

# Red-Black Trees

A **Red-Black Tree** is a self-balancing binary search tree where each node has an extra bit for color (RED or BLACK). The tree uses these colors to ensure the tree remains approximately balanced during insertions and deletions.

---

## Table of Contents
1. [Red-Black Tree Properties](#1-red-black-tree-properties)
2. [Why Red-Black vs AVL?](#2-why-red-black-vs-avl)
3. [Tree Structure Visualization](#3-tree-structure-visualization)
4. [Implementation](#4-implementation)
5. [Insert with Color Fixing](#5-insert-with-color-fixing)
6. [Rotations with Color Changes](#6-rotations-with-color-changes)
7. [Examples](#7-examples)
8. [Complexity Analysis](#8-complexity-analysis)

---

## 1. Red-Black Tree Properties

A Red-Black Tree must satisfy **5 essential properties**:

| # | Property | Description |
|---|----------|-------------|
| 1 | **Color Property** | Every node is either **RED** or **BLACK** |
| 2 | **Root Property** | The root is always **BLACK** |
| 3 | **Leaf Property** | Every leaf (NIL/null) is **BLACK** |
| 4 | **Red Property** | If a node is **RED**, both its children must be **BLACK** (no red-red parent-child) |
| 5 | **Black-Height Property** | Every path from root to any leaf has the **same number of BLACK nodes** |

### Visual Summary:

```
Property 1: Every node is RED or BLACK
┌─────────────────────────────────────┐
│   (B)20         Each node has a     │
│   /   \         color attribute     │
│ (R)10 (R)30                         │
└─────────────────────────────────────┘

Property 2: Root is always BLACK
┌─────────────────────────────────────┐
│   (B)20  ← Root MUST be BLACK       │
│   /   \                             │
│ (R)10 (R)30                         │
└─────────────────────────────────────┘

Property 3: Leaves (NIL) are BLACK
┌─────────────────────────────────────┐
│       (B)20                         │
│       /   \                         │
│    (R)10  (R)30                     │
│    /  \    /  \                     │
│  NIL  NIL NIL NIL  ← All BLACK      │
└─────────────────────────────────────┘

Property 4: RED node → BLACK children
┌─────────────────────────────────────┐
│ VALID:          INVALID:            │
│   (R)10           (R)10             │
│   /   \           /   \             │
│ (B)5  (B)15    (R)5  (B)15  ✗       │
│                 ↑ RED-RED!          │
└─────────────────────────────────────┘

Property 5: Same BLACK count on all paths
┌─────────────────────────────────────┐
│         (B)20                       │
│         /   \                       │
│      (R)10  (B)30                   │
│      /   \                          │
│   (B)5  (B)15                       │
│                                     │
│ Path 20→10→5→NIL:   BLACK=2 (20,5)  │
│ Path 20→10→15→NIL:  BLACK=2 (20,15) │
│ Path 20→30→NIL:     BLACK=2 (20,30) │
└─────────────────────────────────────┘
```

### Why These Properties Matter

These 5 properties **guarantee** that:
- The longest path from root to leaf is **at most twice** the shortest path
- The tree height is **O(log n)**
- Search, insert, and delete all run in **O(log n)** time

---

## 2. Why Red-Black vs AVL?

Both AVL and Red-Black trees are self-balancing BSTs, but they make different trade-offs:

### Comparison Table

| Aspect | AVL Tree | Red-Black Tree |
|--------|----------|----------------|
| **Balance Factor** | Height difference ≤ 1 | Color-based (looser) |
| **Height Bound** | ≤ 1.44 log(n) | ≤ 2 log(n) |
| **Search Speed** | Faster (shorter height) | Slightly slower |
| **Insert/Delete** | More rotations needed | Fewer rotations |
| **Memory** | +1 int (balance factor) | +1 bit (color) |
| **Best For** | Read-heavy workloads | Write-heavy workloads |
| **Real-World Use** | Databases indexes | Java TreeMap, C++ std::map |

### Visual Comparison

```
Same 7 nodes in AVL vs Red-Black:

AVL Tree (strictly balanced):      Red-Black Tree (loosely balanced):

        4                                    4(B)
       / \                                  /    \
      2   6                              2(R)    6(R)
     / \ / \                            /  \    /  \
    1  3 5  7                        1(B) 3(B) 5(B) 7(B)

Height: 3                            Height: 3 (can be up to 4)
Max rotations per insert: 2          Max rotations per insert: 2
But AVL rotates more often!          Recoloring often sufficient
```

### When to Choose Which?

```
Choose AVL when:
├── Lookups >> Insertions/Deletions
├── Need fastest possible search
└── Example: Dictionary lookup, database indexes

Choose Red-Black when:
├── Many insertions and deletions
├── Need good overall performance
└── Example: Language runtime libraries (Java, C++, Linux kernel)
```

---

## 3. Tree Structure Visualization

### Red-Black Tree Structure

Using **(R)** for RED nodes and **(B)** for BLACK nodes:

```
                    (B)13
                   /     \
               (R)8      (R)17
               /  \      /   \
            (B)1  (B)11 (B)15 (B)25
              \             \    \
             (R)6         (R)22 (R)27
```

### Verifying Properties:

```
✓ Property 1: All nodes are R or B
✓ Property 2: Root (13) is BLACK
✓ Property 3: All NIL leaves are BLACK (implicit)
✓ Property 4: RED nodes (8,17,6,22,27) have BLACK children
✓ Property 5: Black-height = 2 on all paths

Black-height verification:
  13→8→1→NIL      = 2 black (13, 1)
  13→8→1→6→NIL    = 2 black (13, 1)
  13→8→11→NIL     = 2 black (13, 11)
  13→17→15→NIL    = 2 black (13, 15)
  13→17→25→22→NIL = 2 black (13, 25)
  13→17→25→27→NIL = 2 black (13, 25)
```

---

## 4. Implementation

Let's implement a Red-Black Tree with full insertion support.

```python
# Color constants
RED = True
BLACK = False


class RBNode:
    """Node for Red-Black Tree with color attribute."""
    
    def __init__(self, data, color=RED):
        self.data = data
        self.color = color  # New nodes are RED by default
        self.parent = None
        self.left = None
        self.right = None
    
    def __repr__(self):
        color_str = "R" if self.color == RED else "B"
        return f"({color_str}){self.data}"
    
    def is_red(self):
        return self.color == RED
    
    def is_black(self):
        return self.color == BLACK
```

```python
class RedBlackTree:
    """Red-Black Tree implementation with insert and fix-up."""
    
    def __init__(self):
        # NIL sentinel node (always BLACK)
        self.NIL = RBNode(data=None, color=BLACK)
        self.root = self.NIL
    
    # ==================== Helper Methods ====================
    
    def is_nil(self, node):
        """Check if node is NIL sentinel."""
        return node == self.NIL or node is None
    
    def is_root(self, node):
        """Check if node is the root."""
        return node.parent is None or node.parent == self.NIL
    
    def is_left_child(self, node):
        """Check if node is a left child."""
        return node.parent and node == node.parent.left
    
    def is_right_child(self, node):
        """Check if node is a right child."""
        return node.parent and node == node.parent.right
    
    def get_uncle(self, node):
        """Get the uncle of a node (parent's sibling)."""
        if node.parent is None or node.parent.parent is None:
            return None
        grandparent = node.parent.parent
        if self.is_left_child(node.parent):
            return grandparent.right
        else:
            return grandparent.left
    
    def get_grandparent(self, node):
        """Get the grandparent of a node."""
        if node.parent:
            return node.parent.parent
        return None
    
    # ==================== Rotation Methods ====================
    
    def rotate_left(self, x):
        """
        Left rotation around node x.
        
              x                    y
             / \                  / \
            a   y      -->       x   c
               / \              / \
              b   c            a   b
        """
        y = x.right
        
        # Turn y's left subtree into x's right subtree
        x.right = y.left
        if not self.is_nil(y.left):
            y.left.parent = x
        
        # Link x's parent to y
        y.parent = x.parent
        if x.parent is None:
            self.root = y
        elif self.is_left_child(x):
            x.parent.left = y
        else:
            x.parent.right = y
        
        # Put x on y's left
        y.left = x
        x.parent = y
        
        return y
    
    def rotate_right(self, y):
        """
        Right rotation around node y.
        
              y                  x
             / \                / \
            x   c     -->      a   y
           / \                    / \
          a   b                  b   c
        """
        x = y.left
        
        # Turn x's right subtree into y's left subtree
        y.left = x.right
        if not self.is_nil(x.right):
            x.right.parent = y
        
        # Link y's parent to x
        x.parent = y.parent
        if y.parent is None:
            self.root = x
        elif self.is_right_child(y):
            y.parent.right = x
        else:
            y.parent.left = x
        
        # Put y on x's right
        x.right = y
        y.parent = x
        
        return x
    
    # ==================== Insert Methods ====================
    
    def insert(self, data):
        """Insert a new value into the Red-Black Tree."""
        # Create new RED node
        new_node = RBNode(data, color=RED)
        new_node.left = self.NIL
        new_node.right = self.NIL
        
        # Standard BST insert
        parent = None
        current = self.root
        
        while not self.is_nil(current):
            parent = current
            if data < current.data:
                current = current.left
            elif data > current.data:
                current = current.right
            else:
                return  # Duplicate, don't insert
        
        new_node.parent = parent
        
        if parent is None:
            self.root = new_node
        elif data < parent.data:
            parent.left = new_node
        else:
            parent.right = new_node
        
        # Fix Red-Black properties
        self._fix_insert(new_node)
    
    def _fix_insert(self, node):
        """
        Fix Red-Black Tree properties after insertion.
        Handles all three cases.
        """
        while node != self.root and node.parent.is_red():
            uncle = self.get_uncle(node)
            grandparent = self.get_grandparent(node)
            
            if self.is_left_child(node.parent):
                # Parent is left child of grandparent
                if uncle and uncle.is_red():
                    # Case 1: Uncle is RED → Recolor
                    node.parent.color = BLACK
                    uncle.color = BLACK
                    grandparent.color = RED
                    node = grandparent  # Move up to check grandparent
                else:
                    # Uncle is BLACK
                    if self.is_right_child(node):
                        # Case 2: Node is right child → Left rotate parent
                        node = node.parent
                        self.rotate_left(node)
                    # Case 3: Node is left child → Right rotate grandparent
                    node.parent.color = BLACK
                    grandparent.color = RED
                    self.rotate_right(grandparent)
            else:
                # Parent is right child of grandparent (mirror cases)
                if uncle and uncle.is_red():
                    # Case 1: Uncle is RED → Recolor
                    node.parent.color = BLACK
                    uncle.color = BLACK
                    grandparent.color = RED
                    node = grandparent
                else:
                    # Uncle is BLACK
                    if self.is_left_child(node):
                        # Case 2: Node is left child → Right rotate parent
                        node = node.parent
                        self.rotate_right(node)
                    # Case 3: Node is right child → Left rotate grandparent
                    node.parent.color = BLACK
                    grandparent.color = RED
                    self.rotate_left(grandparent)
        
        # Ensure root is always BLACK
        self.root.color = BLACK
    
    # ==================== Search Method ====================
    
    def search(self, data):
        """Search for a value in the tree."""
        current = self.root
        while not self.is_nil(current):
            if data == current.data:
                return current
            elif data < current.data:
                current = current.left
            else:
                current = current.right
        return None
    
    # ==================== Traversal Methods ====================
    
    def inorder(self):
        """Return inorder traversal as list."""
        result = []
        self._inorder_helper(self.root, result)
        return result
    
    def _inorder_helper(self, node, result):
        if not self.is_nil(node):
            self._inorder_helper(node.left, result)
            result.append(node)
            self._inorder_helper(node.right, result)
    
    # ==================== Verification Methods ====================
    
    def verify_properties(self):
        """Verify all Red-Black Tree properties."""
        if self.is_nil(self.root):
            return True, "Empty tree is valid"
        
        errors = []
        
        # Property 2: Root is BLACK
        if self.root.is_red():
            errors.append("Property 2 violated: Root is RED")
        
        # Property 4: No red-red & Property 5: Black-height
        black_heights = []
        self._check_properties(self.root, 0, black_heights, errors)
        
        # Check all black heights are equal
        if len(set(black_heights)) > 1:
            errors.append(f"Property 5 violated: Unequal black heights {black_heights}")
        
        if errors:
            return False, errors
        return True, f"Valid Red-Black Tree (black-height: {black_heights[0] if black_heights else 0})"
    
    def _check_properties(self, node, black_count, black_heights, errors):
        """Recursively check properties 4 and 5."""
        if self.is_nil(node):
            black_heights.append(black_count)
            return
        
        # Property 4: Check for red-red violation
        if node.is_red():
            if (not self.is_nil(node.left) and node.left.is_red()):
                errors.append(f"Property 4 violated: RED-RED at {node} → {node.left}")
            if (not self.is_nil(node.right) and node.right.is_red()):
                errors.append(f"Property 4 violated: RED-RED at {node} → {node.right}")
        
        # Count black nodes
        new_count = black_count + (1 if node.is_black() else 0)
        self._check_properties(node.left, new_count, black_heights, errors)
        self._check_properties(node.right, new_count, black_heights, errors)
    
    def black_height(self):
        """Calculate the black-height of the tree."""
        count = 0
        node = self.root
        while not self.is_nil(node):
            if node.is_black():
                count += 1
            node = node.left
        return count
    
    # ==================== Display Methods ====================
    
    def print_tree(self):
        """Print the tree structure with colors."""
        if self.is_nil(self.root):
            print("<Empty Tree>")
            return
        self._print_tree_helper(self.root, "", True)
    
    def _print_tree_helper(self, node, prefix, is_last):
        """Helper for printing tree structure."""
        if self.is_nil(node):
            return
        
        connector = "└── " if is_last else "├── "
        color = "R" if node.is_red() else "B"
        print(f"{prefix}{connector}({color}){node.data}")
        
        new_prefix = prefix + ("    " if is_last else "│   ")
        
        # Check if we have children to print
        has_left = not self.is_nil(node.left)
        has_right = not self.is_nil(node.right)
        
        if has_left or has_right:
            if has_right:
                self._print_tree_helper(node.right, new_prefix, False)
            if has_left:
                self._print_tree_helper(node.left, new_prefix, True)
```
