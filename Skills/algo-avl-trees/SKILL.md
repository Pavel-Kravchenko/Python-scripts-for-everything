---
name: algo-avl-trees
description: "Named after inventors **A**delson-**V**elsky and **L**andis (1962) - the first self-balancing BST."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/05_Tree_Structures/02_avl_trees.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# AVL Trees

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/05_Tree_Structures/02_avl_trees.ipynb`*

# AVL Trees

**Self-Balancing Binary Search Trees**

Named after inventors **A**delson-**V**elsky and **L**andis (1962) - the first self-balancing BST.

---

## 1. Why Self-Balancing Trees?

### The BST Worst Case Problem

A regular Binary Search Tree can degenerate into a **linked list** when elements are inserted in sorted order:

```
Insert: 1, 2, 3, 4, 5

    (1)                 Ideal BST:      (3)
      \                               /     \
      (2)                          (1)       (4)
        \                            \         \
        (3)                          (2)       (5)
          \
          (4)           Height: 2 (balanced)
            \           Operations: O(log n)
            (5)

Height: 4 (degenerate)
Operations: O(n)
```

### Performance Comparison

| Tree Type | Best Case | Average Case | Worst Case |
|-----------|-----------|--------------|------------|
| Regular BST | O(log n) | O(log n) | **O(n)** |
| AVL Tree | O(log n) | O(log n) | **O(log n)** |

**AVL trees guarantee O(log n) performance by maintaining balance after every insertion and deletion.**

## 2. The AVL Property

### Balance Factor

For any node in an AVL tree:

$$\text{Balance Factor} = \text{height(left subtree)} - \text{height(right subtree)}$$

### AVL Invariant

**For every node: $|\text{balance factor}| \leq 1$**

Valid balance factors: **-1, 0, +1**

```
Balance Factor Visualization:

        (8) bf=-1          bf = height(left) - height(right)
       /   \                  = 1 - 2 = -1
    (3)     (10) bf=-1
   bf=0        \
               (14) bf=0

Heights (using h(leaf)=1, h(None)=0):
- Node 14: leaf         → height = 1,  bf = 0 - 0 = 0
- Node 10: right child  → height = 2,  bf = 0 - 1 = -1
- Node 3:  leaf         → height = 1,  bf = 0 - 0 = 0
- Node 8:  root         → height = 3,  bf = h(left) - h(right) = 1 - 2 = -1

All |bf| ≤ 1 → valid AVL tree
```

### Why Balance Factor ≤ 1 Guarantees O(log n) Height

Let $N_h$ = minimum nodes in AVL tree of height $h$

- $N_0 = 1$ (single node)
- $N_1 = 2$ (root + one child)
- $N_h = N_{h-1} + N_{h-2} + 1$ (Fibonacci-like recurrence)

This gives: $N_h > \phi^h$ where $\phi = \frac{1+\sqrt{5}}{2} \approx 1.618$

Therefore: $h < 1.44 \cdot \log_2(n+2)$

**Height is always O(log n)!**

## 3. Rotations

Rotations are the fundamental operations that restore balance. They:
- Maintain BST ordering property
- Execute in **O(1)** time
- Reduce height of unbalanced subtree

---

### 3.1 Right Rotation (LL Case)

**When to use:** Node has `balance_factor = +2` AND left child has `balance_factor ≥ 0`

```
         z                                y
        / \                             /   \
       y   T4      Right Rotate        x     z
      / \          ─────────────►     / \   / \
     x   T3                          T1 T2 T3 T4
    / \
   T1  T2

Steps:
1. y becomes the new root
2. z becomes y's right child
3. y's original right subtree (T3) becomes z's left subtree

BST Property preserved: T1 < x < T2 < y < T3 < z < T4
```

---

### 3.2 Left Rotation (RR Case)

**When to use:** Node has `balance_factor = -2` AND right child has `balance_factor ≤ 0`

```
     z                                   y
    / \                                /   \
   T1  y         Left Rotate          z     x
      / \        ─────────────►      / \   / \
     T2  x                          T1 T2 T3 T4
        / \
       T3  T4

Steps:
1. y becomes the new root
2. z becomes y's left child
3. y's original left subtree (T2) becomes z's right subtree

BST Property preserved: T1 < z < T2 < y < T3 < x < T4
```

---

### 3.3 Left-Right Rotation (LR Case)

**When to use:** Node has `balance_factor = +2` AND left child has `balance_factor = -1`

```
       z                    z                        x
      / \                  / \                     /   \
     y   T4               x   T4                  y     z
    / \      Left        / \        Right        / \   / \
   T1  x     Rotate     y   T3      Rotate      T1 T2 T3 T4
      / \    on y      / \          on z
     T2  T3           T1  T2

Step 1: Left rotate the left child (y)
Step 2: Right rotate the unbalanced node (z)

BST Property preserved: T1 < y < T2 < x < T3 < z < T4
```

---

### 3.4 Right-Left Rotation (RL Case)

**When to use:** Node has `balance_factor = -2` AND right child has `balance_factor = +1`

```
     z                   z                          x
    / \                 / \                       /   \
   T1  y               T1  x                     z     y
      / \    Right        / \       Left        / \   / \
     x   T4  Rotate      T2  y      Rotate     T1 T2 T3 T4
    / \      on y           / \     on z
   T2  T3                  T3  T4

Step 1: Right rotate the right child (y)
Step 2: Left rotate the unbalanced node (z)

BST Property preserved: T1 < z < T2 < x < T3 < y < T4
```

### Rotation Decision Table

| Node bf | Child bf | Case | Rotation |
|---------|----------|------|----------|
| +2 | +1 or 0 | LL (Left-Left) | Right rotation |
| +2 | -1 | LR (Left-Right) | Left-Right rotation |
| -2 | -1 or 0 | RR (Right-Right) | Left rotation |
| -2 | +1 | RL (Right-Left) | Right-Left rotation |

**Memory trick:** The rotation is named opposite to the imbalance direction!

## 4. Insertion with Rebalancing

### Algorithm

1. **Standard BST Insert:** Insert the new node at proper position
2. **Update Heights:** Walk back up, updating heights of ancestors
3. **Check Balance:** At each ancestor, compute balance factor
4. **Rebalance if needed:** Apply appropriate rotation(s)

### Insertion Example: Insert 30, 20, 10

```
Step 1: Insert 30          Step 2: Insert 20         Step 3: Insert 10

    (30) bf=0                 (30) bf=+1                (30) bf=+2  ← UNBALANCED!
                              /                        /
                           (20) bf=0                (20) bf=+1
                                                    /
                                                 (10) bf=0

Detected: LL Case (bf=+2 at 30, bf=+1 at 20)
Action: Right Rotation on node 30

After Right Rotation:

        (20) bf=0
       /    \
    (10)    (30)
    bf=0    bf=0

Tree is now balanced!
```

### Another Example: Insert 30, 10, 20 (LR Case)

```
Step 1: Insert 30          Step 2: Insert 10         Step 3: Insert 20

    (30) bf=0                 (30) bf=+1                (30) bf=+2  ← UNBALANCED!
                              /                        /
                           (10) bf=0                (10) bf=-1
                                                       \
                                                       (20) bf=0

Detected: LR Case (bf=+2 at 30, bf=-1 at 10)
Action: Left-Right Rotation

Step A: Left rotate node 10      Step B: Right rotate node 30

        (30) bf=+2                      (20) bf=0
        /                              /    \
     (20) bf=+1                     (10)    (30)
     /                              bf=0    bf=0
  (10) bf=0

Tree is now balanced!
```

## 5. Deletion with Rebalancing

### Algorithm

1. **Standard BST Delete:** Find and remove the node
   - Leaf: Simply remove
   - One child: Replace with child
   - Two children: Replace with inorder successor/predecessor, then delete that
2. **Update Heights:** Walk back up from deletion point
3. **Check Balance:** At each ancestor, compute balance factor
4. **Rebalance if needed:** May need **multiple rotations** (unlike insertion!)

### Deletion Example

```
Initial Tree:                  After deleting 40:

        (30) bf=0                    (30) bf=+2  ← UNBALANCED!
       /    \                       /
    (20)    (40)                 (20) bf=0
   /  \      bf=0               /  \
(10)  (25)                   (10)  (25)
bf=0  bf=0                   bf=0  bf=0

bf at 30: left_height=2, right_height=0 → bf = +2
bf at 20: left_height=1, right_height=1 → bf = 0

Case: LL (bf=+2 at 30, bf=0 at 20)
Action: Right Rotation on 30

After Right Rotation:

        (20) bf=0
       /    \
    (10)    (30)
    bf=0   /
        (25)
        bf=0

Tree is now balanced!
```

### Key Difference from Insertion

- **Insertion:** At most ONE rotation fixes the tree
- **Deletion:** May need rotations at MULTIPLE ancestors (up to O(log n) rotations)

## 6. Complexity Analysis

| Operation | Time Complexity | Space Complexity |
|-----------|-----------------|------------------|
| Search | O(log n) | O(1) iterative, O(log n) recursive |
| Insert | O(log n) | O(log n)* |
| Delete | O(log n) | O(log n)* |
| Single Rotation | O(1) | O(1) |
| Get Min/Max | O(log n) | O(1) |
| Get Height | O(1) | O(1) |
| Get Balance | O(1) | O(1) |

*Stack space for recursion

### Space Overhead

Each node stores:
- Value
- Left pointer
- Right pointer  
- **Height** (additional field for AVL)

Total space: **O(n)**

---

## Implementation

Based on the AVL implementation from `origin/Algorithms_HW5.ipynb`

```python
from typing import Optional, List, Any


class AVLNode:
    """A node in an AVL tree.
    
    Attributes:
        value: The data stored in the node
        left: Reference to left child
        right: Reference to right child
        height: Height of the node (leaf = 1)
    """
    
    def __init__(self, value: Any) -> None:
        self.value = value
        self.left: Optional['AVLNode'] = None
        self.right: Optional['AVLNode'] = None
        self.height: int = 1
    
    def __repr__(self) -> str:
        return f"AVLNode({self.value})"
```

```python
class AVLTree:
    """AVL Tree - Self-balancing Binary Search Tree.
    
    Maintains balance factor |bf| ≤ 1 for all nodes,
    guaranteeing O(log n) operations.
    
    Reference: origin/Algorithms_HW5.ipynb (AVL_Tree class)
    """
    
    def __init__(self) -> None:
        self.root: Optional[AVLNode] = None
    
    # ==================== Helper Methods ====================
    
    def get_height(self, node: Optional[AVLNode]) -> int:
        """Get height of a node. Empty node has height 0.
        
        Args:
            node: The node to get height of
            
        Returns:
            Height of the node (0 if None)
        """
        if node is None:
            return 0
        return node.height
    
    def get_balance(self, node: Optional[AVLNode]) -> int:
        """Get balance factor of a node.
        
        Balance Factor = height(left) - height(right)
        
        Args:
            node: The node to calculate balance factor for
            
        Returns:
            Balance factor (-1, 0, or 1 for balanced nodes)
        """
        if node is None:
            return 0
        return self.get_height(node.left) - self.get_height(node.right)
    
    def _update_height(self, node: AVLNode) -> None:
        """Update height of a node based on children's heights."""
        node.height = 1 + max(
            self.get_height(node.left),
            self.get_height(node.right)
        )
    
    # ==================== Rotations ====================
    
    def rotate_right(self, z: AVLNode) -> AVLNode:
        """Perform right rotation on node z.
        
             z                 y
            / \               / \
           y   T4    →       x   z
          / \               / \ / \
         x   T3            T1 T2 T3 T4
        / \
       T1  T2
       
        Args:
            z: The unbalanced node to rotate
            
        Returns:
            The new root of this subtree (y)
        """
        y = z.left
        T3 = y.right
        
        # Perform rotation
        y.right = z
        z.left = T3
        
        # Update heights (z first, then y since y is now parent)
        self._update_height(z)
        self._update_height(y)
        
        return y
    
    def rotate_left(self, z: AVLNode) -> AVLNode:
        """Perform left rotation on node z.
        
           z                    y
          / \                  / \
         T1  y       →        z   x
            / \              / \ / \
           T2  x            T1 T2 T3 T4
              / \
             T3  T4
             
        Args:
            z: The unbalanced node to rotate
            
        Returns:
            The new root of this subtree (y)
        """
        y = z.right
        T2 = y.left
        
        # Perform rotation
        y.left = z
        z.right = T2
        
        # Update heights
        self._update_height(z)
        self._update_height(y)
        
        return y
    
    # ==================== Insert ====================
    
    def insert(self, value: Any) -> None:
        """Insert a value into the AVL tree.
        
        Args:
            value: The value to insert
        """
        self.root = self._insert(self.root, value)
    
    def _insert(self, node: Optional[AVLNode], value: Any) -> AVLNode:
        """Recursive insert with rebalancing.
        
        Algorithm:
        1. Standard BST insert
        2. Update height of ancestor
        3. Get balance factor
        4. Rebalance if needed (4 cases)
        
        Args:
            node: Current node in recursion
            value: Value to insert
            
        Returns:
            Root of the (possibly rebalanced) subtree
        """
        # Step 1: Standard BST insert
        if node is None:
            return AVLNode(value)
        
        if value < node.value:
            node.left = self._insert(node.left, value)
        elif value > node.value:
            node.right = self._insert(node.right, value)
        else:
            # Duplicate values not allowed
            return node
        
        # Step 2: Update height of this ancestor
        self._update_height(node)
        
        # Step 3: Get balance factor
        balance = self.get_balance(node)
        
        # Step 4: Rebalance if needed (4 cases)
        
        # Case 1: Left-Left (LL)
        # bf = +2 and new node in left subtree of left child
        if balance > 1 and value < node.left.value:
            return self.rotate_right(node)
        
        # Case 2: Right-Right (RR)
        # bf = -2 and new node in right subtree of right child
        if balance < -1 and value > node.right.value:
            return self.rotate_left(node)
        
        # Case 3: Left-Right (LR)
        # bf = +2 and new node in right subtree of left child
        if balance > 1 and value > node.left.value:
            node.left = self.rotate_left(node.left)
            return self.rotate_right(node)
        
        # Case 4: Right-Left (RL)
        # bf = -2 and new node in left subtree of right child
        if balance < -1 and value < node.right.value:
            node.right = self.rotate_right(node.right)
            return self.rotate_left(node)
        
        return node
    
    # ==================== Delete ====================
    
    def delete(self, value: Any) -> None:
        """Delete a value from the AVL tree.
        
        Args:
            value: The value to delete
        """
        self.root = self._delete(self.root, value)
    
    def _get_min_node(self, node: AVLNode) -> AVLNode:
        """Get the node with minimum value in subtree."""
        current = node
        while current.left is not None:
            current = current.left
        return current
    
    def _delete(self, node: Optional[AVLNode], value: Any) -> Optional[AVLNode]:
        """Recursive delete with rebalancing.
        
        Algorithm:
        1. Standard BST delete
        2. Update height
        3. Get balance factor
        4. Rebalance if needed
        
        Note: Unlike insertion, deletion may require
        multiple rotations up the tree.
        
        Args:
            node: Current node in recursion
            value: Value to delete
            
        Returns:
            Root of the (possibly rebalanced) subtree
        """
        # Step 1: Standard BST delete
        if node is None:
            return None
        
        if value < node.value:
            node.left = self._delete(node.left, value)
        elif value > node.value:
            node.right = self._delete(node.right, value)
        else:
            # Found the node to delete
            
            # Case 1 & 2: Node with one or no child
            if node.left is None:
                return node.right
            elif node.right is None:
                return node.left
            
            # Case 3: Node with two children
            # Get inorder successor (smallest in right subtree)
            successor = self._get_min_node(node.right)
            node.value = successor.value
            node.right = self._delete(node.right, successor.value)
        
        # If tree had only one node
        if node is None:
            return None
        
        # Step 2: Update height
        self._update_height(node)
        
        # Step 3: Get balance factor
        balance = self.get_balance(node)
        
        # Step 4: Rebalance if needed
        
        # Left-Left Case
        if balance > 1 and self.get_balance(node.left) >= 0:
            return self.rotate_right(node)
        
        # Left-Right Case
        if balance > 1 and self.get_balance(node.left) < 0:
            node.left = self.rotate_left(node.left)
            return self.rotate_right(node)
        
        # Right-Right Case
        if balance < -1 and self.get_balance(node.right) <= 0:
            return self.rotate_left(node)
        
        # Right-Left Case
        if balance < -1 and self.get_balance(node.right) > 0:
            node.right = self.rotate_right(node.right)
            return self.rotate_left(node)
        
        return node
    
    # ==================== Search ====================
    
    def search(self, value: Any) -> Optional[AVLNode]:
        """Search for a value in the tree.
        
        Args:
            value: The value to search for
            
        Returns:
            The node containing the value, or None if not found
        """
        node = self.root
        while node is not None:
            if value == node.value:
                return node
            elif value < node.value:
                node = node.left
            else:
                node = node.right
        return None
    
    # ==================== Traversals ====================
    
    def inorder(self) -> List[Any]:
        """Return inorder traversal (sorted order)."""
        result = []
        self._inorder(self.root, result)
        return result
    
    def _inorder(self, node: Optional[AVLNode], result: List[Any]) -> None:
        if node:
            self._inorder(node.left, result)
            result.append(node.value)
            self._inorder(node.right, result)
    
    def preorder(self) -> List[Any]:
        """Return preorder traversal."""
        result = []
        self._preorder(self.root, result)
        return result
    
    def _preorder(self, node: Optional[AVLNode], result: List[Any]) -> None:
        if node:
            result.append(node.value)
            self._preorder(node.left, result)
            self._preorder(node.right, result)
    
    # ==================== Visualization ====================
    
    def print_tree(self) -> None:
        """Print tree structure with balance factors."""
        if self.root is None:
            print("<empty tree>")
            return
        self._print_tree(self.root, "", True)
    
    def _print_tree(self, node: Optional[AVLNode], prefix: str, is_last: bool) -> None:
        """Recursive tree printing helper."""
        if node is not None:
            print(prefix + ("└── " if is_last else "├── ") + 
                  f"{node.value} (h={node.height}, bf={self.get_balance(node)})")
            children = []
            if node.left:
                children.append((node.left, "L"))
            if node.right:
                children.append((node.right, "R"))
            
            for i, (child, side) in enumerate(children):
                is_last_child = (i == len(children) - 1)
                new_prefix = prefix + ("    " if is_last else "│   ")
                self._print_tree(child, new_prefix, is_last_child)
    
    def __repr__(self) -> str:
        return f"AVLTree(inorder={self.inorder()})"
```
