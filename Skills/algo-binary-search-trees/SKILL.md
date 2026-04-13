---
name: algo-binary-search-trees
description: "A comprehensive guide to Binary Search Trees: definition, properties, operations, and implementation."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/05_Tree_Structures/01_binary_search_trees.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Binary Search Trees (BST)

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/05_Tree_Structures/01_binary_search_trees.ipynb`*


A comprehensive guide to Binary Search Trees: definition, properties, operations, and implementation.

**Reference:** Original implementation from `origin/Algorithms_HW5.ipynb`

## Table of Contents

1. [BST Definition & Properties](#1-bst-definition--properties)
2. [Node and Tree Implementation](#2-node-and-tree-implementation)
3. [Insert Operation](#3-insert-operation)
4. [Search Operation](#4-search-operation)
5. [Min/Max Operations](#5-minmax-operations)
6. [Successor/Predecessor](#6-successorpredecessor)
7. [Tree Traversals](#7-tree-traversals)
8. [Delete Operation](#8-delete-operation)
9. [Complexity Analysis](#9-complexity-analysis)
10. [Complete Examples](#10-complete-examples)

---

## 1. BST Definition & Properties

### What is a Binary Search Tree?

A **Binary Search Tree (BST)** is a binary tree data structure where each node satisfies the **BST property**:

> For every node `N`:
> - All values in the **left subtree** are **less than** `N.value`
> - All values in the **right subtree** are **greater than** `N.value`

### Visual Representation

```python
        8
       / \
      3   10
     / \    \
    1   6    14
       / \   /
      4   7 13

BST Property Verification:
├── Node 8:  left subtree {1,3,4,6,7} < 8 < {10,13,14} right subtree  ✓
├── Node 3:  left subtree {1} < 3 < {4,6,7} right subtree             ✓
├── Node 10: left subtree {} < 10 < {13,14} right subtree             ✓
├── Node 6:  left subtree {4} < 6 < {7} right subtree                 ✓
├── Node 14: left subtree {13} < 14 < {} right subtree                ✓
└── Leaf nodes (1, 4, 7, 13): No children to verify                   ✓
```python

### Why Use a BST?

The BST property enables **efficient searching**:

1. **Divide and Conquer**: At each node, we eliminate half of the remaining nodes
2. **Ordered Structure**: Inorder traversal yields sorted sequence
3. **Dynamic Operations**: Supports insert/delete while maintaining order

### Key Properties

| Property | Description |
|----------|-------------|
| Ordering | left < node < right |
| Uniqueness | Typically no duplicate keys |
| Inorder Traversal | Produces sorted sequence |
| Height | Best: O(log n), Worst: O(n) |

---

## 2. Node and Tree Implementation

### Node Class

Each node stores:
- `value`: The key/data
- `left`: Reference to left child
- `right`: Reference to right child
- `parent`: Reference to parent node (useful for deletion)

```python
from typing import Optional, Any, List, Generator
from collections import deque


class Node:
    """A node in a Binary Search Tree.
    
    Attributes:
        value: The key stored in this node.
        left: Reference to the left child (smaller values).
        right: Reference to the right child (larger values).
        parent: Reference to the parent node.
    """
    
    def __init__(self, value: Any) -> None:
        """Initialize a new node with the given value.
        
        Args:
            value: The key to store in this node.
        """
        self.value = value
        self.left: Optional['Node'] = None
        self.right: Optional['Node'] = None
        self.parent: Optional['Node'] = None
    
    def __repr__(self) -> str:
        """Return string representation of the node."""
        return f"Node({self.value})"
    
    def is_leaf(self) -> bool:
        """Check if node is a leaf (no children)."""
        return self.left is None and self.right is None
    
    def has_one_child(self) -> bool:
        """Check if node has exactly one child."""
        return (self.left is None) != (self.right is None)
    
    def has_two_children(self) -> bool:
        """Check if node has two children."""
        return self.left is not None and self.right is not None
```python

```python
class BinarySearchTree:
    """Binary Search Tree implementation with full operations.
    
    Supports: insert, search, delete, min/max, successor/predecessor,
    and all standard traversals.
    
    Attributes:
        root: The root node of the tree.
    """
    
    def __init__(self) -> None:
        """Initialize an empty BST."""
        self.root: Optional[Node] = None
    
    def is_empty(self) -> bool:
        """Check if the tree is empty."""
        return self.root is None
    
    # ==================== INSERT ====================
    
    def insert(self, value: Any) -> Node:
        """Insert a value into the BST (iterative to avoid stack overflow).

        Args:
            value: The value to insert.

        Returns:
            The newly created node.

        Note:
            Duplicate values are not inserted (returns existing node).
        """
        if self.root is None:
            self.root = Node(value)
            return self.root

        node = self.root
        while True:
            if value < node.value:
                if node.left is None:
                    node.left = Node(value)
                    node.left.parent = node
                    return node.left
                node = node.left
            elif value > node.value:
                if node.right is None:
                    node.right = Node(value)
                    node.right.parent = node
                    return node.right
                node = node.right
            else:
                return node  # Duplicate, return existing node
    
    # ==================== SEARCH ====================
    
    def search(self, value: Any) -> Optional[Node]:
        """Search for a value in the BST.
        
        Args:
            value: The value to search for.
            
        Returns:
            The node containing the value, or None if not found.
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
    
    def contains(self, value: Any) -> bool:
        """Check if value exists in the BST."""
        return self.search(value) is not None
    
    # ==================== MIN / MAX ====================
    
    def find_min(self, node: Optional[Node] = None) -> Optional[Node]:
        """Find the node with minimum value.
        
        Args:
            node: Starting node (defaults to root).
            
        Returns:
            Node with minimum value, or None if tree is empty.
        """
        if node is None:
            node = self.root
        if node is None:
            return None
        while node.left is not None:
            node = node.left
        return node
    
    def find_max(self, node: Optional[Node] = None) -> Optional[Node]:
        """Find the node with maximum value.
        
        Args:
            node: Starting node (defaults to root).
            
        Returns:
            Node with maximum value, or None if tree is empty.
        """
        if node is None:
            node = self.root
        if node is None:
            return None
        while node.right is not None:
            node = node.right
        return node
    
    # ==================== SUCCESSOR / PREDECESSOR ====================
    
    def successor(self, node: Node) -> Optional[Node]:
        """Find the inorder successor of a node.
        
        The successor is the node with the smallest value greater than node.value.
        
        Args:
            node: The node to find successor for.
            
        Returns:
            The successor node, or None if node has the maximum value.
        """
        # Case 1: Node has right subtree -> min of right subtree
        if node.right is not None:
            return self.find_min(node.right)
        
        # Case 2: Go up until we find a node that is a left child
        parent = node.parent
        while parent is not None and node == parent.right:
            node = parent
            parent = parent.parent
        return parent
    
    def predecessor(self, node: Node) -> Optional[Node]:
        """Find the inorder predecessor of a node.
        
        The predecessor is the node with the largest value smaller than node.value.
        
        Args:
            node: The node to find predecessor for.
            
        Returns:
            The predecessor node, or None if node has the minimum value.
        """
        # Case 1: Node has left subtree -> max of left subtree
        if node.left is not None:
            return self.find_max(node.left)
        
        # Case 2: Go up until we find a node that is a right child
        parent = node.parent
        while parent is not None and node == parent.left:
            node = parent
            parent = parent.parent
        return parent
    
    # ==================== TRAVERSALS ====================
    
    def inorder(self, node: Optional[Node] = None) -> List[Any]:
        """Inorder traversal (Left-Node-Right) - produces sorted output."""
        if node is None:
            node = self.root
        result = []
        self._inorder_recursive(node, result)
        return result
    
    def _inorder_recursive(self, node: Optional[Node], result: List[Any]) -> None:
        if node is not None:
            self._inorder_recursive(node.left, result)
            result.append(node.value)
            self._inorder_recursive(node.right, result)
    
    def preorder(self, node: Optional[Node] = None) -> List[Any]:
        """Preorder traversal (Node-Left-Right)."""
        if node is None:
            node = self.root
        result = []
        self._preorder_recursive(node, result)
        return result
    
    def _preorder_recursive(self, node: Optional[Node], result: List[Any]) -> None:
        if node is not None:
            result.append(node.value)
            self._preorder_recursive(node.left, result)
            self._preorder_recursive(node.right, result)
    
    def postorder(self, node: Optional[Node] = None) -> List[Any]:
        """Postorder traversal (Left-Right-Node)."""
        if node is None:
            node = self.root
        result = []
        self._postorder_recursive(node, result)
        return result
    
    def _postorder_recursive(self, node: Optional[Node], result: List[Any]) -> None:
        if node is not None:
            self._postorder_recursive(node.left, result)
            self._postorder_recursive(node.right, result)
            result.append(node.value)
    
    def level_order(self) -> List[Any]:
        """Level-order traversal (BFS) - visits nodes level by level."""
        if self.root is None:
            return []
        result = []
        queue = deque([self.root])
        while queue:
            node = queue.popleft()
            result.append(node.value)
            if node.left:
                queue.append(node.left)
            if node.right:
                queue.append(node.right)
        return result
    
    # ==================== DELETE ====================
    
    def delete(self, value: Any) -> bool:
        """Delete a value from the BST.
        
        Handles three cases:
        1. Leaf node: Simply remove
        2. One child: Replace with child
        3. Two children: Replace with inorder successor
        
        Args:
            value: The value to delete.
            
        Returns:
            True if deleted, False if value not found.
        """
        node = self.search(value)
        if node is None:
            return False
        self._delete_node(node)
        return True
    
    def _delete_node(self, node: Node) -> None:
        """Internal method to delete a specific node."""
        # Case 1: Leaf node (no children)
        if node.is_leaf():
            self._replace_node(node, None)
        
        # Case 2: One child
        elif node.has_one_child():
            child = node.left if node.left else node.right
            self._replace_node(node, child)
        
        # Case 3: Two children
        else:
            # Find inorder successor (minimum in right subtree)
            successor = self.find_min(node.right)
            # Copy successor's value to current node
            node.value = successor.value
            # Delete the successor (which has at most one child)
            self._delete_node(successor)
    
    def _replace_node(self, node: Node, new_node: Optional[Node]) -> None:
        """Replace a node with another node (or None)."""
        if node.parent is None:
            # Node is root
            self.root = new_node
        elif node == node.parent.left:
            node.parent.left = new_node
        else:
            node.parent.right = new_node
        
        if new_node is not None:
            new_node.parent = node.parent
    
    # ==================== UTILITY ====================
    
    def height(self, node: Optional[Node] = None) -> int:
        """Calculate the height of the tree (or subtree rooted at node).

        Uses an iterative BFS to avoid hitting Python's recursion limit on
        degenerate (sorted-input) trees.

        Returns -1 for an empty tree, 0 for a single node.
        """
        from collections import deque as _deque
        if node is None:
            node = self.root
        if node is None:
            return -1
        # BFS level-by-level; height = number of levels - 1
        queue = _deque([node])
        h = -1
        while queue:
            h += 1
            for _ in range(len(queue)):
                n = queue.popleft()
                if n.left:
                    queue.append(n.left)
                if n.right:
                    queue.append(n.right)
        return h

    def size(self, node: Optional[Node] = None) -> int:
        """Count the number of nodes in the tree (or subtree rooted at node).

        Uses an iterative traversal to avoid hitting Python's recursion limit.
        """
        if node is None:
            node = self.root
        if node is None:
            return 0
        count = 0
        stack = [node]
        while stack:
            n = stack.pop()
            count += 1
            if n.left:
                stack.append(n.left)
            if n.right:
                stack.append(n.right)
        return count
    
    def clear(self) -> None:
        """Remove all nodes from the tree."""
        self.root = None
    
    def __repr__(self) -> str:
        """Return string representation (inorder traversal)."""
        return f"BST({self.inorder()})"
    
    def __contains__(self, value: Any) -> bool:
        """Support 'in' operator."""
        return self.contains(value)
    
    def __len__(self) -> int:
        """Support len() function."""
        return self.size()
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
