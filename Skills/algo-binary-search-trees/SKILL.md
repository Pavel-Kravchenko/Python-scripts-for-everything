---
name: algo-binary-search-trees
description: "BST operations, complexity, and a clean implementation with parent pointers supporting all standard operations."
tool_type: python
primary_tool: Python
---

# Binary Search Trees (BST)

**Invariant**: left subtree < node < right subtree (no duplicates).

## Complexity

| Operation | Average | Worst (degenerate) |
|-----------|---------|-------------------|
| Search / Insert / Delete | O(log n) | O(n) |
| Min / Max | O(log n) | O(n) |
| Successor / Predecessor | O(log n) | O(n) |
| Inorder traversal | O(n) | O(n) |

Worst case occurs on sorted input — use AVL or Red-Black tree to guarantee O(log n).

## Implementation

```python
from collections import deque
from typing import Any

class Node:
    def __init__(self, value: Any) -> None:
        self.value = value
        self.left = self.right = self.parent = None


class BST:
    def __init__(self) -> None:
        self.root: Node | None = None

    # --- Insert (iterative) ---
    def insert(self, value: Any) -> Node:
        if not self.root:
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
                return node  # duplicate

    # --- Search (iterative) ---
    def search(self, value: Any) -> Node | None:
        node = self.root
        while node:
            if value == node.value:
                return node
            node = node.left if value < node.value else node.right
        return None

    # --- Min / Max ---
    def find_min(self, node: Node | None = None) -> Node | None:
        node = node or self.root
        if not node:
            return None
        while node.left:
            node = node.left
        return node

    def find_max(self, node: Node | None = None) -> Node | None:
        node = node or self.root
        if not node:
            return None
        while node.right:
            node = node.right
        return node

    # --- Successor / Predecessor ---
    def successor(self, node: Node) -> Node | None:
        if node.right:
            return self.find_min(node.right)
        p = node.parent
        while p and node == p.right:
            node, p = p, p.parent
        return p

    def predecessor(self, node: Node) -> Node | None:
        if node.left:
            return self.find_max(node.left)
        p = node.parent
        while p and node == p.left:
            node, p = p, p.parent
        return p

    # --- Delete ---
    def delete(self, value: Any) -> bool:
        node = self.search(value)
        if not node:
            return False
        self._delete_node(node)
        return True

    def _delete_node(self, node: Node) -> None:
        if node.left and node.right:
            succ = self.find_min(node.right)
            node.value = succ.value
            self._delete_node(succ)
        else:
            child = node.left or node.right
            self._transplant(node, child)

    def _transplant(self, u: Node, v: Node | None) -> None:
        if not u.parent:
            self.root = v
        elif u == u.parent.left:
            u.parent.left = v
        else:
            u.parent.right = v
        if v:
            v.parent = u.parent

    # --- Traversals ---
    def inorder(self) -> list[Any]:
        result: list[Any] = []
        def _rec(n):
            if n:
                _rec(n.left); result.append(n.value); _rec(n.right)
        _rec(self.root)
        return result

    def level_order(self) -> list[Any]:
        if not self.root:
            return []
        result, q = [], deque([self.root])
        while q:
            n = q.popleft()
            result.append(n.value)
            if n.left: q.append(n.left)
            if n.right: q.append(n.right)
        return result

    def __contains__(self, value: Any) -> bool:
        return self.search(value) is not None
```

## Delete: Three Cases

| Node type | Action |
|-----------|--------|
| Leaf | Remove directly |
| One child | Replace node with that child |
| Two children | Copy inorder successor's value, delete successor (which has at most one child) |

## Pitfalls

- **Sorted input degenerates to O(n)**: inserting [1,2,3,4,5] creates a linked list; use a self-balancing tree (AVL, Red-Black) for production use.
- **Inorder successor during deletion**: after copying the successor's value to the deleted node, you must delete the successor node — not the original node again.
- **Parent pointer consistency**: when implementing with parent pointers, update `v.parent` in `_transplant` only if `v` is not None.
- **Recursive height on deep trees**: Python's default recursion limit (~1000) is hit on degenerate trees; use iterative BFS for height calculation.
