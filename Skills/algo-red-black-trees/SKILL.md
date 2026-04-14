---
name: algo-red-black-trees
description: "Red-black tree — self-balancing BST with O(log n) operations, 5 invariants, insert fix-up with rotations and recoloring"
tool_type: python
primary_tool: Python
---

# Red-Black Trees

Self-balancing BST where each node is RED or BLACK. Guarantees O(log n) search/insert/delete.

## 5 Invariants

| # | Property |
|---|----------|
| 1 | Every node is RED or BLACK |
| 2 | Root is BLACK |
| 3 | Every leaf (NIL) is BLACK |
| 4 | RED node has only BLACK children (no red-red) |
| 5 | All root-to-leaf paths have equal BLACK node count (black-height) |

These guarantee: longest path <= 2x shortest path, so height is O(log n).

## Red-Black vs AVL

| Aspect | AVL | Red-Black |
|--------|-----|-----------|
| Height bound | 1.44 log n | 2 log n |
| Search | Faster (shorter) | Slightly slower |
| Insert/Delete | More rotations | Fewer rotations (recolor often suffices) |
| Best for | Read-heavy | Write-heavy |
| Used in | Database indexes | Java TreeMap, C++ std::map, Linux CFS |

## Implementation

```python
RED, BLACK = True, False

class RBNode:
    def __init__(self, data, color=RED):
        self.data = data
        self.color = color  # new nodes are RED
        self.parent = self.left = self.right = None

class RedBlackTree:
    def __init__(self):
        self.NIL = RBNode(data=None, color=BLACK)
        self.root = self.NIL

    def rotate_left(self, x):
        y = x.right
        x.right = y.left
        if y.left != self.NIL: y.left.parent = x
        y.parent = x.parent
        if x.parent is None: self.root = y
        elif x == x.parent.left: x.parent.left = y
        else: x.parent.right = y
        y.left = x; x.parent = y

    def rotate_right(self, y):
        x = y.left
        y.left = x.right
        if x.right != self.NIL: x.right.parent = y
        x.parent = y.parent
        if y.parent is None: self.root = x
        elif y == y.parent.right: y.parent.right = x
        else: y.parent.left = x
        x.right = y; y.parent = x

    def insert(self, data):
        node = RBNode(data, color=RED)
        node.left = node.right = self.NIL
        parent = None
        current = self.root
        while current != self.NIL:
            parent = current
            if data < current.data: current = current.left
            elif data > current.data: current = current.right
            else: return  # no duplicates
        node.parent = parent
        if parent is None: self.root = node
        elif data < parent.data: parent.left = node
        else: parent.right = node
        self._fix_insert(node)

    def _fix_insert(self, node):
        while node != self.root and node.parent.color == RED:
            gp = node.parent.parent
            if node.parent == gp.left:
                uncle = gp.right
                if uncle.color == RED:  # Case 1: recolor
                    node.parent.color = uncle.color = BLACK
                    gp.color = RED; node = gp
                else:
                    if node == node.parent.right:  # Case 2: LR -> rotate left
                        node = node.parent; self.rotate_left(node)
                    node.parent.color = BLACK  # Case 3: LL -> rotate right
                    gp.color = RED; self.rotate_right(gp)
            else:  # mirror
                uncle = gp.left
                if uncle.color == RED:
                    node.parent.color = uncle.color = BLACK
                    gp.color = RED; node = gp
                else:
                    if node == node.parent.left:
                        node = node.parent; self.rotate_right(node)
                    node.parent.color = BLACK
                    gp.color = RED; self.rotate_left(gp)
        self.root.color = BLACK

    def search(self, data):
        current = self.root
        while current != self.NIL:
            if data == current.data: return current
            current = current.left if data < current.data else current.right
        return None
```

## Insert Fix-Up Cases (parent is left child of grandparent)

| Case | Condition | Action |
|---|---|---|
| 1 | Uncle is RED | Recolor parent+uncle BLACK, grandparent RED, move up |
| 2 | Uncle BLACK, node is right child | Left-rotate parent (reduces to case 3) |
| 3 | Uncle BLACK, node is left child | Right-rotate grandparent, recolor |

Mirror cases apply when parent is right child.

## Pitfalls

- NIL sentinel must be BLACK and shared — don't create new NIL nodes
- After insert fix-up, always force root to BLACK
- Deletion fix-up is significantly more complex (4 cases + mirrors) — omitted here but follows same rotation/recolor pattern
