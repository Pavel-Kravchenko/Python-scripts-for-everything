---
name: algo-avl-trees
description: "Self-balancing BST (Adelson-Velsky & Landis, 1962) guaranteeing O(log n) operations via rotation-based rebalancing."
tool_type: python
primary_tool: Python
---

# AVL Trees

## Key Invariant
Balance factor (bf) = height(left) - height(right). Every node must have |bf| <= 1.
Height is always <= 1.44 * log2(n), guaranteeing O(log n) operations.

## Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Search / Insert / Delete | O(log n) | O(log n) stack |
| Single rotation | O(1) | O(1) |
| Get height / balance | O(1) | O(1) |

## Rotation Decision Table

| Node bf | Child bf | Case | Fix |
|---------|----------|------|-----|
| +2 | >= 0 | LL | Right rotation |
| +2 | -1 | LR | Left on child, then Right |
| -2 | <= 0 | RR | Left rotation |
| -2 | +1 | RL | Right on child, then Left |

Insertion needs at most **1** rotation. Deletion may need **O(log n)** rotations up the tree.

## Implementation

```python
class AVLNode:
    def __init__(self, value):
        self.value = value
        self.left = self.right = None
        self.height = 1


class AVLTree:
    def get_height(self, node):
        return node.height if node else 0

    def get_balance(self, node):
        return self.get_height(node.left) - self.get_height(node.right) if node else 0

    def _update_height(self, node):
        node.height = 1 + max(self.get_height(node.left), self.get_height(node.right))

    def rotate_right(self, z):
        y, T3 = z.left, z.left.right
        y.right, z.left = z, T3
        self._update_height(z)
        self._update_height(y)
        return y

    def rotate_left(self, z):
        y, T2 = z.right, z.right.left
        y.left, z.right = z, T2
        self._update_height(z)
        self._update_height(y)
        return y

    def _rebalance(self, node, value=None, balance=None):
        """Apply the correct rotation(s). Pass value for insert, balance for delete."""
        bf = self.get_balance(node)
        # LL
        if bf > 1 and (value is not None and value < node.left.value or
                       value is None and self.get_balance(node.left) >= 0):
            return self.rotate_right(node)
        # RR
        if bf < -1 and (value is not None and value > node.right.value or
                        value is None and self.get_balance(node.right) <= 0):
            return self.rotate_left(node)
        # LR
        if bf > 1:
            node.left = self.rotate_left(node.left)
            return self.rotate_right(node)
        # RL
        if bf < -1:
            node.right = self.rotate_right(node.right)
            return self.rotate_left(node)
        return node

    def insert(self, value):
        self.root = self._insert(self.root, value)

    def _insert(self, node, value):
        if not node:
            return AVLNode(value)
        if value < node.value:
            node.left = self._insert(node.left, value)
        elif value > node.value:
            node.right = self._insert(node.right, value)
        else:
            return node  # no duplicates
        self._update_height(node)
        bf = self.get_balance(node)
        if bf > 1 and value < node.left.value:
            return self.rotate_right(node)
        if bf < -1 and value > node.right.value:
            return self.rotate_left(node)
        if bf > 1:
            node.left = self.rotate_left(node.left)
            return self.rotate_right(node)
        if bf < -1:
            node.right = self.rotate_right(node.right)
            return self.rotate_left(node)
        return node

    def delete(self, value):
        self.root = self._delete(self.root, value)

    def _delete(self, node, value):
        if not node:
            return None
        if value < node.value:
            node.left = self._delete(node.left, value)
        elif value > node.value:
            node.right = self._delete(node.right, value)
        else:
            if not node.left:
                return node.right
            if not node.right:
                return node.left
            # two children: replace with inorder successor
            succ = node.right
            while succ.left:
                succ = succ.left
            node.value = succ.value
            node.right = self._delete(node.right, succ.value)
        self._update_height(node)
        bf = self.get_balance(node)
        if bf > 1 and self.get_balance(node.left) >= 0:
            return self.rotate_right(node)
        if bf > 1:
            node.left = self.rotate_left(node.left)
            return self.rotate_right(node)
        if bf < -1 and self.get_balance(node.right) <= 0:
            return self.rotate_left(node)
        if bf < -1:
            node.right = self.rotate_right(node.right)
            return self.rotate_left(node)
        return node

    def __init__(self):
        self.root = None
```

## Pitfalls

- **Update height before checking balance**: height must be refreshed on the way back up the recursion, before computing bf.
- **Deletion may cascade**: unlike insertion (at most one rotation), deletion can require rotations at every ancestor up to the root.
- **Child bf=0 on deletion**: when a node has bf=+2 and its left child has bf=0, that is an LL case (right rotation) — not LR. This differs from insertion where bf=0 never triggers a rotation.
- **Height of None is 0, height of leaf is 1**: mixing up these conventions causes off-by-one errors in balance factor calculation.
- **Do not store balance factor directly**: store height and compute bf on the fly; storing bf directly requires careful bookkeeping during rotations.
