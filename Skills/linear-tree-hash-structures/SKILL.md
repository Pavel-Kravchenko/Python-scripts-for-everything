---
name: linear-tree-hash-structures
description: Linked lists, stacks, queues, dynamic arrays, BST/AVL/Red-Black trees, hash tables with collision resolution, and Bloom filters
primary_tool: Python
---

# Data Structures: Linear, Trees  Hash-Based

## When to Use

| Need | Structure |
|------|-----------|
| O(1) insert/delete at head | Singly linked list |
| O(1) insert/delete at both ends | Doubly linked list |
| LIFO (undo, DFS, bracket matching) | Stack |
| FIFO (BFS, print queue, sliding window) | Queue / deque |
| O(1) amortized append, random access | Dynamic array |
| Sorted order, range queries | BST / AVL |
| Many writes + sorted order | Red-Black tree |
| O(1) avg lookup/insert/delete | Hash table |
| Fast "definitely absent" membership test | Bloom filter |

## Quick Reference: Complexity Table

| Structure | Access | Search | Insert (head/tail) | Insert (mid) | Delete | Space |
|-----------|--------|--------|--------------------|--------------|--------|-------|
| Singly LL | O(n) | O(n) | O(1) / O(1)* | O(n) | O(n)† | O(n) |
| Doubly LL | O(n) | O(n) | O(1) / O(1) | O(n) | O(1)‡ | O(n) |
| Dynamic array | O(1) | O(n) | O(n) / O(1) amort | O(n) | O(n) | O(n) |
| Stack | O(1) top | — | O(1) push | — | O(1) pop | O(n) |
| Queue | O(1) front | — | O(1) enqueue | — | O(1) dequeue | O(n) |
| BST (avg) | O(log n) | O(log n) | O(log n) | — | O(log n) | O(n) |
| BST (worst) | O(n) | O(n) | O(n) | — | O(n) | O(n) |
| AVL | O(log n) | O(log n) | O(log n) | — | O(log n) | O(n) |
| Red-Black | O(log n) | O(log n) | O(log n) | — | O(log n) | O(n) |
| Hash table | O(1) avg | O(1) avg | O(1) avg | — | O(1) avg | O(n) |
| Bloom filter | — | O(k) | O(k) insert | — | N/A | O(m) bits |

\* requires tail pointer; † O(1) with reference to node; ‡ given node reference

## Key Patterns

### Linked Lists
- **Two-pointer**: Floyd's cycle detection (`slow` +1, `fast` +2); find middle; remove nth from end
- **Dummy node**: simplifies edge cases for merge, delete
- **Reverse in-place**: track `prev=None`, `curr=head`, iterate `next=curr.next; curr.next=prev; prev=curr; curr=next`
- **Merge sorted lists**: compare heads, attach smaller, advance pointer

### Stacks
- Bracket/Newick validation: push opening, pop+match on closing
- Postfix evaluation: push operands; on operator, pop two, push result
- Iterative DFS / recursion elimination: push children to explicit stack
- Sliding window max (deque): maintain decreasing deque of indices

### Dynamic Arrays
- Growth factor 2x → O(1) amortized append; total cost ≤ 2n for n appends
- Shrink at ¼ capacity to prevent thrashing
- **Avoid `np.append` in loop** — O(n²); pre-allocate NumPy array when size known
- Python `list` uses ~1.125x growth (memory-efficient CPython formula)

### BST
- **Inorder traversal** yields sorted sequence
- **Delete with two children**: replace with inorder successor (leftmost of right subtree)
- Worst case on sorted input → O(n) height (degenerates to linked list)
- **Treap** (tree + random priority heap) randomises height, expected O(log n)

### AVL Trees
- Balance factor = `height(left) − height(right)`; invariant: |bf| ≤ 1
- Height ≤ 1.44 log₂(n); guaranteed O(log n)
- Rotation decision:

| Node bf | Child bf | Case | Fix |
|---------|----------|------|-----|
| +2 | +1 or 0 | LL | Right rotate |
| +2 | −1 | LR | Left rotate child, right rotate node |
| −2 | −1 or 0 | RR | Left rotate |
| −2 | +1 | RL | Right rotate child, left rotate node |

### Red-Black Trees
- 5 properties: (1) every node RED/BLACK; (2) root BLACK; (3) NIL leaves BLACK; (4) RED node → both children BLACK; (5) all root→leaf paths have equal black-height
- Height ≤ 2 log₂(n+1); insert fix-up: recolor (uncle RED) or 1–2 rotations (uncle BLACK)
- Fewer rotations than AVL on insert/delete → preferred for write-heavy workloads
- Used by: Java `TreeMap`/`TreeSet`, C++ `std::map`/`std::set`, Linux CFS scheduler

### Hash Tables
- Load factor α = n/m; rehash when α > 0.75 (chaining) or > 0.5 (open addressing)
- **Chaining**: bucket = linked list; cache-unfriendly; supports α > 1
- **Linear probing**: probe `(h+i) % m`; cache-friendly; clustering degrades at high α
- **Double hashing**: `(h1(k) + i·h2(k)) % m`; best distribution
- Tombstone marker required for open-addressing deletes (don't break probe chains)

### Bloom Filters
- m bits, k hash functions; false-positive probability: `P ≈ (1 − e^(−kn/m))^k`
- Optimal k = `(m/n) · ln 2 ≈ 0.693 · (m/n)`
- Optimal m = `−n · ln(p) / (ln 2)²`
- ~10 bits/element → ~1% FP; ~16 bits/element → ~0.05% FP
- No false negatives; no deletion (use counting Bloom filter for deletion)

## Code Templates

### Singly Linked List (minimal)
```python
class Node:
    def __init__(self, data, next=None): self.data = data; self.next = next

class SLL:
    def __init__(self): self.head = self.tail = None; self._len = 0
    def add(self, data):               # O(1) append
        n = Node(data)
        if self.tail: self.tail.next = n
        else: self.head = n
        self.tail = n; self._len += 1
    def add_head(self, data):          # O(1) prepend
        n = Node(data, self.head)
        self.head = n
        if not self.tail: self.tail = n
        self._len += 1
    def reverse(self):                 # O(n) in-place
        prev, curr = None, self.head
        self.tail = self.head
        while curr:
            curr.next, prev, curr = prev, curr, curr.next
        self.head = prev
```

### Floyd's Cycle Detection
```python
def has_cycle(head):
    slow = fast = head
    while fast and fast.next:
        slow, fast = slow.next, fast.next.next
        if slow is fast: return True
    return False
```

### Stack (array-backed, Python)
```python
class Stack:
    def __init__(self): self._items = []
    def push(self, x): self._items.append(x)
    def pop(self):
        if not self._items: raise IndexError("pop from empty stack")
        return self._items.pop()
    def peek(self): return self._items[-1]
    def __len__(self): return len(self._items)
```

### Queue (linked-list, O(1) both ends)
```python
class Queue:
    def __init__(self): self._front = self._rear = None; self._size = 0
    def enqueue(self, val):
        n = QueueNode(val)
        if self._rear: self._rear.next = n
        else: self._front = n
        self._rear = n; self._size += 1
    def dequeue(self):
        if not self._front: raise IndexError("dequeue from empty queue")
        val = self._front.value; self._front = self._front.next
        if not self._front: self._rear = None
        self._size -= 1; return val
```

### BST Insert + Search (iterative)
```python
def insert(root, val):
    if root is None: return Node(val)
    node = root
    while True:
        if val < node.val:
            if node.left is None: node.left = Node(val); break
            node = node.left
        elif val > node.val:
            if node.right is None: node.right = Node(val); break
            node = node.right
        else: break
    return root

def search(root, val):
    while root:
        if val == root.val: return root
        root = root.left if val < root.val else root.right
    return None
```

### AVL Rotations (core)
```python
def rotate_right(z):        # LL case
    y, z.left = z.left, z.left.right
    y.right = z
    z.height = 1 + max(height(z.left), height(z.right))
    y.height = 1 + max(height(y.left), height(y.right))
    return y

def rotate_left(z):         # RR case
    y, z.right = z.right, z.right.left
    y.left = z
    z.height = 1 + max(height(z.left), height(z.right))
    y.height = 1 + max(height(y.left), height(y.right))
    return y
```

### Hash Table: Chaining (minimal)
```python
class HashTable:
    def __init__(self, size=7):
        self.buckets = [[] for _ in range(size)]
        self.size = size
    def _h(self, k): return hash(k) % self.size
    def put(self, k, v):
        for item in self.buckets[self._h(k)]:
            if item[0] == k: item[1] = v; return
        self.buckets[self._h(k)].append([k, v])
    def get(self, k):
        for item in self.buckets[self._h(k)]:
            if item[0] == k: return item[1]
        return None
```

### Bloom Filter
```python
import math, hashlib

class BloomFilter:
    def __init__(self, n, fp=0.01):
        self.m = int(-n * math.log(fp) / math.log(2)**2)
        self.k = max(1, int(self.m / n * math.log(2)))
        self.bits = bytearray(self.m)
    def _hashes(self, item):
        h1 = int(hashlib.md5(str(item).encode()).hexdigest(), 16)
        h2 = int(hashlib.sha1(str(item).encode()).hexdigest(), 16)
        return [(h1 + i * h2) % self.m for i in range(self.k)]
    def add(self, item):
        for i in self._hashes(item): self.bits[i] = 1
    def __contains__(self, item):
        return all(self.bits[i] for i in self._hashes(item))
```

## Pitfalls

- **Linked list**: forgetting to update `tail` pointer on head insertion or delete-last
- **Stack/Queue**: off-by-one when using array; need to handle empty-check before pop/dequeue
- **Dynamic array**: using `np.append` in a loop is O(n²); always pre-allocate NumPy arrays
- **BST**: sorted input → O(n) height; always use AVL/RB or shuffle input for production
- **AVL**: update height bottom-up after rotation; update child first, then new parent
- **Red-Black**: new nodes always inserted RED; root must be recolored BLACK after fix-up
- **Hash table open addressing**: must use tombstones on delete or search breaks; table full with α≥1
- **Bloom filter**: cannot delete (bit shared by multiple items); choose m/n ≥ 10 for <1% FP

## Bioinformatics Connections

| Application | Structure | Why |
|-------------|-----------|-----|
| K-mer counting (assembly, variant calling) | Hash table | O(1) count per k-mer; key = k-mer string |
| Known variant lookup (VCF filter) | Bloom filter | 10× memory savings over hash set; FP tolerable |
| Read deduplication (PCR artifacts) | Hash set / Bloom filter | Exact: hash set; approximate: Bloom filter |
| Genome interval store (BED, VCF) | Red-Black tree / BST | Sorted by position; range queries O(log n) |
| Phylogenetic trees (Newick parsing) | Stack | Validate balanced parentheses; build tree recursively |
| Suffix array construction | Dynamic array | Amortised append for suffix list |
| Sequence alignment traceback | Stack | Store path; pop to reconstruct alignment |
| LRU cache for BLAST hits | Doubly LL + hash map | O(1) move-to-front and lookup |
| Chromosome position index | AVL tree | Sorted coordinates; balanced even for sequential inserts |

## Related Skills

- `numpy-pandas-wrangling` — NumPy array internals, pre-allocation patterns
- `python-collections-regex` — Python `collections.deque`, `heapq`, `sortedcontainers`
- `data-visualization-bio` — tree visualisation with `ete3` / `dendropy`
- `complexity-sorting-searching` — Big-O reference for structure operation costs
