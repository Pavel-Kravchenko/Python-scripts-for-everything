---
name: algo-stacks-queues
description: "Stack (LIFO) and Queue (FIFO) — array and linked-list implementations, O(1) operations, common applications"
tool_type: python
primary_tool: Python
---

# Stacks and Queues

## Stack (LIFO)

All operations O(1). Use for: DFS, backtracking, bracket matching, undo/redo, recursion-to-iteration.

### Array-Based Stack

```python
class ArrayStack:
    def __init__(self):
        self._items = []

    def push(self, item):  self._items.append(item)
    def pop(self):
        if not self._items: raise IndexError("Pop from empty stack")
        return self._items.pop()
    def peek(self):
        if not self._items: raise IndexError("Peek from empty stack")
        return self._items[-1]
    def is_empty(self): return len(self._items) == 0
    def __len__(self):  return len(self._items)
    def __bool__(self): return len(self._items) > 0
```

### Linked-List Stack

Truly O(1) push (no amortized resize), but extra pointer memory per node.

```python
class StackNode:
    def __init__(self, value, next_node=None):
        self.value = value
        self.next = next_node

class LinkedStack:
    def __init__(self):
        self._top = None
        self._size = 0

    def push(self, item):
        self._top = StackNode(item, self._top)
        self._size += 1

    def pop(self):
        if not self._top: raise IndexError("Pop from empty stack")
        value = self._top.value
        self._top = self._top.next
        self._size -= 1
        return value

    def peek(self):
        if not self._top: raise IndexError("Peek from empty stack")
        return self._top.value

    def is_empty(self): return self._top is None
    def __len__(self):  return self._size
```

## Queue (FIFO)

All operations O(1). Use for: BFS, task scheduling, message buffers, simulation.

```python
class QueueNode:
    def __init__(self, value, next_node=None):
        self.value = value
        self.next = next_node

class LinkedQueue:
    def __init__(self):
        self._front = self._rear = None
        self._size = 0

    def enqueue(self, item):
        node = QueueNode(item)
        if self._rear is None:
            self._front = self._rear = node
        else:
            self._rear.next = node
            self._rear = node
        self._size += 1

    def dequeue(self):
        if not self._front: raise IndexError("Dequeue from empty queue")
        value = self._front.value
        self._front = self._front.next
        if self._front is None: self._rear = None
        self._size -= 1
        return value

    def peek(self):
        if not self._front: raise IndexError("Peek from empty queue")
        return self._front.value

    def is_empty(self): return self._front is None
    def __len__(self):  return self._size
```

## Python stdlib equivalents

- Stack: just use `list` (append/pop)
- Queue: `collections.deque` (appendleft/append + popleft/pop, all O(1))
- Priority queue: `heapq` or `queue.PriorityQueue`

## Pitfalls

- Using `list.pop(0)` as a queue is O(n) — use `collections.deque.popleft()` instead
- Array-based stack `push` is O(1) amortized but O(n) worst case during resize
- Queue with a plain list requires shifting elements on dequeue — always use deque or linked list
