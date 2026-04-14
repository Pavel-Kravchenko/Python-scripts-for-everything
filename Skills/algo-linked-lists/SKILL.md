---
name: algo-linked-lists
description: "Singly linked list — full implementation with head/tail pointers, insert/delete/search/reverse, complexity table"
tool_type: python
primary_tool: Python
---

# Linked Lists

## Array vs Linked List

| Feature | Array | Linked List |
|---------|-------|-------------|
| Memory | Contiguous | Scattered |
| Random Access | O(1) | O(n) |
| Insert at head | O(n) | O(1) |
| Cache performance | Excellent | Poor |

## Time Complexity

| Operation | Time | Notes |
|-----------|------|-------|
| Access by index | O(n) | Must traverse |
| Insert at head | O(1) | |
| Insert at tail | O(1)* | *With tail pointer |
| Delete at head | O(1) | |
| Delete at tail | O(n) | Must find second-to-last |
| Search | O(n) | |

## Implementation

```python
from typing import Optional, Any, Iterator

class Node:
    def __init__(self, data: Any, next_node: Optional['Node'] = None):
        self.data = data
        self.next = next_node

class SinglyLinkedList:
    def __init__(self):
        self.head: Optional[Node] = None
        self.tail: Optional[Node] = None
        self._length: int = 0

    def __len__(self): return self._length
    def __iter__(self) -> Iterator[Any]:
        current = self.head
        while current:
            yield current.data
            current = current.next
    def __str__(self):
        return " -> ".join(f"[{d}]" for d in self) + " -> null" if self.head else "empty list"
    def is_empty(self): return self.head is None

    # --- Insertion ---
    def add_at_head(self, data):
        new_node = Node(data, self.head)
        self.head = new_node
        if self.tail is None:
            self.tail = new_node
        self._length += 1

    def add_at_tail(self, data):
        new_node = Node(data)
        if self.tail is None:
            self.head = self.tail = new_node
        else:
            self.tail.next = new_node
            self.tail = new_node
        self._length += 1

    def insert_at_index(self, index, data):
        if index < 0 or index > self._length: return False
        if index == 0: self.add_at_head(data); return True
        if index == self._length: self.add_at_tail(data); return True
        current = self.head
        for _ in range(index - 1):
            current = current.next
        current.next = Node(data, current.next)
        self._length += 1
        return True

    # --- Deletion ---
    def delete_at_head(self):
        if not self.head: return None
        data = self.head.data
        self.head = self.head.next
        if not self.head: self.tail = None
        self._length -= 1
        return data

    def delete_by_value(self, data):
        if not self.head: return False
        if self.head.data == data:
            self.delete_at_head(); return True
        current = self.head
        while current.next and current.next.data != data:
            current = current.next
        if not current.next: return False
        if current.next == self.tail: self.tail = current
        current.next = current.next.next
        self._length -= 1
        return True

    # --- Search ---
    def search(self, data):
        current, index = self.head, 0
        while current:
            if current.data == data: return index
            current = current.next
            index += 1
        return -1

    # --- Reverse ---
    def reverse(self):
        self.tail = self.head
        prev, current = None, self.head
        while current:
            next_node = current.next
            current.next = prev
            prev = current
            current = next_node
        self.head = prev
```

## Pitfalls

- Forgetting to update `tail` when deleting the last node or when reversing
- Delete at tail is O(n) in singly linked lists — use doubly linked list if frequent tail deletion is needed
- Memory overhead: each node carries a pointer (~8 bytes on 64-bit) in addition to data
