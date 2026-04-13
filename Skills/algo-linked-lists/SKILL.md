---
name: algo-linked-lists
description: "A comprehensive guide to singly and doubly linked lists with implementations, visualizations, and practical examples."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/04_Linear_Data_Structures/01_linked_lists.ipynb"
---

# Linked Lists

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/04_Linear_Data_Structures/01_linked_lists.ipynb`*

# Linked Lists

A comprehensive guide to singly and doubly linked lists with implementations, visualizations, and practical examples.

---

## Table of Contents

1. [Introduction to Linked Lists](#1-introduction-to-linked-lists)
2. [Singly Linked List](#2-singly-linked-list)
3. [Doubly Linked List](#3-doubly-linked-list)
4. [Common Interview Problems](#4-common-interview-problems)
5. [Summary](#5-summary)

---

## 1. Introduction to Linked Lists

### What is a Linked List?

A **linked list** is a linear data structure where elements are stored in nodes, and each node points to the next node in the sequence. Unlike arrays, linked list elements are not stored in contiguous memory locations.

### Why Use Linked Lists?

| Advantage | Description |
|-----------|-------------|
| **Dynamic Size** | Can grow or shrink during runtime |
| **Efficient Insertion/Deletion** | O(1) at known positions (no shifting required) |
| **Memory Efficiency** | Allocates memory as needed |
| **No Wasted Space** | No pre-allocated unused slots |

### Comparison: Linked Lists vs Arrays

| Feature | Array | Linked List |
|---------|-------|-------------|
| **Memory Layout** | Contiguous | Non-contiguous |
| **Size** | Fixed (static) or costly resize | Dynamic |
| **Random Access** | O(1) | O(n) |
| **Insertion at Beginning** | O(n) | O(1) |
| **Memory Overhead** | None | Extra pointer(s) per node |
| **Cache Performance** | Excellent (locality) | Poor (scattered) |

### Memory Representation

```
Array (contiguous memory):
┌─────┬─────┬─────┬─────┬─────┐
│  5  │  3  │  8  │  1  │  9  │   Memory addresses: 100, 104, 108, 112, 116
└─────┴─────┴─────┴─────┴─────┘

Linked List (scattered memory):
┌─────┐     ┌─────┐     ┌─────┐     ┌─────┐
│  5  │────→│  3  │────→│  8  │────→│  1  │────→ null
└─────┘     └─────┘     └─────┘     └─────┘
addr:200    addr:508    addr:124    addr:336
```

---

## 2. Singly Linked List

### Theory

A **singly linked list** is the simplest form of linked list where each node contains:
- **Data**: The value stored in the node
- **Next**: A pointer/reference to the next node

The list maintains a **head** pointer to the first node. The last node's `next` pointer is `null`.

### Visual Structure

```
head
  ↓
┌───┬───┐   ┌───┬───┐   ┌───┬───┐   ┌───┬───┐
│ 5 │ ●─┼──→│ 3 │ ●─┼──→│ 8 │ ●─┼──→│ 1 │ ○ │
└───┴───┘   └───┴───┘   └───┴───┘   └───┴───┘
 data next   data next   data next   data null
```

### Time Complexity

| Operation | Time Complexity | Notes |
|-----------|-----------------|-------|
| Access by index | O(n) | Must traverse from head |
| Insert at head | O(1) | Direct pointer update |
| Insert at tail | O(n) or O(1)* | *O(1) if tail pointer maintained |
| Insert at position | O(n) | Must find position first |
| Delete at head | O(1) | Direct pointer update |
| Delete at tail | O(n) | Must find second-to-last |
| Delete by value | O(n) | Must search for value |
| Search | O(n) | Linear traversal |
| Get length | O(1) or O(n) | O(1) if length cached |

```python
from typing import Optional, Any, Iterator


class Node:
    """
    A node in a singly linked list.
    
    Attributes:
        data: The value stored in this node
        next: Reference to the next node (None if this is the last node)
    """
    
    def __init__(self, data: Any, next_node: Optional['Node'] = None) -> None:
        self.data = data
        self.next = next_node
    
    def __repr__(self) -> str:
        return f"Node({self.data})"


class SinglyLinkedList:
    """
    A singly linked list implementation with head and tail pointers.
    
    Maintains both head and tail pointers for O(1) insertion at both ends.
    Also caches length for O(1) size queries.
    
    Attributes:
        head: Reference to the first node
        tail: Reference to the last node
        _length: Number of nodes in the list
    """
    
    def __init__(self) -> None:
        """Initialize an empty linked list."""
        self.head: Optional[Node] = None
        self.tail: Optional[Node] = None
        self._length: int = 0
    
    def __len__(self) -> int:
        """Return the number of nodes in the list. O(1)"""
        return self._length
    
    def __iter__(self) -> Iterator[Any]:
        """Iterate over values in the list."""
        current = self.head
        while current is not None:
            yield current.data
            current = current.next
    
    def __str__(self) -> str:
        """Return string representation: [5] -> [3] -> [8] -> null"""
        if self.head is None:
            return "empty list"
        
        nodes = []
        current = self.head
        while current is not None:
            nodes.append(f"[{current.data}]")
            current = current.next
        nodes.append("null")
        return " -> ".join(nodes)
    
    def is_empty(self) -> bool:
        """Check if the list is empty. O(1)"""
        return self.head is None
    
    # ==================== INSERTION OPERATIONS ====================
    
    def add_at_head(self, data: Any) -> None:
        """
        Insert a new node at the beginning of the list. O(1)
        
        Visual:
            Before: head -> [3] -> [8] -> null
            After:  head -> [5] -> [3] -> [8] -> null
        """
        new_node = Node(data, self.head)
        self.head = new_node
        
        # If list was empty, new node is also the tail
        if self.tail is None:
            self.tail = new_node
        
        self._length += 1
    
    def add_at_tail(self, data: Any) -> None:
        """
        Insert a new node at the end of the list. O(1) with tail pointer.
        
        Visual:
            Before: head -> [5] -> [3] -> null
            After:  head -> [5] -> [3] -> [8] -> null
        """
        new_node = Node(data)
        
        if self.tail is None:
            # List is empty
            self.head = self.tail = new_node
        else:
            self.tail.next = new_node
            self.tail = new_node
        
        self._length += 1
    
    def add(self, data: Any) -> None:
        """Alias for add_at_tail - adds element to the end. O(1)"""
        self.add_at_tail(data)
    
    def insert_at_index(self, index: int, data: Any) -> bool:
        """
        Insert a new node at the specified index. O(n)
        
        Args:
            index: Position to insert (0-based)
            data: Value to insert
            
        Returns:
            True if successful, False if index out of bounds
        """
        if index < 0 or index > self._length:
            return False
        
        if index == 0:
            self.add_at_head(data)
            return True
        
        if index == self._length:
            self.add_at_tail(data)
            return True
        
        # Find node at position (index - 1)
        current = self.head
        for _ in range(index - 1):
            current = current.next
        
        # Insert new node after current
        new_node = Node(data, current.next)
        current.next = new_node
        self._length += 1
        return True
    
    # ==================== DELETION OPERATIONS ====================
    
    def delete_at_head(self) -> Optional[Any]:
        """
        Remove and return the first node's data. O(1)
        
        Visual:
            Before: head -> [5] -> [3] -> [8] -> null
            After:  head -> [3] -> [8] -> null
            Returns: 5
        """
        if self.head is None:
            return None
        
        data = self.head.data
        self.head = self.head.next
        
        # If list becomes empty, update tail
        if self.head is None:
            self.tail = None
        
        self._length -= 1
        return data
    
    def delete_at_tail(self) -> Optional[Any]:
        """
        Remove and return the last node's data. O(n)
        
        Note: O(n) because we need to find the second-to-last node
        to update its next pointer.
        """
        if self.head is None:
            return None
        
        # Single node case
        if self.head == self.tail:
            data = self.head.data
            self.head = self.tail = None
            self._length -= 1
            return data
        
        # Find second-to-last node
        current = self.head
        while current.next != self.tail:
            current = current.next
        
        data = self.tail.data
        current.next = None
        self.tail = current
        self._length -= 1
        return data
    
    def delete_by_value(self, data: Any) -> bool:
        """
        Remove the first node with the specified value. O(n)
        
        Visual:
            Delete 3:
            Before: [5] -> [3] -> [8] -> null
            After:  [5] ---------> [8] -> null
        
        Returns:
            True if node was found and deleted, False otherwise
        """
        if self.head is None:
            return False
        
        # Special case: deleting head
        if self.head.data == data:
            self.delete_at_head()
            return True
        
        # Find the node before the one to delete
        current = self.head
        while current.next is not None and current.next.data != data:
            current = current.next
        
        # Node not found
        if current.next is None:
            return False
        
        # Delete the node
        if current.next == self.tail:
            self.tail = current
        
        current.next = current.next.next
        self._length -= 1
        return True
    
    def delete_at_index(self, index: int) -> Optional[Any]:
        """
        Remove and return the node at the specified index. O(n)
        
        Returns:
            The deleted value, or None if index out of bounds
        """
        if index < 0 or index >= self._length:
            return None
        
        if index == 0:
            return self.delete_at_head()
        
        # Find node at position (index - 1)
        current = self.head
        for _ in range(index - 1):
            current = current.next
        
        data = current.next.data
        
        # Update tail if needed
        if current.next == self.tail:
            self.tail = current
        
        current.next = current.next.next
        self._length -= 1
        return data
    
    # ==================== SEARCH OPERATIONS ====================
    
    def search(self, data: Any) -> int:
        """
        Find the index of the first node with the specified value. O(n)
        
        Returns:
            Index of the node (0-based), or -1 if not found
        """
        current = self.head
        index = 0
        
        while current is not None:
            if current.data == data:
                return index
            current = current.next
            index += 1
        
        return -1
    
    def contains(self, data: Any) -> bool:
        """Check if the list contains the specified value. O(n)"""
        return self.search(data) != -1
    
    def get_at_index(self, index: int) -> Optional[Any]:
        """
        Get the value at the specified index. O(n)
        
        Returns:
            The value at the index, or None if index out of bounds
        """
        if index < 0 or index >= self._length:
            return None
        
        current = self.head
        for _ in range(index):
            current = current.next
        
        return current.data
    
    # ==================== UTILITY OPERATIONS ====================
    
    def reverse(self) -> None:
        """
        Reverse the linked list in place. O(n)
        
        Algorithm:
            Before: [1] -> [2] -> [3] -> null
            
            Step 1: null <- [1]    [2] -> [3] -> null
            Step 2: null <- [1] <- [2]    [3] -> null  
            Step 3: null <- [1] <- [2] <- [3]
            
            After:  [3] -> [2] -> [1] -> null
        """
        # Update tail to current head before reversing
        self.tail = self.head
        
        prev = None
        current = self.head
        
        while current is not None:
            next_node = current.next  # Save next
            current.next = prev       # Reverse pointer
            prev = current            # Move prev forward
            current = next_node       # Move current forward
        
        self.head = prev
    
    def to_list(self) -> list:
        """Convert linked list to Python list. O(n)"""
        return list(self)
    
    def clear(self) -> None:
        """Remove all nodes from the list. O(1)"""
        self.head = None
        self.tail = None
        self._length = 0
```
