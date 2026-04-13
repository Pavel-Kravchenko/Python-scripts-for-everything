---
name: algo-stacks-queues
description: "**Fundamental Linear Data Structures for Managing Collections**"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/04_Linear_Data_Structures/02_stacks_queues.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Stacks and Queues

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/04_Linear_Data_Structures/02_stacks_queues.ipynb`*

# Stacks and Queues

**Fundamental Linear Data Structures for Managing Collections**

This notebook covers the essential stack and queue data structures, their implementations, and practical applications.

## Table of Contents
1. [Stack (LIFO)](#1-stack-lifo)
   - Array-based implementation
   - Linked-list-based implementation
2. [Queue (FIFO)](#2-queue-fifo)
3. [Queue Using Two Stacks](#3-queue-using-two-stacks)
4. [Applications](#4-applications)
   - Bracket matching
   - Postfix expression evaluation
   - Iterative QuickSort
5. [Complexity Summary](#5-complexity-summary)

---

## 1. Stack (LIFO)

### What is a Stack?

A **stack** is a linear data structure that follows the **Last-In-First-Out (LIFO)** principle. The last element added to the stack is the first one to be removed.

### Real-World Analogies

- **Stack of plates**: You add plates to the top and remove from the top
- **Browser history**: Back button goes to the most recently visited page
- **Undo functionality**: Most recent action is undone first
- **Function call stack**: Most recently called function returns first

### When to Use?

- **Backtracking algorithms** (maze solving, parsing)
- **Expression evaluation** (postfix, prefix notation)
- **Undo/redo mechanisms**
- **Converting recursion to iteration**
- **Depth-First Search (DFS)**
- **Bracket matching / syntax validation**

### Core Operations

| Operation | Description |
|-----------|-------------|
| `push(item)` | Add item to the top |
| `pop()` | Remove and return top item |
| `peek()` / `top()` | Return top item without removing |
| `is_empty()` | Check if stack is empty |
| `size()` | Return number of elements |

### ASCII Art Visualization

```
         push(5)     push(3)     pop()       push(7)
           ↓           ↓          ↑            ↓
         ┌───┐       ┌───┐      ┌───┐       ┌───┐
         │   │       │ 3 │      │   │       │ 7 │
         ├───┤       ├───┤      ├───┤       ├───┤
         │ 5 │       │ 5 │      │ 5 │       │ 5 │
         └───┘       └───┘      └───┘       └───┘
        top=0       top=1      top=0       top=1
                              returns 3
```

**Internal array representation:**
```
Array: [5, 3, _, _, _]     After pop(): [5, _, _, _, _]
        ↑  ↑                             ↑
      [0][1]                           [0]
        top=1                          top=0
```

### 1.1 Array-Based Stack Implementation

```python
from typing import Any, Optional, List


class ArrayStack:
    """
    Stack implementation using a dynamic array (Python list).
    
    The stack follows LIFO (Last-In-First-Out) principle.
    Elements are added and removed from the end of the array
    for O(1) amortized time complexity.
    
    Attributes:
        _items: Internal list storing stack elements
    
    Time Complexity:
        - push: O(1) amortized
        - pop: O(1)
        - peek: O(1)
        - is_empty: O(1)
    
    Space Complexity: O(n) where n is the number of elements
    """
    
    def __init__(self) -> None:
        """Initialize an empty stack."""
        self._items: List[Any] = []
    
    def push(self, item: Any) -> None:
        """
        Add an item to the top of the stack.
        
        Args:
            item: The element to be added to the stack
        """
        self._items.append(item)
    
    def pop(self) -> Any:
        """
        Remove and return the top item from the stack.
        
        Returns:
            The top element of the stack
        
        Raises:
            IndexError: If the stack is empty
        """
        if self.is_empty():
            raise IndexError("Pop from empty stack")
        return self._items.pop()
    
    def peek(self) -> Any:
        """
        Return the top item without removing it.
        
        Returns:
            The top element of the stack
        
        Raises:
            IndexError: If the stack is empty
        """
        if self.is_empty():
            raise IndexError("Peek from empty stack")
        return self._items[-1]
    
    def is_empty(self) -> bool:
        """Return True if the stack is empty."""
        return len(self._items) == 0
    
    def size(self) -> int:
        """Return the number of elements in the stack."""
        return len(self._items)
    
    def __len__(self) -> int:
        """Return the number of elements (enables len(stack))."""
        return len(self._items)
    
    def __bool__(self) -> bool:
        """Return True if stack is non-empty (enables if stack:)."""
        return len(self._items) > 0
    
    def __str__(self) -> str:
        """Return string representation of the stack."""
        if not self._items:
            return "Stack: [] (empty)"
        items_str = " ".join(str(item) for item in self._items)
        return f"Stack: [{items_str}] ← top"
    
    def __repr__(self) -> str:
        """Return detailed string representation."""
        return f"ArrayStack({self._items})"
```

```python
# Demonstration: Array-based Stack
print("=" * 50)
print("Array-Based Stack Demonstration")
print("=" * 50)

stack = ArrayStack()
print(f"\nInitial state: {stack}")
print(f"Is empty? {stack.is_empty()}")

# Push operations
for value in [10, 20, 30, 40]:
    stack.push(value)
    print(f"push({value}): {stack}")

print(f"\nPeek: {stack.peek()}")
print(f"Size: {stack.size()}")

# Pop operations
print("\nPopping elements:")
while not stack.is_empty():
    value = stack.pop()
    print(f"pop() → {value}, {stack}")
```

### 1.2 Linked-List-Based Stack Implementation

**Advantages over array-based:**
- No need to resize (truly O(1) push, not amortized)
- No wasted space from pre-allocated capacity

**Disadvantages:**
- Extra memory for node pointers
- No random access (not needed for stacks)
- Slightly worse cache locality

```python
class StackNode:
    """
    A node in the linked-list-based stack.
    
    Attributes:
        value: The data stored in this node
        next: Reference to the next node (towards bottom of stack)
    """
    
    def __init__(self, value: Any, next_node: Optional['StackNode'] = None) -> None:
        self.value = value
        self.next = next_node


class LinkedStack:
    """
    Stack implementation using a singly linked list.
    
    Elements are added and removed from the head of the list,
    which represents the top of the stack.
    
    Attributes:
        _top: Reference to the top node of the stack
        _size: Number of elements in the stack
    
    Time Complexity:
        - push: O(1)
        - pop: O(1)
        - peek: O(1)
        - is_empty: O(1)
    
    Space Complexity: O(n) where n is the number of elements
    """
    
    def __init__(self) -> None:
        """Initialize an empty stack."""
        self._top: Optional[StackNode] = None
        self._size: int = 0
    
    def push(self, item: Any) -> None:
        """
        Add an item to the top of the stack.
        
        Creates a new node and makes it the new top.
        
        Args:
            item: The element to be added
        """
        new_node = StackNode(item, self._top)
        self._top = new_node
        self._size += 1
    
    def pop(self) -> Any:
        """
        Remove and return the top item from the stack.
        
        Returns:
            The top element
        
        Raises:
            IndexError: If the stack is empty
        """
        if self.is_empty():
            raise IndexError("Pop from empty stack")
        value = self._top.value
        self._top = self._top.next
        self._size -= 1
        return value
    
    def peek(self) -> Any:
        """
        Return the top item without removing it.
        
        Returns:
            The top element
        
        Raises:
            IndexError: If the stack is empty
        """
        if self.is_empty():
            raise IndexError("Peek from empty stack")
        return self._top.value
    
    def is_empty(self) -> bool:
        """Return True if the stack is empty."""
        return self._top is None
    
    def size(self) -> int:
        """Return the number of elements in the stack."""
        return self._size
    
    def __len__(self) -> int:
        """Return the number of elements."""
        return self._size
    
    def __bool__(self) -> bool:
        """Return True if stack is non-empty."""
        return self._top is not None
    
    def __str__(self) -> str:
        """Return string representation of the stack."""
        if self.is_empty():
            return "Stack: [] (empty)"
        
        elements = []
        current = self._top
        while current:
            elements.append(str(current.value))
            current = current.next
        
        # Reverse to show bottom-to-top order
        return f"Stack: [{' '.join(reversed(elements))}] ← top"
```

```python
# Demonstration: Linked-list-based Stack
print("=" * 50)
print("Linked-List-Based Stack Demonstration")
print("=" * 50)

stack = LinkedStack()
print(f"\nInitial state: {stack}")

# Push operations
for value in ['A', 'B', 'C', 'D']:
    stack.push(value)
    print(f"push('{value}'): {stack}")

print(f"\nPeek: {stack.peek()}")
print(f"Size: {stack.size()}")

# Pop operations
print("\nPopping elements:")
while stack:
    value = stack.pop()
    print(f"pop() → '{value}', {stack}")
```

---

## 2. Queue (FIFO)

### What is a Queue?

A **queue** is a linear data structure that follows the **First-In-First-Out (FIFO)** principle. The first element added is the first one to be removed.

### Real-World Analogies

- **Line at a store**: First person in line gets served first
- **Print queue**: Documents print in the order they were sent
- **Traffic at intersection**: Cars pass in arrival order
- **Customer service calls**: Calls answered in order received

### When to Use?

- **Breadth-First Search (BFS)**
- **Task scheduling** (process scheduling, print jobs)
- **Message buffers** (producer-consumer problems)
- **Simulation systems** (event-driven simulation)
- **Caching** (LRU cache eviction)

### Core Operations

| Operation | Description |
|-----------|-------------|
| `enqueue(item)` | Add item to the rear |
| `dequeue()` | Remove and return front item |
| `front()` / `peek()` | Return front item without removing |
| `is_empty()` | Check if queue is empty |
| `size()` | Return number of elements |

### ASCII Art Visualization

```
enqueue(5)  enqueue(3)  enqueue(8)   dequeue()    enqueue(2)
    ↓           ↓           ↓           ↑            ↓
  ┌───┐      ┌───┬───┐   ┌───┬───┬───┐  ┌───┬───┐   ┌───┬───┬───┐
  │ 5 │      │ 5 │ 3 │   │ 5 │ 3 │ 8 │  │ 3 │ 8 │   │ 3 │ 8 │ 2 │
  └───┘      └───┴───┘   └───┴───┴───┘  └───┴───┘   └───┴───┴───┘
  front       front       front         front        front
  rear        rear        rear          rear         rear
                                       returns 5
```

**Linked-list representation:**
```
front                                  rear
  ↓                                     ↓
┌───┐    ┌───┐    ┌───┐    ┌───┐
│ 5 │───►│ 3 │───►│ 8 │───►│ 2 │───► None
└───┘    └───┘    └───┘    └───┘

dequeue(): remove from front (5)
enqueue(): add at rear
```

```python
class QueueNode:
    """
    A node in the linked-list-based queue.
    
    Attributes:
        value: The data stored in this node
        next: Reference to the next node (towards rear)
    """
    
    def __init__(self, value: Any, next_node: Optional['QueueNode'] = None) -> None:
        self.value = value
        self.next = next_node


class LinkedQueue:
    """
    Queue implementation using a singly linked list.
    
    Maintains front and rear pointers for O(1) operations.
    Elements are added at the rear and removed from the front.
    
    Attributes:
        _front: Reference to the front node (for dequeue)
        _rear: Reference to the rear node (for enqueue)
        _size: Number of elements in the queue
    
    Time Complexity:
        - enqueue: O(1)
        - dequeue: O(1)
        - peek: O(1)
        - is_empty: O(1)
    
    Space Complexity: O(n) where n is the number of elements
    """
    
    def __init__(self) -> None:
        """Initialize an empty queue."""
        self._front: Optional[QueueNode] = None
        self._rear: Optional[QueueNode] = None
        self._size: int = 0
    
    def enqueue(self, item: Any) -> None:
        """
        Add an item to the rear of the queue.
        
        Args:
            item: The element to be added
        """
        new_node = QueueNode(item)
        
        if self.is_empty():
            self._front = self._rear = new_node
        else:
            self._rear.next = new_node
            self._rear = new_node
        
        self._size += 1
    
    def dequeue(self) -> Any:
        """
        Remove and return the front item from the queue.
        
        Returns:
            The front element
        
        Raises:
            IndexError: If the queue is empty
        """
        if self.is_empty():
            raise IndexError("Dequeue from empty queue")
        
        value = self._front.value
        self._front = self._front.next
        self._size -= 1
        
        # If queue becomes empty, update rear as well
        if self._front is None:
            self._rear = None
        
        return value
    
    def peek(self) -> Any:
        """
        Return the front item without removing it.
        
        Returns:
            The front element
        
        Raises:
            IndexError: If the queue is empty
        """
        if self.is_empty():
            raise IndexError("Peek from empty queue")
        return self._front.value
    
    def is_empty(self) -> bool:
        """Return True if the queue is empty."""
        return self._front is None
    
    def size(self) -> int:
        """Return the number of elements in the queue."""
        return self._size
    
    def __len__(self) -> int:
        """Return the number of elements."""
        return self._size
    
    def __bool__(self) -> bool:
        """Return True if queue is non-empty."""
        return self._front is not None
    
    def __str__(self) -> str:
        """Return string representation of the queue."""
        if self.is_empty():
            return "Queue: [] (empty)"
        
        elements = []
        current = self._front
        while current:
            elements.append(str(current.value))
            current = current.next
        
        return f"Queue: front → [{' '.join(elements)}] ← rear"
```
