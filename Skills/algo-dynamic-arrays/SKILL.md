---
name: algo-dynamic-arrays
description: "A comprehensive guide to dynamic arrays: how they work, why we need them, and their performance characteristics."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/04_Linear_Data_Structures/03_dynamic_arrays.ipynb"
---

# Dynamic Arrays

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/04_Linear_Data_Structures/03_dynamic_arrays.ipynb`*

# Dynamic Arrays

A comprehensive guide to dynamic arrays: how they work, why we need them, and their performance characteristics.

**Topics Covered:**
1. Static vs Dynamic Arrays
2. Dynamic Array Implementation
3. Amortized Analysis
4. Python's `list` Under the Hood

---

## 1. Static vs Dynamic Arrays

### Why Can't We Just Use Fixed-Size Arrays?

**Static arrays** have a fixed size determined at creation time. This leads to two problems:

1. **Wasted Memory**: If we allocate too much space, we waste memory
2. **Overflow**: If we allocate too little, we can't store all our data

```
Static Array Problem:

Case 1: Over-allocation (waste)
┌───┬───┬───┬───┬───┬───┬───┬───┐
│ A │ B │ C │   │   │   │   │   │  ← 5 slots wasted!
└───┴───┴───┴───┴───┴───┴───┴───┘

Case 2: Under-allocation (overflow)
┌───┬───┬───┬───┐
│ A │ B │ C │ D │  ← FULL! Can't add 'E'!
└───┴───┴───┴───┘
```

**Dynamic arrays** solve this by automatically resizing when needed.

### Real-World Analogy

Think of it like a parking lot:
- **Static**: Fixed number of spaces. When full, cars must go elsewhere.
- **Dynamic**: When lot fills up, you automatically build a bigger lot and move all cars there.

---

## 2. How Does Automatic Resizing Work?

When a dynamic array runs out of space:

1. **Allocate** a new, larger array
2. **Copy** all existing elements to the new array
3. **Free** the old array (in languages with manual memory management)
4. **Insert** the new element

```
Array Resizing Visualization:

Initial: capacity=4, size=4
┌───┬───┬───┬───┐
│ 1 │ 2 │ 3 │ 4 │  ← FULL!
└───┴───┴───┴───┘

append(5) triggers resize:

Step 1: Allocate new array (2x capacity)
┌───┬───┬───┬───┬───┬───┬───┬───┐
│   │   │   │   │   │   │   │   │  capacity=8
└───┴───┴───┴───┴───┴───┴───┴───┘

Step 2: Copy elements
┌───┬───┬───┬───┬───┬───┬───┬───┐
│ 1 │ 2 │ 3 │ 4 │   │   │   │   │
└───┴───┴───┴───┴───┴───┴───┴───┘

Step 3: Add new element
┌───┬───┬───┬───┬───┬───┬───┬───┐
│ 1 │ 2 │ 3 │ 4 │ 5 │   │   │   │  size=5
└───┴───┴───┴───┴───┴───┴───┴───┘
```

### Growth Factor Trade-offs

| Growth Factor | Pros | Cons |
|--------------|------|------|
| **2x** (doubling) | Fewer resizes, simpler math | More wasted space (up to 50%) |
| **1.5x** | Less wasted space (~33%) | More frequent resizes |
| **1.25x** | Even less waste | Many more resizes |

**Common choices:**
- Python `list`: ~1.125x (complex formula, very space-efficient)
- Java `ArrayList`: 1.5x
- C++ `vector`: 2x (GCC) or 1.5x (MSVC)

---

## 3. Dynamic Array Implementation

### Original Implementation (from Algorithms_HW7)

The basic `DynamicMassive` class uses numpy arrays as the underlying storage:

```python
import numpy as np

# Original DynamicMassive from Algorithms_HW7
class DynamicMassive:
    """Basic dynamic array using numpy for storage."""
    
    def __init__(self):
        self.count = 0
        self.size = 1
        self.A = self._make_array(self.size)
    
    def __len__(self):
        return self.count
    
    def __getitem__(self, k):
        if k < 0 or k >= self.count:
            raise IndexError('Index out of bounds!')
        return self.A[k]
    
    def append(self, elem):
        if self.count == self.size:
            # Double the capacity
            B = self._make_array(2 * self.size)
            for i in range(self.count):
                B[i] = self.A[i]
            self.A = B
            self.size = 2 * self.size
        
        self.A[self.count] = elem
        self.count += 1
    
    def _make_array(self, new_cap):
        return np.empty(new_cap, dtype=int, order='C')
```

### Enhanced Implementation with Full Features

Here's a complete implementation with all standard operations:

```python
from typing import Any, Iterator, Optional
import ctypes


class DynamicArray:
    """
    A dynamic array implementation that automatically resizes.
    
    Resize Strategy:
    - Grow: When full, double the capacity (2x growth factor)
    - Shrink: When 1/4 full, halve the capacity (to avoid thrashing)
    
    The 2x growth factor ensures O(1) amortized append operations.
    The 1/4 shrink threshold prevents resize thrashing at boundaries.
    
    Attributes:
        _size: Number of elements currently stored
        _capacity: Total slots available in underlying array
        _array: Raw ctypes array for storage
    """
    
    GROWTH_FACTOR: float = 2.0
    SHRINK_THRESHOLD: float = 0.25  # Shrink when size <= capacity * threshold
    MIN_CAPACITY: int = 1
    
    def __init__(self, initial_capacity: int = 1) -> None:
        """
        Initialize an empty dynamic array.
        
        Args:
            initial_capacity: Starting capacity (default: 1)
        """
        self._size: int = 0
        self._capacity: int = max(initial_capacity, self.MIN_CAPACITY)
        self._array: ctypes.Array = self._make_array(self._capacity)
    
    def __len__(self) -> int:
        """Return the number of elements in the array."""
        return self._size
    
    def __getitem__(self, index: int) -> Any:
        """
        Get element at index. Supports negative indexing.
        
        Args:
            index: Position to access (-size <= index < size)
            
        Returns:
            Element at the specified index
            
        Raises:
            IndexError: If index is out of bounds
        """
        if index < 0:
            index += self._size
        if not 0 <= index < self._size:
            raise IndexError(f'Index {index} out of bounds for size {self._size}')
        return self._array[index]
    
    def __setitem__(self, index: int, value: Any) -> None:
        """
        Set element at index. Supports negative indexing.
        
        Args:
            index: Position to modify (-size <= index < size)
            value: New value to store
            
        Raises:
            IndexError: If index is out of bounds
        """
        if index < 0:
            index += self._size
        if not 0 <= index < self._size:
            raise IndexError(f'Index {index} out of bounds for size {self._size}')
        self._array[index] = value
    
    def __iter__(self) -> Iterator[Any]:
        """Iterate over elements in the array."""
        for i in range(self._size):
            yield self._array[i]
    
    def __repr__(self) -> str:
        """String representation showing contents."""
        elements = [str(self._array[i]) for i in range(self._size)]
        return f"DynamicArray([{', '.join(elements)}])"
    
    @property
    def capacity(self) -> int:
        """Current capacity of the underlying array."""
        return self._capacity
    
    def append(self, value: Any) -> None:
        """
        Add element to the end of the array.
        
        If the array is full, it will be resized to double its capacity.
        This operation is O(1) amortized, O(n) worst case.
        
        Args:
            value: Element to append
        """
        if self._size == self._capacity:
            self._resize(int(self._capacity * self.GROWTH_FACTOR))
        
        self._array[self._size] = value
        self._size += 1
    
    def insert(self, index: int, value: Any) -> None:
        """
        Insert element at specified index.
        
        Elements at and after the index are shifted right.
        This operation is O(n) due to shifting.
        
        Args:
            index: Position to insert at (0 <= index <= size)
            value: Element to insert
            
        Raises:
            IndexError: If index is out of bounds
        """
        if index < 0:
            index += self._size
        if not 0 <= index <= self._size:
            raise IndexError(f'Insert index {index} out of bounds')
        
        if self._size == self._capacity:
            self._resize(int(self._capacity * self.GROWTH_FACTOR))
        
        # Shift elements right
        for i in range(self._size, index, -1):
            self._array[i] = self._array[i - 1]
        
        self._array[index] = value
        self._size += 1
    
    def delete(self, index: int) -> Any:
        """
        Remove and return element at specified index.
        
        Elements after the index are shifted left.
        If size drops below 1/4 capacity, array is shrunk.
        This operation is O(n) due to shifting.
        
        Args:
            index: Position to delete (-size <= index < size)
            
        Returns:
            The removed element
            
        Raises:
            IndexError: If index is out of bounds or array is empty
        """
        if self._size == 0:
            raise IndexError('Cannot delete from empty array')
        
        if index < 0:
            index += self._size
        if not 0 <= index < self._size:
            raise IndexError(f'Delete index {index} out of bounds')
        
        value = self._array[index]
        
        # Shift elements left
        for i in range(index, self._size - 1):
            self._array[i] = self._array[i + 1]
        
        self._size -= 1
        
        # Shrink if too empty (but maintain minimum capacity)
        if (self._size <= self._capacity * self.SHRINK_THRESHOLD 
            and self._capacity > self.MIN_CAPACITY):
            self._resize(max(self._capacity // 2, self.MIN_CAPACITY))
        
        return value
    
    def pop(self) -> Any:
        """
        Remove and return the last element.
        
        Returns:
            The removed element
            
        Raises:
            IndexError: If array is empty
        """
        return self.delete(self._size - 1)
    
    def search(self, value: Any) -> int:
        """
        Find the index of the first occurrence of value.
        
        Args:
            value: Element to search for
            
        Returns:
            Index of the element, or -1 if not found
        """
        for i in range(self._size):
            if self._array[i] == value:
                return i
        return -1
    
    def _resize(self, new_capacity: int) -> None:
        """
        Resize the underlying array to new_capacity.
        
        This is a private method that:
        1. Creates a new array with the specified capacity
        2. Copies all existing elements
        3. Replaces the old array
        
        Time Complexity: O(n) where n is current size
        
        Args:
            new_capacity: The new capacity (must be >= size)
        """
        new_array = self._make_array(new_capacity)
        
        for i in range(self._size):
            new_array[i] = self._array[i]
        
        self._array = new_array
        self._capacity = new_capacity
    
    def _make_array(self, capacity: int) -> ctypes.Array:
        """
        Create a new raw array of Python objects.
        
        Uses ctypes to create a low-level array that can hold
        references to any Python object.
        
        Args:
            capacity: Size of the array to create
            
        Returns:
            A ctypes array of py_object
        """
        return (capacity * ctypes.py_object)()
```

---

## 4. Complexity Analysis

### Time & Space Complexity Table

| Operation | Worst Case | Amortized | Space |
|-----------|------------|-----------|-------|
| Access `[i]` | O(1) | O(1) | - |
| Append | O(n) | **O(1)** | O(n) |
| Insert at i | O(n) | O(n) | - |
| Delete at i | O(n) | O(n) | - |
| Search | O(n) | O(n) | - |
| Pop (last) | O(n)* | **O(1)** | - |

*Pop is O(n) worst case only when shrinking occurs

---

## 5. Amortized Analysis: Why Append is O(1)

### The Accounting Method

Let's trace what happens when we append 9 elements:

```
Amortized Analysis of append() with 2x growth factor:

Op  Size  Capacity  Cost  Cumulative  Notes
─────────────────────────────────────────────────────────
1   1     1         1     1           
2   2     2         2     3           resize: copy 1 + insert 1
3   3     4         2     5           resize: copy 2 + insert 1
4   4     4         1     6           
5   5     8         5     11          resize: copy 4 + insert 1
6   6     8         1     12          
7   7     8         1     13          
8   8     8         1     14          
9   9     16        9     23          resize: copy 8 + insert 1

Total cost: 23 for 9 operations
Amortized: 23/9 ≈ 2.56 = O(1) per operation!
```

### Why Does This Work?

The key insight: **We resize infrequently enough that the expensive copies are "paid for" by the cheap inserts.**

```
Cost breakdown for n appends:

Resize copies: 1 + 2 + 4 + 8 + ... + n/2 = n - 1 copies
Regular inserts: n inserts (cost 1 each)

Total cost ≤ n + (n - 1) = 2n - 1
Amortized cost = (2n - 1) / n ≈ 2 = O(1)
```
