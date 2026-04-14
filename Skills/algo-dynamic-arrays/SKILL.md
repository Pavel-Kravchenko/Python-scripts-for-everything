---
name: algo-dynamic-arrays
description: "Dynamic arrays: amortized O(1) append via doubling, O(n) insert/delete, and Python list internals."
tool_type: python
primary_tool: NumPy
---

# Dynamic Arrays

## Complexity

| Operation | Worst | Amortized |
|-----------|-------|-----------|
| Access `[i]` | O(1) | O(1) |
| Append | O(n) | **O(1)** |
| Insert at i | O(n) | O(n) |
| Delete at i | O(n) | O(n) |
| Pop (last) | O(n)* | **O(1)** |

*O(n) only when a shrink occurs.

## Growth Factor Trade-offs

| Factor | Memory waste | Resize frequency |
|--------|-------------|-----------------|
| 2x | up to 50% | least |
| 1.5x | up to 33% | moderate |
| Python list | ~12% (complex formula) | most efficient |

Python `list` uses ~1.125x growth. Java `ArrayList` = 1.5x. C++ `vector` = 2x.

## Amortized Analysis (2x growth)

Total cost of n appends = n inserts + (1 + 2 + 4 + ... + n/2) copies = n + (n-1) ≈ 2n.
Amortized cost per append = 2n/n = **O(1)**.

## Shrink Threshold

Shrink at 1/4 capacity (not 1/2) to avoid thrashing at the boundary when alternating append/delete.

## Implementation

```python
import ctypes

class DynamicArray:
    GROWTH = 2.0
    SHRINK_AT = 0.25

    def __init__(self):
        self._size = 0
        self._capacity = 1
        self._data = (1 * ctypes.py_object)()

    def __len__(self) -> int:
        return self._size

    def __getitem__(self, index: int):
        if index < 0:
            index += self._size
        if not 0 <= index < self._size:
            raise IndexError(f"index {index} out of range")
        return self._data[index]

    def __setitem__(self, index: int, value) -> None:
        if index < 0:
            index += self._size
        if not 0 <= index < self._size:
            raise IndexError(f"index {index} out of range")
        self._data[index] = value

    def append(self, value) -> None:
        if self._size == self._capacity:
            self._resize(int(self._capacity * self.GROWTH))
        self._data[self._size] = value
        self._size += 1

    def insert(self, index: int, value) -> None:
        if self._size == self._capacity:
            self._resize(int(self._capacity * self.GROWTH))
        for i in range(self._size, index, -1):
            self._data[i] = self._data[i - 1]
        self._data[index] = value
        self._size += 1

    def pop(self, index: int = -1):
        if self._size == 0:
            raise IndexError("pop from empty array")
        if index < 0:
            index += self._size
        value = self._data[index]
        for i in range(index, self._size - 1):
            self._data[i] = self._data[i + 1]
        self._size -= 1
        if self._size <= self._capacity * self.SHRINK_AT and self._capacity > 1:
            self._resize(max(1, self._capacity // 2))
        return value

    def _resize(self, new_cap: int) -> None:
        new_data = (new_cap * ctypes.py_object)()
        for i in range(self._size):
            new_data[i] = self._data[i]
        self._data = new_data
        self._capacity = new_cap
```

## Python `list` vs `numpy.ndarray`

| Aspect | `list` | `ndarray` |
|--------|--------|-----------|
| Element type | Any Python object | Fixed dtype (e.g. float64) |
| Memory | Pointers + objects | Contiguous raw bytes |
| Append | O(1) amortized | Not supported (copy required) |
| Vectorized math | No | Yes (SIMD) |

Use `ndarray` for numeric data; use `list` when elements are mixed types or you need frequent appends.

## Pitfalls

- **Shrink at 1/4, not 1/2**: shrinking at 1/2 capacity causes O(n) cost on alternating append/pop at the boundary.
- **`ctypes.py_object` array raises `ValueError` on uninitialized access**: only read slots `< self._size`.
- **`numpy.empty` for the original `DynamicMassive`**: `np.empty` leaves memory uninitialized — never read beyond `self.count` or you get garbage values.
- **Python `list.insert(0, x)` is O(n)**: every element shifts right; use `collections.deque` if you need O(1) front inserts.
