"""
Dynamic array (resizable array) implementation.

Backed by a fixed-capacity NumPy integer array.  When the array is full,
a new array of double the capacity is allocated and all elements are copied.
This gives amortized O(1) append.
"""
from __future__ import annotations

import numpy


class DynamicArray:
    """Resizable array of integers backed by a NumPy buffer.

    The internal buffer doubles in size whenever capacity is exhausted.
    """

    def __init__(self):
        self.count = 0     # number of elements currently stored
        self._size = 1     # current allocated capacity
        self._A = self._make_array(self._size)

    def __len__(self) -> int:
        return self.count

    def capacity(self) -> int:
        """Return the current allocated capacity."""
        return self._size

    def __getitem__(self, k: int) -> int:
        if k < 0 or k >= self.count:
            raise IndexError("Index out of range")
        return self._A[k]

    def append(self, elem: int) -> None:
        """Append elem to the end of the array, resizing if necessary."""
        if self.count == self._size:
            # Buffer is full — allocate a new one with double capacity
            new_buf = self._make_array(2 * self._size)
            for i in range(self.count):
                new_buf[i] = self._A[i]
            self._A = new_buf
            self._size = 2 * self._size

        self._A[self.count] = elem
        self.count += 1

    @staticmethod
    def _make_array(capacity: int):
        """Allocate a raw NumPy integer array of the given capacity."""
        return numpy.empty(capacity, dtype=int, order="C")


if __name__ == "__main__":
    da = DynamicArray()

    # Append elements and observe capacity doublings
    for i in range(20):
        da.append(i * 10)
        print(f"len={len(da):2d}  capacity={da.capacity():4d}  last={da[len(da)-1]}")

    print(f"\nElement at index 5: {da[5]}")
