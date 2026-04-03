"""
Hash table with separate chaining using a singly-linked list.

Each bucket in the table holds a LinkedList of (key, value) pairs.
Collisions are resolved by appending to the chain at the target bucket.
"""
from __future__ import annotations

import random
import string
from random import randint


# ---------------------------------------------------------------------------
# Singly-linked list used as chain storage
# ---------------------------------------------------------------------------

class _ListNode:
    """Node in the singly-linked list."""

    def __init__(self, value=None, next: "_ListNode | None" = None):
        self.value = value  # stored as [key, data]
        self.next = next


class LinkedList:
    """Minimal singly-linked list supporting append and key-based lookup."""

    def __init__(self):
        self.first: _ListNode | None = None
        self.last: _ListNode | None = None
        self.length = 0

    def add(self, elem) -> None:
        """Append elem to the end of the list."""
        self.length += 1
        if self.first is None:
            self.last = self.first = _ListNode(elem)
        else:
            self.last.next = self.last = _ListNode(elem)

    def search_by_key(self, key):
        """Return the [key, value] pair whose first element matches key,
        or -1 if not found.
        """
        if self.first is None:
            return -1
        curr = self.first
        while curr is not None:
            if curr.value[0] == key:
                return curr.value
            curr = curr.next
        return -1

    def search_by_index(self, index: int):
        """Return the value stored at the given 0-based index, or -1."""
        if self.first is None:
            return -1
        curr = self.first
        count = 0
        while curr is not None:
            if count == index:
                return curr.value
            curr = curr.next
            count += 1
        return -1

    def delete_by_value(self, data) -> int:
        """Remove the first node whose value equals data.
        Returns -1 if the list is empty or data is not found.
        """
        if self.first is None:
            return -1
        curr = self.first
        prev = None
        while curr.value != data and curr.next is not None:
            prev = curr
            curr = curr.next
        if curr.value == data:
            if curr is self.first:
                self.first = curr.next if curr.next else None
            else:
                prev.next = curr.next if curr.next else None
                if curr.next is None:
                    self.last = prev
            return 0
        return -1

    def reverse(self) -> None:
        """Reverse the list in place."""
        curr = self.first
        prev = None
        while curr is not None:
            nxt = curr.next
            curr.next = prev
            prev = curr
            curr = nxt
        self.first = prev

    def __str__(self) -> str:
        if self.first is None:
            return "LinkedList []"
        curr = self.first
        parts = [str(curr.value)]
        while curr.next is not None:
            curr = curr.next
            parts.append(str(curr.value))
        return "LinkedList [ " + " ".join(parts) + " ]"


# ---------------------------------------------------------------------------
# Hash table with chaining
# ---------------------------------------------------------------------------

class HashTable:
    """Hash table using separate chaining for collision resolution.

    Each bucket is a LinkedList of [key, value] pairs.

    Args:
        size: Number of buckets.
    """

    def __init__(self, size: int = 100):
        self.size = size
        self.cells: list[LinkedList] = [LinkedList() for _ in range(self.size)]

    def _bucket(self, key) -> int:
        return hash(key) % self.size

    def put(self, key, value) -> None:
        """Insert key-value pair.  Duplicate keys result in a second entry
        in the chain (original student behaviour).
        """
        bucket_idx = self._bucket(key)
        bucket = self.cells[bucket_idx]
        if bucket.first is None:
            bucket.add([key, value])
        else:
            # Original: only skip if the first node's value equals the new data
            if bucket.first.value == value:
                return
            bucket.add([key, value])

    def get(self, key):
        """Return the [key, value] pair for the given key, or -1 if absent."""
        bucket_idx = self._bucket(key)
        return self.cells[bucket_idx].search_by_key(key)

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        self.put(key, value)


if __name__ == "__main__":
    ht = HashTable(size=16)

    # Simple key-value insertions
    ht["ATGC"] = "adenine-rich"
    ht["GCTA"] = "gc-rich"
    ht["TTTT"] = "poly-T"

    print(f"get 'ATGC': {ht['ATGC']}")
    print(f"get 'GCTA': {ht['GCTA']}")
    print(f"get 'AAAA': {ht['AAAA']}")  # -1 (not inserted)

    # Stress test
    big = HashTable(size=200)
    test_keys = [randint(0, 100) for _ in range(100)]
    for k in test_keys:
        w = "".join(random.choice(string.ascii_letters) for _ in range(5))
        big[k] = w
    last_key = test_keys[-1]
    print(f"Last inserted key {last_key} -> {big[last_key]}")
