"""
Hash table with open addressing (linear probing).

Collision resolution: when a slot is occupied, the table scans forward
linearly (wrapping around) until an empty slot or the same key is found.

Two parallel arrays are used:
- cells: stores the keys
- data:  stores the associated values
"""
from __future__ import annotations

import random
import string
from random import randint


class HashTable:
    """Hash table using linear-probing open addressing.

    Args:
        size: Number of slots in the table.
    """

    def __init__(self, size: int = 100):
        self.size = size
        self.cells: list = [None] * self.size
        self.data: list = [None] * self.size

    def put(self, key, value) -> None:
        """Insert or update the key-value pair."""
        slot = hash(key) % self.size

        if self.cells[slot] is None:
            # Slot is empty — insert directly
            self.cells[slot] = key
            self.data[slot] = value
        elif self.cells[slot] == key:
            # Key already present — update value
            self.data[slot] = value
        else:
            # Collision: scan forward with linear probing
            next_slot = (slot + 1) % self.size
            steps = 0
            while (self.cells[next_slot] is not None
                   and self.cells[next_slot] != key
                   and steps < 2 * self.size):
                next_slot = (next_slot + 1) % self.size
                steps += 1

            if self.cells[next_slot] is None:
                self.cells[next_slot] = key
                self.data[next_slot] = value
            else:
                # Key found in probing chain — update value
                self.data[next_slot] = value

    def get(self, key):
        """Return the value for key, or None if the key is not present."""
        start = hash(key) % self.size
        position = start
        found = False
        stop = False
        value = None

        while self.cells[position] is not None and not found and not stop:
            if self.cells[position] == key:
                found = True
                value = self.data[position]
            else:
                position = (position + 1) % self.size
                if position == start:
                    # Wrapped all the way around — key not present
                    stop = True

        return value

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        self.put(key, value)


if __name__ == "__main__":
    ht = HashTable(size=16)

    # Insert some DNA k-mer frequencies
    kmers = ["ATGC", "GCTA", "TTTT", "AAAA", "CCGG", "ATGC"]
    for kmer in kmers:
        existing = ht.get(kmer)
        ht[kmer] = (existing or 0) + 1

    print(f"'ATGC' count: {ht['ATGC']}")   # 2 (inserted twice)
    print(f"'GCTA' count: {ht['GCTA']}")   # 1
    print(f"'TATA' count: {ht['TATA']}")   # None (not inserted)

    # Stress test: insert random key-value pairs
    big = HashTable(size=200)
    test_keys = [randint(0, 100) for _ in range(100)]
    for k in test_keys:
        w = "".join(random.choice(string.ascii_letters) for _ in range(5))
        big[k] = w
    last_key = test_keys[-1]
    print(f"Last inserted key {last_key} -> {big[last_key]}")
