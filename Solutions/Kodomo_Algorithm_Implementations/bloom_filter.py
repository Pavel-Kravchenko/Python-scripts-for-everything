"""
Bloom filter implementation.

This is a simplified Bloom filter that uses Python's built-in hash() function
as the single hash function and maps it into a fixed-size cell array.

Note: A production Bloom filter uses multiple independent hash functions to
reduce false positives.  This version uses one hash, which gives a higher
false-positive rate but demonstrates the core idea: O(1) insert and O(1)
membership test with no false negatives.
"""
from __future__ import annotations


class BloomFilter:
    """Single-hash Bloom filter backed by a fixed-size array.

    Args:
        capacity: Number of cells in the backing array.  Larger values
                  reduce the false-positive rate.
    """

    def __init__(self, capacity: int = 1000):
        self.size = capacity
        self.cells: list[str | None] = [None] * self.size

    def _slot(self, key: str) -> int:
        """Map key to a cell index."""
        return hash(key) % self.size

    def add(self, key: str) -> None:
        """Insert key into the filter."""
        self.cells[self._slot(key)] = key

    def might_contain(self, key: str) -> bool:
        """Return True if key was probably inserted, False if definitely not.

        A False result is definitive (no false negatives).
        A True result may be a false positive.
        """
        return self.cells[self._slot(key)] is not None


def bloom_filter_check(words: list[str], word: str, capacity_multiplier: int = 100) -> bool:
    """Insert all words into a filter sized len(words) * capacity_multiplier
    and check whether word might be present.

    Mirrors the original student function signature.
    """
    bf = BloomFilter(capacity=len(words) * capacity_multiplier)
    for w in words:
        bf.add(w)
    return bf.might_contain(word)


if __name__ == "__main__":
    dictionary = ["ATGC", "GCTA", "TTTT", "AAAA", "CCGG"]

    bf = BloomFilter(capacity=500)
    for seq in dictionary:
        bf.add(seq)

    print(f"'ATGC' in filter: {bf.might_contain('ATGC')}")   # True (inserted)
    print(f"'GCTA' in filter: {bf.might_contain('GCTA')}")   # True (inserted)
    print(f"'TATA' in filter: {bf.might_contain('TATA')}")   # False (not inserted)
    print(f"'ACGT' in filter: {bf.might_contain('ACGT')}")   # False (not inserted)
