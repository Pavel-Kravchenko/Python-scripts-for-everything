---
name: algo-hash-tables-bloom
description: "Hash tables (chaining vs open addressing) and Bloom filters: complexity, trade-offs, and implementation patterns."
tool_type: python
primary_tool: Python
---

# Hash Tables and Bloom Filters

## Complexity

| Operation | Average | Worst |
|-----------|---------|-------|
| Insert / Search / Delete | O(1) | O(n) |
| Space | O(n) | O(n) |

Worst case occurs when all keys collide (pathological hash function or adversarial input).

## Hash Functions

```python
def mod_hash(key: int, size: int) -> int:
    return key % size                      # size should be prime

def poly_hash(key: str, size: int) -> int:
    h = 0
    for ch in key:
        h = h * 31 + ord(ch)
    return h % size

def mult_hash(key: int, size: int) -> int:
    A = 0.6180339887  # (sqrt(5)-1)/2
    return int(size * ((key * A) % 1))
```

## Collision Resolution

### Chaining (Separate Chaining)
Each bucket holds a linked list. Load factor can exceed 1.0.

```python
class HashTableChaining:
    def __init__(self, size: int = 7):
        self.size = size
        self.buckets: list[list] = [[] for _ in range(size)]
        self.count = 0

    def _hash(self, key) -> int:
        return hash(key) % self.size

    def put(self, key, value) -> None:
        i = self._hash(key)
        for j, (k, _) in enumerate(self.buckets[i]):
            if k == key:
                self.buckets[i][j] = (key, value)
                return
        self.buckets[i].append((key, value))
        self.count += 1

    def get(self, key):
        for k, v in self.buckets[self._hash(key)]:
            if k == key:
                return v
        return None

    def delete(self, key) -> bool:
        i = self._hash(key)
        for j, (k, _) in enumerate(self.buckets[i]):
            if k == key:
                self.buckets[i].pop(j)
                self.count -= 1
                return True
        return False

    def load_factor(self) -> float:
        return self.count / self.size
```

### Open Addressing Probe Sequences

```python
# Linear probing
h(key, i) = (h(key) + i) % m

# Quadratic probing
h(key, i) = (h(key) + i*i) % m

# Double hashing (best distribution)
h(key, i) = (h1(key) + i * h2(key)) % m
# h2 must never return 0: h2(k) = 1 + (k % (m-1))
```

### Chaining vs Open Addressing

| Aspect | Chaining | Open Addressing |
|--------|----------|----------------|
| Load factor | Can exceed 1.0 | Must stay < 1.0 |
| Cache | Poor (pointer chasing) | Better (contiguous) |
| Delete | Simple | Needs tombstones |
| Clustering | None | Primary/secondary |

## Load Factor & Rehashing

- Rehash when load factor > 0.7 (chaining) or > 0.5 (open addressing)
- Double the table size (use next prime) and reinsert all keys

## Bloom Filter

Probabilistic set: **no false negatives, possible false positives**.

```python
import mmh3
from bitarray import bitarray

class BloomFilter:
    def __init__(self, capacity: int, error_rate: float = 0.01):
        self.size = self._optimal_size(capacity, error_rate)
        self.hash_count = self._optimal_hashes(self.size, capacity)
        self.bits = bitarray(self.size)
        self.bits.setall(0)

    def _optimal_size(self, n, p):
        import math
        return int(-n * math.log(p) / (math.log(2) ** 2))

    def _optimal_hashes(self, m, n):
        import math
        return max(1, int((m / n) * math.log(2)))

    def add(self, item: str) -> None:
        for seed in range(self.hash_count):
            self.bits[mmh3.hash(item, seed) % self.size] = 1

    def __contains__(self, item: str) -> bool:
        return all(
            self.bits[mmh3.hash(item, seed) % self.size]
            for seed in range(self.hash_count)
        )
```

## Pitfalls

- **Table size should be prime**: reduces clustering in modular hash functions.
- **Open addressing requires tombstones on delete**: simply clearing a slot breaks search chains for keys inserted after a collision at that slot.
- **Never resize a Bloom filter**: it cannot be resized without rebuilding from scratch; pre-size using the capacity formula.
- **Python `hash()` is randomized per-process**: do not persist Python's `hash()` across runs; use a stable hash (e.g., `mmh3`, `hashlib`) for on-disk or cross-process use.
- **Load factor threshold matters**: at load 0.9 with linear probing, expected probe length exceeds 5; keep load < 0.7.
