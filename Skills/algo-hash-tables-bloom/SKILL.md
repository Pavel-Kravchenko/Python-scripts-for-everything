---
name: algo-hash-tables-bloom
description: "This notebook covers hash-based data structures that provide **O(1) average-case** operations for insert, search, and delete."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/06_Hash_Based_Structures/01_hash_tables_bloom.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Hash Tables and Bloom Filters

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/06_Hash_Based_Structures/01_hash_tables_bloom.ipynb`*


This notebook covers hash-based data structures that provide **O(1) average-case** operations for insert, search, and delete.

## Table of Contents
1. [Hash Functions](#1-hash-functions)
2. [Hash Tables - Direct Addressing](#2-hash-tables---direct-addressing)
3. [Collision Resolution](#3-collision-resolution)
4. [Load Factor & Rehashing](#4-load-factor--rehashing)
5. [Bloom Filters](#5-bloom-filters)

---

## 1. Hash Functions

### What is a Hash Function?

A **hash function** maps data of arbitrary size to fixed-size values (hash codes). It converts keys into array indices for fast lookups.

```python
key ──────► h(key) ──────► index
"hello"      hash()        42
"world"      hash()        17
12345        hash()        85
```python

### Properties of a Good Hash Function

| Property | Description |
|----------|-------------|
| **Deterministic** | Same input always produces same output |
| **Uniform Distribution** | Spreads keys evenly across the table |
| **Fast Computation** | O(1) time to compute hash value |
| **Minimizes Collisions** | Different keys rarely map to same index |

### Why O(1) Average Lookup?

With a good hash function and proper table size:
1. Computing `h(key)` takes O(1) time
2. Direct array access at index `h(key)` takes O(1) time
3. If load factor is kept low, collisions are rare → O(1) on average

**Worst case O(n)** occurs when all keys hash to the same index (pathological case).

```python
# Simple Hash Function Examples

def simple_mod_hash(key: int, table_size: int) -> int:
    """
    Basic modulo hash function for integers.
    
    Args:
        key: Integer key to hash
        table_size: Size of the hash table (preferably prime)
    
    Returns:
        Hash index in range [0, table_size-1]
    """
    return key % table_size


def string_hash(key: str, table_size: int) -> int:
    """
    Polynomial rolling hash for strings.
    
    Uses base 31 (common choice for lowercase strings).
    h(s) = s[0] + s[1]*31 + s[2]*31^2 + ... + s[n]*31^n
    
    Args:
        key: String key to hash
        table_size: Size of the hash table
    
    Returns:
        Hash index in range [0, table_size-1]
    """
    hash_value = 0
    base = 31
    for char in key:
        hash_value = hash_value * base + ord(char)
    return hash_value % table_size


def multiplication_hash(key: int, table_size: int) -> int:
    """
    Multiplication method hash (Knuth's suggestion).
    
    h(k) = floor(m * (k * A mod 1))
    where A = (√5 - 1) / 2 ≈ 0.6180339887 (golden ratio conjugate)
    
    Args:
        key: Integer key to hash
        table_size: Size of the hash table (m)
    
    Returns:
        Hash index in range [0, table_size-1]
    """
    A = 0.6180339887  # (sqrt(5) - 1) / 2
    return int(table_size * ((key * A) % 1))


# Demonstrate hash functions
print("=== Hash Function Examples ===")
print(f"\nTable size: 7 (prime)")
print(f"\nModulo hash:")
for key in [15, 22, 8, 1, 29, 36]:
    print(f"  h({key}) = {key} % 7 = {simple_mod_hash(key, 7)}")

print(f"\nString hash (table size 10):")
for word in ["hello", "world", "hash", "table"]:
    print(f"  h('{word}') = {string_hash(word, 10)}")
```python

---

## 2. Hash Tables - Direct Addressing

### Concept

A hash table uses an array where elements are stored at indices determined by a hash function.

```python
Direct Addressing (no collisions):

Keys: {15, 3, 22, 8}
h(key) = key % 10

Index:  0    1    2    3    4    5    6    7    8    9
      ┌────┬────┬────┬────┬────┬────┬────┬────┬────┬────┐
      │    │    │ 22 │  3 │    │ 15 │    │    │  8 │    │
      └────┴────┴────┴────┴────┴────┴────┴────┴────┴────┘
              ↑    ↑         ↑              ↑
         22%10=2  3%10=3   15%10=5       8%10=8
```python

### Complexity Analysis

| Operation | Average Case | Worst Case | Notes |
|-----------|--------------|------------|-------|
| Insert | O(1) | O(n) | Worst when all keys collide |
| Search | O(1) | O(n) | Worst when searching collision chain |
| Delete | O(1) | O(n) | Must handle collision structure |
| Space | O(n) | O(n) | n = number of stored elements |

```python
class SimpleHashTable:
    """
    Simple hash table without collision handling.
    
    Demonstrates basic hash table concept - will overwrite on collision.
    NOT suitable for production use.
    """
    
    def __init__(self, size: int = 10):
        """
        Initialize hash table with given size.
        
        Args:
            size: Number of slots in the table
        """
        self.size = size
        self.keys = [None] * size
        self.values = [None] * size
    
    def _hash(self, key) -> int:
        """Compute hash index for a key."""
        return hash(key) % self.size
    
    def put(self, key, value):
        """Insert or update a key-value pair."""
        index = self._hash(key)
        self.keys[index] = key
        self.values[index] = value
    
    def get(self, key):
        """Retrieve value by key. Returns None if not found."""
        index = self._hash(key)
        if self.keys[index] == key:
            return self.values[index]
        return None
    
    def display(self):
        """Print the hash table contents."""
        print("Index | Key    | Value")
        print("-" * 25)
        for i in range(self.size):
            key_str = str(self.keys[i]) if self.keys[i] is not None else "empty"
            val_str = str(self.values[i]) if self.values[i] is not None else "-"
            print(f"  {i}   | {key_str:6} | {val_str}")


# Demonstrate simple hash table
print("=== Simple Hash Table Demo ===")
ht = SimpleHashTable(7)

# Insert some values
ht.put(15, "apple")
ht.put(3, "banana")
ht.put(8, "cherry")

print("\nAfter inserting 15, 3, 8:")
ht.display()

print(f"\nget(15) = {ht.get(15)}")
print(f"get(3) = {ht.get(3)}")
print(f"get(99) = {ht.get(99)}")
```python

---

## 3. Collision Resolution

A **collision** occurs when two different keys hash to the same index. There are two main strategies to handle collisions:

### 3.1 Chaining (Closed Addressing)

Each slot contains a linked list of all elements that hash to that index.

```python
Hash function: h(key) = key % 7

Insert: 15, 22, 8, 1, 29, 36

Calculations:
  15 % 7 = 1    22 % 7 = 1    8 % 7 = 1
   1 % 7 = 1    29 % 7 = 1   36 % 7 = 1

All hash to index 1! (pathological case for demonstration)

Index:  0    1    2    3    4    5    6
      ┌────┬────┬────┬────┬────┬────┬────┐
      │    │ ●  │    │    │    │    │    │
      └────┴─│──┴────┴────┴────┴────┴────┘
             │
             ↓
           [15] → [22] → [8] → [1] → [29] → [36] → null

Insert/Search: Walk the chain at index h(key)
```python

### 3.2 Open Addressing

All elements are stored in the array itself. On collision, probe for next available slot.

#### Linear Probing

```python
h(key, i) = (h(key) + i) % m

Insert 15, 22, 29 (all hash to 1)

Step 1: h(15)=1, slot 1 empty → insert at 1
Index:  0    1    2    3    4    5    6
      ┌────┬────┬────┬────┬────┬────┬────┐
      │    │ 15 │    │    │    │    │    │
      └────┴────┴────┴────┴────┴────┴────┘

Step 2: h(22)=1, collision! probe i=1 → slot 2 empty → insert at 2
Index:  0    1    2    3    4    5    6
      ┌────┬────┬────┬────┬────┬────┬────┐
      │    │ 15 │ 22 │    │    │    │    │
      └────┴────┴────┴────┴────┴────┴────┘
            ↑────↑
            collision, move right

Step 3: h(29)=1 → try 2 (full) → try 3 → insert at 3
Index:  0    1    2    3    4    5    6
      ┌────┬────┬────┬────┬────┬────┬────┐
      │    │ 15 │ 22 │ 29 │    │    │    │
      └────┴────┴────┴────┴────┴────┴────┘
            ↑────↑────↑
            probe sequence
```python

#### Quadratic Probing

```python
h(key, i) = (h(key) + c₁*i + c₂*i²) % m

Typically: h(key, i) = (h(key) + i²) % m

Probes: +0, +1, +4, +9, +16, ...
Reduces primary clustering but may not visit all slots
```python

#### Double Hashing

```python
h(key, i) = (h₁(key) + i * h₂(key)) % m

h₁(key) = key % m
h₂(key) = 1 + (key % (m-1))  # must never return 0

Best distribution but requires two hash functions
```python

### Comparison: Chaining vs Open Addressing

| Aspect | Chaining | Open Addressing |
|--------|----------|----------------|
| **Memory** | Extra for pointers | No extra (in-place) |
| **Cache** | Poor (pointer chasing) | Better (contiguous) |
| **Load factor** | Can exceed 1.0 | Must stay < 1.0 |
| **Delete** | Easy | Needs tombstones |
| **Clustering** | No clustering | Primary/secondary clustering |

```python
# Linked List Node for Chaining
class Node:
    """
    Node for linked list used in chaining collision resolution.
    
    Attributes:
        key: The key stored in this node
        value: The value associated with the key
        next: Reference to the next node in the chain
    """
    
    def __init__(self, key, value):
        self.key = key
        self.value = value
        self.next = None


class HashTableChaining:
    """
    Hash table with separate chaining collision resolution.
    
    Each bucket contains a linked list of key-value pairs that
    hash to the same index.
    
    Attributes:
        size: Number of buckets in the table
        buckets: Array of linked list heads (one per bucket)
        count: Number of key-value pairs stored
    """
    
    def __init__(self, size: int = 7):
        """
        Initialize hash table with given number of buckets.
        
        Args:
            size: Number of buckets (preferably prime for better distribution)
        """
        self.size = size
        self.buckets = [None] * size
        self.count = 0
    
    def _hash(self, key) -> int:
        """Compute bucket index for a key."""
        return hash(key) % self.size
    
    def put(self, key, value) -> None:
        """
        Insert or update a key-value pair.
        
        Args:
            key: Key to insert
            value: Value associated with the key
        """
        index = self._hash(key)
        
        # Search for existing key in chain
        current = self.buckets[index]
        while current:
            if current.key == key:
                current.value = value  # Update existing
                return
            current = current.next
        
        # Key not found, insert at head of chain
        new_node = Node(key, value)
        new_node.next = self.buckets[index]
        self.buckets[index] = new_node
        self.count += 1
    
    def get(self, key):
        """
        Retrieve value by key.
        
        Args:
            key: Key to search for
        
        Returns:
            Value associated with key, or None if not found
        """
        index = self._hash(key)
        current = self.buckets[index]
        
        while current:
            if current.key == key:
                return current.value
            current = current.next
        
        return None
    
    def delete(self, key) -> bool:
        """
        Remove a key-value pair from the table.
        
        Args:
            key: Key to remove
        
        Returns:
            True if key was found and removed, False otherwise
        """
        index = self._hash(key)
        current = self.buckets[index]
        prev = None
        
        while current:
            if current.key == key:
                if prev:
                    prev.next = current.next
                else:
                    self.buckets[index] = current.next
                self.count -= 1
                return True
            prev = current
            current = current.next
        
        return False
    
    def load_factor(self) -> float:
        """Return current load factor (elements / buckets)."""
        return self.count / self.size
    
    def display(self):
        """Print visual representation of the hash table."""
        print(f"Hash Table (size={self.size}, count={self.count}, load={self.load_factor():.2f})")
        print("=" * 50)
        for i in range(self.size):
            chain = []
            current = self.buckets[i]
            while current:
                chain.append(f"({current.key}: {current.value})")
                current = current.next
            chain_str = " → ".join(chain) if chain else "empty"
            print(f"[{i}]: {chain_str}")
    
    def __getitem__(self, key):
        return self.get(key)
    
    def __setitem__(self, key, value):
        self.put(key, value)


# Demonstrate chaining
print("=== Hash Table with Chaining ===")
ht = HashTableChaining(7)

# Insert values that will cause collisions
data = [(15, "apple"), (22, "banana"), (8, "cherry"), 
        (1, "date"), (29, "elderberry"), (36, "fig")]

print("\nInserting keys and their hash values (mod 7):")
for key, value in data:
    print(f"  {key} % 7 = {key % 7}")
    ht.put(key, value)

print("\nResulting hash table:")
ht.display()

print(f"\nLookup examples:")
print(f"  get(22) = {ht.get(22)}")
print(f"  get(36) = {ht.get(36)}")
print(f"  get(100) = {ht.get(100)}")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
