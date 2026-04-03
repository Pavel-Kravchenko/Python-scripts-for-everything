# 06. Hash-Based Structures

**Tier 4: Algorithms & Data Structures**

Hash tables deliver O(1) average-case lookup by trading space for speed -- a bargain that underlies k-mer counting, read deduplication, and every Python `dict`. This module also covers Bloom filters, a probabilistic structure that answers membership queries in constant space, used in contamination screening pipelines.

## Topics Covered

- Hash functions: properties, polynomial rolling hash, avalanche effect
- Collision resolution: separate chaining, open addressing (linear probing, quadratic probing, double hashing)
- Load factor and rehashing thresholds
- Bloom filters: bit array, multiple hash functions, false-positive rate
- Counting Bloom filters and their limitations
- Applications: set membership, deduplication, frequency estimation

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_hash_tables_bloom.ipynb](01_hash_tables_bloom.ipynb) | Hash functions, collision resolution strategies, load factors, Bloom filters (20 cells) |

## Bioinformatics Connections

- k-mer counting with Jellyfish uses a hash table internally; understanding load factor and hash quality explains its memory usage. See [NGS Fundamentals](../../Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/01_ngs_fundamentals.ipynb).
- Read deduplication (e.g., Picard MarkDuplicates) hashes read start positions to find identical reads in O(n) time.
- Contamination screening tools use Bloom filters to quickly reject reads that cannot belong to the target organism.

## Prerequisites

- [01. Complexity Analysis](../01_Complexity_Analysis/README.md) -- amortized analysis of rehashing

---

[<< 05. Tree Structures](../05_Tree_Structures/README.md) | [Tier 4 Overview](../README.md) | [07. String Algorithms >>](../07_String_Algorithms/README.md)
