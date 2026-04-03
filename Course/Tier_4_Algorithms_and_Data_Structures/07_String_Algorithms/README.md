# 07. String Algorithms

**Tier 4: Algorithms & Data Structures**

Pattern matching algorithms that power sequence search tools. The progression from O(nm) brute force to O(n+m) KMP and O(n) rolling-hash methods mirrors the engineering history of BLAST and other alignment tools. The DFA notebook shows how regular-expression engines are compiled -- the same mechanism behind motif scanning.

## Topics Covered

- Brute-force (naive) pattern matching and its O(nm) worst case
- Sliding window technique
- Knuth-Morris-Pratt (KMP): failure function construction, O(n+m) search
- Rabin-Karp: rolling hash, fingerprinting, multiple-pattern search
- Deterministic finite automaton (DFA) construction for a pattern
- DFA matching: O(n) search after O(m|Σ|) preprocessing
- Relationship between these algorithms and the seed-and-extend strategy in BLAST

## Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [01_naive_pattern_matching.ipynb](01_naive_pattern_matching.ipynb) | Brute-force matching, sliding window |
| 2 | [02_kmp_algorithm.ipynb](02_kmp_algorithm.ipynb) | Failure function construction, KMP search |
| 3 | [03_rabin_karp.ipynb](03_rabin_karp.ipynb) | Rolling hash, fingerprinting, multiple patterns |
| 4 | [04_dfa_matching.ipynb](04_dfa_matching.ipynb) | DFA construction, O(n) matching |

## Complexity Reference

| Algorithm | Preprocessing | Search | Space |
|-----------|---------------|--------|-------|
| Naive | O(1) | O(nm) | O(1) |
| KMP | O(m) | O(n) | O(m) |
| Rabin-Karp | O(m) | O(n) average, O(nm) worst | O(1) |
| DFA | O(m\|Σ\|) | O(n) | O(m\|Σ\|) |

n = text length, m = pattern length, |Σ| = alphabet size

## Bioinformatics Connections

- BLAST seed-and-extend heuristic uses a hash-based seed lookup (Rabin-Karp ideas) followed by ungapped extension (similar to naive matching on a short window). See [BLAST Searching](../../Tier_2_Core_Bioinformatics/04_BLAST_Searching/01_blast_searching.ipynb).
- Motif scanning across a genome (e.g., finding all TATA box occurrences) is a direct application of the DFA and KMP algorithms from this module.

## Prerequisites

- [01. Complexity Analysis](../01_Complexity_Analysis/README.md)
- [06. Hash-Based Structures](../06_Hash_Based_Structures/README.md) -- rolling hash builds on hash function concepts

---

[<< 06. Hash-Based Structures](../06_Hash_Based_Structures/README.md) | [Tier 4 Overview](../README.md) | [08. Advanced String Structures >>](../08_Advanced_String_Structures/README.md)
