---
name: python-bio-sets
description: "Sets and Counter for bioinformatics: gene list operations (intersection/union/difference), k-mer deduplication, and unique element tracking. Bio-specific patterns only."
tool_type: python
primary_tool: Python
---

# Sets in Bioinformatics

## Pitfalls

- **Sets are unordered:** `my_set[0]` raises `TypeError`. Convert to `sorted(my_set)` when you need ordering.
- **Set elements must be hashable:** Lists and dicts cannot be set members. Use `frozenset` for a set of sets.
- **Empty set:** `{}` creates an empty dict. Use `set()`.
- **`in` on a set is O(1):** Converting a list to a set before repeated membership tests is a significant speedup for large gene lists.

## Set Operations Reference

| Operation | Syntax | Meaning |
|---|---|---|
| Union | `A \| B` | All elements in either set |
| Intersection | `A & B` | Elements in both sets |
| Difference | `A - B` | Elements in A but not B |
| Symmetric diff | `A ^ B` | Elements in exactly one set |
| Subset test | `A <= B` | All of A is in B |

## Key Patterns

### Gene list comparisons
```python
tumor_genes   = {'TP53', 'BRCA1', 'EGFR', 'KRAS', 'MYC'}
pathway_genes = {'EGFR', 'KRAS', 'PIK3CA', 'AKT1', 'MYC'}

in_pathway   = tumor_genes & pathway_genes   # {'EGFR', 'KRAS', 'MYC'}
all_genes    = tumor_genes | pathway_genes
tumor_only   = tumor_genes - pathway_genes
```

### Unique k-mers from a sequence
```python
seq = "ATGCGATCGATCGATCGATCG"
unique_3mers = set(seq[i:i+3] for i in range(len(seq) - 2))
```

### Fast membership lookup (gene universe filter)
```python
# O(1) per lookup vs O(n) for list
known_oncogenes = set(["TP53", "KRAS", "MYC", "EGFR", ...])
hits = [g for g in candidate_list if g in known_oncogenes]
```

### Deduplicate while preserving some information
```python
# Count unique nucleotides in a sequence
unique_nt = set(sequence)   # {'A', 'T', 'G', 'C'} — or fewer if sequence has ambiguities
has_ambiguous = bool(unique_nt - {'A', 'T', 'G', 'C'})
```

### Counter for nucleotide/k-mer frequencies
```python
from collections import Counter

nt_counts = Counter(sequence)
gc_pct = (nt_counts['G'] + nt_counts['C']) / len(sequence) * 100

kmers = [sequence[i:i+3] for i in range(len(sequence) - 2)]
kmer_counts = Counter(kmers)
top5 = kmer_counts.most_common(5)

# Compare profiles between two sequences
kmers_a = Counter(seq_a[i:i+3] for i in range(len(seq_a) - 2))
kmers_b = Counter(seq_b[i:i+3] for i in range(len(seq_b) - 2))
shared   = kmers_a & kmers_b   # min counts
combined = kmers_a + kmers_b   # sum counts
```
