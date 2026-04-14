---
name: python-bio-comprehensions
description: Python comprehensions and generators for bioinformatics — list/dict/set comprehensions, generator expressions for large sequence files, bio-specific patterns.
tool_type: python
primary_tool: Python
---

# Python Comprehensions for Bioinformatics

## Pitfalls

- **Nested comprehension loop order**: In a flat nested comprehension, outer loop comes first — `[expr for outer in outer_list for inner in inner_list]`. Use `[[expr for inner in ...] for outer in ...]` for list-of-lists.
- **`if` filter vs `if`/`else` transform**: `if` after `for` is a filter (removes items). `if`/`else` in the expression is a conditional transform (keeps all items, maps differently).
- **Generator expressions are single-pass**: Once consumed (`sum(...)`, `list(...)`, etc.), the generator is exhausted. Use a list comprehension when you need multiple passes.
- **Dict comprehension with `.count()` is O(n²)**: `{kmer: seq.count(kmer) for kmer in set(kmers)}` recomputes count for every unique key. Use `collections.Counter` instead.
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based.
- **Deep vs shallow copy**: `list.copy()` and `[:]` only copy the top level. Use `copy.deepcopy()` for nested structures.

## Quick Reference

```python
# Filter: keeps only G and C
[nt for nt in dna if nt in 'GC']

# Transform: relabels every base
["purine" if nt in "AG" else "pyrimidine" for nt in dna]

# Filter AND transform: keep long seqs, return (name, gc)
[(name, gc(seq)) for name, seq in fasta.items() if len(seq) > 100]

# Dict comprehension
{name: gc(seq) for name, seq in fasta.items()}

# Set comprehension
{aa for aa in protein}

# Generator (lazy, memory-efficient)
avg_gc = sum(gc(s) for s in sequences) / len(sequences)
```

## Bio-Specific Patterns

```python
# GC content
def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq)

# Codon splitting
codons = [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]

# Translation (up to first stop)
protein = ''.join(
    aa for aa in (GENETIC_CODE.get(c, 'X') for c in codons)
    if aa != '*'
)

# Reverse complement
rc_map = {'A':'T','T':'A','G':'C','C':'G'}
rev_comp = ''.join(rc_map[n] for n in seq[::-1])

# Sliding window GC profile
window = 100
gc_profile = [(i, (seq[i:i+window].count('G') + seq[i:i+window].count('C')) / window)
              for i in range(len(seq) - window + 1)]

# Three reading frames (list of lists)
frames = [[dna[i:i+3] for i in range(f, len(dna)-2, 3)] for f in range(3)]

# Flatten all codons across frames
all_codons = [codon for frame in frames for codon in frame]

# All k-mers of length k
from itertools import product
all_kmers = [''.join(c) for c in product('ATGC', repeat=k)]

# Unique k-mer spectrum
spectrum = sorted({seq[i:i+k] for i in range(len(seq) - k + 1)})

# Codons encoding a given amino acid
leu_codons = {codon for codon, aa in GENETIC_CODE.items() if aa == 'L'}

# Filter FASTA by GC and length
filtered = {name: seq for name, seq in fasta.items()
            if len(seq) > 100 and 0.4 <= gc_content(seq) <= 0.6}

# Sequences with GC > 55% — generator, no intermediate list
high_gc_count = sum(1 for seq in sequences if gc_content(seq) > 0.55)

# Gene-to-pathway reverse map
pathways = {'DNA_repair': {'BRCA1','BRCA2','ATM'}, 'Cell_cycle': {'TP53','RB1','CDK4'}}
gene_to_pathways = {gene: [p for p, genes in pathways.items() if gene in genes]
                    for gene in {g for genes in pathways.values() for g in genes}}

# any/all with generators
any_pure_gc  = any(set(s) <= {'G', 'C'} for s in sequences)
all_valid    = all(set(s) <= {'A', 'T', 'G', 'C'} for s in sequences)
```

## Memory: Generator vs List

```python
import sys
n = 100_000
list_comp = [x**2 for x in range(n)]   # ~800 KB in memory
gen_expr  = (x**2 for x in range(n))   # ~200 bytes
# Use generator when: processing genome-scale data, streaming files, passing to sum/max/any/all
# Use list when: you need multiple passes, indexing, or len()
```
