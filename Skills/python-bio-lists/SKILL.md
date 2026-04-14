---
name: python-bio-lists
description: List and tuple patterns for bioinformatics — codon splitting, gene coordinate records, named tuples, sorting by GC/length, and copy pitfalls
tool_type: python
primary_tool: Python
---

# Lists and Tuples for Bioinformatics

## Pitfalls

- **Assignment does not copy:** `alias = original` makes both point to the same list. Use `original.copy()` or `original[:]` for an independent copy.
- **`append()` vs `extend()`:** `genes.append(["BRCA1", "TP53"])` adds a list as one element; `genes.extend(["BRCA1", "TP53"])` adds two elements.
- **`sort()` returns `None`:** it modifies in place. `sorted()` returns a new list.
- **Tuples are immutable and hashable:** use for fixed records (coordinates, SNPs) and as dict keys.

## Bio Recipes

### Split CDS into Codons

```python
codons = [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]
```

### Sort Sequences by GC Content

```python
def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq)

by_gc = sorted(sequences, key=gc_content)
```

### Gene Coordinate Tuples

```python
genes = [
    ("BRCA1", 43044295, 43125483),
    ("TP53",  7661779,  7687538),
    ("EGFR",  55019017, 55207337),
]
longest = max(genes, key=lambda g: g[2] - g[1])
```

### Named Tuples for Readable Records

```python
from collections import namedtuple

Gene = namedtuple('Gene', ['name', 'chromosome', 'start', 'end', 'strand'])
brca1 = Gene('BRCA1', 'chr17', 43044295, 43125483, '-')
length = brca1.end - brca1.start   # access by name, not index
```

### When to Use Tuples vs Lists

| Tuples | Lists |
|--------|-------|
| Fixed records (gene coords, SNPs) | Collections that grow/shrink |
| Dict keys (must be hashable) | Data you will sort or reorder |
| Multiple return values | Intermediate computation results |
