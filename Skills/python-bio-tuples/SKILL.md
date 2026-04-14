---
name: python-bio-tuples
description: "Lists and tuples for bioinformatics: codon extraction, genomic coordinate records, named tuples, and copy semantics. Bio-specific patterns only."
tool_type: python
primary_tool: Python
---

# Lists and Tuples in Bioinformatics

## Pitfalls

- **Assignment does not copy:** `alias = original` — both names point to the same list. Use `original[:]` or `original.copy()` for an independent copy.
- **`append()` vs `extend()`:** `genes.append(["BRCA1", "TP53"])` adds a nested list (one element); `genes.extend(["BRCA1", "TP53"])` adds two elements.
- **`sort()` vs `sorted()`:** `my_list.sort()` modifies in place, returns `None`. `sorted(my_list)` leaves original unchanged.
- **Tuples are immutable:** Cannot append, remove, or change elements. Use tuples for fixed records (coordinates, SNPs) and lists for mutable collections.
- **Mutable default arguments:** Never `def f(x=[])` — use `def f(x=None)` and assign inside.
- **Off-by-one:** Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based.
- **Deep vs shallow copy:** `list.copy()` only copies the top level; nested structures require `copy.deepcopy()`.

## Key Patterns

### Codon extraction from CDS
```python
cds = "ATGGCCGATCGATAG"
codons = [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]
```

### Sort by biological property
```python
sequences = ["ATATATAT", "GCGCGCGC", "ATGCATGC"]
gc = lambda s: (s.count('G') + s.count('C')) / len(s)
by_gc = sorted(sequences, key=gc)

gene_data = [("BRCA1", 81189), ("TP53", 19149), ("EGFR", 188307)]
longest = max(gene_data, key=lambda g: g[1])
```

### Tuple unpacking for genomic records
```python
genes = [
    ("BRCA1", 43044295, 43125483),
    ("TP53",  7661779,  7687538),
]
for name, start, end in genes:
    print(f"{name}: {end - start:,} bp")

# Star unpacking
header = ("BRCA1", "breast cancer type 1", "chr17", 43044295, 43125483, "-")
gene_id, description, *location_info = header
```

### Named tuples for readable records
```python
from collections import namedtuple

Gene = namedtuple('Gene', ['name', 'chromosome', 'start', 'end', 'strand'])
brca1 = Gene('BRCA1', 'chr17', 43044295, 43125483, '-')

print(f"{brca1.name} on {brca1.chromosome}:{brca1.start}-{brca1.end}")
by_length = sorted([brca1, tp53], key=lambda g: g.end - g.start, reverse=True)
```

### Safe copy pattern
```python
original = ["ATG", "GCC", "GAT"]
copy = original[:]      # or original.copy()
copy.append("TAG")      # original unchanged
```

## Decision Table: Tuple vs List

| Use case | Type |
|---|---|
| Fixed genomic record (name, chr, start, end) | Tuple |
| Dictionary key (coordinates) | Tuple (hashable) |
| Return multiple values from function | Tuple |
| Collection that grows during pipeline | List |
| Sequences to sort or reorder | List |
| Intermediate computation results | List |
