---
name: python-bio-dictionaries
description: "Dictionaries and sets for bioinformatics: codon tables, nucleotide frequency, defaultdict grouping, Counter k-mers, and set operations on gene lists."
tool_type: python
primary_tool: Python
---

# Dictionaries and Sets in Bioinformatics

## Pitfalls

- **`KeyError` vs `get()`:** Direct access `d["NNN"]` raises `KeyError` on unknown codons. Use `GENETIC_CODE.get(codon, 'X')` to handle ambiguous bases.
- **Keys must be hashable:** Lists and dicts cannot be dict keys. Tuples like `("chr17", 43044295)` work; lists do not.
- **Iterating and modifying simultaneously:** Changing dict size during iteration raises `RuntimeError`. Collect changes separately, then apply.
- **`in` checks keys, not values:** `"ATG" in codon_table` is True if `"ATG"` is a key. For values: `"Met" in codon_table.values()`.
- **Sets are unordered:** `my_set[0]` raises `TypeError`. Use `sorted(my_set)` to get ordered elements.
- **Empty set:** `{}` creates an empty dict, not a set. Use `set()`.
- **Set elements must be hashable:** Use `frozenset` when you need a set of sets.

## Key Patterns

### Genetic code translation
```python
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def translate(dna):
    protein = []
    for i in range(0, len(dna) - 2, 3):
        aa = GENETIC_CODE.get(dna[i:i+3], 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)
```

### Nucleotide frequency
```python
freq = {}
for nt in sequence:
    freq[nt] = freq.get(nt, 0) + 1
gc_pct = (freq.get('G', 0) + freq.get('C', 0)) / len(sequence) * 100
```

### defaultdict for grouping genes by chromosome
```python
from collections import defaultdict

by_chrom = defaultdict(list)
for gene, chrom in gene_locations:
    by_chrom[chrom].append(gene)
```

### Counter for k-mer profiles
```python
from collections import Counter

nt_counts = Counter(sequence)
nt_counts.most_common(2)

k = 3
kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
kmer_counts = Counter(kmers)

# Counter arithmetic — compare k-mer profiles
kmers_a = Counter(seq_a[i:i+3] for i in range(len(seq_a) - 2))
kmers_b = Counter(seq_b[i:i+3] for i in range(len(seq_b) - 2))
shared   = kmers_a & kmers_b   # min of each count
combined = kmers_a + kmers_b   # sum of each count
```

### Set operations on gene lists
```python
tumor_genes   = {'TP53', 'BRCA1', 'EGFR', 'KRAS'}
pathway_genes = {'EGFR', 'KRAS', 'PIK3CA', 'AKT1'}

shared   = tumor_genes & pathway_genes   # intersection
all_genes = tumor_genes | pathway_genes  # union
tumor_only = tumor_genes - pathway_genes # difference
unique_to_either = tumor_genes ^ pathway_genes  # symmetric difference
```

### Nested dict as annotation database
```python
gene_db = {
    "BRCA1": {
        "chromosome": "17",
        "coordinates": (43044295, 43125483),
        "strand": "-",
        "go_terms": ["DNA repair", "tumor suppression"],
    },
}
start, end = gene_db["BRCA1"]["coordinates"]
```
