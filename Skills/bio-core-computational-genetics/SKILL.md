---
name: bio-core-computational-genetics
description: - Build the standard genetic code programmatically and translate DNA to protein - Understand codon degeneracy and compute codon usage statistics (RSCU, CAI) - Perform virtual restriction enzyme digest
tool_type: python
primary_tool: NumPy
---

# Computational Genetics

## Key Concepts

- **Degeneracy vs. synonymous**: Four-fold degenerate sites (e.g., Pro: CCN) are effectively neutral; two-fold degenerate sites have one constrained position. All synonymous codons encode the same AA but not all positions are equally degenerate.
- **CAI reference set**: CAI measures codon adaptation relative to a reference set of highly expressed genes. Wrong organism or non-HEG reference gives misleading values; use ribosomal protein or metabolic genes.
- **Restriction enzyme naming**: Type II named from organism (EcoRI = *E. coli* RY13). Roman numeral distinguishes isolates. Type IIS (BsaI) cuts outside recognition sequence — essential for Golden Gate.
- **Three-point cross**: Double-crossover class is smallest. Middle gene identified by comparing parental vs. double-crossover classes. Map distance = (single COs + 2×double COs) / total × 100 cM.
- **HWE violations**: Departure signals population structure, selection, inbreeding, genotyping error, or CNV. SNPs failing HWE in controls (p < 1e-6) in GWAS are likely genotyping errors.

## Environment check (run this first)

```python
from collections import defaultdict, Counter
import numpy as np
import matplotlib.pyplot as plt
import random

bases = 'TCAG'
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = {}
for i, first in enumerate(bases):
    for j, second in enumerate(bases):
        for k, third in enumerate(bases):
            codon = first + second + third
            codon_table[codon] = amino_acids[i*16 + j*4 + k]

stop_codons = [c for c, a in codon_table.items() if a == '*']
print(f"Standard genetic code: {len(codon_table)} codons, stop: {stop_codons}")
print("Bacteria alt starts: GTG, TTG (still incorporate Met)")
```

## Codon Degeneracy

| Fold | Description | Example |
|------|-------------|--------|
| **1-fold** | Every change changes the AA | AUG (Met) |
| **2-fold** | Two alternatives at 3rd position | AAA/AAG → Lys |
| **4-fold** | All four bases at 3rd position same AA | GCN → Ala |

Third position is most often degenerate — synonymous substitutions accumulate faster than non-synonymous.

```python
aa_to_codons = defaultdict(list)
for codon, aa in codon_table.items():
    aa_to_codons[aa].append(codon)

degeneracy_counts = defaultdict(int)
for aa, codons in aa_to_codons.items():
    if aa != '*':
        degeneracy_counts[len(codons)] += 1

for fold in sorted(degeneracy_counts):
    print(f'  {fold}-fold: {degeneracy_counts[fold]} amino acids')
```

Stop codons: `TAA` (ochre), `TAG` (amber), `TGA` (opal). In some organisms reassigned: `TGA` → Sec (selenoproteins), `TAG` → Pyl (methanogenic archaea).

### Exercise 1.1: translate() with frame support

```python
def translate_frame(dna_sequence, codon_table, frame=0):
    return translate(dna_sequence[frame:], codon_table)

test_seq = 'TAATGCCCGAATTTGCCTAAATGGGCAAATAG'
for frame in range(3):
    protein = translate_frame(test_seq, codon_table, frame)
    print(f'Frame {frame}: {protein if protein else "(no ORF found)"}')
```

## Codon Usage Bias

Arises from: tRNA abundance (preferred codons → faster elongation), GC content, translational selection (HEGs show stronger bias), mutational pressure. Practical: heterologous expression benefits from codon optimisation to match host.

```python
def count_codons(coding_sequence):
    seq = coding_sequence.upper()
    counts = Counter()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3 and 'N' not in codon:
            counts[codon] += 1
    return counts
```

### RSCU (Relative Synonymous Codon Usage)

$$\text{RSCU}_{ij} = \frac{X_{ij}}{\frac{1}{n_i}\sum_{j=1}^{n_i} X_{ij}}$$

- RSCU = 1.0: equal usage; > 1.0: over-represented; < 1.0: disfavoured

```python
def compute_rscu(codon_counts, codon_table):
    aa_groups = defaultdict(list)
    for codon, aa in codon_table.items():
        if aa != '*':
            aa_groups[aa].append(codon)
    rscu = {}
    for aa, synonyms in aa_groups.items():
        total = sum(codon_counts.get(c, 0) for c in synonyms)
        n = len(synonyms)
        expected = total / n if n > 0 else 0
        for codon in synonyms:
            rscu[codon] = (codon_counts.get(codon, 0) / expected) if expected > 0 else 0.0
    return rscu
```

### CAI (Codon Adaptation Index)

Relative adaptiveness: $w_{ij} = \text{RSCU}_{ij} / \text{RSCU}_{i,\max}$

Geometric mean over all codons in gene:

$$\text{CAI} = \left(\prod_{k=1}^{L} w_k\right)^{1/L}$$

CAI 0 (poorly adapted) → 1.0 (perfectly adapted to reference).

## Pitfalls

- **Coordinate systems**: BED 0-based half-open; VCF/GFF 1-based inclusive — mixing causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) for thousands of features
