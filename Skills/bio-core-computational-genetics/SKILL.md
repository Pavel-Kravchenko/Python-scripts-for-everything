---
name: bio-core-computational-genetics
description: - Build the standard genetic code programmatically and translate DNA to protein - Understand codon degeneracy and compute codon usage statistics (RSCU, CAI) - Perform virtual restriction enzyme digest
tool_type: python
primary_tool: NumPy
---

# Computational Genetics

## Complicated moments explained

- **Degeneracy vs. synonymous**: All synonymous codons encode the same amino acid, but not all positions in a codon are equally degenerate. Four-fold degenerate sites (e.g., the third position of Pro: CCA, CCG, CCT, CCC all encode Pro) are effectively neutral markers. Two-fold degenerate sites (e.g., first two Leu codons: CTN = four Leu codons, but TTA/TTG are two more) have constrained positions.
- **CAI reference set matters**: CAI measures codon usage adaptation relative to a reference set of highly expressed genes. Using the wrong reference set (wrong organism, or genes that aren't actually highly expressed) gives misleading CAI values. For a new organism, use ribosomal protein genes or metabolic genes as the reference.
- **Restriction enzyme naming conventions**: Type II enzymes are named for the organism (first three letters = genus + species initial + strain): EcoRI = *Escherichia coli* RY13, HindIII = *Haemophilus influenzae* Rd. The Roman numeral distinguishes enzymes from the same organism. Type IIS enzymes (like BsaI) cut outside their recognition sequence — essential for Golden Gate assembly.
- **Three-point cross logic**: The double-crossover class has the smallest frequency. The gene in the middle is identified by comparing parental and double-crossover classes — the middle gene changes relative to both flanking markers in double crossovers. Map distance = (single COs + 2×double COs) / total × 100 cM.
- **HWE violations are informative**: Significant departure from Hardy-Weinberg expectations signals population structure, selection, inbreeding, genotyping error, or copy number variation. In GWAS, SNPs that fail HWE in controls (p < 1e-6) are often flagged as potential genotyping errors.

## Environment check (run this first)

```python
# Environment check
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import random

print("Imports ready.")

# Quick genetic code demo
bases = 'TCAG'
amino_acids = (
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
)
codon_table = {}
for i, first in enumerate(bases):
    for j, second in enumerate(bases):
        for k, third in enumerate(bases):
            codon = first + second + third
            aa = amino_acids[i*16 + j*4 + k]
            codon_table[codon] = aa

stop_codons = [c for c, a in codon_table.items() if a == '*']
start_codons = ['ATG']
print(f"Standard genetic code: {len(codon_table)} codons -> 20 amino acids + 3 stop codons ({stop_codons})")
print(f"Start codon(s): {start_codons} (bacteria also use GTG, TTG)")

# Degeneracy
from collections import Counter
aa_counts = Counter(codon_table.values())
aa_counts.pop('*')
by_fold = defaultdict(list)
for aa, count in aa_counts.items():
    by_fold[count].append(aa)
for fold in sorted(by_fold):
    print(f"  {fold}-fold degenerate: {', '.join(sorted(by_fold[fold]))}")
```

### Codon Degeneracy

Not all codon positions contribute equally to amino acid identity:

| Fold | Description | Example |
|------|-------------|--------|
| **1-fold** | Every change at that position changes the amino acid | AUG (Met) — only codon |
| **2-fold** | Two alternatives at the 3rd position encode the same aa | AAA/AAG → Lys |
| **4-fold** | All four bases at the 3rd position encode the same aa | GCN → Ala |

The third codon position is most often degenerate — this is why synonymous substitutions
accumulate faster than non-synonymous ones.

```python
from collections import defaultdict

# Group codons by amino acid
aa_to_codons = defaultdict(list)
for codon, aa in codon_table.items():
    aa_to_codons[aa].append(codon)

# Classify degeneracy of the 3rd position for each synonymous family
degeneracy_counts = defaultdict(int)  # fold -> number of amino acid families
for aa, codons in aa_to_codons.items():
    if aa == '*':
        continue
    fold = len(codons)
    degeneracy_counts[fold] += 1

print('Codon family sizes (fold-degeneracy) for the 20 amino acids:')
for fold in sorted(degeneracy_counts):
    print(f'  {fold}-fold: {degeneracy_counts[fold]} amino acids')

print()
print('4-fold degenerate amino acids (any nucleotide at position 3):')
for aa, codons in sorted(aa_to_codons.items()):
    if len(codons) == 4 and aa != '*':
        print(f'  {aa}: {codons}')
```

### Start and Stop Codons

- **Standard start codon**: `ATG` (Met) — universal across domains of life
- **Alternative starts**: `TTG`, `GTG`, `CTG` used in bacteria (still incorporate Met)
- **Stop codons**: `TAA` (ochre), `TAG` (amber), `TGA` (opal/umber)

In some organisms stop codons are reassigned (e.g., `TGA` → Sec in selenoproteins,
`TAG` → Pyl in methanogenic archaea).

```python
# Demonstrate alternative start codons in bacteria
alt_starts = ['ATG', 'TTG', 'GTG', 'CTG']
stop_codons = [c for c, a in codon_table.items() if a == '*']

print('Start codons and their standard amino acid translation:')
for codon in alt_starts:
    aa = codon_table[codon]
    note = ' <- canonical' if codon == 'ATG' else ' <- alt start (bacteria)'
    print(f'  {codon} -> {aa}{note}')

print()
print('Stop codons (traditional names):')
stop_names = {'TAA': 'ochre', 'TAG': 'amber', 'TGA': 'opal'}
for codon in stop_codons:
    print(f'  {codon} ({stop_names[codon]})')
```

### Exercise 1.1: Implement `translate()` with Frame Support

Extend the `translate()` function to accept a `frame` parameter (0, 1, or 2)
that shifts the reading frame before searching for the start codon.
Then translate the sequence below in all three forward frames.


```python
# Solution
def translate_frame(dna_sequence, codon_table, frame=0):
    """Translate from a given reading frame (0, 1, or 2)."""
    return translate(dna_sequence[frame:], codon_table)


test_seq = 'TAATGCCCGAATTTGCCTAAATGGGCAAATAG'
print(f'Sequence: {test_seq} (len={len(test_seq)})')
print()
for frame in range(3):
    protein = translate_frame(test_seq, codon_table, frame)
    print(f'Frame {frame}: {protein if protein else "(no ORF found)"}')
```


## Codon Usage Bias

Even though synonymous codons encode the same amino acid, organisms do not use them equally.
This **codon usage bias** arises from:

- **tRNA abundance**: abundant tRNAs match preferred codons → faster elongation
- **GC content**: genomic GC bias pushes toward GC-rich or AT-rich codons
- **Translational selection**: highly expressed genes show stronger bias
- **Mutational pressure**: background mutation rates shape the baseline

Practical consequence: heterologous protein expression benefits from **codon optimisation**
to match the host organism's preferences.

```python
from collections import Counter

def count_codons(coding_sequence):
    """Count triplet codons in a coding sequence (must start at codon position 0)."""
    seq = coding_sequence.upper()
    counts = Counter()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3 and 'N' not in codon:
            counts[codon] += 1
    return counts


# Synthetic E. coli-like coding sequence (200 codons)
import random
random.seed(42)

# E. coli strongly prefers certain Leu codons: CTG  others
ecoli_codon_weights = {
    # Leu
    'CTG': 0.50, 'CTT': 0.10, 'CTC': 0.10, 'CTA': 0.04, 'TTA': 0.13, 'TTG': 0.13,
    # Ala
    'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33,
    # Val
    'GTG': 0.37, 'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.15,
    # Gly
    'GGC': 0.34, 'GGT': 0.35, 'GGG': 0.15, 'GGA': 0.11, 'GGN': 0.05,
}

# Build a realistic synthetic gene for E. coli
ecoli_synonymous_families = {
    'Leu': ['CTG', 'CTT', 'CTC', 'CTA', 'TTA', 'TTG'],
    'Ala': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Arg': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'Gly': ['GGT', 'GGC', 'GGA', 'GGG'],
    'Pro': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Thr': ['ACT', 'ACC', 'ACA', 'ACG'],
    'Val': ['GTT', 'GTC', 'GTA', 'GTG'],
    'Ser': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
}

ecoli_aa_weights = {
    'Leu': [0.50, 0.10, 0.10, 0.04, 0.13, 0.13],
    'Ala': [0.18, 0.26, 0.23, 0.33],
    'Arg': [0.36, 0.36, 0.07, 0.11, 0.07, 0.04],
    'Gly': [0.35, 0.34, 0.11, 0.15, 0.05],
    'Pro': [0.16, 0.10, 0.20, 0.55],
    'Thr': [0.19, 0.40, 0.17, 0.25],
    'Val': [0.28, 0.20, 0.15, 0.37],
    'Ser': [0.15, 0.15, 0.14, 0.15, 0.15, 0.26],
}

def build_synthetic_gene(aa_families, aa_weights, n_codons=200):
    aas = list(aa_families.keys())
    gene_codons = []
    for _ in range(n_codons):
        aa = random.choice(aas)
        codons = aa_families[aa]
        weights = aa_weights[aa][:len(codons)]
        chosen = random.choices(codons, weights=weights)[0]
        gene_codons.append(chosen)
    return 'ATG' + ''.join(gene_codons) + 'TAA'

ecoli_gene = build_synthetic_gene(ecoli_synonymous_families, ecoli_aa_weights)
ecoli_counts = count_codons(ecoli_gene)

print(f'E. coli synthetic gene length: {len(ecoli_gene)} nt ({len(ecoli_gene)//3} codons)')
print('Top 10 most frequent codons:')
for codon, count in ecoli_counts.most_common(10):
    aa = codon_table.get(codon, '?')
    print(f'  {codon} ({aa}): {count}')
```

### Relative Synonymous Codon Usage (RSCU)

RSCU normalises codon counts by the expected frequency if all synonymous codons
were used equally:

$$\text{RSCU}_{ij} = \frac{X_{ij}}{\bar{X}_i} = \frac{X_{ij}}{\frac{1}{n_i}\sum_{j=1}^{n_i} X_{ij}}$$

where $X_{ij}$ is the count of the $j$-th codon for amino acid $i$, and $n_i$ is the
number of synonymous codons for that amino acid.

- **RSCU = 1.0**: codon used as expected under equal usage
- **RSCU > 1.0**: over-represented (preferred)
- **RSCU < 1.0**: under-represented (disfavoured)

```python
def compute_rscu(codon_counts, codon_table):
    """Compute RSCU for each codon given observed counts."""
    # Group codons by amino acid
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
            observed = codon_counts.get(codon, 0)
            rscu[codon] = (observed / expected) if expected > 0 else 0.0
    return rscu


ecoli_rscu = compute_rscu(ecoli_counts, codon_table)

# Show RSCU for Leucine codons
leu_codons = ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']
print('RSCU for Leucine (Leu) codons in E. coli synthetic gene:')
print(f'{"Codon":8} {"Count":8} {"RSCU":8}')
print('-' * 26)
for codon in leu_codons:
    count = ecoli_counts.get(codon, 0)
    rscu_val = ecoli_rscu.get(codon, 0)
    bar = '*' * int(rscu_val * 5)
    print(f'{codon:8} {count:8} {rscu_val:8.2f}  {bar}')
```

### Codon Adaptation Index (CAI)

The **CAI** (Sharp & Li, 1987) measures how well a gene's codon usage matches the organism's
highly expressed reference gene set. It uses the **Relative Adaptiveness** ($w_{ij}$):

$$w_{ij} = \frac{\text{RSCU}_{ij}}{\text{RSCU}_{i,\max}}$$

The CAI for a gene of $L$ codons is the geometric mean of the $w$ values:

$$\text{CAI} = \left(\prod_{k=1}^{L} w_k\right)^{1/L}$$

CAI ranges from 0 (poorly adapted) to 1.0 (perfectly adapted to the reference).

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
