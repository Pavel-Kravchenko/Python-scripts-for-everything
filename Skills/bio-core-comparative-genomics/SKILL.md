---
name: bio-core-comparative-genomics
description: "- Understand the goals and rationale of comparative genomics - Implement dot plots from scratch and interpret visual patterns (inversions, duplications, translocations) - Understand synteny analysis a"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/12_Comparative_Genomics/01_comparative_genomics.ipynb"
---

# Comparative Genomics

*Source: Course notebook `Tier_2_Core_Bioinformatics/12_Comparative_Genomics/01_comparative_genomics.ipynb`*

# Comparative Genomics

---

## Learning Objectives

By the end of this notebook you will be able to:

- Understand the goals and rationale of comparative genomics
- Implement dot plots from scratch and interpret visual patterns (inversions, duplications, translocations)
- Understand synteny analysis and conserved gene order
- Distinguish orthologs from paralogs and know how to find them
- Navigate genome browsers (UCSC, Ensembl, NCBI)
- Understand whole-genome alignment concepts
- Apply these tools to real biological questions (human vs. mouse, E. coli strain comparison)

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section

## Complicated moments explained

- **Ortholog vs. paralog in comparative genomics**: Orthologs (genes that diverged via speciation) tend to retain the same function — they are useful for functional inference across species. Paralogs (genes that diverged via duplication) may gain new functions. The most reliable ortholog inference combines sequence similarity with syntenic context — genes that are similar AND in the same genomic neighborhood are more likely to be genuine orthologs.
- **Dot plot diagonals encode structural information**: A continuous diagonal from bottom-left to top-right = conserved linear segment. A diagonal running bottom-right to top-left = conserved but inverted (inversion). Repeated diagonals = tandem repeats. An off-main-diagonal block = transposition. The slope and length of diagonals encode the nature of the genomic relationship.
- **Percent identity thresholds for orthologs are heuristic**: Commonly used thresholds (>30% identity for bacterial orthologs, >70% for strain-level comparison) are empirical rules of thumb. The correct approach is to use synteny, phylogeny, and functional data together. For bacterial species, the species-level boundary is typically defined as <5% ANI (Average Nucleotide Identity) divergence.
- **Core genome vs. pan-genome**: The core genome consists of genes present in all strains of a species. The pan-genome = core + accessory (some strains) + unique (single strains). Bacterial pan-genomes are 'open' (they keep growing as more strains are sequenced); eukaryotic pan-genomes are 'closed'.
- **MUMmer vs. LASTZ vs. minimap2**: For large-genome alignment, sequence length and divergence guide tool choice. MUMmer/nucmer is fast and good for closely related genomes (>90% ANI). LASTZ handles >70% identity. minimap2 is the modern choice for any pairwise comparison, especially with long-read data.

## Environment check (run this first)

```python
# Environment check
import numpy as np
import matplotlib.pyplot as plt
import shutil

print("Imports ready.")

# Check for comparative genomics tools
tools = ['mummer', 'nucmer', 'lastz', 'minimap2', 'mafft']
for tool in tools:
    found = shutil.which(tool)
    status = f"found at {found}" if found else "NOT FOUND"
    print(f"  {tool}: {status}")
    
print("\nNote: Most demonstrations in this notebook use only Python (numpy/matplotlib)")
print("and do not require external tools.")
```

```python
def filtered_dotplot(seq1, seq2, window=5, stringency=0.6):
    """
    Dot plot with sliding window filter.

    Parameters
    ----------
    seq1, seq2 : str
        Sequences to compare.
    window : int
        Window size.
    stringency : float
        Fraction of positions in the window that must match (0.0 to 1.0).

    Returns
    -------
    2D numpy array with match scores.
    """
    threshold = int(window * stringency)
    matrix = np.zeros((len(seq2), len(seq1)))

    for i in range(len(seq2) - window + 1):
        for j in range(len(seq1) - window + 1):
            matches = sum(1 for k in range(window) if seq1[j + k] == seq2[i + k])
            if matches >= threshold:
                ci = i + window // 2
                cj = j + window // 2
                matrix[ci, cj] = matches / window

    return matrix


# Compare different window sizes
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, w in zip(axes, [1, 3, 5]):
    if w == 1:
        mat = simple_dotplot(seq1, seq2).astype(float)
    else:
        mat = filtered_dotplot(seq1, seq2, window=w, stringency=0.6)
    ax.imshow(mat, cmap='Blues', aspect='auto', interpolation='nearest')
    ax.set_title(f'Window = {w}')
    ax.set_xlabel('Sequence 1')
    ax.set_ylabel('Sequence 2')

plt.suptitle('Effect of Window Size on Dot Plot Noise', y=1.02)
plt.tight_layout()
plt.show()

print("Larger windows remove noise but may also lose short regions of similarity.")
```

### 2.3 Detecting Genomic Rearrangements

Dot plots are powerful for visualizing structural changes. Let us construct synthetic sequences with known rearrangements and see how they appear.

```python
np.random.seed(42)

def random_dna(length):
    """Generate a random DNA sequence."""
    return ''.join(np.random.choice(list('ACGT'), size=length))


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(c, c) for c in reversed(seq))


# Build a reference genome of ~600 bp with distinct segments
seg_A = random_dna(150)  # Segment A
seg_B = random_dna(150)  # Segment B
seg_C = random_dna(150)  # Segment C
seg_D = random_dna(150)  # Segment D

reference = seg_A + seg_B + seg_C + seg_D  # A-B-C-D

# Rearranged genome: A - C_inverted - B - D  (inversion of B-C region)
rearranged = seg_A + reverse_complement(seg_C) + seg_B + seg_D

print(f"Reference genome:  A -- B -- C -- D  ({len(reference)} bp)")
print(f"Rearranged genome: A -- C'-- B -- D  ({len(rearranged)} bp)")
print(f"  (C' = reverse complement of C)")
print(f"\nThis represents an inversion of segment C and a translocation of B and C.")
```

```python
def dna_dotplot(seq1, seq2, window=11, stringency=0.7, check_revcomp=True):
    """
    Create a dot plot for DNA sequences, optionally detecting both
    forward and reverse-complement matches.

    Returns (forward_matrix, revcomp_matrix) -- both are 2D arrays.
    """
    threshold = int(window * stringency)
    fwd = np.zeros((len(seq2), len(seq1)))
    rev = np.zeros((len(seq2), len(seq1)))

    seq2_rc = reverse_complement(seq2) if check_revcomp else ''

    for i in range(len(seq2) - window + 1):
        for j in range(len(seq1) - window + 1):
            # Forward match
            fwd_matches = sum(1 for k in range(window) if seq1[j+k] == seq2[i+k])
            if fwd_matches >= threshold:
                fwd[i + window//2, j + window//2] = fwd_matches / window

            # Reverse complement match
            if check_revcomp:
                rev_i = len(seq2) - 1 - i
                rev_matches = sum(1 for k in range(window) if seq1[j+k] == seq2_rc[i+k])
                if rev_matches >= threshold:
                    rev[rev_i - window//2, j + window//2] = rev_matches / window

    return fwd, rev


fwd_mat, rev_mat = dna_dotplot(reference, rearranged, window=11, stringency=0.7)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Forward matches
axes[0].imshow(fwd_mat, cmap='Blues', aspect='auto', interpolation='nearest')
axes[0].set_title('Forward matches')
axes[0].set_xlabel('Reference (A-B-C-D)')
axes[0].set_ylabel('Rearranged (A-C\'-B-D)')

# Reverse complement matches
axes[1].imshow(rev_mat, cmap='Reds', aspect='auto', interpolation='nearest')
axes[1].set_title('Reverse complement matches')
axes[1].set_xlabel('Reference (A-B-C-D)')
axes[1].set_ylabel('Rearranged (A-C\'-B-D)')

# Combined
combined = np.zeros((*fwd_mat.shape, 3))
combined[:, :, 2] = np.clip(fwd_mat, 0, 1)      # Blue = forward
combined[:, :, 0] = np.clip(rev_mat, 0, 1)       # Red = reverse complement
axes[2].imshow(combined, aspect='auto', interpolation='nearest')
axes[2].set_title('Combined (blue=forward, red=inverted)')
axes[2].set_xlabel('Reference (A-B-C-D)')
axes[2].set_ylabel('Rearranged (A-C\'-B-D)')

# Mark segment boundaries
for ax in axes:
    for pos in [150, 300, 450]:
        ax.axvline(pos, color='grey', linewidth=0.5, linestyle='--')
        ax.axhline(pos, color='grey', linewidth=0.5, linestyle='--')

plt.tight_layout()
plt.show()

print("Interpretation:")
print("  - Segment A maps to A (diagonal in top-left)")
print("  - Segment C' (inverted C) shows as an anti-diagonal (red)")
print("  - Segment B has shifted position (off-diagonal forward match)")
print("  - Segment D maps to D (diagonal in bottom-right)")
```

### 2.4 Self-Comparison Dot Plots

Comparing a sequence against itself reveals **internal repeats** and **palindromic sequences** (important for DNA regulatory elements).

```python
def self_dotplot(sequence, window=7, stringency=0.8):
    """
    Self-comparison dot plot with the main diagonal removed.
    Off-diagonal signals reveal internal repeats.
    """
    matrix = filtered_dotplot(sequence, sequence, window, stringency)

    # Remove trivial main diagonal
    for i in range(min(matrix.shape)):
        lo = max(0, i - window)
        hi = min(matrix.shape[1], i + window + 1)
        matrix[i, lo:hi] = 0

    return matrix


# Create a protein with two repeated domains
domain = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDAT"
linker = "GGSGGSGGSGGS"
protein_with_repeats = domain + linker + domain + linker + domain

matrix = self_dotplot(protein_with_repeats, window=5, stringency=0.7)

fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(matrix, cmap='Reds', aspect='equal', interpolation='nearest')
ax.set_xlabel('Sequence position')
ax.set_ylabel('Sequence position')
ax.set_title('Self-comparison Dot Plot\n(off-diagonal = internal repeats)')

# Mark domain boundaries
for pos in [len(domain), len(domain) + len(linker),
            2 * len(domain) + len(linker), 2 * len(domain) + 2 * len(linker)]:
    ax.axvline(pos, color='grey', linewidth=0.5, linestyle='--')
    ax.axhline(pos, color='grey', linewidth=0.5, linestyle='--')

plt.tight_layout()
plt.show()

print(f"Protein length: {len(protein_with_repeats)} aa")
print(f"Domain length: {len(domain)} aa")
print(f"The parallel off-diagonal lines show 3 copies of the repeated domain.")
```

---

## 3. Synteny Analysis

**Synteny** refers to the conservation of gene order and orientation between chromosomes of different species. Micro-synteny looks at a few neighboring genes; macro-synteny considers entire chromosomes.

```
Human chr 17:       --[TP53]--[WRAP53]--[EFNB3]--[DLGAP1]--[ATP1B2]--
                        |        |         |         |          |
Mouse chr 11:       --[Trp53]--[Wrap53]--[Efnb3]--[Dlgap1]--[Atp1b2]--

=> The gene order is conserved: this region is syntenic.


Human chr 2:        --[HOXD13]--[HOXD12]--[HOXD11]--[HOXD10]--
                        |          |          |          |
Mouse chr 2:        --[Hoxd13]--[Hoxd12]--[Hoxd11]--[Hoxd10]--

=> Perfect synteny in the HOX cluster -- strong evolutionary constraint.
```

Synteny breaks indicate **rearrangements** (inversions, translocations, fusions, fissions) that occurred after species diverged.

```python
# Synteny analysis with simulated gene maps

# Species A: chromosome with genes in order
species_a_genes = [
    {'name': 'GeneA', 'start': 100, 'end': 500, 'strand': '+'},
    {'name': 'GeneB', 'start': 800, 'end': 1200, 'strand': '+'},
    {'name': 'GeneC', 'start': 1500, 'end': 2000, 'strand': '+'},
    {'name': 'GeneD', 'start': 2200, 'end': 2700, 'strand': '-'},
    {'name': 'GeneE', 'start': 3000, 'end': 3500, 'strand': '+'},
    {'name': 'GeneF', 'start': 3800, 'end': 4200, 'strand': '+'},
    {'name': 'GeneG', 'start': 4500, 'end': 5000, 'strand': '-'},
    {'name': 'GeneH', 'start': 5300, 'end': 5800, 'strand': '+'},
]

# Species B: same genes but with an inversion of D-E-F region
species_b_genes = [
    {'name': 'GeneA', 'start': 200, 'end': 600, 'strand': '+'},
    {'name': 'GeneB', 'start': 900, 'end': 1300, 'strand': '+'},
    {'name': 'GeneC', 'start': 1600, 'end': 2100, 'strand': '+'},
    # Inverted region: F-E-D (reversed order, flipped strands)
    {'name': 'GeneF', 'start': 2300, 'end': 2700, 'strand': '-'},
    {'name': 'GeneE', 'start': 3000, 'end': 3500, 'strand': '-'},
    {'name': 'GeneD', 'start': 3800, 'end': 4300, 'strand': '+'},
    # Back to normal
    {'name': 'GeneG', 'start': 4600, 'end': 5100, 'strand': '-'},
    {'name': 'GeneH', 'start': 5400, 'end': 5900, 'strand': '+'},
]


def find_ortholog_pairs(genes_a, genes_b):
    """Match genes by name (simulating ortholog identification)."""
    name_to_b = {g['name']: g for g in genes_b}
    pairs = []
    for ga in genes_a:
        gb = name_to_b.get(ga['name'])
        if gb:
            pairs.append((ga, gb))
    return pairs


def detect_synteny_blocks(pairs):
    """
    Detect synteny blocks: consecutive runs of genes in the same
    relative order and orientation.
    """
    blocks = []
    current_block = [pairs[0]]

    for i in range(1, len(pairs)):
        prev_a, prev_b = pairs[i - 1]
        curr_a, curr_b = pairs[i]

        # Check if order is maintained and strand relationship is consistent
        same_order = (curr_b['start'] > prev_b['start'])
        same_strand_rel = ((curr_a['strand'] == curr_b['strand']) ==
                           (prev_a['strand'] == prev_b['strand']))

        if same_order and same_strand_rel:
            current_block.append(pairs[i])
        else:
            blocks.append(current_block)
            current_block = [pairs[i]]

    blocks.append(current_block)
    return blocks


pairs = find_ortholog_pairs(species_a_genes, species_b_genes)
blocks = detect_synteny_blocks(pairs)

print("Synteny Analysis")
print("=" * 65)
print(f"Ortholog pairs found: {len(pairs)}")
print(f"Synteny blocks: {len(blocks)}")
print()

for i, block in enumerate(blocks, 1):
    gene_names = [ga['name'] for ga, _ in block]
    gene_names_b = [gb['name'] for _, gb in block]
    # Check if block is inverted
    inverted = (block[0][0]['strand'] != block[0][1]['strand'])
    status = 'INVERTED' if inverted else 'COLLINEAR'
    print(f"Block {i} ({status}):")
    print(f"  Species A: {' -> '.join(gene_names)}")
    print(f"  Species B: {' -> '.join(gene_names_b)}")
```
