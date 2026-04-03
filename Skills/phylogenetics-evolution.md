---
name: phylogenetics-evolution
description: Phylogenetic tree construction (NJ, UPGMA, ML), bootstrap analysis, comparative genomics, synteny, and ortholog detection
---

# Phylogenetics & Comparative Genomics

## When to Use
- Inferring evolutionary relationships from aligned sequences
- Comparing gene families across species (orthologs/paralogs)
- Detecting genomic rearrangements, conserved synteny, or duplications
- Assessing statistical support for tree clades (bootstrap)
- Whole-genome alignment and pan-genome analysis

## Quick Reference

### Tree Anatomy
- **Leaf / terminal node**: observed taxon (sequence or species)
- **Internal node**: inferred ancestor
- **Branch length**: evolutionary distance (substitutions/site) or time
- **Clade**: a node and all its descendants
- **Root**: common ancestor of all taxa; absent in unrooted trees
- **Bipartition / split**: an internal branch divides leaves into two sets
- **MRCA**: Most Recent Common Ancestor

### Method Comparison
| Method | Type | Molecular clock? | Complexity | Best For |
|--------|------|-----------------|------------|----------|
| UPGMA | Distance | Required | O(n³) | Guide trees for MSA |
| Neighbor-Joining | Distance | No | O(n³) | Quick exploratory analysis |
| Max Parsimony | Character | No | NP-hard | Few taxa, morphological data |
| Max Likelihood | Statistical | Optional | Heuristic | Rigorous hypothesis testing |
| Bayesian (MCMC) | Statistical | Optional | Very slow | Divergence time estimation |

### Bootstrap Thresholds
- ≥95%: strong support
- 70–94%: moderate support
- <70%: weak, treat clade with caution

### Distance Models
- **p-distance**: raw fraction of differing sites; underestimates at high divergence
- **JC69**: corrects for back-mutations, assumes equal base frequencies; `d = -0.75 * ln(1 - 4p/3)`
- **Kimura 2P**: separate rates for transitions (Ts) and transversions (Tv)
- **Rule**: use a correction model whenever p-distance > 0.1

## Key Patterns

### UPGMA vs Neighbor-Joining
- UPGMA assumes a **molecular clock** → produces ultrametric trees (all tips equidistant from root); fails when rates vary across lineages
- NJ uses the Q-matrix to correct for rate variation → preferred for most real data
- Both are O(n³) and fast enough for hundreds of taxa

### Outgroup Rooting
- NJ trees are unrooted by default
- Root by placing the outgroup (a known distant taxon) as sister to all ingroup taxa
- Without an outgroup, mid-point rooting (longest branch) is a fallback

### Robinson-Foulds Distance
- Counts splits present in one tree but absent in the other
- RF = 0: identical topology; RF = 2(n-3): maximally different (n taxa)
- Use to compare recovered tree vs. true tree in simulation studies

### Dot Plot Patterns
| Visual pattern | Biological meaning |
|---------------|-------------------|
| Continuous diagonal | Collinear homology |
| Parallel diagonals | Tandem/dispersed duplications |
| Anti-diagonal line | Inversion |
| Shifted diagonal | Translocation |
| Off-diagonal cluster (self-dot) | Internal repeat |

### Ortholog vs Paralog
- **Ortholog**: same gene in different species separated by speciation → usually same function
- **Paralog**: copies in the same or different species separated by gene duplication → often diverge in function
- **In-paralog**: duplication after speciation
- **Out-paralog**: duplication before speciation (both lineages carry both copies)
- Annotation transfer between species is reliable only for 1:1 orthologs

### Pan-genome (Bacteria)
- **Core genome**: genes in all strains; housekeeping functions
- **Accessory genome**: genes in some strains; often on genomic islands, plasmids
- **Unique genes**: strain-specific; phage insertions, novel virulence factors

## Code Templates

### BioPython: distance matrix → UPGMA/NJ tree
```python
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import io

# From a MultipleSeqAlignment object
calculator = DistanceCalculator('identity')   # or 'blosum62' for proteins
dm = calculator.get_distance(alignment)

constructor = DistanceTreeConstructor()
upgma_tree = constructor.upgma(dm)
nj_tree    = constructor.nj(dm)

Phylo.draw_ascii(nj_tree)
```

### BioPython: read/write/inspect trees
```python
# Parse Newick
tree = Phylo.read(io.StringIO("((Human:0.1,Chimp:0.12):0.2,(Mouse:0.25,Rat:0.23):0.15);"), 'newick')

# Write
buf = io.StringIO()
Phylo.write(tree, buf, 'newick')

# Inspect
tree.count_terminals()          # number of leaves
tree.total_branch_length()
tree.is_bifurcating()
[t.name for t in tree.get_terminals()]
[n.branch_length for n in tree.get_nonterminals()]

# Distances and MRCA
tree.distance('Human', 'Chimp')
mrca = tree.common_ancestor('Human', 'Chimp')
[t.name for t in mrca.get_terminals()]

# Rooting
tree.root_with_outgroup('Mouse')   # in-place

# Draw with matplotlib
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 6))
Phylo.draw(tree, axes=ax, do_show=False)
```

### JC69 distance
```python
import numpy as np

def jukes_cantor_dist(seq1, seq2):
    diffs = sum(a != b for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    compared = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    p = diffs / compared
    if p >= 0.75:
        return float('inf')   # saturation
    return -0.75 * np.log(1 - 4 * p / 3)
```

### Bootstrap analysis
```python
import random

def bootstrap_replicate(alignment):
    """Sample alignment columns with replacement."""
    L = len(alignment[0])
    cols = [random.randint(0, L - 1) for _ in range(L)]
    return [''.join(seq[c] for c in cols) for seq in alignment]

# Pattern: run N replicates, count how often each clade appears
clade_counts = {}
for _ in range(100):
    boot_aln = bootstrap_replicate(seqs)
    boot_tree = build_tree(boot_aln, names)   # your NJ/UPGMA function
    for clade in get_clades(boot_tree):
        clade_counts[clade] = clade_counts.get(clade, 0) + 1
# support = clade_counts[clade] / n_replicates * 100
```

### Robinson-Foulds distance
```python
def get_splits(tree):
    """Return set of frozensets, each a bipartition of leaf names."""
    all_leaves = frozenset(t.name for t in tree.get_terminals())
    splits = set()
    for clade in tree.get_nonterminals():
        below = frozenset(t.name for t in clade.get_terminals())
        above = all_leaves - below
        if len(below) > 1 and len(above) > 1:
            splits.add(frozenset([below, above]))
    return splits

def robinson_foulds(tree1, tree2):
    s1, s2 = get_splits(tree1), get_splits(tree2)
    return len(s1.symmetric_difference(s2))
```

### Dot plot (filtered)
```python
import numpy as np

def filtered_dotplot(seq1, seq2, window=11, stringency=0.7):
    threshold = int(window * stringency)
    matrix = np.zeros((len(seq2), len(seq1)), dtype=np.int8)
    for i in range(len(seq1) - window + 1):
        for j in range(len(seq2) - window + 1):
            matches = sum(seq1[i+k] == seq2[j+k] for k in range(window))
            if matches >= threshold:
                matrix[j + window//2, i + window//2] = 1
    return matrix
# For DNA: also run on reverse complement to detect inversions
```

### Ortholog detection (BBH – bidirectional best hits)
```python
# Blast seq_a vs seq_b and seq_b vs seq_a
# A and B are orthologs if A is B's best hit AND B is A's best hit
def bbh(blast_a_vs_b, blast_b_vs_a):
    best_ab = {q: hits[0] for q, hits in blast_a_vs_b.items()}
    best_ba = {q: hits[0] for q, hits in blast_b_vs_a.items()}
    orthologs = []
    for gene_a, gene_b in best_ab.items():
        if best_ba.get(gene_b) == gene_a:
            orthologs.append((gene_a, gene_b))
    return orthologs
```

## Common Pitfalls

- **UPGMA on rate-varying data**: produces wrong topology when lineage rates differ; use NJ or ML instead
- **Long branch attraction (LBA)**: fast-evolving lineages cluster together artifactually; use ML with good models or break long branches with more taxa
- **Ignoring multiple substitutions**: p-distance saturates above ~0.75; always apply a correction
- **Tip order ≠ evolutionary distance**: rotating branches around a node is cosmetic; read clade membership, not left-to-right order
- **Unrooted ≠ rooted**: NJ/parsimony trees are unrooted; adding an outgroup is required before interpreting directionality of evolution
- **Gene tree ≠ species tree**: HGT, incomplete lineage sorting, and gene duplication mean a gene tree can differ from the species tree
- **Bootstrap inflates with long alignments**: high bootstrap can occur even with systematic errors; ML bootstrap ≥70% is the typical publication threshold
- **Dot plot noise**: use window filtering (window ≥9, stringency ≥0.7 for DNA); single-character plots are unreliable
- **Paralog annotation transfer**: only 1:1 orthologs reliably share function; in/out-paralogs may have diverged

## Related Skills
- `sequence-alignment` — MSA quality directly determines tree accuracy
- `biopython-databases` — distance calculations, substitution models, sequence retrieval
- `ngs-variant-calling` — whole-genome data for comparative analyses
