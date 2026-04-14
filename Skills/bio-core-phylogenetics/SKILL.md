---
name: bio-core-phylogenetics
description: "Phylogenetics with BioPython: distance models (p-distance, JC69, K2P), UPGMA vs NJ tree construction, Newick parsing, and bootstrap interpretation."
tool_type: python
primary_tool: NumPy
---

# Phylogenetics: Reconstructing Evolutionary History

## Pitfalls

- **Tip order does not imply closeness:** Only shared branching points (nodes) indicate relationships. Rotating branches around any internal node produces an equivalent tree. Always look at which node is shared, not which tips are adjacent.
- **UPGMA assumes a molecular clock:** Produces rooted, ultrametric trees (all tips equidistant from root). This assumption is frequently violated. Neighbor-Joining (NJ) does not assume a clock — prefer NJ for most datasets.
- **Branch lengths have units:** In sequence-based trees, branch lengths = expected substitutions per site. They are not time units unless the tree is explicitly time-calibrated.
- **Bootstrap values are NOT probabilities:** Bootstrap 95 means 95% of replicates recovered that clade — not 95% probability of being correct. Thresholds: ≥70 for NJ/parsimony, ≥80 for ML.
- **Model selection matters for ML:** Default to IQ-TREE's built-in model selection (`-m TEST`) rather than hardcoding JC69. GTR+Γ+I is common for nucleotide data; LG or WAG for protein.
- **p-distance underestimates true distance:** Multiple hits at the same site are invisible. JC69 correction becomes undefined at p ≥ 0.75 (saturation).

## Distance Models

| Model | Formula | Assumption |
|---|---|---|
| p-distance | diffs / compared_sites | No correction |
| JC69 | −3/4 · ln(1 − 4p/3) | Equal rates, equal base freq |
| K2P | −½·ln(1−2S−V) − ¼·ln(1−2V) | Distinguishes transitions (S) from transversions (V) |

```python
import numpy as np

def p_distance(seq1, seq2):
    pairs = [(a, b) for a, b in zip(seq1, seq2) if a != '-' and b != '-']
    diffs = sum(1 for a, b in pairs if a != b)
    return diffs / len(pairs) if pairs else 0

def jukes_cantor(seq1, seq2):
    p = p_distance(seq1, seq2)
    if p >= 0.75:
        return float('inf')
    return -0.75 * np.log(1 - 4 * p / 3)

def kimura_k2p(seq1, seq2):
    purines, pyrimidines = set('AG'), set('CT')
    pairs = [(a, b) for a, b in zip(seq1, seq2) if a != '-' and b != '-']
    n = len(pairs)
    S = sum(1 for a, b in pairs if a != b and
            ((a in purines and b in purines) or (a in pyrimidines and b in pyrimidines))) / n
    V = sum(1 for a, b in pairs if a != b and
            not ((a in purines and b in purines) or (a in pyrimidines and b in pyrimidines))) / n
    t1, t2 = 1 - 2*S - V, 1 - 2*V
    if t1 <= 0 or t2 <= 0:
        return float('inf')
    return -0.5 * np.log(t1) - 0.25 * np.log(t2)
```

## BioPython Tree Construction

```python
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import io

# From an alignment
aln = AlignIO.read("alignment.fasta", "fasta")
calc = DistanceCalculator('identity')   # or 'blosum62' for protein
dm = calc.get_distance(aln)

constructor = DistanceTreeConstructor()
nj_tree   = constructor.nj(dm)      # Neighbor-Joining (preferred)
upgma_tree = constructor.upgma(dm)  # UPGMA (assumes clock)

Phylo.draw_ascii(nj_tree)
nj_tree.distance('Human', 'Chimp')  # pairwise distance on tree
```

## Newick Format

```python
import io
from Bio import Phylo

newick = "((Human:0.01,Chimp:0.012):0.02,(Mouse:0.25,Rat:0.23):0.3,Zebrafish:0.6);"
tree = Phylo.read(io.StringIO(newick), "newick")

# Navigate
tree.get_terminals()          # leaf nodes
tree.get_nonterminals()       # internal nodes
tree.distance('Human', 'Chimp')
tree.common_ancestor('Human', 'Chimp')
```

## Assembler/Tool Reference

| Tool | Algorithm | Use case |
|---|---|---|
| BioPython `DistanceTreeConstructor` | NJ, UPGMA | Quick trees from distance matrix |
| IQ-TREE2 | ML + model selection | Publication-quality ML trees |
| RAxML-NG | ML | Large datasets, fast |
| MrBayes | Bayesian | Posterior probability support |
| FastTree | Approximate ML | Very large alignments |

```bash
# IQ-TREE with automatic model selection and 1000 ultrafast bootstraps
iqtree2 -s alignment.fasta -m TEST -B 1000 -T AUTO
```
