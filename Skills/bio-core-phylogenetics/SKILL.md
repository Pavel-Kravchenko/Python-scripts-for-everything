---
name: bio-core-phylogenetics
description: "Phylogenetic trees are hypotheses about the evolutionary relationships among organisms, genes, or other biological entities. This notebook covers the theory, algorithms, and practical tools for buildi"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/06_Phylogenetics/01_phylogenetics.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Phylogenetics: Reconstructing Evolutionary History

*Source: Course notebook `Tier_2_Core_Bioinformatics/06_Phylogenetics/01_phylogenetics.ipynb`*


Phylogenetic trees are hypotheses about the evolutionary relationships among organisms, genes, or other biological entities. This notebook covers the theory, algorithms, and practical tools for building and interpreting phylogenetic trees.

---

## Learning Objectives

- Master phylogenetic tree terminology (nodes, branches, clades, rooted vs unrooted)
- Implement distance-based methods from scratch: UPGMA and Neighbor-Joining
- Understand distance matrices and models of sequence evolution (Jukes-Cantor, Kimura)
- Grasp character-based and probabilistic methods (Parsimony, Maximum Likelihood, Bayesian)
- Perform and interpret bootstrap analysis
- Build and visualize trees with BioPython (Bio.Phylo)
- Interpret phylogenetic trees in biological context

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **Tip order does not imply closeness**: A common mistake is reading left-to-right tip order as indicating evolutionary closeness. Only *shared branching points* (nodes) indicate relationships. Rotating branches around any internal node produces an equivalent tree — the same topology, different visual arrangement. Always look at which node is shared, not which tips are adjacent.
- **UPGMA assumes a molecular clock**: UPGMA produces a rooted, ultrametric tree (all tips equidistant from the root) by assuming all lineages evolve at the same rate. This assumption is often violated. Neighbor-Joining (NJ) does not assume a molecular clock and is preferred for most datasets.
- **Branch lengths have units**: In phylogenies inferred from sequences, branch lengths typically represent *expected substitutions per site* — the average number of base/amino acid changes per position along that branch. They are not time units unless the tree is explicitly time-calibrated.
- **Bootstrap values are NOT probabilities**: A bootstrap value of 95 does not mean there is a 95% probability the clade is correct. It means 95% of bootstrap replicates recovered that clade. Values ≥70 are generally considered well-supported for NJ/parsimony; for maximum likelihood, values ≥80 are more stringent.
- **Model of evolution matters for ML**: For nucleotide sequences, models range from simple (JC69: equal base frequencies, equal rates) to complex (GTR+Γ+I: unequal rates, gamma-distributed rate variation, invariant sites). Use ModelTest/IQ-TREE's built-in model selection rather than always defaulting to JC69.

## Environment check (run this first)

```python
# Environment check
import numpy as np
import matplotlib.pyplot as plt
import io
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

print("BioPython Phylo modules loaded.")

# Quick tree parsing demo
newick = "((Human:0.01,Chimp:0.012):0.02,(Mouse:0.25,Rat:0.23):0.3,Zebrafish:0.6);"
tree = Phylo.read(io.StringIO(newick), "newick")
print(f"\nExample primate tree:")
Phylo.draw_ascii(tree)
print(f"Terminals: {[c.name for c in tree.get_terminals()]}")
print(f"Human-Chimp distance: {tree.distance('Human','Chimp'):.3f}")
print(f"Human-Mouse distance: {tree.distance('Human','Mouse'):.3f}")
print("\nProceed to Section 1.")
```python

```python
import re

def parse_newick(newick):
    """
    Parse a Newick string into a PhyloNode tree.
    Handles names, branch lengths, and nested clades.
    """
    newick = newick.strip().rstrip(';')
    
    def _parse(s, pos=0):
        node = PhyloNode()
        children = []
        
        if s[pos] == '(':
            pos += 1  # skip '('
            while True:
                child, pos = _parse(s, pos)
                children.append(child)
                if pos < len(s) and s[pos] == ',':
                    pos += 1
                elif pos < len(s) and s[pos] == ')':
                    pos += 1
                    break
                else:
                    break
            for child in children:
                node.add_child(child)
        
        # Parse name and branch length
        name = ''
        while pos < len(s) and s[pos] not in ',):;':
            name += s[pos]
            pos += 1
        
        if ':' in name:
            parts = name.split(':')
            node.name = parts[0] if parts[0] else None
            node.branch_length = float(parts[1]) if parts[1] else 0.0
        elif name:
            node.name = name
        
        return node, pos
    
    tree, _ = _parse(newick, 0)
    return tree


def print_tree_ascii(node, prefix="", is_last=True, is_root=True):
    """
    Print an ASCII representation of the tree.
    """
    connector = "" if is_root else ("\u2514\u2500 " if is_last else "\u251c\u2500 ")
    name = node.name or ""
    bl = f":{node.branch_length:.4f}" if node.branch_length > 0 else ""
    print(f"{prefix}{connector}{name}{bl}")
    
    child_prefix = prefix + ("   " if is_last or is_root else "\u2502  ")
    for i, child in enumerate(node.children):
        print_tree_ascii(child, child_prefix, i == len(node.children) - 1, False)


# Parse and display
newick_str = "((Human:0.1,Chimp:0.12):0.2,(Mouse:0.25,Rat:0.23):0.15);"
tree = parse_newick(newick_str)

print(f"Newick: {newick_str}")
print(f"\nTree structure:")
print_tree_ascii(tree)
```python

---

## 3. Distance Matrices and Models of Evolution

### 3.1 From Sequences to Distances

The simplest distance between two aligned sequences is the **p-distance** -- the proportion of sites that differ:

$$p = \frac{\text{number of differences}}{\text{number of compared sites}}$$

However, p-distance **underestimates** the true evolutionary distance because it ignores **multiple substitutions** at the same site (back-mutations, parallel mutations).

### 3.2 The Jukes-Cantor Model (JC69)

The simplest substitution model assumes:
- All substitutions are equally likely
- All nucleotides have equal frequency (25% each)

The JC69-corrected distance:

$$d = -\frac{3}{4} \ln\left(1 - \frac{4p}{3}\right)$$

This correction accounts for unobserved substitutions. Note that $d$ is undefined when $p \geq 0.75$ (sequences are so divergent that the model breaks down).

### 3.3 The Kimura Two-Parameter Model (K2P)

A more realistic model that distinguishes **transitions** (purine-to-purine: A<->G, or pyrimidine-to-pyrimidine: C<->T) from **transversions** (purine-to-pyrimidine and vice versa):

$$d = -\frac{1}{2} \ln(1 - 2S - V) - \frac{1}{4} \ln(1 - 2V)$$

Where $S$ = proportion of transitions and $V$ = proportion of transversions.

```python
import numpy as np
import matplotlib.pyplot as plt

def p_distance(seq1, seq2):
    """Proportion of differing sites (ignoring gaps)."""
    diffs = 0
    compared = 0
    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            compared += 1
            if a != b:
                diffs += 1
    return diffs / compared if compared > 0 else 0


def jukes_cantor_distance(seq1, seq2):
    """JC69-corrected distance."""
    p = p_distance(seq1, seq2)
    if p >= 0.75:
        return float('inf')  # Sequences too divergent
    return -0.75 * np.log(1 - 4 * p / 3)


def kimura_distance(seq1, seq2):
    """Kimura 2-parameter corrected distance."""
    transitions = 0
    transversions = 0
    compared = 0
    
    purines = set('AG')
    pyrimidines = set('CT')
    
    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            compared += 1
            if a != b:
                if (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
                    transitions += 1
                else:
                    transversions += 1
    
    if compared == 0:
        return 0
    
    S = transitions / compared
    V = transversions / compared
    
    term1 = 1 - 2 * S - V
    term2 = 1 - 2 * V
    
    if term1 <= 0 or term2 <= 0:
        return float('inf')
    
    return -0.5 * np.log(term1) - 0.25 * np.log(term2)


# Compare distance metrics
seq_a = "ATCGATCGATCGATCGATCG"
seq_b = "ATCGTTCGATCGATCGATCG"  # 1 transition
seq_c = "TTCGATCGATCGATCAATCG"  # 1 transversion + 1 transversion

print("Distance comparison:")
print(f"{'Pair':<15} {'p-distance':>12} {'JC69':>12} {'Kimura':>12}")
print("-" * 55)

for name, s1, s2 in [("A vs B", seq_a, seq_b), ("A vs C", seq_a, seq_c), ("B vs C", seq_b, seq_c)]:
    pd = p_distance(s1, s2)
    jc = jukes_cantor_distance(s1, s2)
    k2 = kimura_distance(s1, s2)
    print(f"{name:<15} {pd:>12.4f} {jc:>12.4f} {k2:>12.4f}")
```python

```python
# Visualize how corrections diverge from p-distance
p_values = np.arange(0.01, 0.74, 0.01)
jc_values = [-0.75 * np.log(1 - 4 * p / 3) for p in p_values]

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(p_values, p_values, 'b-', linewidth=2, label='p-distance (observed)')
ax.plot(p_values, jc_values, 'r-', linewidth=2, label='JC69-corrected (estimated actual)')
ax.plot([0, 0.75], [0, 0.75], 'k:', alpha=0.3)

ax.set_xlabel('Observed p-distance', fontsize=12)
ax.set_ylabel('Estimated evolutionary distance', fontsize=12)
ax.set_title('Jukes-Cantor Correction: Why p-distance Underestimates', fontsize=13)
ax.legend(fontsize=11)
ax.set_xlim(0, 0.75)
ax.set_ylim(0, 2.0)
ax.grid(True, alpha=0.3)

# Annotate the divergence
ax.annotate('Multiple substitutions\nbecome significant',
            xy=(0.5, 0.5), xytext=(0.3, 1.2),
            arrowprops=dict(arrowstyle='->', color='gray'),
            fontsize=10, color='gray')

plt.tight_layout()
plt.show()

print("Key insight: at p=0.30, JC69 estimates d=0.38 -- 27% more than observed.")
print("At p=0.50, JC69 estimates d=0.82 -- 64% more than observed.")
print("At p>=0.75, saturation -- random expectation for 4 equally likely nucleotides.")
```python

```python
def build_distance_matrix(sequences, names, method='jukes-cantor'):
    """
    Calculate pairwise distance matrix from aligned sequences.
    """
    n = len(sequences)
    matrix = np.zeros((n, n))
    
    distance_fn = {
        'p-distance': p_distance,
        'jukes-cantor': jukes_cantor_distance,
        'kimura': kimura_distance,
    }[method]
    
    for i in range(n):
        for j in range(i + 1, n):
            d = distance_fn(sequences[i], sequences[j])
            matrix[i, j] = matrix[j, i] = d
    
    return matrix


def print_distance_matrix(matrix, names):
    """Pretty-print a distance matrix."""
    n = len(names)
    max_name = max(len(name) for name in names)
    
    header = ' ' * (max_name + 2)
    for name in names:
        header += f"{name:>8}"
    print(header)
    
    for i in range(n):
        row = f"{names[i]:>{max_name}}  "
        for j in range(n):
            row += f"{matrix[i, j]:>8.4f}"
        print(row)


# Example with primate and rodent sequences
aligned_seqs = [
    "ATCGATCGATCGATCGATCG",  # Species A
    "ATCGTTCGATCGATCGATCG",  # Species B (1 change from A)
    "TTCGATCGATCGATCGATCG",  # Species C (1 change from A)
    "TTCGTTCGATCGATCGATCG",  # Species D (2 changes from A)
    "AACCGGTTAACCGGTTAACC",  # Species E (very different)
]
seq_names = ["Sp_A", "Sp_B", "Sp_C", "Sp_D", "Sp_E"]

dm = build_distance_matrix(aligned_seqs, seq_names, method='jukes-cantor')
print("Distance matrix (Jukes-Cantor):")
print_distance_matrix(dm, seq_names)
```python

```python
def plot_distance_matrix(matrix, names, title="Distance Matrix"):
    """Visualize distance matrix as a heatmap."""
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(matrix, cmap='YlOrRd', vmin=0)
    
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=10)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=10)
    
    for i in range(len(names)):
        for j in range(len(names)):
            color = 'white' if matrix[i, j] > matrix.max() * 0.6 else 'black'
            ax.text(j, i, f"{matrix[i, j]:.3f}", ha='center', va='center',
                    fontsize=9, color=color)
    
    plt.colorbar(im, label='Distance')
    ax.set_title(title, fontsize=13)
    plt.tight_layout()
    plt.show()


plot_distance_matrix(dm, seq_names, title="JC69 Distance Matrix")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
