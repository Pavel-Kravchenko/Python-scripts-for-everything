---
name: bio-core-multiple-sequence-alignment
description: "Multiple Sequence Alignment (MSA) is one of the most fundamental techniques in bioinformatics. While pairwise alignment compares two sequences, MSA simultaneously aligns three or more sequences to rev"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/05_Multiple_Sequence_Alignment/01_multiple_sequence_alignment.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Multiple Sequence Alignment: From Theory to Practice

*Source: Course notebook `Tier_2_Core_Bioinformatics/05_Multiple_Sequence_Alignment/01_multiple_sequence_alignment.ipynb`*

# Multiple Sequence Alignment: From Theory to Practice

Multiple Sequence Alignment (MSA) is one of the most fundamental techniques in bioinformatics. While pairwise alignment compares two sequences, MSA simultaneously aligns three or more sequences to reveal conserved regions, functional domains, and evolutionary relationships that pairwise comparisons would miss.

---

## Learning Objectives

- Understand why MSA is essential and what biological questions it answers
- Learn the progressive alignment strategy and its algorithmic foundations
- Score MSA quality using sum-of-pairs and column-based metrics
- Use major MSA tools (ClustalW, MUSCLE, MAFFT, T-Coffee)
- Read, write, and manipulate alignments with BioPython AlignIO
- Compute consensus sequences and conservation scores
- Build and interpret sequence logos
- Visualize alignments programmatically

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **Progressive alignment errors propagate**: In progressive alignment, once two sequences are aligned incorrectly in an early step, that error is locked in for all subsequent alignments. This is the main weakness of ClustalW and MUSCLE default mode. MAFFT's iterative refinement and T-Coffee's library approach help mitigate this.
- **Tool selection by dataset size**: ClustalW is slow for large datasets (guide tree computation is O(N²)); MUSCLE is fast and accurate up to ~1,000 sequences; MAFFT with `--auto` scales to tens of thousands. For >100,000 sequences, use MAFFT `--parttree` or Clustal Omega. T-Coffee produces highest quality but is too slow beyond ~500 sequences.
- **Gap treatment is critical**: Columns with >50% gaps are often unreliable and should be masked before phylogenetic analysis. Tools like trimAl and Gblocks remove poorly aligned columns.
- **Sum-of-Pairs (SP) score vs. Column Score (CS)**: SP score rewards each correctly aligned pair of residues; CS rewards only columns where all pairs are correct. For benchmarking, CS is more stringent. Real alignments are evaluated against reference alignments from databases like BAliBASE and HOMSTRAD.
- **Amino acid vs. nucleotide alignment**: For protein-coding genes, align protein sequences and then back-translate to codons (codon-aware alignment). Direct nucleotide alignment of coding regions is often misleading because synonymous substitutions saturate at the third codon position.

## Environment check (run this first)

```python
# Environment check
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment, PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import io, subprocess, shutil

print("Core imports ready.")

# Check which external MSA tools are installed
for tool in ['mafft', 'muscle', 'clustalw', 'clustalo']:
    found = shutil.which(tool)
    status = f"found at {found}" if found else "NOT FOUND (install via conda/brew)"
    print(f"  {tool}: {status}")

print("\nNote: if no external MSA tools are found, the notebook will use")
print("a built-in progressive alignment implementation for demonstrations.")
```

```python
def exact_msa_three_sequences(seq1, seq2, seq3, match=1, mismatch=-1, gap=-2):
    """
    Exact MSA for three short sequences via 3D dynamic programming.
    Only practical for very short sequences (demonstration purposes).
    """
    n1, n2, n3 = len(seq1), len(seq2), len(seq3)
    
    def score(a, b):
        if a == '-' or b == '-':
            return gap
        return match if a == b else mismatch
    
    # Initialize 3D DP matrix
    dp = np.full((n1 + 1, n2 + 1, n3 + 1), -np.inf)
    dp[0, 0, 0] = 0
    
    # Fill borders
    for i in range(1, n1 + 1):
        dp[i, 0, 0] = dp[i-1, 0, 0] + 2 * gap
    for j in range(1, n2 + 1):
        dp[0, j, 0] = dp[0, j-1, 0] + 2 * gap
    for k in range(1, n3 + 1):
        dp[0, 0, k] = dp[0, 0, k-1] + 2 * gap
    
    # Fill edges
    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            dp[i, j, 0] = max(
                dp[i-1, j-1, 0] + score(seq1[i-1], seq2[j-1]) + 2 * gap,
                dp[i-1, j, 0] + 2 * gap,
                dp[i, j-1, 0] + 2 * gap
            )
    for i in range(1, n1 + 1):
        for k in range(1, n3 + 1):
            dp[i, 0, k] = max(
                dp[i-1, 0, k-1] + score(seq1[i-1], seq3[k-1]) + 2 * gap,
                dp[i-1, 0, k] + 2 * gap,
                dp[i, 0, k-1] + 2 * gap
            )
    for j in range(1, n2 + 1):
        for k in range(1, n3 + 1):
            dp[0, j, k] = max(
                dp[0, j-1, k-1] + score(seq2[j-1], seq3[k-1]) + 2 * gap,
                dp[0, j-1, k] + 2 * gap,
                dp[0, j, k-1] + 2 * gap
            )
    
    # Fill 3D matrix
    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            for k in range(1, n3 + 1):
                dp[i, j, k] = max(
                    dp[i-1, j-1, k-1] + score(seq1[i-1], seq2[j-1]) + score(seq1[i-1], seq3[k-1]) + score(seq2[j-1], seq3[k-1]),
                    dp[i, j-1, k-1] + score(seq2[j-1], seq3[k-1]) + gap,
                    dp[i-1, j, k-1] + score(seq1[i-1], seq3[k-1]) + gap,
                    dp[i-1, j-1, k] + score(seq1[i-1], seq2[j-1]) + gap,
                    dp[i, j, k-1] + 2 * gap,
                    dp[i, j-1, k] + 2 * gap,
                    dp[i-1, j, k] + 2 * gap
                )
    
    print(f"3D DP matrix size: {n1+1} x {n2+1} x {n3+1} = {(n1+1)*(n2+1)*(n3+1)} cells")
    print(f"Optimal score: {dp[n1, n2, n3]}")
    return dp[n1, n2, n3]

# Demo with very short sequences
s1 = "ATCG"
s2 = "ATGG"
s3 = "ACCG"

optimal = exact_msa_three_sequences(s1, s2, s3)
print(f"\nFor 3 sequences of length ~4: manageable")
print(f"For 3 sequences of length 100: {101**3:,} cells")
print(f"For 10 sequences of length 100: {101**10:.2e} cells -- impossible!")
```

---

## 3. Progressive Alignment

Since exact MSA is intractable, practical tools use **heuristic** approaches. The most widely used is **progressive alignment**, pioneered by Feng and Doolittle (1987) and implemented in ClustalW.

### The Progressive Strategy

```
Progressive Alignment Pipeline
================================

Step 1: All-vs-all pairwise alignment
  s1 vs s2, s1 vs s3, ..., s4 vs s5
       |
       v
Step 2: Build distance matrix from pairwise scores
       |                s1    s2    s3    s4    s5
       |          s1  [0.00  0.23  0.45  0.67  0.71]
       |          s2  [0.23  0.00  0.38  0.62  0.65]
       |          ...
       v
Step 3: Build guide tree (UPGMA or NJ)
       |          +----- s1
       |       +--|
       |       |  +----- s3
       |    +--|
       |    |  +-------- s2
       |  --|
       |    |  +-------- s4
       |    +--|
       |       +-------- s5
       v
Step 4: Align sequences following guide tree order
  - Align s1 + s3  (most similar pair)
  - Align (s1,s3) + s2  (alignment to sequence)
  - Align s4 + s5
  - Align (s1,s3,s2) + (s4,s5)  (alignment to alignment)
       |
       v
Step 5: Final MSA
```

### The Chicken-and-Egg Problem

We need a tree to guide the alignment, but we need an alignment to build a good tree. Progressive alignment solves this by:
1. Building a **rough** guide tree from pairwise distances (no MSA needed)
2. Using this guide tree to produce an MSA
3. Optionally iterating: use the MSA to build a better tree, then re-align

```python
def kmer_distance(seq1, seq2, k=3):
    """
    Alignment-free distance based on shared k-mers.
    Returns fraction of k-mers unique to either sequence.
    """
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
    
    if not kmers1 and not kmers2:
        return 0.0
    
    shared = kmers1 & kmers2
    total = kmers1 | kmers2
    
    return 1.0 - len(shared) / len(total)


def build_distance_matrix(sequences, names, k=3):
    """
    Build a distance matrix from sequences using k-mer distance.
    """
    n = len(sequences)
    matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i + 1, n):
            d = kmer_distance(sequences[i], sequences[j], k)
            matrix[i, j] = matrix[j, i] = d
    
    return matrix


# Example sequences (globin-like)
sequences = {
    'HBA_HUMAN':  'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH',
    'HBA_MOUSE':  'MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSH',
    'HBB_HUMAN':  'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLST',
    'HBB_MOUSE':  'MVHLTDAEKAAVNGLWGKVNSDEVGGEALGRLLVVYPWTQRYFDSFGDLSS',
    'MYG_HUMAN':  'MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLK',
}

names = list(sequences.keys())
seqs = list(sequences.values())

dm = build_distance_matrix(seqs, names)

print("K-mer distance matrix (k=3):")
print(f"{'':>12}", end='')
for name in names:
    print(f"{name:>12}", end='')
print()

for i, name in enumerate(names):
    print(f"{name:>12}", end='')
    for j in range(len(names)):
        print(f"{dm[i,j]:>12.3f}", end='')
    print()
```

### 3.1 Building the Guide Tree with UPGMA

UPGMA (Unweighted Pair Group Method with Arithmetic Mean) is the simplest hierarchical clustering method. At each step:

1. Find the closest pair of clusters
2. Merge them into a new cluster
3. Update distances as the arithmetic mean of all pairwise distances
4. Repeat until one cluster remains

UPGMA assumes a molecular clock (equal rates of evolution), which is generally false -- but it is fast and sufficient for building a rough guide tree.

```python
def upgma(dist_matrix, names):
    """
    UPGMA clustering. Returns a Newick-format tree string.
    """
    n = len(names)
    D = dist_matrix.copy()
    node_names = list(names)
    # Track cluster sizes for weighted averaging
    cluster_sizes = [1] * n
    
    while len(node_names) > 1:
        m = len(node_names)
        
        # Find minimum distance (off-diagonal)
        min_dist = np.inf
        min_i, min_j = 0, 1
        for i in range(m):
            for j in range(i + 1, m):
                if D[i, j] < min_dist:
                    min_dist = D[i, j]
                    min_i, min_j = i, j
        
        i, j = min_i, min_j
        branch_length = min_dist / 2
        
        # Create new node
        new_name = f"({node_names[i]}:{branch_length:.4f},{node_names[j]}:{branch_length:.4f})"
        new_size = cluster_sizes[i] + cluster_sizes[j]
        
        # Compute new distances (arithmetic mean)
        new_D = np.zeros((m - 1, m - 1))
        new_names = []
        new_sizes = []
        idx_map = []
        
        for k in range(m):
            if k != i and k != j:
                new_names.append(node_names[k])
                new_sizes.append(cluster_sizes[k])
                idx_map.append(k)
        new_names.append(new_name)
        new_sizes.append(new_size)
        
        # Copy existing distances
        for a, ka in enumerate(idx_map):
            for b, kb in enumerate(idx_map):
                new_D[a, b] = D[ka, kb]
        
        # Compute distances to new merged cluster
        for a, ka in enumerate(idx_map):
            d_new = (cluster_sizes[i] * D[ka, i] + cluster_sizes[j] * D[ka, j]) / new_size
            new_D[a, -1] = new_D[-1, a] = d_new
        
        D = new_D
        node_names = new_names
        cluster_sizes = new_sizes
    
    return node_names[0] + ";"


guide_tree = upgma(dm, names)
print("Guide tree (Newick):")
print(guide_tree)
```

### 3.2 Aligning Alignments: The Profile Concept

In progressive alignment, we don't just align single sequences. We also need to:
- Align a sequence to an existing alignment
- Align two alignments together

This is done by converting each alignment column into a **profile** -- a vector of character frequencies. The substitution score between two profile columns is the **average** of all pairwise scores:

$$
S(\text{col}_1, \text{col}_2) = \frac{\sum_{a \in \text{col}_1} \sum_{b \in \text{col}_2} s(a, b)}{|\text{col}_1| \times |\text{col}_2|}
$$

For example, aligning column `[A, C]` against column `[T, G]` with a substitution matrix $m$:

$$
S = \frac{m[A][T] + m[A][G] + m[C][T] + m[C][G]}{2 \times 2}
$$
