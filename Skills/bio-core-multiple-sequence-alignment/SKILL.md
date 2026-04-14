---
name: bio-core-multiple-sequence-alignment
description: Progressive MSA algorithms, tool selection by dataset size, guide tree construction, and profile alignment
tool_type: python
primary_tool: NumPy
---

# Multiple Sequence Alignment

## Pitfalls

- **Progressive alignment errors propagate**: Early misalignments are locked in. MAFFT iterative refinement and T-Coffee's library approach mitigate this.
- **Tool selection by dataset size**: ClustalW is O(N^2), slow beyond ~200 seqs. MUSCLE handles ~1K seqs. MAFFT `--auto` scales to tens of thousands; `--parttree` for >100K. T-Coffee is highest quality but too slow beyond ~500.
- **Gap treatment**: Columns with >50% gaps are unreliable -- mask with trimAl or Gblocks before phylogenetic analysis.
- **SP score vs Column Score**: SP rewards each correct pair; CS requires all pairs correct in a column. CS is more stringent for benchmarking.
- **Protein-coding genes**: Align protein sequences, then back-translate to codons. Direct nucleotide alignment is misleading due to synonymous substitution saturation at 3rd codon position.
- **Coordinate systems**: BED = 0-based half-open; VCF/GFF = 1-based inclusive.

## Exact MSA: 3D DP (Demo Only)

```python
import numpy as np

def exact_msa_three_sequences(seq1, seq2, seq3, match=1, mismatch=-1, gap=-2):
    """3D DP for three short sequences. O(n^3) -- impractical beyond ~10 residues."""
    n1, n2, n3 = len(seq1), len(seq2), len(seq3)
    dp = np.full((n1+1, n2+1, n3+1), -np.inf)
    dp[0, 0, 0] = 0
    def score(a, b):
        return gap if a == '-' or b == '-' else (match if a == b else mismatch)
    for i in range(1, n1+1): dp[i,0,0] = dp[i-1,0,0] + 2*gap
    for j in range(1, n2+1): dp[0,j,0] = dp[0,j-1,0] + 2*gap
    for k in range(1, n3+1): dp[0,0,k] = dp[0,0,k-1] + 2*gap
    # Fill -- 7 transitions per cell (all combos of match/gap across 3 seqs)
    for i in range(1, n1+1):
        for j in range(1, n2+1):
            for k in range(1, n3+1):
                dp[i,j,k] = max(
                    dp[i-1,j-1,k-1] + score(seq1[i-1],seq2[j-1]) + score(seq1[i-1],seq3[k-1]) + score(seq2[j-1],seq3[k-1]),
                    dp[i,j-1,k-1] + score(seq2[j-1],seq3[k-1]) + gap,
                    dp[i-1,j,k-1] + score(seq1[i-1],seq3[k-1]) + gap,
                    dp[i-1,j-1,k] + score(seq1[i-1],seq2[j-1]) + gap,
                    dp[i,j,k-1] + 2*gap, dp[i,j-1,k] + 2*gap, dp[i-1,j,k] + 2*gap)
    return dp[n1, n2, n3]
# 3 seqs len 100 = 101^3 ~ 1M cells; 10 seqs len 100 = 101^10 -- impossible
```

## Progressive Alignment Pipeline

1. All-vs-all pairwise distances (k-mer or alignment-based)
2. Build guide tree (UPGMA or NJ)
3. Align sequences following tree order (profiles merge at each node)

## K-mer Distance and UPGMA

```python
def kmer_distance(seq1, seq2, k=3):
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1)-k+1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2)-k+1))
    if not kmers1 and not kmers2: return 0.0
    return 1.0 - len(kmers1 & kmers2) / len(kmers1 | kmers2)

def build_distance_matrix(sequences, k=3):
    n = len(sequences)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = kmer_distance(sequences[i], sequences[j], k)
            matrix[i,j] = matrix[j,i] = d
    return matrix

def upgma(dist_matrix, names):
    """UPGMA clustering -> Newick string. Assumes molecular clock (often false)."""
    n = len(names)
    D = dist_matrix.copy()
    node_names = list(names)
    cluster_sizes = [1] * n
    while len(node_names) > 1:
        m = len(node_names)
        min_dist, mi, mj = np.inf, 0, 1
        for i in range(m):
            for j in range(i+1, m):
                if D[i,j] < min_dist: min_dist, mi, mj = D[i,j], i, j
        bl = min_dist / 2
        new_name = f"({node_names[mi]}:{bl:.4f},{node_names[mj]}:{bl:.4f})"
        new_size = cluster_sizes[mi] + cluster_sizes[mj]
        idx_map = [k for k in range(m) if k != mi and k != mj]
        new_D = np.zeros((len(idx_map)+1, len(idx_map)+1))
        for a, ka in enumerate(idx_map):
            for b, kb in enumerate(idx_map):
                new_D[a,b] = D[ka,kb]
            new_D[a,-1] = new_D[-1,a] = (cluster_sizes[mi]*D[ka,mi] + cluster_sizes[mj]*D[ka,mj]) / new_size
        node_names = [node_names[k] for k in idx_map] + [new_name]
        cluster_sizes = [cluster_sizes[k] for k in idx_map] + [new_size]
        D = new_D
    return node_names[0] + ";"
```

## Profile Alignment Scoring

When aligning two sub-alignments, score each column pair as the average of all pairwise substitution scores:

S(col1, col2) = sum(s(a,b) for a in col1 for b in col2) / (|col1| * |col2|)
