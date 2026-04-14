---
name: bio-core-pairwise-sequence-alignment
description: Pairwise Sequence Alignment with NumPy
tool_type: python
primary_tool: NumPy
---

# Pairwise Sequence Alignment

## Key Concepts

- **Global vs. local**: Needleman-Wunsch aligns end-to-end — use for full-length orthologs of similar length. Smith-Waterman finds the best-matching subregion — use for shared domains or short query vs. long subject. BLAST uses heuristic local alignment.
- **BLOSUM numbering**: Higher number = more conserved sequences used to build it. BLOSUM80 for close sequences; BLOSUM45 for distant; BLOSUM62 is BLAST default. **PAM is opposite**: higher PAM = more divergent (PAM250 ≈ 250 mutations per 100 sites).
- **Affine gap penalties**: Linear adds same cost per position. Affine charges large opening (e.g. −10) + small extension (e.g. −0.5). Biologically realistic — insertions/deletions occur in runs, not scattered single positions.
- **Score comparability**: A score of 50 for 100 residues ≠ score of 50 for 20 residues. Use percent identity and bit scores for comparison.

## Environment Check

```python
import numpy as np
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

seq1 = Seq("MVLSPADKTNVKAAWGKVGAHAG")  # HBA N-term
seq2 = Seq("MVHLTPEEKSAVTALWGKVNVDE")  # HBB N-term
best = list(aligner.align(seq1, seq2))[0]
print(best); print(f"Score: {best.score:.1f}")
```

## Dot Plots

```python
def windowed_dot_plot(seq1, seq2, window=3, threshold=2, title="Windowed Dot Plot"):
    s1, s2 = str(seq1), str(seq2)
    rows = len(s2) - window + 1
    cols = len(s1) - window + 1
    matrix = np.zeros((rows, cols), dtype=int)
    for i in range(rows):
        for j in range(cols):
            if sum(1 for k in range(window) if s2[i+k] == s1[j+k]) >= threshold:
                matrix[i, j] = 1
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(matrix, cmap="Greys", aspect="auto", interpolation="none")
    ax.set_xlabel("Sequence 1"); ax.set_ylabel("Sequence 2"); ax.set_title(title)
    plt.tight_layout(); plt.show()
```

Parallel diagonals in self-comparison reveal internal repeats.

## Substitution Matrices

For DNA, +1/−1 match/mismatch is often sufficient. For proteins, substitution matrices encode log-odds of observed substitution frequencies:

$$S(a, b) = \log_2 \frac{q_{ab}}{p_a \cdot p_b}$$

Positive = conservative substitution; negative = disruptive; diagonal always positive.

**Amino acid conservative groups:**

| Group | Amino acids |
|-------|-------------|
| Small nonpolar | A, V, I, L, M |
| Aromatic | F, Y, W |
| Polar uncharged | S, T, N, Q |
| Positively charged | K, R, H |
| Negatively charged | D, E |
| Special | C, P, G |

**PAM vs BLOSUM:**

| Feature | PAM | BLOSUM |
|---|---|---|
| Developed by | Dayhoff (1978) | Henikoff & Henikoff (1992) |
| Basis | Closely related, extrapolated | Conserved blocks, diverse proteins |
| Numbering | Higher = more distant (PAM250) | Higher = more similar (BLOSUM80) |
| Default (BLAST) | PAM250 | **BLOSUM62** |

PAM1 = 1% divergence, extrapolated by matrix exponentiation. BLOSUM directly computed from conserved blocks at given % identity threshold.

```python
blosum62 = substitution_matrices.load("BLOSUM62")
pairs = [("L","I","hydrophobic"), ("D","E","neg. charged"), ("L","D","disruptive")]
for a, b, reason in pairs:
    print(f"{a}<->{b}: {blosum62[a,b]:+d}  {reason}")
```

## Gap Penalties

**Linear:** $\text{cost}(k) = -d \times k$

**Affine:** $\text{cost}(k) = -d - (k-1) \times e$ where $e < d$

| Scheme | Typical values |
|---|---|
| Linear | $d = 8$ |
| Affine | $d = 10$, $e = 0.5$ |

Affine makes opening expensive but extension cheap — one long gap is biologically more plausible than scattered single-position gaps (replication slippage).

## Needleman-Wunsch (Global Alignment)

Recurrence relation:

$$F(i, j) = \max \begin{cases} F(i-1, j-1) + s(x_j, y_i) \\ F(i-1, j) - d \\ F(i, j-1) - d \end{cases}$$

Initialize first row/column with cumulative gap penalties. Traceback from bottom-right corner.

- **Time:** $O(m \times n)$
- **Space:** $O(m \times n)$ for full traceback matrix

## Pitfalls

- **Coordinate systems**: BED 0-based half-open; VCF/GFF 1-based inclusive — mixing causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) for thousands of features
