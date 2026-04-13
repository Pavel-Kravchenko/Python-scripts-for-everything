---
name: bio-core-pairwise-sequence-alignment
description: "1. Explain **why** biologists align sequences and what an alignment represents biologically 2. Construct and interpret **dot plots** as a visual alignment tool 3. Describe substitution matrices (PAM, "
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/03_Pairwise_Sequence_Alignment/01_pairwise_sequence_alignment.ipynb"
---

# Pairwise Sequence Alignment

*Source: Course notebook `Tier_2_Core_Bioinformatics/03_Pairwise_Sequence_Alignment/01_pairwise_sequence_alignment.ipynb`*

# Pairwise Sequence Alignment

**Tier 2 -- Core Bioinformatics**

---

## Learning Objectives

By the end of this notebook you will be able to:

1. Explain **why** biologists align sequences and what an alignment represents biologically
2. Construct and interpret **dot plots** as a visual alignment tool
3. Describe substitution matrices (PAM, BLOSUM) and gap penalty schemes (linear, affine)
4. Implement the **Needleman-Wunsch** (global) algorithm from scratch with dynamic programming
5. Implement the **Smith-Waterman** (local) algorithm from scratch
6. Use BioPython's alignment tools for practical work
7. Interpret alignment statistics: scores, E-values, and percent identity

---

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **Global vs. local alignment**: Needleman-Wunsch forces an end-to-end alignment of both sequences — good for comparing full-length orthologs of similar length. Smith-Waterman finds the best-matching subregion — good for finding shared domains in otherwise divergent proteins, or locating a short query within a long subject. BLAST uses a heuristic local alignment.
- **BLOSUM numbering is counterintuitive**: Higher BLOSUM number = more conserved sequences were used to build the matrix. BLOSUM80 is for closely related sequences; BLOSUM45 is for distantly related ones. For default BLAST searches, BLOSUM62 is standard. **PAM is the opposite**: higher PAM number = more divergent (PAM250 = ~250 mutations per 100 sites).
- **Affine vs. linear gap penalties**: Linear gap penalty adds the same cost for each gap position. Affine gap penalty charges a large opening cost (e.g., -10) plus a small extension cost (e.g., -0.5 per additional position). Affine is biologically realistic because insertions/deletions tend to occur in continuous runs, not scattered single positions.
- **Alignment score is not directly comparable**: A score of 50 for a 100-residue alignment is not equivalent to a score of 50 for a 20-residue alignment. Use percent identity and bit scores for comparison.

## Environment check (run this first)

```python
# Environment check
import numpy as np
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

# Quick alignment demo
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

seq1 = Seq("MVLSPADKTNVKAAWGKVGAHAG")  # human hemoglobin alpha N-term
seq2 = Seq("MVHLTPEEKSAVTALWGKVNVDE")  # human hemoglobin beta N-term

alignments = list(aligner.align(seq1, seq2))
best = alignments[0]
print("Global alignment of HBA vs HBB (N-terminal region):")
print(best)
print(f"Score: {best.score:.1f}")
print("\nAll imports ready. Proceed to Section 1.")
```

The strong diagonal indicates high similarity between the two hemoglobin sequences. The few off-diagonal dots are noise from individual matching residues.

```python
def windowed_dot_plot(seq1, seq2, window=3, threshold=2, title="Windowed Dot Plot"):
    """Dot plot with a sliding window to reduce noise.
    
    A dot is placed at (i,j) if at least `threshold` of the `window`
    residues starting at (i,j) match.
    """
    s1, s2 = str(seq1), str(seq2)
    rows = len(s2) - window + 1
    cols = len(s1) - window + 1
    matrix = np.zeros((rows, cols), dtype=int)
    for i in range(rows):
        for j in range(cols):
            matches = sum(1 for k in range(window) if s2[i + k] == s1[j + k])
            if matches >= threshold:
                matrix[i, j] = 1
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(matrix, cmap="Greys", aspect="auto", interpolation="none")
    ax.set_xlabel("Sequence 1")
    ax.set_ylabel("Sequence 2")
    ax.set_title(title)
    plt.tight_layout()
    plt.show()


windowed_dot_plot(hba_human, hba_mouse, window=5, threshold=4,
                  title="Windowed Dot Plot (w=5, t=4)")
```

The windowed version dramatically reduces noise: only regions where several consecutive residues match are shown. This makes the alignment signal much clearer.

### Dot Plots Reveal Structural Features

Let us compare a sequence to itself to reveal internal repeats.

```python
# A sequence with an internal repeat
repeat_seq = "ACGTACGTACGTNNNNACGTACGTACGT"
windowed_dot_plot(repeat_seq, repeat_seq, window=4, threshold=4,
                  title="Self Dot Plot -- Internal Repeat")
```

The parallel diagonals in a self-comparison reveal repeated elements within the sequence.

---

## 3. Scoring Alignments

To decide which alignment is best, we need a **scoring system** with three components:

1. **Match/mismatch scores** -- how much to reward identical positions and penalise mismatches
2. **Substitution matrices** -- position-specific scores for amino acid pairs
3. **Gap penalties** -- the cost of introducing gaps (indels)

### 3.1 Simple Match/Mismatch Scoring

For DNA, a simple scheme often suffices:
- Match: +1 (or +2)
- Mismatch: -1
- Gap: -2

This works because there are only 4 nucleotides and transitions vs. transversions can be handled separately if needed.

### 3.2 Substitution Matrices for Proteins

For proteins, not all substitutions are equal. Replacing leucine (L) with isoleucine (I) -- both hydrophobic, branched-chain amino acids -- is far more likely in evolution than replacing leucine with aspartate (D), a negatively charged residue.

**Substitution matrices** encode the log-odds of observed substitution frequencies:

$$S(a, b) = \log_2 \frac{q_{ab}}{p_a \cdot p_b}$$

where $q_{ab}$ is the observed frequency of substitution and $p_a, p_b$ are background amino acid frequencies.

- **Positive score**: substitution observed more often than expected (conservative)
- **Negative score**: substitution observed less often than expected (disruptive)
- **Zero**: observed at expected frequency

### 3.2 Amino Acid Substitution Matrices

For protein sequences, a simple +1/-1 match/mismatch scoring is too crude because:
- Some amino acid substitutions are chemically conservative (e.g., Leu↔Ile, Asp↔Glu) and should be penalized less
- Some are radical (e.g., Trp↔Pro, charged↔hydrophobic) and should be penalized more

**Substitution matrices** encode the biological "cost" of replacing one amino acid with another, derived from observed frequencies in real protein families.

**Key properties of a good substitution matrix:**
- Diagonal values are positive (matching yourself is rewarded)
- Conservative substitutions (similar physicochemical properties) have positive or small negative scores
- Radical substitutions have large negative scores
- The matrix is symmetric: score(A→B) = score(B→A)

**Amino acid property groups** that tend to substitute conservatively:

| Group | Amino acids |
|-------|-------------|
| Small nonpolar | Ala (A), Val (V), Ile (I), Leu (L), Met (M) |
| Aromatic | Phe (F), Tyr (Y), Trp (W) |
| Polar uncharged | Ser (S), Thr (T), Asn (N), Gln (Q) |
| Positively charged | Lys (K), Arg (R), His (H) |
| Negatively charged | Asp (D), Glu (E) |
| Special | Cys (C), Pro (P), Gly (G) |

The two main matrix families — PAM and BLOSUM — differ in how they were derived (see Section 3.3).

### 3.3 PAM vs BLOSUM Matrices

| Feature | PAM | BLOSUM |
|---|---|---|
| **Full name** | Point Accepted Mutation | BLOcks SUbstitution Matrix |
| **Developed by** | Margaret Dayhoff (1978) | Henikoff & Henikoff (1992) |
| **Basis** | Closely related sequences, extrapolated | Conserved blocks across diverse proteins |
| **Numbering** | Higher number = more distant (PAM250) | Higher number = more similar (BLOSUM80) |
| **For close sequences** | PAM40 | BLOSUM80 |
| **For distant sequences** | PAM250 | BLOSUM45 |
| **Most commonly used** | PAM250 | **BLOSUM62** (default in BLAST) |

**Key insight**: PAM1 represents 1% divergence and is extrapolated by matrix exponentiation to higher evolutionary distances. BLOSUM matrices are directly computed from blocks of conserved sequences at the specified percent identity threshold.

```python
from Bio.Align import substitution_matrices

# Load BLOSUM62 -- the most widely used substitution matrix
blosum62 = substitution_matrices.load("BLOSUM62")

# Display a focused subset: hydrophobic amino acids
aa_subset = "AILMFYWV"
print("BLOSUM62 -- Hydrophobic amino acids:")
print("     ", "  ".join(f"{aa:>3}" for aa in aa_subset))
for aa1 in aa_subset:
    row = [f"{blosum62[aa1, aa2]:>3}" for aa2 in aa_subset]
    print(f"  {aa1}  ", "  ".join(row))
```

```python
# Compare a few informative pairs
pairs = [
    ("L", "I", "Both hydrophobic, branched-chain"),
    ("F", "Y", "Both aromatic"),
    ("D", "E", "Both negatively charged"),
    ("K", "R", "Both positively charged"),
    ("S", "T", "Both small, hydroxyl-containing"),
    ("L", "D", "Hydrophobic vs charged -- disruptive"),
    ("W", "G", "Large aromatic vs smallest amino acid"),
]

print(f"{'Pair':<8} {'Score':>6}  Reason")
print("-" * 55)
for a, b, reason in pairs:
    score = blosum62[a, b]
    print(f"{a} <-> {b}  {score:>+4}    {reason}")
```

```python
# Visualize the full BLOSUM62 matrix as a heatmap
standard_aa = "ARNDCQEGHILKMFPSTWYV"
matrix_data = np.zeros((20, 20))
for i, a in enumerate(standard_aa):
    for j, b in enumerate(standard_aa):
        matrix_data[i, j] = blosum62[a, b]

fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(matrix_data, cmap="RdBu_r", vmin=-4, vmax=11)
ax.set_xticks(range(20))
ax.set_yticks(range(20))
ax.set_xticklabels(list(standard_aa))
ax.set_yticklabels(list(standard_aa))
ax.set_title("BLOSUM62 Substitution Matrix")
plt.colorbar(im, label="Score")
plt.tight_layout()
plt.show()
```

**Observations from the heatmap:**
- The diagonal (self-substitutions) is always positive: identical residues are rewarded
- Tryptophan (W) has the highest self-score (11) -- it is the rarest amino acid, so seeing it conserved is very informative
- Clusters of positive off-diagonal scores reflect chemical similarity groups (e.g., I/L/M/V, D/E, K/R)

### 3.4 Gap Penalties

Gaps represent insertions or deletions (indels). Two common penalty schemes:

**Linear gap penalty:**
$$\text{gap\_cost}(k) = -d \times k$$

where $d$ is the penalty per gap position and $k$ is the gap length. Every gap position costs the same.

**Affine gap penalty:**
$$\text{gap\_cost}(k) = -d - (k - 1) \times e$$

where $d$ is the **gap opening** penalty and $e$ is the **gap extension** penalty ($e < d$).

| Scheme | Advantages | Typical values |
|---|---|---|
| Linear | Simpler, fewer parameters | $d = 8$ |
| Affine | More biologically realistic | $d = 10$, $e = 0.5$ |

**Why affine?** In biology, a single mutational event can insert or delete several residues at once (e.g., replication slippage). A single long gap is biologically more plausible than many scattered single-residue gaps. Affine penalties encourage this by making gap *opening* expensive but gap *extension* cheap.

```python
# Compare linear vs affine gap costs
gap_lengths = np.arange(1, 16)

d_linear = 8
linear_cost = d_linear * gap_lengths

d_open, e_extend = 10, 0.5
affine_cost = d_open + (gap_lengths - 1) * e_extend

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(gap_lengths, linear_cost, "o-", label=f"Linear (d={d_linear})")
ax.plot(gap_lengths, affine_cost, "s-", label=f"Affine (d={d_open}, e={e_extend})")
ax.set_xlabel("Gap length (residues)")
ax.set_ylabel("Total penalty")
ax.set_title("Linear vs Affine Gap Penalty")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print("With affine penalties, short gaps are penalized MORE than linear,")
print("but long gaps are penalized LESS -- favouring fewer, longer gaps.")
```

---

## 4. Global Alignment: The Needleman-Wunsch Algorithm

Published in 1970 by Saul Needleman and Christian Wunsch, this was one of the first applications of **dynamic programming** to biology.

**Global alignment** aligns two sequences end-to-end in their entirety. It is appropriate when the sequences are of similar length and suspected to be homologous over their full length.

### Algorithm Overview

Given sequences $X$ of length $m$ and $Y$ of length $n$:

1. **Create** a score matrix $F$ of size $(n+1) \times (m+1)$ and a traceback matrix $T$
2. **Initialize** the first row and column with cumulative gap penalties
3. **Fill** each cell using the recurrence relation
4. **Traceback** from the bottom-right corner to recover the alignment

### Recurrence Relation

$$F(i, j) = \max \begin{cases} F(i-1, j-1) + s(x_j, y_i) & \text{(match/mismatch)} \\ F(i-1, j) - d & \text{(gap in sequence X)} \\ F(i, j-1) - d & \text{(gap in sequence Y)} \end{cases}$$

where $s(x_j, y_i)$ is the substitution score and $d$ is the gap penalty.

### Time and Space Complexity

- **Time**: $O(m \times n)$ -- we fill every cell in the matrix
- **Space**: $O(m \times n)$ -- we store the full matrix for traceback
