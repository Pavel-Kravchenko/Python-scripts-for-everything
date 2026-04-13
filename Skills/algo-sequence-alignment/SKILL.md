---
name: algo-sequence-alignment
description: "1. Implement Needleman-Wunsch global alignment 2. Implement Smith-Waterman local alignment 3. Understand scoring matrices and gap penalties"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/04_sequence_alignment.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# 🧬 Sequence Alignment with Dynamic Programming

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/04_sequence_alignment.ipynb`*

# 🧬 Sequence Alignment with Dynamic Programming

## Learning Objectives

1. Implement Needleman-Wunsch global alignment
2. Implement Smith-Waterman local alignment
3. Understand scoring matrices and gap penalties

---

This is the **most important** DP application in bioinformatics!

## Why Sequence Alignment Matters

Living organisms accumulate mutations over time: single nucleotide substitutions, insertions, and deletions all cause related sequences to diverge from a common ancestor. Despite this divergence, evolutionary relationships leave a recognizable signal — two sequences that share a common origin will tend to match far more often than random chance would predict. Alignment is the computational procedure that recovers this signal by arranging the two sequences so that corresponding positions are placed side by side, introducing gap characters where one lineage inserted or deleted residues relative to the other.

This single idea underpins much of modern bioinformatics. BLAST searches a database by efficiently finding high-scoring local alignments between a query and millions of stored sequences. Phylogenetic tree-building methods require a multiple sequence alignment before inferring branching order. AlphaFold and related structure-prediction pipelines use evolutionary co-variation signals extracted from large alignments to constrain 3-D models. In each case the dynamic programming algorithms developed here — Needleman-Wunsch and Smith-Waterman — are either used directly or form the conceptual foundation for faster heuristic variants.

For the practical BioPython implementation of pairwise alignment, see [Tier 2 Module 3](../../Tier_2_Core_Bioinformatics/03_Pairwise_Sequence_Alignment/01_pairwise_sequence_alignment.ipynb).

## Scoring Systems

The simplest scoring scheme awards a fixed bonus for a match and a fixed penalty for a mismatch. This works for DNA where all four bases are treated symmetrically, but proteins require something more nuanced. Not all amino acid substitutions are equally likely: a leucine-to-isoleucine swap is biochemically conservative and occurs frequently, while a tryptophan-to-glycine swap is disruptive and rare.

**Substitution matrices** encode observed substitution frequencies derived from curated alignments:

- **BLOSUM62** (BLOcks SUbstitution Matrix) is built from blocks of ungapped alignments with less than 62% identity. Each entry is a log-odds score reflecting how often that pair appears in real alignments versus random chance.
- **PAM matrices** (Point Accepted Mutations) model a Markov chain of evolutionary change; PAM250 corresponds to roughly 250 accepted point mutations per 100 residues.

For DNA work, a simple match/mismatch matrix is usually sufficient, though **transition/transversion** weighting adds realism — transitions (purine↔purine: A↔G, pyrimidine↔pyrimidine: C↔T) occur roughly twice as often as transversions and so deserve a less severe penalty.

A minimal 4×4 DNA scoring matrix:

|   | A  | T  | G  | C  |
|---|----|----|----|----|
| A | +2 | -1 | -1 | -1 |
| T | -1 | +2 | -1 | -1 |
| G | -1 | -1 | +2 | -1 |
| C | -1 | -1 | -1 | +2 |

```python
import numpy as np

# Simple DNA scoring matrix — symmetric, match=+2, any mismatch=-1
DNA_SCORE = {
    ('A','A'): 2, ('A','T'):-1, ('A','G'):-1, ('A','C'):-1,
    ('T','A'):-1, ('T','T'): 2, ('T','G'):-1, ('T','C'):-1,
    ('G','A'):-1, ('G','T'):-1, ('G','G'): 2, ('G','C'):-1,
    ('C','A'):-1, ('C','T'):-1, ('C','G'):-1, ('C','C'): 2,
}
# Transitions (A<->G, C<->T) could be scored differently for realism
# e.g. DNA_SCORE[('A','G')] = DNA_SCORE[('G','A')] = -0.5

bases = ['A', 'T', 'G', 'C']
print("DNA scoring matrix:")
print(f"{'':4}", end="")
for b in bases:
    print(f"{b:>5}", end="")
print()
for r in bases:
    print(f"{r:4}", end="")
    for c in bases:
        print(f"{DNA_SCORE[(r,c)]:>5}", end="")
    print()
```

## DP Table Fill: Step by Step

Both Needleman-Wunsch and Smith-Waterman fill a 2-D matrix of size `(len(seq1)+1) x (len(seq2)+1)`. The algorithm proceeds in three phases:

**1. Initialization**

The first row and first column represent alignments where one sequence is consumed against nothing but gaps.
- For Needleman-Wunsch: `dp[i][0] = i * gap_penalty`, `dp[0][j] = j * gap_penalty`.
- For Smith-Waterman: the entire border stays at 0 (a local alignment can start anywhere).

**2. Cell-by-cell fill**

For each interior cell `dp[i][j]`, we consider three predecessors:

| Move | From | Meaning |
|------|------|---------|
| Diagonal | `dp[i-1][j-1]` | Align `seq1[i-1]` with `seq2[j-1]` (match or mismatch) |
| Up | `dp[i-1][j]` | Consume `seq1[i-1]`, insert gap in seq2 |
| Left | `dp[i][j-1]` | Insert gap in seq1, consume `seq2[j-1]` |

We take the maximum of these three options (and 0 for Smith-Waterman).

**3. Traceback**

Starting from the optimal cell (bottom-right for NW, maximum cell for SW), we follow the path of choices backwards until we reach the origin (or a zero cell for SW), reconstructing the aligned sequences.

```python
import numpy as np

def build_nw_dp(seq1, seq2, match=1, mismatch=-1, gap=-2):
    # Build only the DP matrix (no traceback) for display purposes
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(m + 1):
        dp[i][0] = i * gap
    for j in range(n + 1):
        dp[0][j] = j * gap
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(dp[i-1][j-1] + s, dp[i-1][j] + gap, dp[i][j-1] + gap)
    return dp

def print_dp_table(seq1, seq2, dp):
    # Header row: blank + gap column label + seq2 characters
    col_w = 5
    header = f"{'':>{col_w}}{'':>{col_w}}" + "".join(f"{c:>{col_w}}" for c in seq2)
    print(header)
    # Row 0 is the gap-initialised row (no seq1 character)
    row0 = f"{'':>{col_w}}" + "".join(f"{dp[0][j]:>{col_w}}" for j in range(len(seq2) + 1))
    print(row0)
    # Remaining rows: labelled with seq1 characters
    for i in range(1, len(seq1) + 1):
        row = f"{seq1[i-1]:>{col_w}}" + "".join(f"{dp[i][j]:>{col_w}}" for j in range(len(seq2) + 1))
        print(row)

# Demonstrate on a short pair so the table fits on screen
s1, s2 = "AGTC", "AATC"
dp = build_nw_dp(s1, s2)
print(f"NW DP table for seq1='{s1}', seq2='{s2}' (match=1, mismatch=-1, gap=-2):\n")
print_dp_table(s1, s2, dp)
```

## Traceback Visualization

Once the DP table is filled, the optimal alignment is recovered by following a path of arrows backwards from the goal cell. Each arrow type corresponds to one alignment operation:

- **Diagonal (\\)** — came from `dp[i-1][j-1]`: the two characters at positions `i` and `j` were aligned (match or mismatch).
- **Up (|)** — came from `dp[i-1][j]`: character `seq1[i-1]` was aligned against a gap in seq2 (gap inserted in seq2).
- **Left (--)** — came from `dp[i][j-1]`: a gap was inserted in seq1 and character `seq2[j-1]` was consumed.

When multiple predecessors achieve the same score (ties), any valid arrow can be chosen — different choices produce different optimal alignments with identical scores. The function below reconstructs the arrows and prints them alongside the scores so you can trace the path by eye.

```python
import numpy as np

def visualize_traceback(seq1, seq2, match=1, mismatch=-1, gap=-2):
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(m + 1):
        dp[i][0] = i * gap
    for j in range(n + 1):
        dp[0][j] = j * gap
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(dp[i-1][j-1] + s, dp[i-1][j] + gap, dp[i][j-1] + gap)

    # Determine arrow for each cell: D=diagonal, U=up, L=left
    arrows = [['.' for _ in range(n + 1)] for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1] + s:
                arrows[i][j] = 'D'   # diagonal
            elif dp[i][j] == dp[i-1][j] + gap:
                arrows[i][j] = 'U'   # up
            else:
                arrows[i][j] = 'L'   # left

    # Mark the traceback path
    on_path = [[False] * (n + 1) for _ in range(m + 1)]
    i, j = m, n
    while i > 0 or j > 0:
        on_path[i][j] = True
        a = arrows[i][j]
        if i == 0:
            j -= 1
        elif j == 0:
            i -= 1
        elif a == 'D':
            i -= 1; j -= 1
        elif a == 'U':
            i -= 1
        else:
            j -= 1
    on_path[0][0] = True

    arrow_char = {'D': '\\', 'U': '|', 'L': '-', '.': ' '}
    col_w = 6

    # Header
    print(f"{'':>{col_w}}{'':>{col_w}}", end="")
    for c in seq2:
        print(f"{c:>{col_w}}", end="")
    print()

    for i in range(m + 1):
        label = seq1[i-1] if i > 0 else ' '
        print(f"{label:>{col_w}}", end="")
        for j in range(n + 1):
            cell_val = dp[i][j]
            arr = arrow_char[arrows[i][j]]
            marker = '*' if on_path[i][j] else ' '
            print(f"{marker}{arr}{cell_val:>3} ", end="")
        print()

    print("\nLegend: * = traceback path  D=\\ diagonal  U=| up  L=- left")

visualize_traceback("AGTC", "AATC")
```

## Global vs Local Alignment

Before looking at the implementations, it is worth being precise about when to use each algorithm.

**Global alignment (Needleman-Wunsch)** aligns two sequences end-to-end. Every residue in both sequences must be accounted for — either matched/mismatched against a residue in the other sequence, or paired with a gap. Use global alignment when:
- The two sequences are roughly the same length.
- You expect them to be homologs across their full extent (e.g., comparing two orthologs of a protein from different species).
- You want to measure overall divergence for phylogenetic analysis.

**Local alignment (Smith-Waterman)** finds the highest-scoring sub-region of similarity between two sequences. The algorithm is free to ignore flanking regions that do not contribute positively to the score. Use local alignment when:
- One sequence is much longer than the other (e.g., mapping a 150 bp read against a chromosome).
- You are searching for a conserved functional domain embedded in otherwise diverged proteins.
- The sequences share a common motif but have different architectures overall (e.g., multi-domain proteins that share one domain).

The practical difference comes down to initialisation and the floor on cell values: NW permits negative scores in the matrix and always traces back from the corner, whereas SW floors each cell at zero and traces back from the global maximum.

```python
import numpy as np

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Global sequence alignment.
    Returns: (score, aligned_seq1, aligned_seq2)
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize scoring matrix
    dp = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(m + 1):
        dp[i][0] = i * gap
    for j in range(n + 1):
        dp[0][j] = j * gap
    
    # Fill matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(
                dp[i-1][j-1] + score,  # Match/mismatch
                dp[i-1][j] + gap,       # Gap in seq2
                dp[i][j-1] + gap        # Gap in seq1
            )
    
    # Traceback
    align1, align2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1] + score:
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1
    
    return dp[m][n], ''.join(reversed(align1)), ''.join(reversed(align2))

seq1 = "ATCGTAC"
seq2 = "ATGTTAT"
score, a1, a2 = needleman_wunsch(seq1, seq2)

print(f"Global Alignment (Needleman-Wunsch)")
print(f"Score: {score}")
print(f"Seq1: {a1}")
print(f"Seq2: {a2}")
```
