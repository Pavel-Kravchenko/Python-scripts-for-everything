---
name: algo-sequence-alignment
description: "Needleman-Wunsch global and Smith-Waterman local sequence alignment — DP table fill, traceback, scoring matrices"
tool_type: python
primary_tool: NumPy
---

# Sequence Alignment with Dynamic Programming

## Global vs Local

| | Needleman-Wunsch (Global) | Smith-Waterman (Local) |
|---|---|---|
| Aligns | Full sequences end-to-end | Best matching sub-region |
| Init | First row/col = cumulative gap penalty | First row/col = 0 |
| Cell floor | No floor (can go negative) | Floor at 0 |
| Traceback from | Bottom-right corner | Maximum cell in matrix |
| Use when | Same-length homologs, full divergence | Domain search, read mapping |

## DP Table Fill

For each cell dp[i][j], take max of three predecessors:

| Move | Source | Meaning |
|---|---|---|
| Diagonal | dp[i-1][j-1] + score(seq1[i-1], seq2[j-1]) | Match/mismatch |
| Up | dp[i-1][j] + gap | Gap in seq2 |
| Left | dp[i][j-1] + gap | Gap in seq1 |

## Scoring

Simple DNA matrix: match=+2, mismatch=-1. For proteins, use substitution matrices (BLOSUM62, PAM250) that reflect observed evolutionary substitution frequencies.

```python
DNA_SCORE = {(a,b): 2 if a==b else -1 for a in 'ATGC' for b in 'ATGC'}
```

## Needleman-Wunsch Implementation

```python
import numpy as np

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """Global alignment. Returns (score, aligned_seq1, aligned_seq2)."""
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(m + 1): dp[i][0] = i * gap
    for j in range(n + 1): dp[0][j] = j * gap

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(dp[i-1][j-1] + s, dp[i-1][j] + gap, dp[i][j-1] + gap)

    # Traceback
    align1, align2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1] + s:
                align1.append(seq1[i-1]); align2.append(seq2[j-1])
                i -= 1; j -= 1; continue
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            align1.append(seq1[i-1]); align2.append('-'); i -= 1
        else:
            align1.append('-'); align2.append(seq2[j-1]); j -= 1

    return dp[m][n], ''.join(reversed(align1)), ''.join(reversed(align2))
```

## Smith-Waterman (Local) — Key Differences

```python
# Initialization: all zeros (no cumulative gap penalty)
# Cell computation: max(diagonal, up, left, 0)  <-- floor at 0
# Traceback: start from global maximum cell, stop when hitting 0
```

## Complexity

O(mn) time and space for both algorithms. Space can be reduced to O(min(m,n)) with rolling arrays (loses traceback).

## Pitfalls

- Traceback ties (multiple predecessors with same score) produce different equally-optimal alignments
- Affine gap penalties (open + extend) are more realistic than linear gaps but require 3 DP matrices
- NW forces alignment of entire sequences — if one is much longer, use SW or semi-global alignment
- BLOSUM62 is default for protein BLAST; PAM250 is for more distant homologs
