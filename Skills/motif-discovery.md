---
name: motif-discovery
description: Quantitative motif analysis — PWM/PPM construction, information content, KDIC scoring, Fisher enrichment, TomTom matching, BED/FASTA pipeline design
---

## When to Use

Use this skill when:
- Building Position Weight Matrices (PWMs) from aligned binding sites
- Computing information content and KDIC scores for motif quality assessment
- Generating IUPAC consensus sequences from PPMs
- Scoring genomic sequences against a PWM
- Performing motif enrichment tests (Fisher's exact + Benjamini-Hochberg)
- Matching motifs to JASPAR / HOCOMOCO databases (TomTom concept)

## Quick Reference

| Concept | Formula | Notes |
|---|---|---|
| PPM (Position Probability Matrix) | PPM[i,k] = (count[i,k] + pseudocount) / (N + 4·pseudocount) | Shape: L×4 (length × ACGT) |
| PWM (log-odds) | PWM[i,k] = log2(PPM[i,k] / bg[k]) | bg = background (0.25 each) |
| IC per position | IC[i] = 2 + Σ_k PPM[i,k]·log2(PPM[i,k]) | bits; max = 2 |
| Total IC | sum(IC) | bits; higher = more specific |
| KDIC | mean(IC) / 2 | normalized to [0,1] |
| IUPAC code | degenerate base for positions with IC < 1 | R=AG, Y=CT, N=any, etc. |
| Sequence score | Σ_i PWM[i, seq[i]] | higher = better match |
| Fisher enrichment | Fisher's exact test on 2×2 table: (hits in peaks, hits in bg) | p-value + FDR |

## Key Patterns

**Pattern 1: Build PPM/PWM from aligned sequences**
```python
import numpy as np

ALPHABET = "ACGT"
BASE_TO_IDX = {b: i for i, b in enumerate(ALPHABET)}

def make_ppm(seqs, pseudocount=0.1):
    """seqs: list of equal-length strings (ACGT). Returns L×4 PPM."""
    L = len(seqs[0])
    counts = np.zeros((L, 4))
    for s in seqs:
        for i, b in enumerate(s.upper()):
            if b in BASE_TO_IDX:
                counts[i, BASE_TO_IDX[b]] += 1
    return (counts + pseudocount) / (len(seqs) + 4 * pseudocount)

def make_pwm(ppm, bg=0.25):
    return np.log2(ppm / bg)
```

**Pattern 2: Information content and KDIC**
```python
def information_content(ppm):
    """Returns per-position IC array (bits) and total IC."""
    ic = 2 + np.sum(ppm * np.log2(ppm + 1e-9), axis=1)
    return ic, ic.sum()

def kdic(ppm):
    ic, _ = information_content(ppm)
    return ic.mean() / 2   # normalized to [0, 1]
```

**Pattern 3: Score a sequence against a PWM**
```python
def score_sequence(seq, pwm):
    """Returns PWM log-odds score for a sequence string."""
    score = 0.0
    for i, b in enumerate(seq.upper()):
        idx = BASE_TO_IDX.get(b, -1)
        if idx >= 0:
            score += pwm[i, idx]
    return score
```

**Pattern 4: Score distribution (exact for short motifs)**
```python
def exact_score_distribution(pwm, bg=0.25):
    """Enumerate all L-mers and return (scores, probabilities). Only for L ≤ 10."""
    L = pwm.shape[0]
    assert L <= 10, "Too slow for L > 10; use Monte Carlo"
    from itertools import product
    scores, probs = [], []
    for bases in product(range(4), repeat=L):
        s = sum(pwm[i, b] for i, b in enumerate(bases))
        p = np.prod([bg ** 1] * L)  # uniform bg
        scores.append(s); probs.append(p)
    return np.array(scores), np.array(probs)
```

**Pattern 5: Fisher enrichment with BH correction**
```python
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def motif_enrichment(hits_peak, total_peak, hits_bg, total_bg):
    table = [[hits_peak, total_peak - hits_peak],
             [hits_bg,   total_bg  - hits_bg]]
    _, p = fisher_exact(table, alternative="greater")
    return p

# BH correction for multiple motifs
pvals = [motif_enrichment(...) for motif in motifs]
_, qvals, _, _ = multipletests(pvals, method="fdr_bh")
```

## Code Templates

**Template 1: Full motif analysis pipeline**
```python
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

ALPHABET = "ACGT"
BASE_TO_IDX = {b: i for i, b in enumerate(ALPHABET)}
IUPAC = {
    frozenset("A"): "A", frozenset("C"): "C",
    frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y",
    frozenset("GC"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D",
    frozenset("ACT"): "H", frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}

def consensus(ppm, ic_threshold=1.0):
    cons = []
    for pos in ppm:
        dominant = frozenset(ALPHABET[i] for i, p in enumerate(pos) if p == pos.max())
        ic = 2 + sum(p * np.log2(p + 1e-9) for p in pos)
        if ic < ic_threshold:
            cons.append("N")
        else:
            cons.append(IUPAC.get(dominant, "N"))
    return "".join(cons)
```

**Template 2: Sequence logo data (heights)**
```python
def logo_heights(ppm):
    """Return (L, 4) matrix of letter heights for sequence logo."""
    ic, _ = information_content(ppm)
    return ppm * ic[:, None]  # scale each column by its IC
```

## Common Pitfalls

- **PPM shape confusion:** always use L×4 (rows=positions, cols=ACGT); transpose errors cause silent wrong scores
- **Pseudocount required:** zero counts produce -inf in log2; use pseudocount ≥ 0.01
- **IC < 0 at degenerate positions:** possible with pseudocounts when N is small; clip to 0
- **Monte Carlo convergence:** always report sample size and CI; n=100k is minimum for reliable p-values at p < 0.001
- **IUPAC N vs gap:** N = any base; `-` = alignment gap; keep them separate
- **Fisher vs chi-squared:** use Fisher's exact when expected cell counts < 5; chi-squared otherwise

## Related Skills

- `structural-bioinformatics` — PWM/PROSITE motifs in protein context
- `hic-analysis` — pileup analysis uses motif sites (CTCF) as feature set
- `ngs-variant-calling` — ChIP-seq peaks are the typical input for motif enrichment
