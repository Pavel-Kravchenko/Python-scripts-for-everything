---
name: bio-core-motif-discovery
description: "*Prerequisites: Modules 1–14. Module 10 (Sequence Motifs and Domains) for PWM background.*"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/15_Motif_Discovery/15_motif_discovery.ipynb"
---

# Motif Discovery

*Source: Course notebook `Tier_2_Core_Bioinformatics/15_Motif_Discovery/15_motif_discovery.ipynb`*

# Motif Discovery

**Tier 2 — Core Bioinformatics | Module 15**

*Prerequisites: Modules 1–14. Module 10 (Sequence Motifs and Domains) for PWM background.*

---

This module covers quantitative motif analysis: building Position Weight Matrices from known sites, scoring and scanning sequences, computing information content, performing statistical enrichment testing, and the principles of de novo motif discovery.

**By the end of this module you will be able to:**

1. Build a PPM/PWM from a set of aligned binding sites
2. Calculate per-position information content and KDIC
3. Score sequences against a PWM and select a significance threshold
4. Test motif enrichment between foreground and background sequences using Fisher's exact test
5. Explain how MEME's Expectation-Maximization discovers motifs de novo
6. Interpret TomTom output for matching discovered motifs to known TFs

## Why this notebook matters

Motif discovery is the bridge between functional genomics experiments and mechanistic understanding. When you run ChIP-seq for a transcription factor and get 10,000 peaks, the first question is: *what sequence does this factor actually recognize?* De novo motif discovery answers this. When you find a novel regulatory element, motif enrichment analysis tells you which TFs might bind it. Understanding how PWMs work, how MEME's EM algorithm discovers motifs, and how to calibrate score thresholds is essential for any regulatory genomics analysis.

## Complicated moments explained

- **PFM → PPM → PWM pipeline**: PFM = raw counts; PPM = frequencies (add pseudocount before dividing); PWM = log-odds = log2(PPM_ij / background_i). The PWM is the only version suitable for scoring — using raw probabilities directly would make very rare bases dominate logarithmically.
- **Pseudocount choice matters**: Too small a pseudocount and zero counts cause -∞ scores; too large and it drowns out the real signal. For a dataset with N sequences, a common rule is pseudocount = sqrt(N)/N or simply 0.1-0.5. The goal is to regularize without distorting well-sampled positions.
- **Background model**: The PWM log-odds is relative to a background model. Using equal (0.25) background is simplest; using a 2nd-order Markov model (dinucleotide frequencies matching the input sequences' composition) is more powerful for GC-rich genomes.
- **MEME EM algorithm**: MEME treats binding site positions as hidden variables and iterates between: E-step (estimate which window in each sequence is the best site, given current motif model) and M-step (update the motif model using weighted contributions from all windows). It converges to a local optimum — run MEME multiple times or use STREME for large datasets.
- **Threshold selection**: There is no universal threshold. Options: (1) p-value from score distribution (compute for all 4^L sequences); (2) fraction of maximum score (e.g., 80%); (3) empirical calibration from ChIP-seq data. Always report the threshold used and the resulting false-positive rate.

## Background: From Binding Sites to Matrices

Transcription factors recognize short degenerate DNA sequences (~6–20 bp). The same TF can bind many sequences because some positions tolerate multiple bases.

**Representation levels:**

| Matrix | Content | Use |
|---|---|---|
| PFM (Position Frequency Matrix) | Raw counts per base per position | Input data |
| PPM (Position Probability Matrix) | Frequencies (+ pseudocounts) | Normalized |
| PWM (Position Weight Matrix) | Log-odds vs background | Scoring |

**Information content (IC)** measures how much a position deviates from random:
- IC = 0 bits → completely degenerate (any base equally likely)
- IC = 2 bits → fully conserved (only one base allowed)
- IC per position: IC[i] = 2 + Σ_k p[i,k] · log₂(p[i,k])

**KDIC (Kullback-Discrimination Information Criterion):** mean IC across all positions, normalized to [0, 1] by dividing by 2. Useful for ranking motifs by specificity.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from itertools import product
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")

ALPHABET = "ACGT"
BASE_TO_IDX = {b: i for i, b in enumerate(ALPHABET)}
COLORS = {"A": "#2ecc71", "C": "#3498db", "G": "#f39c12", "T": "#e74c3c"}

plt.rcParams.update({"figure.dpi": 120, "axes.spines.top": False, "axes.spines.right": False})
print("Libraries loaded.")
```

## 1. Building a PPM and PWM

We start with a set of aligned sequences — typically ChIP-seq peak summits or SELEX-derived binding sites. We will use a synthetic set of CTCF-like sequences (CCGCGNGGNGGCAG consensus) for demonstration.

```python
# Synthetic CTCF-like binding sites (based on public JASPAR MA0139.1 consensus)
rng = np.random.default_rng(42)

def sample_ctcf_like(n=200):
    """Generate synthetic CTCF-like 19-mer sequences."""
    consensus_weights = [
        [0.1, 0.7, 0.1, 0.1],  # pos 0: C dominant
        [0.1, 0.7, 0.1, 0.1],  # pos 1: C dominant
        [0.1, 0.1, 0.7, 0.1],  # pos 2: G
        [0.1, 0.7, 0.1, 0.1],  # pos 3: C
        [0.1, 0.1, 0.7, 0.1],  # pos 4: G
        [0.25]*4,               # pos 5: N
        [0.1, 0.1, 0.7, 0.1],  # pos 6: G
        [0.1, 0.1, 0.7, 0.1],  # pos 7: G
        [0.25]*4,               # pos 8: N
        [0.1, 0.1, 0.7, 0.1],  # pos 9: G
        [0.1, 0.1, 0.7, 0.1],  # pos 10: G
        [0.1, 0.7, 0.1, 0.1],  # pos 11: C
        [1.0, 0.0, 0.0, 0.0],  # pos 12: A
        [0.1, 0.1, 0.7, 0.1],  # pos 13: G
    ]
    seqs = []
    for _ in range(n):
        seq = "".join(
            rng.choice(list(ALPHABET), p=w)
            for w in consensus_weights
        )
        seqs.append(seq)
    return seqs

seqs = sample_ctcf_like(200)
print(f"Generated {len(seqs)} sequences, length {len(seqs[0])}")
print("Sample:", seqs[:3])

def make_ppm(seqs, pseudocount=0.1):
    L = len(seqs[0])
    counts = np.zeros((L, 4))
    for s in seqs:
        for i, b in enumerate(s.upper()):
            if b in BASE_TO_IDX:
                counts[i, BASE_TO_IDX[b]] += 1
    ppm = (counts + pseudocount) / (len(seqs) + 4 * pseudocount)
    return ppm

def make_pwm(ppm, bg=0.25):
    return np.log2(ppm / bg)

ppm = make_ppm(seqs)
pwm = make_pwm(ppm)

print(f"\nPPM shape: {ppm.shape}  (L\u00d74: {ppm.shape[0]} positions \u00d7 4 bases)")
print(f"PPM row sums (should \u2248 1): {ppm.sum(axis=1).round(4)}")
print(f"\nPWM (first 5 positions):\n{pd.DataFrame(pwm[:5], columns=list(ALPHABET)).round(3)}")
```

## 2. Information Content and KDIC

IC per position quantifies how much that position constrains base identity. A fully conserved position (IC = 2 bits) contributes maximally; a degenerate position (IC ≈ 0) contributes nothing to discrimination.

**KDIC** summarizes the whole motif with one number: mean IC / 2, ranging from 0 (no information) to 1 (fully conserved). This is useful for comparing motifs of different lengths.

```python
def information_content(ppm):
    ic = 2 + np.sum(ppm * np.log2(np.clip(ppm, 1e-9, 1)), axis=1)
    ic = np.clip(ic, 0, 2)
    return ic

def kdic(ppm):
    return information_content(ppm).mean() / 2

ic = information_content(ppm)
k = kdic(ppm)

print(f"IC per position (bits): {ic.round(3)}")
print(f"Total IC: {ic.sum():.2f} bits")
print(f"KDIC: {k:.4f}  (higher = more specific motif)")

fig, axes = plt.subplots(1, 2, figsize=(14, 4))

# IC bar plot
axes[0].bar(range(len(ic)), ic, color=[
    "#e74c3c" if v > 1.5 else "#f39c12" if v > 0.8 else "#95a5a6"
    for v in ic
])
axes[0].axhline(1.0, color="gray", ls="--", lw=0.8, label="IC = 1 bit")
axes[0].set_xlabel("Position"); axes[0].set_ylabel("IC (bits)")
axes[0].set_title(f"Information content per position\nTotal IC = {ic.sum():.1f} bits | KDIC = {k:.3f}")
axes[0].legend()

# Sequence logo (letter heights)
logo_h = ppm * ic[:, None]  # L\u00d74 heights
x = np.arange(len(ic))
bottom = np.zeros(len(ic))
for j, base in enumerate(ALPHABET):
    axes[1].bar(x, logo_h[:, j], bottom=bottom, color=COLORS[base], label=base, width=0.8)
    bottom += logo_h[:, j]
axes[1].set_xlabel("Position"); axes[1].set_ylabel("Height (bits)")
axes[1].set_title("Sequence logo (information-content heights)")
axes[1].legend(frameon=False, ncol=4)

plt.tight_layout(); plt.show()
```

## 3. IUPAC Consensus Sequence

The IUPAC alphabet has degenerate codes for positions with multiple likely bases. We assign a code based on which bases exceed a frequency threshold, and mark low-IC positions as `N`.

```python
IUPAC_MAP = {
    frozenset("A"): "A", frozenset("C"): "C",
    frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y",
    frozenset("GC"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D",
    frozenset("ACT"): "H", frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}

def iupac_consensus(ppm, ic_threshold=1.0, freq_threshold=0.25):
    ic = information_content(ppm)
    cons = []
    for pos_idx, (row, ic_val) in enumerate(zip(ppm, ic)):
        if ic_val < ic_threshold:
            cons.append("N")
        else:
            dominant = frozenset(ALPHABET[j] for j, p in enumerate(row) if p >= freq_threshold)
            cons.append(IUPAC_MAP.get(dominant, "N"))
    return "".join(cons)

consensus_seq = iupac_consensus(ppm)
print(f"IUPAC consensus: {consensus_seq}")
print(f"Length: {len(consensus_seq)}")

# Compare with max-IC base at each position
max_base = "".join(ALPHABET[np.argmax(ppm[i])] for i in range(len(ppm)))
print(f"Max-probability: {max_base}")
```

## 4. Scoring Sequences Against a PWM

Given a PWM and a query sequence, the score is the sum of log-odds values at each position. Higher scores = better match.

For genome-wide scanning: slide the motif along every position and collect scores above a threshold. The threshold is typically set at the percentile of the score distribution.

```python
def score_sequence(seq, pwm):
    score = 0.0
    for i, b in enumerate(seq.upper()):
        j = BASE_TO_IDX.get(b, -1)
        if j >= 0 and i < len(pwm):
            score += pwm[i, j]
    return score

def scan_sequence(genome_seq, pwm, rc=True):
    """Scan a sequence and its reverse complement. Returns list of (pos, score, strand)."""
    L = len(pwm)
    hits = []
    rc_map = str.maketrans("ACGT", "TGCA")

    def _scan(seq, strand):
        for i in range(len(seq) - L + 1):
            s = score_sequence(seq[i:i+L], pwm)
            hits.append((i, s, strand))

    _scan(genome_seq, "+")
    if rc:
        _scan(genome_seq[::-1].translate(rc_map), "-")
    return hits

# Test
test_seqs = seqs[:10] + ["AAAAAAAAAAAAA" + "A" * (len(seqs[0]) - 13)]  # motif + random
scores = [score_sequence(s, pwm) for s in test_seqs]
print("Motif sequence scores:", [f"{s:.2f}" for s in scores[:5]])
print("Random sequence score:", f"{scores[-1]:.2f}")

# Score distribution visualization
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist([score_sequence(s, pwm) for s in seqs], bins=30, color="steelblue", alpha=0.8, label="Motif seqs")
rnd = ["".join(rng.choice(list(ALPHABET), size=len(seqs[0]))) for _ in range(500)]
ax.hist([score_sequence(s, pwm) for s in rnd], bins=30, color="salmon", alpha=0.8, label="Random seqs")
ax.set_xlabel("PWM score"); ax.set_ylabel("Count")
ax.set_title("Score distributions: motif sites vs random sequences")
ax.legend(); plt.tight_layout(); plt.show()
```

## 5. Score Distribution and Threshold Selection

To call a hit, we need a significance threshold. There are two approaches:

1. **Exact enumeration** (motif length ≤ 10): compute scores for all 4^L possible sequences, weighted by background probability
2. **Monte Carlo** (motif length > 10): sample random sequences under the background model and compute empirical p-values

The threshold is usually the score at a p-value of 10⁻⁴ to 10⁻⁶ relative to the background distribution.

```python
L = len(pwm)
if L <= 10:
    print(f"Motif length = {L} \u2264 10 \u2192 exact enumeration")
    all_scores = []
    for bases in product(range(4), repeat=L):
        s = sum(pwm[i, b] for i, b in enumerate(bases))
        all_scores.append(s)
    all_scores = np.array(all_scores)
    threshold_1e4 = np.percentile(all_scores, 99.99)  # top 0.01%
    print(f"Threshold (p < 1e-4): {threshold_1e4:.2f}")
else:
    print(f"Motif length = {L} > 10 \u2192 Monte Carlo (n=50,000)")
    mc_scores = np.array([
        score_sequence("".join(rng.choice(list(ALPHABET), size=L)), pwm)
        for _ in range(50_000)
    ])
    threshold_1e4 = np.percentile(mc_scores, 99.99)
    print(f"Threshold (p < 1e-4): {threshold_1e4:.2f}")
    all_scores = mc_scores

fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(all_scores, bins=60, color="steelblue", alpha=0.8)
ax.axvline(threshold_1e4, color="red", lw=2, label=f"p < 1e-4 threshold ({threshold_1e4:.1f})")
ax.set_xlabel("PWM score"); ax.set_ylabel("Count")
ax.set_title(f"Background score distribution (L={L})")
ax.legend(); plt.tight_layout(); plt.show()
```
