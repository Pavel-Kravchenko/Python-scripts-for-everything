---
name: bio-core-motif-discovery
description: PWM construction, scoring, threshold selection, and motif scanning for transcription factor binding sites
tool_type: python
primary_tool: Matplotlib
---

# Motif Discovery

## Pitfalls

- **PFM -> PPM -> PWM pipeline**: PFM = raw counts; PPM = frequencies (add pseudocount first); PWM = log2(PPM/background). Only PWM is suitable for scoring.
- **Pseudocount choice**: Too small -> log(0) = -inf; too large -> drowns signal. Common: 0.1-0.5 or sqrt(N)/N.
- **Background model**: Equal (0.25) is simplest; dinucleotide Markov model is better for GC-rich genomes.
- **MEME EM**: Converges to local optimum. Run multiple times or use STREME for large datasets.
- **Threshold selection**: No universal threshold. Options: p-value from score distribution, fraction of max score (e.g. 80%), or empirical calibration from ChIP-seq. Always report threshold and FPR.
- **Coordinate systems**: BED = 0-based half-open; VCF/GFF = 1-based inclusive.

## Matrix Representations

| Matrix | Content | Use |
|---|---|---|
| PFM | Raw counts per base per position | Input data |
| PPM | Frequencies (+ pseudocounts) | Normalized |
| PWM | Log-odds vs background | Scoring |

**Information content**: IC[i] = 2 + sum(p[i,k] * log2(p[i,k])). IC=0 bits -> degenerate; IC=2 bits -> fully conserved.

**KDIC**: mean IC / 2, normalized to [0,1]. Useful for ranking motifs by specificity.

## Core Functions

```python
import numpy as np

ALPHABET = "ACGT"
BASE_TO_IDX = {b: i for i, b in enumerate(ALPHABET)}

def make_ppm(seqs, pseudocount=0.1):
    L = len(seqs[0])
    counts = np.zeros((L, 4))
    for s in seqs:
        for i, b in enumerate(s.upper()):
            if b in BASE_TO_IDX:
                counts[i, BASE_TO_IDX[b]] += 1
    return (counts + pseudocount) / (len(seqs) + 4 * pseudocount)

def make_pwm(ppm, bg=0.25):
    return np.log2(ppm / bg)

def information_content(ppm):
    ic = 2 + np.sum(ppm * np.log2(np.clip(ppm, 1e-9, 1)), axis=1)
    return np.clip(ic, 0, 2)

def kdic(ppm):
    return information_content(ppm).mean() / 2
```

## IUPAC Consensus

```python
IUPAC_MAP = {
    frozenset("A"): "A", frozenset("C"): "C", frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M", frozenset("CGT"): "B", frozenset("AGT"): "D",
    frozenset("ACT"): "H", frozenset("ACG"): "V", frozenset("ACGT"): "N",
}

def iupac_consensus(ppm, ic_threshold=1.0, freq_threshold=0.25):
    ic = information_content(ppm)
    cons = []
    for row, ic_val in zip(ppm, ic):
        if ic_val < ic_threshold:
            cons.append("N")
        else:
            dominant = frozenset(ALPHABET[j] for j, p in enumerate(row) if p >= freq_threshold)
            cons.append(IUPAC_MAP.get(dominant, "N"))
    return "".join(cons)
```

## Scoring and Scanning

```python
def score_sequence(seq, pwm):
    return sum(pwm[i, BASE_TO_IDX[b]] for i, b in enumerate(seq.upper())
               if b in BASE_TO_IDX and i < len(pwm))

def scan_sequence(genome_seq, pwm, rc=True):
    """Scan sequence and reverse complement. Returns [(pos, score, strand)]."""
    L = len(pwm)
    rc_map = str.maketrans("ACGT", "TGCA")
    hits = []
    def _scan(seq, strand):
        for i in range(len(seq) - L + 1):
            hits.append((i, score_sequence(seq[i:i+L], pwm), strand))
    _scan(genome_seq, "+")
    if rc:
        _scan(genome_seq[::-1].translate(rc_map), "-")
    return hits
```

## Threshold Selection

- **Motif length <= 10**: exact enumeration over all 4^L sequences
- **Motif length > 10**: Monte Carlo (sample 50K+ random sequences)

```python
from itertools import product

L = len(pwm)
if L <= 10:
    all_scores = np.array([sum(pwm[i, b] for i, b in enumerate(bases))
                           for bases in product(range(4), repeat=L)])
    threshold = np.percentile(all_scores, 99.99)
else:
    rng = np.random.default_rng(42)
    mc_scores = np.array([score_sequence("".join(rng.choice(list(ALPHABET), size=L)), pwm)
                          for _ in range(50_000)])
    threshold = np.percentile(mc_scores, 99.99)
```
