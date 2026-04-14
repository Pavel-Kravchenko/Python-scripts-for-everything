---
name: bio-core-sequence-motifs
description: PWM/PFM/PPM construction, scoring, scanning, information content, sequence logos, and PROSITE pattern conversion.
tool_type: python
primary_tool: NumPy
---

# Sequence Motifs

## PFM → PPM → PWM Pipeline

```python
import numpy as np

BASES = ["A", "C", "G", "T"]

def build_pfm(sequences: list[str]) -> np.ndarray:
    """Raw count matrix (4 x L)."""
    L = len(sequences[0])
    pfm = np.zeros((4, L), dtype=int)
    for seq in sequences:
        for pos, base in enumerate(seq.upper()):
            if base in BASES:
                pfm[BASES.index(base), pos] += 1
    return pfm

def pfm_to_ppm(pfm: np.ndarray, pseudocount: float = 0.1) -> np.ndarray:
    """Divide by sequence count; add pseudocount to avoid log(0)."""
    n_seq = pfm.sum(axis=0)[0]
    return (pfm + pseudocount) / (n_seq + 4 * pseudocount)

def ppm_to_pwm(ppm: np.ndarray, background: float = 0.25) -> np.ndarray:
    """Log-odds scores relative to uniform background."""
    return np.log2(ppm / background)
```

## Scoring and Scanning

```python
def score_sequence(pwm: np.ndarray, sequence: str) -> float:
    """Sum log-odds scores at each position."""
    return sum(pwm[BASES.index(b), i] for i, b in enumerate(sequence.upper()) if b in BASES)

def scan_sequence(pwm: np.ndarray, sequence: str, threshold: float | None = None) -> list[tuple]:
    """Return (pos, subseq, score) hits above threshold, sorted by score desc."""
    motif_len = pwm.shape[1]
    max_score = np.sum(np.max(pwm, axis=0))
    if threshold is None:
        threshold = 0.6 * max_score  # 60% of max is a common default
    hits = []
    for i in range(len(sequence) - motif_len + 1):
        subseq = sequence[i:i + motif_len]
        s = score_sequence(pwm, subseq)
        if s >= threshold:
            hits.append((i, subseq, s))
    return sorted(hits, key=lambda x: -x[2])

def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ATGCatgc", "TACGtacg")
    return seq.translate(comp)[::-1]

def scan_both_strands(pwm, sequence, threshold=None):
    """Scan forward and reverse strand; convert rev-hit positions to fwd coords."""
    motif_len = pwm.shape[1]
    fwd = [(pos, seq, score, "+") for pos, seq, score in scan_sequence(pwm, sequence, threshold)]
    rev = scan_sequence(pwm, reverse_complement(sequence), threshold)
    rev_fwd = [(len(sequence) - pos - motif_len, seq, score, "-") for pos, seq, score in rev]
    return sorted(fwd + rev_fwd, key=lambda x: -x[2])
```

## Information Content

```python
def information_content(ppm: np.ndarray) -> np.ndarray:
    """IC per position (bits). Max = 2 for DNA (fully conserved)."""
    ic = np.zeros(ppm.shape[1])
    for pos in range(ppm.shape[1]):
        p = ppm[:, pos]
        entropy = -np.sum(p * np.log2(p + 1e-12))
        ic[pos] = 2.0 - entropy  # log2(4) - H
    return ic
```

IC interpretation:
- 2 bits — fully conserved position
- 0 bits — random / no conservation
- Total IC (sum over all positions) reflects motif specificity
  - CTCF: ~15–16 bits; TATA box: ~10–12 bits

## PROSITE Pattern → Python Regex

```python
import re

def prosite_to_regex(pattern: str) -> str:
    """Convert PROSITE pattern string to Python regex string.

    Syntax: N-{P}-[ST]-{P}  →  N[^P][ST][^P]
    x = any AA, [ABC] = one of, {P} = not P, (n) = repeat, (n,m) = range
    """
    pattern = pattern.strip(".")
    parts = []
    for elem in pattern.split("-"):
        m = re.match(r'^(.+?)\((\d+)(?:,(\d+))?\)$', elem)
        core = m.group(1) if m else elem
        low, high = (m.group(2), m.group(3)) if m else (None, None)

        if core == "x":
            r = "."
        elif core.startswith("[") and core.endswith("]"):
            r = core
        elif core.startswith("{") and core.endswith("}"):
            r = f"[^{core[1:-1]}]"
        else:
            r = core

        if low:
            r += f"{{{low},{high}}}" if high else f"{{{low}}}"
        parts.append(r)
    return "".join(parts)

# Examples:
# N-glycosylation:   N-{P}-[ST]-{P}  →  N[^P][ST][^P]
# Kinase C site:     [ST]-x-[RK]     →  [ST].[RK]
# Zinc finger C2H2:  C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H
```

## Pitfalls

- **PFM vs PPM vs PWM**: only the PWM (log-odds) is suitable for scoring. Raw PPM values near zero cause extreme sensitivity. PFM counts are just for display.
- **Pseudocounts are mandatory**: a single missing base at any position gives log(0) = -∞ in the PWM. Use 0.5 (or 0.1× background freq) as the pseudocount.
- **Threshold selection is empirical**: no universal correct value. Common approaches: 60–80% of max possible score, or calibrated from ChIP-seq positive set.
- **Always scan both strands**: TF binding sites can be on either strand; single-strand scanning misses ~half of sites.
- **PROSITE is qualitative, PWM is quantitative**: PROSITE patterns are all-or-nothing; for degenerate sites, PWMs give far better sensitivity/specificity.
- **Information content vs raw frequency**: two positions can have identical consensus but very different IC — a 60%/40% split vs 99%/1% split look similar in sequence but differ in motif strength.
