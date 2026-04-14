---
name: bio-applied-regulatory-analysis
description: "Promoter and regulatory sequence analysis: TATA box detection, CpG island scanning, PWM/PFM construction, and TFBS scanning. Reference for computational promoter analysis."
tool_type: python
primary_tool: NumPy
---

# Promoter and Regulatory Sequence Analysis

## Key Regulatory Elements

| Element | Location | Function |
|---------|----------|----------|
| Promoter | −1 to −1000 bp from TSS | Recruits RNA Pol II |
| Enhancer | Distal (kb–Mb away) | Boosts transcription |
| Silencer | Variable | Represses transcription |
| Insulator | Between elements | Blocks enhancer–promoter crosstalk |

- **TSS** = position +1; upstream = negative coordinates
- **TATA box**: consensus TATAAA at ~−30; present in ~10–20% of human genes (TATA-less promoters use Inr, DPE)
- **CpG islands**: near ~70% of human gene promoters; criteria: ≥200 bp, GC ≥50%, CpG O/E ≥0.6

## TATA Box Detection

```python
import re

def find_tata_boxes(sequence, strict=True):
    sequence = sequence.upper()
    pattern = 'TATAAA' if strict else r'TATA[AT]A[AT]'
    return [(m.start(), m.group()) for m in re.finditer(pattern, sequence)]
```

## CpG Island Scanner

```python
def cpg_island_scanner(sequence, window=200, step=50, gc_thresh=0.5, oe_thresh=0.6):
    """Sliding-window CpG island detection. Use step=1 for exact boundaries."""
    sequence = sequence.upper()
    islands = []
    for i in range(0, len(sequence) - window + 1, step):
        win = sequence[i:i + window]
        n_c, n_g = win.count('C'), win.count('G')
        n_cpg = win.count('CG')
        gc = (n_c + n_g) / window
        oe = (n_cpg * window) / (n_c * n_g) if n_c > 0 and n_g > 0 else 0.0
        if gc >= gc_thresh and oe >= oe_thresh:
            islands.append((i, i + window, gc, oe))
    return islands
```

## PWM/PFM Construction and Scanning

```python
import numpy as np

def build_pfm(sites):
    """Build position frequency matrix from aligned binding sites (same length)."""
    length = len(sites[0])
    pfm = {b: [0] * length for b in 'ACGT'}
    for site in sites:
        for i, base in enumerate(site.upper()):
            pfm[base][i] += 1
    return pfm

def pfm_to_pwm(pfm, pseudocount=0.5, bg=None):
    """Convert PFM to log2 odds PWM."""
    bg = bg or {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    n = sum(pfm[b][0] for b in 'ACGT')
    length = len(pfm['A'])
    pwm = {b: [] for b in 'ACGT'}
    for b in 'ACGT':
        for i in range(length):
            freq = (pfm[b][i] + pseudocount) / (n + 4 * pseudocount)
            pwm[b].append(np.log2(freq / bg[b]))
    return pwm

def scan_with_pwm(sequence, pwm, threshold=0.0):
    """Return (pos, subseq, score) tuples above threshold."""
    motif_len = len(pwm['A'])
    hits = []
    sequence = sequence.upper()
    for i in range(len(sequence) - motif_len + 1):
        sub = sequence[i:i + motif_len]
        score = sum(pwm.get(b, {i: 0})[j] if b in pwm else 0
                    for j, b in enumerate(sub))
        if score >= threshold:
            hits.append((i, sub, score))
    return hits

# Cleaner score_sequence helper
def score_sequence(subseq, pwm):
    return sum(pwm[b][i] for i, b in enumerate(subseq.upper()) if b in pwm)
```

## TSS Prediction Signals

| Signal | Peak location | Method |
|--------|--------------|--------|
| TATA box | −30 | Motif scan |
| CpG island | centered on TSS | O/E ratio |
| TFBS density | −200 upstream | PWM scan |
| CAGE signal | +1 | Experimental |

## TF Motif Databases

| Database | Description | URL |
|----------|-------------|-----|
| JASPAR | Open-access, curated | jaspar.elixir.no |
| HOCOMOCO | Human/mouse from ChIP-seq | hocomoco11.autosome.org |
| TRANSFAC | Comprehensive (commercial) | genexplain.com/transfac |

## Pitfalls

- **CpG underrepresentation**: vertebrate bulk genome has CpG O/E ~0.2 due to methylation-driven mutation; islands (O/E ≥0.6) are real anomalies
- **Coordinate systems**: BED = 0-based half-open; VCF/GFF = 1-based inclusive — mixing causes off-by-one errors
- **PWM threshold selection**: score ≥80% of max PWM score is a common heuristic; too low → flood of false positives
- **Pseudocount matters**: without pseudocounts, a single zero count in PFM gives −∞ PWM score for any sequence containing that base at that position
- **Multiple testing**: scanning a genome-wide promoter set requires FDR correction (Benjamini-Hochberg)
