---
name: bio-applied-promoter
description: "Core promoter elements: TATA box, CpG islands, PWM construction, and TFBS scanning. Companion reference card to bio-applied-regulatory-analysis."
tool_type: python
primary_tool: NumPy
---

# Core Promoter Elements

> For full regulatory analysis including GRN context, see `bio-applied-regulatory-analysis`.

## Promoter Architecture

| Element | Position | Notes |
|---------|----------|-------|
| TATA box | −25 to −35 | TATAAA consensus; ~10–20% of human genes |
| Initiator (Inr) | +1 | Present in TATA-less promoters |
| DPE | +28 to +32 | Downstream promoter element; compensates for absent TATA |
| CpG island | −200 to +100 | ~70% of human promoters; mark of constitutive expression |

## TATA Box Finder

```python
import re

def find_tata_boxes(sequence, strict=True):
    """strict=True: exact TATAAA; False: TATAWAW (IUPAC)."""
    pattern = 'TATAAA' if strict else r'TATA[AT]A[AT]'
    return [(m.start(), m.group())
            for m in re.finditer(pattern, sequence.upper())]
```

## CpG Island Scanner

```python
def cpg_island_scanner(sequence, window=200, step=50, gc_thresh=0.5, oe_thresh=0.6):
    """CpG O/E = (N_CpG * L) / (N_C * N_G). Criteria: len>=200, GC>=50%, O/E>=0.6."""
    sequence = sequence.upper()
    islands = []
    for i in range(0, len(sequence) - window + 1, step):
        win = sequence[i:i + window]
        n_c, n_g, n_cpg = win.count('C'), win.count('G'), win.count('CG')
        gc = (n_c + n_g) / window
        oe = (n_cpg * window) / (n_c * n_g) if n_c > 0 and n_g > 0 else 0.0
        if gc >= gc_thresh and oe >= oe_thresh:
            islands.append((i, i + window, round(gc, 3), round(oe, 3)))
    return islands
```

## PWM Construction and Scanning

```python
import numpy as np

def build_pwm(sites, pseudocount=0.5, bg=None):
    """Build log2-odds PWM directly from aligned binding site strings."""
    bg = bg or {b: 0.25 for b in 'ACGT'}
    n, L = len(sites), len(sites[0])
    counts = {b: [0] * L for b in 'ACGT'}
    for site in sites:
        for i, b in enumerate(site.upper()):
            counts[b][i] += 1
    pwm = {}
    for b in 'ACGT':
        pwm[b] = [np.log2((counts[b][i] + pseudocount) /
                           (n + 4 * pseudocount) / bg[b])
                  for i in range(L)]
    return pwm

def scan_pwm(sequence, pwm, threshold=0.0):
    L = len(pwm['A'])
    sequence = sequence.upper()
    return [(i, sequence[i:i+L],
             sum(pwm.get(b, [0]*L)[j] for j, b in enumerate(sequence[i:i+L])))
            for i in range(len(sequence) - L + 1)
            if sum(pwm.get(b, [0]*L)[j] for j, b in enumerate(sequence[i:i+L])) >= threshold]
```

## Pitfalls

- **CpG O/E formula**: denominator is `N_C * N_G`, not `(N_C + N_G)^2 / 4` — use the Gardiner-Garden & Frommer (1987) formula shown above
- **Zero counts in PFM → −∞ PWM**: always use pseudocounts (0.5 is standard)
- **TATA box false positives**: random AT-rich sequences generate TATAAA by chance; context (−30 relative to TSS) matters
- **Coordinate systems**: BED = 0-based half-open; VCF/GFF = 1-based inclusive
- **step size in scanner**: step=1 gives exact island boundaries; step=50 misses islands smaller than the step
