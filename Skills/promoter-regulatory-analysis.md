---
name: promoter-regulatory-analysis
description: Promoter element detection (TATA box, CpG islands, Inr, DPE), PWM/TFBS construction and scanning, TSS prediction concepts, methylation context, regulatory sequence analysis
---

# Promoter & Regulatory Sequence Analysis

## When to Use
- Scanning sequences for core promoter elements (TATA box, Inr, DPE)
- Detecting CpG islands near transcription start sites
- Building Position Weight Matrices (PWMs) from aligned TF binding sites
- Scoring genomic sequences against PWMs for TFBS prediction
- Understanding regulatory architecture upstream of a gene of interest

---

## Quick Reference

### Core Promoter Elements
| Element | Position (rel. to TSS) | Consensus | Notes |
|---------|------------------------|-----------|-------|
| TATA box | −25 to −30 bp | TATAAA (TATAWAW) | ~10–20% human genes; tissue-specific |
| CpG island | Near TSS | GC ≥ 50%, CpG O/E ≥ 0.6, ≥ 200 bp | ~70% human gene promoters |
| Inr | +1 (TSS) | YYANWYY | TATA-less promoters |
| DPE | +28 to +32 | RGWYV | Pairs with Inr in TATA-less promoters |
| BRE | −37 (upstream) | SSRCGCC | TFIIB recognition element |

CpG O/E formula: `(N_CpG × L) / (N_C × N_G)` — values ≥ 0.6 define an island

### Two Promoter Architectures
| Type | Elements | Gene examples |
|------|----------|---------------|
| TATA-containing | TATA box + Inr | Highly regulated, tissue-specific genes |
| CpG-island | CpG island + Inr ± DPE | Housekeeping genes; ubiquitous expression |

### PWM Key Formulas
| Quantity | Formula | Notes |
|----------|---------|-------|
| PFM | raw counts per base per position | Integer matrix, L × 4 |
| PPM | `(count + pc) / (N + 4·pc)` | pc = pseudocount (0.5 typical) |
| PWM (log-odds) | `log₂(PPM[b,i] / bg[b])` | bg = 0.25 for uniform background |
| Position IC | `2 + Σ_b PPM[b,i] · log₂(PPM[b,i])` | bits; 2 = fully conserved; 0 = uniform |
| Total IC | `Σ_i IC[i]` | Higher total IC = more specific motif |
| KDIC | `mean(IC) / 2` | Normalised to [0, 1] |

---

## Key Patterns

### TATA Box Detection
```python
import re

def find_tata_boxes(seq: str, strict: bool = True) -> list[tuple[int, str]]:
    """Find TATA box motifs. strict=True: exact TATAAA; False: TATAWAW."""
    pattern = 'TATAAA' if strict else r'TATA[AT]A[AT]'
    return [(m.start(), m.group()) for m in re.finditer(pattern, seq.upper())]

# Look for TATA boxes in the −40 to −20 window upstream of your TSS
promoter_seq = "GCGATCGTATAAAAGTCGATCG"
hits = find_tata_boxes(promoter_seq, strict=False)
```

### CpG Island Scanner
```python
def cpg_island_scanner(seq: str, window: int = 200, step: int = 1,
                        gc_thresh: float = 0.5,
                        oe_thresh: float = 0.6) -> list[tuple]:
    """Slide a window and return (start, end, gc, oe) for qualifying windows."""
    seq = seq.upper()
    islands = []
    for i in range(0, len(seq) - window + 1, step):
        w = seq[i:i + window]
        nc, ng = w.count('C'), w.count('G')
        ncpg = w.count('CG')
        gc = (nc + ng) / window
        oe = (ncpg * window) / (nc * ng) if nc > 0 and ng > 0 else 0.0
        if gc >= gc_thresh and oe >= oe_thresh:
            islands.append((i, i + window, gc, oe))
    return islands


def merge_cpg_windows(windows: list[tuple], gap: int = 100) -> list[tuple]:
    """Merge overlapping/adjacent CpG island windows."""
    if not windows:
        return []
    merged = [list(windows[0])]
    for start, end, gc, oe in windows[1:]:
        if start <= merged[-1][1] + gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end, gc, oe])
    return [tuple(m) for m in merged]
```

### PWM Construction
```python
import numpy as np

BASES = 'ACGT'
BASE_IDX = {b: i for i, b in enumerate(BASES)}

def build_pfm(sites: list[str]) -> np.ndarray:
    """Build PFM (L × 4) from aligned binding-site strings."""
    L = len(sites[0])
    pfm = np.zeros((L, 4), dtype=int)
    for site in sites:
        for i, b in enumerate(site.upper()):
            if b in BASE_IDX:
                pfm[i, BASE_IDX[b]] += 1
    return pfm

def pfm_to_pwm(pfm: np.ndarray,
               pseudocount: float = 0.5,
               bg: float = 0.25) -> np.ndarray:
    """Convert PFM → PPM → PWM (log-odds, L × 4)."""
    n = pfm.sum(axis=1, keepdims=True)
    ppm = (pfm + pseudocount) / (n + 4 * pseudocount)
    return np.log2(ppm / bg)
```

### PWM Scanning
```python
def scan_pwm(seq: str, pwm: np.ndarray,
             threshold: float = 0.0) -> list[tuple[int, str, float]]:
    """Scan seq with PWM (L × 4). Returns (pos, kmer, score) above threshold."""
    L = pwm.shape[0]
    hits = []
    seq = seq.upper()
    for i in range(len(seq) - L + 1):
        kmer = seq[i:i + L]
        score = sum(pwm[j, BASE_IDX[b]]
                    for j, b in enumerate(kmer) if b in BASE_IDX)
        if score >= threshold:
            hits.append((i, kmer, score))
    return hits
```

---

## Code Templates

### Complete Promoter Scan Pipeline
```python
def scan_promoter(upstream_seq: str,
                  tf_sites: list[str],
                  pwm_threshold: float = 0.0,
                  cpg_window: int = 200) -> dict:
    """
    Full promoter analysis:
      - Find TATA boxes (strict + loose)
      - Detect CpG islands
      - Build PWM from provided TF sites and scan for TFBS
    Returns a summary dict.
    """
    results = {}

    # TATA boxes
    results['tata_strict'] = find_tata_boxes(upstream_seq, strict=True)
    results['tata_loose']  = find_tata_boxes(upstream_seq, strict=False)

    # CpG islands (use step=10 for long sequences)
    step = 1 if len(upstream_seq) < 5000 else 10
    raw_islands = cpg_island_scanner(upstream_seq, window=cpg_window, step=step)
    results['cpg_islands'] = merge_cpg_windows(raw_islands)

    # TFBS scan
    pfm = build_pfm(tf_sites)
    pwm = pfm_to_pwm(pfm)
    results['tfbs_hits'] = scan_pwm(upstream_seq, pwm, threshold=pwm_threshold)

    return results
```

### Information Content per Position
```python
def information_content(ppm: np.ndarray) -> tuple[np.ndarray, float]:
    """Return per-position IC (bits) and total IC. ppm: L × 4."""
    ic = 2 + np.sum(ppm * np.log2(ppm + 1e-9), axis=1)
    return ic, float(ic.sum())

def kdic(ppm: np.ndarray) -> float:
    """KDIC — normalised motif quality score in [0, 1]."""
    ic, _ = information_content(ppm)
    return float(ic.mean() / 2)
```

### IUPAC Consensus from PPM
```python
IUPAC = {
    frozenset('A'): 'A', frozenset('C'): 'C',
    frozenset('G'): 'G', frozenset('T'): 'T',
    frozenset('AG'): 'R', frozenset('CT'): 'Y',
    frozenset('CG'): 'S', frozenset('AT'): 'W',
    frozenset('GT'): 'K', frozenset('AC'): 'M',
    frozenset('CGT'): 'B', frozenset('AGT'): 'D',
    frozenset('ACT'): 'H', frozenset('ACG'): 'V',
    frozenset('ACGT'): 'N',
}

def ppm_to_iupac(ppm: np.ndarray, ic_thresh: float = 1.0) -> str:
    """Generate IUPAC consensus. Degenerate symbol when IC < ic_thresh."""
    ic, _ = information_content(ppm)
    consensus = []
    for i in range(len(ppm)):
        if ic[i] >= ic_thresh:
            consensus.append(BASES[np.argmax(ppm[i])])
        else:
            dominant = {BASES[j] for j in range(4) if ppm[i, j] >= 0.25}
            consensus.append(IUPAC.get(frozenset(dominant), 'N'))
    return ''.join(consensus)
```

---

## Common Pitfalls

- **CpG island detection with step=1 on large sequences** — use `step=10` or `step=50` for sequences >10 kb; then merge overlapping windows.
- **PWM threshold too low** — produces many false positives. Use the 80–90th percentile of scores from random sequence as a baseline cutoff.
- **Forgetting pseudocounts in PWM** — a single missing base in the training set causes −∞ log-odds for any sequence containing that base at that position.
- **Background frequency assumption** — the default 0.25 flat background underestimates scores in AT-rich promoters; use regional GC content for more accurate scoring.
- **TATA box prevalence** — only ~10–20% of human genes have a canonical TATA box; absence of TATA does not mean absence of regulation (CpG island promoters are TATA-less).
- **CpG island O/E ≥ 0.6 alone is insufficient** — always require GC ≥ 50% and length ≥ 200 bp; many GC-rich regions satisfy the ratio alone but are not functional islands.
- **Single PWM hits ≠ functional binding** — combine PWM scanning with open-chromatin data (ATAC-seq, DNase-seq) for meaningful predictions.

---

## Related Skills
- `rnaseq` — RNA-seq differential expression; Module 03 (upstream of promoter analysis)
- `motif-discovery` — PPM/PWM construction, IC/KDIC scoring, Fisher enrichment, TomTom matching
- `tf-footprinting-atac` — ATAC-seq integration for detecting TF occupancy from chromatin accessibility
- `chipseq-epigenomics` — ChIP-seq for validating predicted TFBS experimentally
