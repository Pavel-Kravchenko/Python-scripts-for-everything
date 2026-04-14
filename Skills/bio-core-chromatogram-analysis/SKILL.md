---
name: bio-core-chromatogram-analysis
description: Sanger sequencing chromatogram analysis — reading .ab1 files with BioPython, Phred quality scores, trace visualization, and quality-based trimming
tool_type: python
primary_tool: NumPy
---

# Chromatogram Analysis: Sanger Sequencing

## Phred Quality Score Reference

| Q Score | Error Probability | Accuracy | Interpretation |
|:-------:|:-----------------:|:--------:|---------------|
| 10 | 1/10 | 90% | Poor — avoid |
| 20 | 1/100 | 99% | Minimum acceptable |
| 30 | 1/1,000 | 99.9% | High confidence |
| 40 | 1/10,000 | 99.99% | Excellent |

Formula: Q = -10 × log₁₀(P_error)

Practical thresholds: Q<20 = trim/re-sequence; Q20-29 = acceptable; Q≥30 = high confidence.

## Sanger Read Quality Profile

- Bases 1–25: poor quality (primer artifact, dye blobs) — always trim
- Bases 25–700: high quality plateau (Q>30 typical)
- Bases 700+: signal decay, quality drops sharply — trim tail

## Reading .ab1 Files

```python
from Bio import SeqIO
import numpy as np

record = SeqIO.read("sample.ab1", "abi")

# Basic info
sequence = str(record.seq)
quals = record.letter_annotations["phred_quality"]  # list of ints

# Raw trace data
raw = record.annotations["abif_raw"]
base_order = raw.get("FWO_1", b"GATC").decode()  # e.g. "GATC"

# Extract processed trace channels
channels = {}
for i, base in enumerate(base_order):
    channels[base] = np.array(raw[f"DATA{9 + i}"])   # DATA9..DATA12

peak_locs = list(raw["PLOC1"])   # scan positions of called bases
```

### Key .ab1 Tags

| Tag | Content |
|-----|---------|
| `DATA9`–`DATA12` | Processed trace channels (order from `FWO_1`) |
| `DATA1`–`DATA4` | Raw trace channels |
| `PLOC1` / `PLOC2` | Peak locations (scan numbers) |
| `PBAS1` / `PBAS2` | Called bases (raw / edited) |
| `PCON1` / `PCON2` | Quality values |
| `FWO_1` | Base order for channels (e.g. `b'GATC'`) |

Note: `.ab1` and `.abi` are the same format; `.scf` is an older MRC-LMB format. BioPython reads all via `"abi"`.

## Trace Visualization

```python
import matplotlib.pyplot as plt

TRACE_COLORS = {"A": "green", "C": "blue", "G": "black", "T": "red"}

def plot_chromatogram(channels, peak_locs, sequence, start=0, end=50):
    scan_start = peak_locs[start] - 10
    scan_end = peak_locs[min(end, len(peak_locs)) - 1] + 10
    x = np.arange(scan_start, scan_end)

    fig, ax = plt.subplots(figsize=(14, 3.5))
    for base, color in TRACE_COLORS.items():
        ax.plot(x, channels[base][scan_start:scan_end], color=color, lw=0.8, label=base)

    for i in range(start, min(end, len(peak_locs))):
        base = sequence[i]
        ax.text(peak_locs[i], 0, base, ha="center", fontsize=7,
                fontweight="bold", color=TRACE_COLORS.get(base, "gray"))
    ax.legend(ncol=4, fontsize=8)
    ax.set_xlabel("Scan position"); ax.set_ylabel("Fluorescence")
    plt.tight_layout(); plt.show()
```

## Quality Trimming

```python
def trim_by_quality(sequence, quals, min_q=20, window=10):
    """Trim low-quality ends using sliding window."""
    # Find first window where mean quality >= min_q (5' end)
    start = 0
    for i in range(len(quals) - window + 1):
        if np.mean(quals[i:i+window]) >= min_q:
            start = i
            break
    # Find last window (3' end)
    end = len(quals)
    for i in range(len(quals) - window, -1, -1):
        if np.mean(quals[i:i+window]) >= min_q:
            end = i + window
            break
    return sequence[start:end], quals[start:end]
```

## Pitfalls

- **Double peaks = heterozygosity OR contamination**: two overlapping peaks of comparable height at one position — check surrounding context; genuine het SNPs form a consistent pattern across the read; contamination creates irregular double peaks throughout
- **First 20–30 bases always low quality**: this is normal (primer signal + dye artifacts) — never report variants from this region
- **Why 4 dyes**: each ddNTP (ddA, ddC, ddG, ddT) carries a different fluorescent dye; capillary electrophoresis separates by size, laser detects dye — the instrument reconstructs sequence from size+color
- **`.ab1` stores raw fluorescence, not just sequence**: the base calls are computed by KB Basecaller software; you can re-call from the raw traces with third-party tools
- **`DATA9`–`DATA12` are processed; `DATA1`–`DATA4` are raw**: use processed channels for visualization (they have baseline correction and color deconvolution applied)
- **Sanger max read length ~800 bp**: beyond this, band resolution collapses — design primers to place your region of interest in the Q>30 window (bases ~30–700)
