---
name: ai-science-splicing-models
description: "Splicing Models: SpliceAI and AlphaGenome with NumPy"
tool_type: python
primary_tool: NumPy
---

# Splicing Models: SpliceAI and AlphaGenome

## Splice Signal Refresher

Canonical motifs:
- **Donor (5' splice site)**: usually `GT`
- **Acceptor (3' splice site)**: usually `AG`

## Delta Score Framework

Four SpliceAI delta scores (range 0–1):
- `DS_DG`: donor gain — variant creates a NEW donor site
- `DS_DL`: donor loss — variant weakens an existing donor
- `DS_AG`: acceptor gain — variant creates a NEW acceptor site
- `DS_AL`: acceptor loss — variant weakens an existing acceptor

Thresholds: ≥0.2 for clinical flagging; ≥0.8 high confidence. Gain scores (DS_DG/DS_AG) require evaluating the new exon/intron structure — harder to interpret than losses.

**SpliceAI**: deep residual ConvNet, 10,000 bp window, reports delta scores for 50 nearest splice sites.
**Pangolin**: extends SpliceAI to tissue-specific splicing.
**AlphaGenome**: multi-task model adding expression/chromatin/contact signals alongside splice deltas.

```python
import numpy as np

def donor_score(window: str) -> float:
    # Toy donor score based on central dinucleotide and nearby G-rich context
    center = window[len(window)//2 - 1: len(window)//2 + 1]
    score = 0.7 if center == "GT" else 0.1
    score += 0.03 * window.count("G")
    return min(score, 1.0)

def acceptor_score(window: str) -> float:
    center = window[len(window)//2 - 1: len(window)//2 + 1]
    score = 0.7 if center == "AG" else 0.1
    score += 0.02 * window.count("T")
    return min(score, 1.0)
```

## Compute Splice Delta Scores

```python
def mutate(seq: str, pos: int, alt: str) -> str:
    return seq[:pos] + alt + seq[pos + 1:]

def splice_deltas(seq: str, pos: int, alt: str, w: int = 9):
    half = w // 2
    start = max(0, pos - half)
    end = min(len(seq), pos + half + 1)

    ref_window = seq[start:end]
    alt_seq = mutate(seq, pos, alt)
    alt_window = alt_seq[start:end]

    d_ref = donor_score(ref_window)
    d_alt = donor_score(alt_window)
    a_ref = acceptor_score(ref_window)
    a_alt = acceptor_score(alt_window)

    return {
        "DS_DG": max(0.0, d_alt - d_ref),
        "DS_DL": max(0.0, d_ref - d_alt),
        "DS_AG": max(0.0, a_alt - a_ref),
        "DS_AL": max(0.0, a_ref - a_alt),
    }

example = "CCTGACTGGTGAGTCTCAGGTTAC"
pos = 9
for alt in "ACGT":
    if alt == example[pos]:
        continue
    print(example[pos], ">", alt, splice_deltas(example, pos, alt))
```

## Rank Candidate Variants

```python
candidates = [
    (8, "A"), (8, "C"), (9, "A"), (10, "T"), (14, "G"), (17, "A")
]

ranked = []
for p, alt in candidates:
    if alt == example[p]:
        continue
    ds = splice_deltas(example, p, alt)
    max_ds = max(ds.values())
    ranked.append((p, example[p], alt, max_ds, ds))

ranked.sort(key=lambda x: x[3], reverse=True)

for r in ranked:
    print(f"pos={r[0]} {r[1]}>{r[2]} maxDS={r[3]:.3f} details={r[4]}")
```

## Multi-Modal Integration

```python
def integrated_priority(max_ds: float, expr_delta: float, chrom_delta: float) -> float:
    # weighted integration (toy)
    return 0.6 * max_ds + 0.25 * abs(expr_delta) + 0.15 * abs(chrom_delta)

demo = [
    {"id": "v1", "max_ds": 0.85, "expr": 0.10, "chrom": 0.05},
    {"id": "v2", "max_ds": 0.40, "expr": 0.50, "chrom": 0.30},
    {"id": "v3", "max_ds": 0.65, "expr": 0.20, "chrom": 0.05},
]
for v in demo:
    v["priority"] = integrated_priority(v["max_ds"], v["expr"], v["chrom"])

for v in sorted(demo, key=lambda d: d["priority"], reverse=True):
    print(v)
```

## References

- [Jaganathan et al. 2019, SpliceAI (Cell)](https://doi.org/10.1016/j.cell.2018.12.015)
- [SpliceAI repository](https://github.com/Illumina/SpliceAI)
- [AlphaGenome paper (Nature 2026)](https://www.nature.com/articles/s41586-025-10014-0)
- [AlphaGenome research repository](https://github.com/google-deepmind/alphagenome_research)

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
