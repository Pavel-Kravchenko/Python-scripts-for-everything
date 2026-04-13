---
name: ai-science-splicing-models
description: "**Tier 5 — Modern AI for Science | Module 05 · Notebook 3**"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/03_splicing_models.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Splicing Models: SpliceAI and AlphaGenome

*Source: Course notebook `Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/03_splicing_models.ipynb`*


**Tier 5 — Modern AI for Science | Module 05 · Notebook 3**

*Prerequisites: Notebook 2 (Enformer)*

---

**By the end of this notebook you will be able to:**
1. Explain donor and acceptor splice-site logic
2. Compute toy splice delta scores for SNVs
3. Prioritize splice-disrupting variants from candidate lists
4. Contrast splice-specialized vs multi-task interpretation

## Why this notebook matters

Splicing mutations are responsible for ~10–15% of all pathogenic variants causing human disease. Many of these variants do not disrupt the canonical GT/AG dinucleotides — instead they weaken splice site recognition, create cryptic splice sites, or disrupt splicing regulatory elements. SpliceAI (Jaganathan et al., Cell 2019) was the first deep learning model to predict these effects from raw sequence with near-clinical accuracy. Understanding the delta score framework is essential for rare disease variant triage.

## How to work through this notebook

1. Read the splice signal refresher (Section 1) to establish vocabulary: donor = 5' end of intron (GT), acceptor = 3' end of intron (AG).
2. The delta score framework (Section 2) is the core clinical output — understand all four delta values before the code.
3. Sections 3–4 show prioritization and multi-modal integration workflows.

## Common sticking points

- **SpliceAI's 50 nt context window**: SpliceAI uses a deep residual convolutional network that sees a 10,000 bp window centered on each variant, but the reported delta scores are for the 50 nearest splice sites. The toy model here uses a 9 bp window — much shorter — as an illustration.
- **Four delta scores**: DS_DG (donor gain), DS_DL (donor loss), DS_AG (acceptor gain), DS_AL (acceptor loss). Each ranges 0–1. A threshold of 0.2 is commonly used for clinical flagging; 0.8 is high confidence.
- **Cryptic splice sites**: DS_DG and DS_AG > 0 mean the variant creates a NEW splice site that didn't exist before. These are harder to interpret than losses (which simply weaken an existing site) and require evaluating the new exon/intron structure.
- **SpliceAI vs Pangolin**: Pangolin extends SpliceAI to predict tissue-specific splicing, important for diseases where splicing patterns differ between tissues (e.g., brain vs liver isoforms).

## 1. Splice Signals Refresher

Canonical motifs:
- **Donor (5' splice site)**: usually `GT`
- **Acceptor (3' splice site)**: usually `AG`

SpliceAI-like models learn richer context, but motif disruption is still a useful intuition layer.

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
```python

## 2. Compute Splice Delta Scores (DS)

We mimic SpliceAI-style deltas:
- `DS_DG`: donor gain
- `DS_DL`: donor loss
- `DS_AG`: acceptor gain
- `DS_AL`: acceptor loss

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
```python

## 3. Rank Candidate Variants

A practical heuristic is to rank by `max(DS_*)`.

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
```python

## 4. Splice-Specialized vs Multi-Task Outputs

- **SpliceAI-style**: best for splice disruption prioritization with explicit splice deltas.
- **AlphaGenome-style multi-task**: adds expression/chromatin/contact signals for broader interpretation.

In practice, many teams use both: splice model for high-sensitivity filtering, multi-task model for mechanism context.

```python
def integrated_priority(max_ds: float, expr_delta: float, chrom_delta: float) -> float:
    # weighted integration example (toy)
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
```python

## Summary

- Splice deltas are highly actionable for variant triage.
- `max(DS_*)` is a robust first-pass ranking signal.
- Integrating splice scores with broader multi-omic predictions improves interpretation depth.

## Source-backed Context

- SpliceAI is explicitly designed to annotate variant effects on splicing with donor/acceptor gain/loss-style outputs.
- AlphaGenome presents unified variant-effect prediction across expression, splicing, chromatin, and contact modalities.

## Validated Sources

Checked online during content expansion.

- [Jaganathan et al. 2019, SpliceAI (Cell)](https://doi.org/10.1016/j.cell.2018.12.015)
- [SpliceAI repository](https://github.com/Illumina/SpliceAI)
- [AlphaGenome paper (Nature 2026)](https://www.nature.com/articles/s41586-025-10014-0)
- [AlphaGenome research repository](https://github.com/google-deepmind/alphagenome_research)

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
