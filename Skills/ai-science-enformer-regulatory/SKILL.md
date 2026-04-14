---
name: ai-science-enformer-regulatory
description: "Enformer architecture for regulatory prediction from DNA, in-silico mutagenesis (ISM), and variant prioritization"
tool_type: python
primary_tool: NumPy
---

# Enformer: Regulatory Activity from DNA

## Key Facts

- **Input**: one-hot DNA (196,608 bp × 4). Centered on region of interest
- **Output**: 896 bins × 128 bp (~115 kb predicted region) across ~5,313 tracks (CAGE, DNase, ATAC, ChIP-seq, 100+ cell types)
- **Architecture**: Stem conv → 6 dilated conv blocks → 11 transformer layers (8 heads, relative pos encodings) → pointwise conv
- **Training**: Poisson log-likelihood loss, softplus output, hg38 + mm10
- **Borzoi** extends to 524 kb input, 32 bp resolution, RNA-seq coverage (including splicing patterns)

## Gotchas

- **196 kb context**: outer flanks provide context but are not predicted — only central ~115 kb
- **ISM vs gradient methods**: ISM is simple but O(3L) model calls. DeepLIFT/GradCAM are faster but need model internals
- The 11 transformer layers are what differentiate Enformer from conv-only baselines (Basenji2)

## In-Silico Mutagenesis (ISM)

```python
import numpy as np

def mutate_base(seq: str, pos: int, alt: str) -> str:
    return seq[:pos] + alt + seq[pos + 1:]

def ism_scores(seq: str, pos: int, model_fn):
    baseline = model_fn(seq)
    ref = seq[pos]
    effects = {}
    for alt in "ACGT":
        if alt == ref:
            continue
        effects[alt] = model_fn(mutate_base(seq, pos, alt)) - baseline
    return ref, baseline, effects
```

## Variant Prioritization

Rank by absolute predicted effect on track of interest. Combine with population frequency, GWAS/eQTL evidence, and experimental context.

```python
def score_variant(seq, pos, alt, track, model_fn):
    base = model_fn(seq)[track]
    mut = model_fn(mutate_base(seq, pos, alt))[track]
    return float(mut - base)
```

## Sources

- [Avsec et al. 2021, Enformer (Nature Methods)](https://www.nature.com/articles/s41592-021-01252-x)
- [DeepMind Enformer](https://github.com/google-deepmind/deepmind-research/tree/master/enformer)
- [Borzoi](https://github.com/calico/borzoi)
