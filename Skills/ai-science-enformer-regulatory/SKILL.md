---
name: ai-science-enformer-regulatory
description: "**Tier 5 — Modern AI for Science | Module 05 · Notebook 2**"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/02_enformer_regulatory.ipynb"
---

# Enformer: Predicting Regulatory Activity from DNA

*Source: Course notebook `Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/02_enformer_regulatory.ipynb`*

# Enformer: Predicting Regulatory Activity from DNA

**Tier 5 — Modern AI for Science | Module 05 · Notebook 2**

*Prerequisites: Notebook 1 (Genomic LLMs)*

---

**By the end of this notebook you will be able to:**
1. Explain Enformer-style sequence-to-track prediction
2. Build a toy sequence-to-signal model with motif logic
3. Run in-silico mutagenesis (ISM) to score variant effects
4. Interpret directional effect sizes for variant prioritization

## Why this notebook matters

Most genomic features — gene expression levels, chromatin accessibility, TF binding sites — are determined by the local DNA sequence and its regulatory context. Enformer (Avsec et al., 2021) demonstrated that a single deep-learning model can predict thousands of these tracks simultaneously from 196,608 bp of DNA input. This sequence-to-function paradigm is the foundation for computational variant effect prediction: instead of asking "is this variant in a known functional element?" we ask "does this variant change the predicted regulatory signal?"

## How to work through this notebook

1. Read Section 1 (Enformer architecture) to understand what the model inputs, outputs, and predicts.
2. The toy predictor (Section 2) is intentionally simple — focus on the workflow, not the model.
3. Run the ISM code (Section 3) to understand how single nucleotide changes translate to regulatory signal changes.
4. The variant prioritization pattern (Section 4) is the core practical output of this workflow.

## Common sticking points

- **Enformer's 196,608 bp context**: the input is centered on the target region. Enformer takes the full 196 kb window and predicts 896 bins of 128 bp each (covering ~115 kb of the central region). The outer flanks provide context but are not predicted.
- **Track interpretation**: real Enformer predicts ~5,313 tracks covering CAGE (RNA 5' ends), DNase-seq, ChIP-seq for 100s of TFs and histone marks, and ATAC-seq — across 100+ cell types. The toy model here has only 3 synthetic tracks.
- **ISM vs attribution methods**: ISM is simple (compute output for each possible single mutation) but slow for large sequences. Gradient-based attribution (DeepLIFT, GradCAM) is faster and gives continuous attribution scores, but requires access to the model's internals.
- **Borzoi context**: Borzoi (Linder et al.) uses a longer 524 kb input and predicts RNA-seq coverage at 32 bp resolution, rather than CAGE peaks. It is better suited for predicting read coverage patterns, including splicing and isoform-specific signals.

## 1. Enformer Architecture

Enformer (Avsec et al., 2021, Nature Methods) processes a 196,608 bp DNA window and predicts ~5,313 genomic tracks. Its architecture:

### Input
- One-hot encoded DNA: shape (196608, 4) — each position is a 4D vector (A, C, G, T)
- The sequence is always centered on the region of interest

### Trunk: Deep conv → Transformer
1. **Stem**: two 1D convolution layers with GELU activation reduce the sequence from 196,608 positions to a coarser representation
2. **Tower of conv blocks**: 6 residual conv blocks with exponentially increasing dilation, progressively compressing spatial resolution. After the tower, the representation is 1,536 positions × 768 channels.
3. **Transformer encoder**: 11 multi-head attention blocks (8 heads, relative positional encodings) running on the 1,536 position sequence. This allows long-range interactions across the 196 kb window.
4. **Pointwise convolution**: final 1×1 conv to produce 896 output bins × 3,072 channels

### Output heads
- **Human head**: 5,313 tracks (CAGE, DNase, ATAC, ChIP-seq for TFs and histone marks, across many cell types)
- **Mouse head**: 1,643 tracks
- Each track is a sequence of 896 bins × 128 bp = ~115 kb predicted region

### Training
- Trained on hg38 (human) and mm10 (mouse) with leave-out chromosomes for validation
- Loss: Poisson log-likelihood (suitable for count data)
- Soft-plus output activation to ensure non-negative predictions

### Key design decisions
- **Relative positional encodings in the transformer**: enables the model to generalize to different promoter-enhancer distance patterns
- **The 11 transformer layers** integrate long-range context that the conv tower cannot capture — this is what gives Enformer substantially better performance than conv-only baselines (like Basenji2)

In this notebook we use a toy model with the same input/output structure but no learned weights, to practice the ISM workflow.

```python
import numpy as np

np.random.seed(11)

def random_dna(n: int) -> str:
    return "".join(np.random.choice(list("ACGT"), size=n))

def count_motif(seq: str, motif: str) -> int:
    return sum(1 for i in range(len(seq) - len(motif) + 1) if seq[i:i + len(motif)] == motif)
```

## 2. Toy Multi-Track Regulatory Predictor

We define three synthetic tracks:
- **Track 0 (promoter-like)**: activated by `TATAAA`
- **Track 1 (enhancer-like)**: activated by `GATA`
- **Track 2 (repressor-like)**: decreased by `CGCG` (repressor motif)

This is a didactic stand-in for real Enformer outputs.

```python
def toy_regulatory_model(seq: str) -> np.ndarray:
    promoter = 0.4 + 0.25 * count_motif(seq, "TATAAA")
    enhancer = 0.2 + 0.10 * count_motif(seq, "GATA")
    repressor = 1.0 - 0.15 * count_motif(seq, "CGCG")
    vec = np.array([promoter, enhancer, repressor], dtype=float)
    return np.clip(vec, 0.0, 2.0)

seq = random_dna(300)
scores = toy_regulatory_model(seq)
print("Toy track scores:", np.round(scores, 3))
```

## 3. In-Silico Mutagenesis (ISM)

ISM scores each possible nucleotide substitution at a chosen position by comparing model output before and after mutation.

```python
def mutate_base(seq: str, pos: int, alt: str) -> str:
    return seq[:pos] + alt + seq[pos + 1:]

def ism_scores(seq: str, pos: int, model_fn):
    baseline = model_fn(seq)
    ref = seq[pos]
    effects = {}
    for alt in "ACGT":
        if alt == ref:
            continue
        mut = mutate_base(seq, pos, alt)
        effects[alt] = model_fn(mut) - baseline
    return ref, baseline, effects

test_seq = "A" * 120 + "TATAAA" + "A" * 120
pos = 122  # inside motif
ref, baseline, effects = ism_scores(test_seq, pos, toy_regulatory_model)

print("Reference base:", ref)
print("Baseline:", np.round(baseline, 3))
for alt, delta in effects.items():
    print(f"{ref}>{alt} delta:", np.round(delta, 3))
```

## 4. Variant Prioritization by Effect Size

A common pattern is ranking variants by absolute predicted effect on a track of interest.

```python
def score_variant(seq: str, pos: int, alt: str, track: int = 0) -> float:
    base = toy_regulatory_model(seq)[track]
    mut = toy_regulatory_model(mutate_base(seq, pos, alt))[track]
    return float(mut - base)

candidate_positions = [118, 120, 122, 124, 126]
variants = []
for p in candidate_positions:
    for alt in "ACGT":
        if alt == test_seq[p]:
            continue
        d = score_variant(test_seq, p, alt, track=0)
        variants.append((p, test_seq[p], alt, d))

variants = sorted(variants, key=lambda x: abs(x[3]), reverse=True)
print("Top promoter-impact variants:")
for v in variants[:6]:
    print(f"pos={v[0]} {v[1]}>{v[2]} delta={v[3]:+.3f}")
```

## 5. Interpreting Direction

- Positive delta: predicted increase in track signal
- Negative delta: predicted decrease in track signal
- Larger absolute values indicate stronger predicted regulatory impact

In real studies, these scores are combined with population frequency, GWAS/eQTL evidence, and experimental context.

## Optional: Real Enformer Loading (commented)

```python
# from enformer_pytorch import from_pretrained
# model = from_pretrained('EleutherAI/enformer-official-rough', target_length=-1)
# model.eval()
```

The course notebook keeps default execution lightweight while preserving the interpretation workflow.

## Summary

- Enformer-style models map long DNA windows to many molecular tracks.
- ISM is a practical and intuitive way to estimate variant impact.
- Ranking by absolute delta is a useful first-pass prioritization strategy.

## Source-backed Context

- Enformer is framed as a sequence-to-function model for gene-expression/chromatin regulatory prediction over long windows.
- Borzoi extends sequence-to-signal modeling to RNA-seq coverage at 32 bp resolution with 524 kb input windows.

## Validated Sources

Checked online during content expansion.

- [Avsec et al. 2021, Enformer (Nature Methods)](https://www.nature.com/articles/s41592-021-01252-x)
- [DeepMind Enformer implementation](https://github.com/google-deepmind/deepmind-research/tree/master/enformer)
- [Borzoi repository](https://github.com/calico/borzoi)
