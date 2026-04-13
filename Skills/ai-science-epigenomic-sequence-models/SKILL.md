---
name: ai-science-epigenomic-sequence-models
description: "**Tier 5 — Modern AI for Science | Module 05 · Notebook 4**"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/04_epigenomic_sequence_models.ipynb"
---

# Epigenomic Sequence Models: Borzoi and Epiformer

*Source: Course notebook `Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/04_epigenomic_sequence_models.ipynb`*

# Epigenomic Sequence Models: Borzoi and Epiformer

**Tier 5 — Modern AI for Science | Module 05 · Notebook 4**

*Prerequisites: Notebook 2 (Enformer), Notebook 3 (Splicing Models)*

---

**By the end of this notebook you will be able to:**
1. Distinguish expression-style and accessibility-style outputs
2. Build toy predictors for RNA-seq coverage and chromatin accessibility
3. Compare modality-specific performance across synthetic samples
4. Choose a model family based on assay target and context

## Why this notebook matters

Chromatin accessibility (open vs closed chromatin) and RNA-seq coverage are related but distinct aspects of gene regulation. A variant can disrupt a DNase hypersensitive site without changing RNA levels, or change RNA levels without altering accessibility. Different model families are optimized for different output modalities, and using the wrong model for a task produces poor results. This notebook builds intuition for matching model families to tasks.

## How to work through this notebook

1. Read Section 1 (output modalities) carefully — the distinction between RNA-seq coverage style outputs and accessibility style outputs drives all model selection decisions.
2. The synthetic benchmark (Section 2) demonstrates why you must evaluate models on the modality you care about, not a proxy modality.
3. Section 3 (cross-task mismatch) is a deliberate failure demonstration — run it to see what happens when you use the wrong model.

## Common sticking points

- **Borzoi vs Enformer outputs**: Enformer predicts 128 bp bin averages of ChIP/CAGE/DNase. Borzoi predicts RNA-seq coverage at 32 bp resolution — a much finer scale — enabling it to distinguish splicing patterns, 5'/3' UTR usage, and isoform-specific signals.
- **Epiformer's conservation feature**: Epiformer takes both the DNA sequence AND a conservation track (e.g., PhyloP scores) as input. This means you cannot use Epiformer without a precomputed conservation track for your genome region. This is different from Borzoi, which takes sequence only.
- **Why conservation helps for accessibility**: evolutionarily conserved non-coding regions are enriched for functional regulatory elements. Including conservation as a feature gives the model a strong prior for which regions are regulatory-active.
- **AlphaGenome's unified approach**: AlphaGenome (DeepMind/Google, 2025) predicts expression, splicing, chromatin, and 3D contact maps from the same model. This avoids model-selection decisions but requires more compute.

```python
import numpy as np

np.random.seed(19)

def random_dna(n: int) -> str:
    return "".join(np.random.choice(list("ACGT"), size=n))

def motif_count(seq: str, motif: str) -> int:
    return sum(1 for i in range(len(seq) - len(motif) + 1) if seq[i:i + len(motif)] == motif)
```

## 1. Output Modalities

- **Borzoi-like outputs**: dense RNA-seq coverage style signals
- **Epiformer-like outputs**: chromatin accessibility predictions with sequence + conservation features

Different modalities can disagree for the same sequence, and that disagreement is often biologically informative.

```python
def toy_borzoi_expression(seq: str) -> float:
    # expression-like score from promoter and enhancer motifs
    score = 0.2 + 0.35 * motif_count(seq, "TATAAA") + 0.08 * motif_count(seq, "GATA")
    return float(np.clip(score, 0.0, 3.0))

def toy_epiformer_accessibility(seq: str, conservation: float) -> float:
    # accessibility-like score from open-chromatin motifs + conservation
    score = 0.15 + 0.12 * motif_count(seq, "ATAC") + 0.55 * conservation
    return float(np.clip(score, 0.0, 2.5))

s = random_dna(400)
print("Expression-like:", toy_borzoi_expression(s))
print("Accessibility-like:", toy_epiformer_accessibility(s, conservation=0.62))
```

## 2. Generate Synthetic Benchmark Dataset

We create synthetic samples with sequence and conservation values, then define "observed" targets by adding controlled noise.

```python
N = 120
seqs = [random_dna(500) for _ in range(N)]
cons = np.random.uniform(0.1, 0.95, size=N)

true_expr = np.array([toy_borzoi_expression(s) for s in seqs]) + np.random.normal(0, 0.08, size=N)
true_acc = np.array([toy_epiformer_accessibility(s, c) for s, c in zip(seqs, cons)]) + np.random.normal(0, 0.06, size=N)

pred_expr = np.array([toy_borzoi_expression(s) for s in seqs])
pred_acc = np.array([toy_epiformer_accessibility(s, c) for s, c in zip(seqs, cons)])

def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])

print("Expression correlation:", round(corr(pred_expr, true_expr), 3))
print("Accessibility correlation:", round(corr(pred_acc, true_acc), 3))
```

## 3. Cross-Task Mismatch Demonstration

A model optimized for one modality can perform poorly on another. We illustrate this by trying to predict accessibility from expression scores only.

```python
bad_acc_proxy = (pred_expr - pred_expr.min()) / (pred_expr.max() - pred_expr.min() + 1e-9)
print("Accessibility correlation using wrong modality proxy:", round(corr(bad_acc_proxy, true_acc), 3))
```

## 4. Practical Model Choice Matrix

| Goal | Preferred model family |
|---|---|
| RNA-seq coverage from sequence | Borzoi-like sequence-to-coverage models |
| ATAC/chromatin accessibility prediction | Epiformer-like accessibility models |
| Broad multi-modal variant effects | AlphaGenome-like unified models |
| Splicing triage | SpliceAI-like specialized models |

Use task-specific models for primary scoring and multi-task models for context integration.

## Summary

- Epigenomic modeling is modality-specific: expression and accessibility are related but not interchangeable.
- Sequence-only and sequence+conservation models can complement each other.
- Benchmark each model on the output modality you actually need.

## Source-backed Context

- Epiformer reports sequence+conservation inputs for chromatin-accessibility prediction around ~100 kb context.
- Borzoi and Epiformer target related but distinct modalities (RNA-seq coverage vs accessibility), motivating modality-aware benchmarking.

## Validated Sources

Checked online during content expansion.

- [Borzoi repository and preprint link](https://github.com/calico/borzoi)
- [Epiformer repository](https://github.com/yal054/epiformer)
- [Epiformer-linked Science paper context](https://github.com/yal054/epiformer)
