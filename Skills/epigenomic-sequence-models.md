---
name: epigenomic-sequence-models
description: Borzoi/Epiformer-style sequence modeling for RNA-seq coverage and chromatin accessibility tasks.
---

# epigenomic-sequence-models

## When to Use

Use this skill when:
- Predicting expression-like or accessibility-like outputs from sequence
- Comparing modality-specific model behavior
- Designing task-appropriate epigenomic benchmarks

## Quick Reference

| Model family | Typical target |
|---|---|
| Borzoi-style | RNA-seq coverage / expression-like tracks |
| Epiformer-style | Chromatin accessibility (often with conservation features) |
| AlphaGenome-style | Unified multi-modal outputs |

## Key Patterns

**Pattern 1: modality-specific scoring**
```python
def expr_score(seq):
    return 0.2 + 0.35 * seq.count('TATAAA') + 0.08 * seq.count('GATA')

def acc_score(seq, conservation):
    return 0.15 + 0.12 * seq.count('ATAC') + 0.55 * conservation
```

**Pattern 2: simple benchmark correlation**
```python
import numpy as np

def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])
```

## Code Templates

### Cross-task mismatch check
```python
def modality_mismatch(expr_pred, acc_truth):
    proxy = (expr_pred - expr_pred.min()) / (expr_pred.max() - expr_pred.min() + 1e-9)
    return corr(proxy, acc_truth)
```

### Selection matrix
```python
MODEL_BY_TASK = {
    'rna_seq_coverage': 'Borzoi-style',
    'chromatin_accessibility': 'Epiformer-style',
    'multi_modal_variant_effect': 'AlphaGenome-style',
}
```

## Common Pitfalls

- Evaluating expression models on accessibility targets (or vice versa)
- Ignoring conservation/cell-type context in accessibility tasks
- Assuming one model family dominates every modality

## Related Skills

- `chipseq-epigenomics`
- `enformer-regulatory-prediction`
- `genomic-foundation-models`

