---
name: enformer-regulatory-prediction
description: Enformer-style sequence-to-track prediction, in-silico mutagenesis, and regulatory variant ranking.
---

# enformer-regulatory-prediction

## When to Use

Use this skill when:
- Predicting regulatory track changes from DNA sequence
- Running in-silico mutagenesis (ISM) for SNP effect estimates
- Prioritizing variants by direction and magnitude of predicted impact

## Quick Reference

| Concept | Value |
|---|---|
| Enformer canonical input | 196,608 bp |
| Typical output | Thousands of tracks at binned resolution |
| Variant effect primitive | `prediction(mut) - prediction(ref)` |
| Fast first-pass ranking | sort by `abs(delta)` |

## Key Patterns

**Pattern 1: one-base mutation helper**
```python
def mutate_base(seq, pos, alt):
    return seq[:pos] + alt + seq[pos+1:]
```

**Pattern 2: ISM deltas for all non-reference alleles**
```python
def ism_deltas(seq, pos, model_fn):
    ref = seq[pos]
    baseline = model_fn(seq)
    out = {}
    for alt in 'ACGT':
        if alt == ref:
            continue
        out[alt] = model_fn(mutate_base(seq, pos, alt)) - baseline
    return out
```

## Code Templates

### Minimal toy sequence-to-track model
```python
import numpy as np

def toy_regulatory_model(seq):
    promoter = 0.4 + 0.25 * seq.count('TATAAA')
    enhancer = 0.2 + 0.10 * seq.count('GATA')
    repressor = 1.0 - 0.15 * seq.count('CGCG')
    return np.clip(np.array([promoter, enhancer, repressor]), 0, 2.0)
```

### Optional real model load (commented)
```python
# from enformer_pytorch import from_pretrained
# model = from_pretrained('EleutherAI/enformer-official-rough', target_length=-1)
# model.eval()
```

## Common Pitfalls

- Using wrong window length/alignment around the variant
- Interpreting small deltas without replicate/context filters
- Treating track deltas as causal without orthogonal evidence

## Related Skills

- `genomic-llm-embeddings`
- `splicing-variant-models`
- `genomic-foundation-models`

