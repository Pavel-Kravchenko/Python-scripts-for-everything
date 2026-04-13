---
name: genomic-llm-embeddings
description: DNA tokenization, k-mer baselines, and sequence embeddings for genomic foundation model workflows.
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, transformers 4.38+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# genomic-llm-embeddings

## When to Use

Use this skill when:
- Building genomic sequence embeddings for classification or retrieval
- Choosing between k-mer, character-level, and pretrained tokenization
- Creating lightweight baselines before large-model inference

## Quick Reference

| Strategy | Typical Context | Strength | Limitation |
|---|---:|---|---|
| Character tokens (A/C/G/T/N) | Long | Max resolution | Long sequences/tokens |
| k-mer (k=3..6) | Short/medium | Fast baseline | Loses positional nuance |
| DNABERT-2 | ~512 bp windows | Strong short-window tasks | Limited context |
| Nucleotide Transformer | kb-scale | Good transfer embeddings | Higher memory |
| HyenaDNA | up to 1M bp | Long-range signals | Heavier training/inference |

## Key Patterns

**Pattern 1: k-mer embedding baseline**
```python
from collections import Counter
import numpy as np

def kmer_embedding(seq, vocab, k=3):
    tokens = [seq[i:i+k] for i in range(len(seq)-k+1)]
    cnt = Counter(tokens)
    vec = np.array([cnt[v] for v in vocab], dtype=float)
    return vec / (vec.sum() + 1e-9)
```

**Pattern 2: quick probe with nearest centroid**
```python
def nearest_centroid_predict(X_train, y_train, X_test):
    c0 = X_train[y_train == 0].mean(axis=0)
    c1 = X_train[y_train == 1].mean(axis=0)
    d0 = ((X_test - c0) ** 2).sum(axis=1)
    d1 = ((X_test - c1) ** 2).sum(axis=1)
    return (d1 < d0).astype(int)
```

## Code Templates

### Build 3-mer vocabulary
```python
alphabet = ['A', 'C', 'G', 'T']
vocab3 = [a+b+c for a in alphabet for b in alphabet for c in alphabet]
```

### Optional pretrained embedding load (commented)
```python
# from transformers import AutoTokenizer, AutoModel
# model_name = 'InstaDeepAI/nucleotide-transformer-v2-500m-multi-species'
# tok = AutoTokenizer.from_pretrained(model_name)
# model = AutoModel.from_pretrained(model_name).eval()
```

## Common Pitfalls

- Mixing tokenization schemes between train and inference
- Comparing models without matched sequence windows
- Skipping a simple baseline (hard to detect overfitting/bugs)

## Related Skills

- `genomic-foundation-models`
- `enformer-regulatory-prediction`
- `ml-deep-learning-bio`

