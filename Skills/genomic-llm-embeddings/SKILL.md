---
name: genomic-llm-embeddings
description: DNA tokenization strategies, k-mer baselines, and pretrained sequence embeddings for genomic ML.
tool_type: python
primary_tool: NumPy
---

# genomic-llm-embeddings

## Tokenization Strategy Selection

| Strategy | Context | Strength | Limitation |
|---|---:|---|---|
| Character (A/C/G/T/N) | Long | Max resolution | Long token sequences |
| k-mer (k=3..6) | Short/medium | Fast baseline | Loses positional nuance |
| DNABERT-2 (BPE) | ~512 bp windows | Strong short-window tasks | Limited context |
| Nucleotide Transformer (6-mer) | kb-scale | Good transfer embeddings | Higher memory |
| HyenaDNA | up to 1M bp | Long-range signals | Heavier training/inference |

## Key Patterns

**k-mer embedding baseline**
```python
from collections import Counter
import numpy as np

def kmer_embedding(seq, vocab, k=3):
    tokens = [seq[i:i+k] for i in range(len(seq)-k+1)]
    cnt = Counter(tokens)
    vec = np.array([cnt[v] for v in vocab], dtype=float)
    return vec / (vec.sum() + 1e-9)
```

**Nearest centroid probe (sanity check before deep models)**
```python
def nearest_centroid_predict(X_train, y_train, X_test):
    c0 = X_train[y_train == 0].mean(axis=0)
    c1 = X_train[y_train == 1].mean(axis=0)
    d0 = ((X_test - c0) ** 2).sum(axis=1)
    d1 = ((X_test - c1) ** 2).sum(axis=1)
    return (d1 < d0).astype(int)
```

## Pitfalls

- Mixing tokenization schemes between train and inference silently destroys performance
- Always compare against a simple k-mer baseline -- hard to detect overfitting/bugs without one
- Match sequence window lengths when comparing models

## Related Skills

- `genomic-foundation-models` -- pretrained model details and fine-tuning
- `protein-language-models` -- ESM2 embeddings for protein sequences
