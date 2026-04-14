---
name: ai-science-genomic-llms
description: "Genomic Foundation Models: Nucleotide Transformers, HyenaDNA, and Evo with NumPy"
tool_type: python
primary_tool: NumPy
---

# Genomic Foundation Models: Nucleotide Transformers, HyenaDNA, and Evo

## Tokenization Strategies

- **Character-level** (`A,C,G,T,N`): highest resolution, long sequences
- **k-mer tokens** (e.g., k=6): compressed representation; k=6 → 4096-token vocab, k=8 → 65,536 — use k≤6 for explicit k-mer tokenization
- **BPE/subword**: data-driven token units (used in some genomic LMs)
- **Nucleotide Transformer**: 6-mer tokens, stride=1, 4096-vocab; ~L/6 tokens per sequence — loses single-nucleotide resolution
- **HyenaDNA**: replaces self-attention with Hyena state-space operator; handles up to 1M nucleotides without O(L²) attention cost
- **Evo**: trained on prokaryotic/viral genomes (not human) — use for microbial sequence design, not mammalian regulatory biology

```python
import numpy as np
from collections import Counter

np.random.seed(7)
```

```python
def kmers(seq: str, k: int = 6):
    seq = seq.upper()
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]

sequence = "ATGCGTACGTTAGCGTATCGATCGGATCGA"
print("Sequence length:", len(sequence))
print("First 10 6-mers:", kmers(sequence, 6)[:10])
```

## k-mer Frequency Embedding Baseline

```python
def kmer_embedding(seq: str, vocab: list[str], k: int = 3) -> np.ndarray:
    tokens = kmers(seq, k)
    counts = Counter(tokens)
    vec = np.array([counts[v] for v in vocab], dtype=float)
    return vec / (vec.sum() + 1e-9)

alphabet = ["A", "C", "G", "T"]
vocab_3 = [a + b + c for a in alphabet for b in alphabet for c in alphabet]

example_vec = kmer_embedding("ATGATGATGCCC", vocab_3, k=3)
print("Embedding dimension:", example_vec.shape[0])
print("Non-zero features:", int((example_vec > 0).sum()))
```

## Promoter Probe Task (Nearest-Centroid Classifier)

Synthetic sequences with TATA motif injection; nearest-centroid on k-mer embeddings.

```python
def random_dna(n: int) -> str:
    return "".join(np.random.choice(list("ACGT"), size=n))

def inject_motif(seq: str, motif: str, pos: int) -> str:
    return seq[:pos] + motif + seq[pos + len(motif):]

n_samples = 120
length = 80
motif = "TATAAA"

seqs, labels = [], []
for _ in range(n_samples):
    s = random_dna(length)
    if np.random.rand() < 0.5:
        s = inject_motif(s, motif, pos=20)
        labels.append(1)
    else:
        labels.append(0)
    seqs.append(s)

X = np.stack([kmer_embedding(s, vocab_3, k=3) for s in seqs])
y = np.array(labels)

train_idx = np.arange(0, 90)
test_idx = np.arange(90, n_samples)

X_train, y_train = X[train_idx], y[train_idx]
X_test, y_test = X[test_idx], y[test_idx]

c0 = X_train[y_train == 0].mean(axis=0)
c1 = X_train[y_train == 1].mean(axis=0)

d0 = ((X_test - c0) ** 2).sum(axis=1)
d1 = ((X_test - c1) ** 2).sum(axis=1)
pred = (d1 < d0).astype(int)

acc = (pred == y_test).mean()
print(f"Nearest-centroid probe accuracy: {acc:.3f}")
```

## Long-Range Context

Enhancer-promoter interactions can span 100kb+; short windows miss distal signals. HyenaDNA and similar long-context models address this.

```python
def distal_interaction_label(seq: str) -> int:
    # Toy long-range rule: motif A near start AND motif B near end
    has_a = "GATA" in seq[:120]
    has_b = "CACC" in seq[-120:]
    return int(has_a and has_b)

toy_long = random_dna(1000)
toy_long = inject_motif(toy_long, "GATA", 40)
toy_long = inject_motif(toy_long, "CACC", 920)

print("Long-range label:", distal_interaction_label(toy_long))
print("If truncated to 200 bp, end motif is lost -> label can flip.")
```

## Model Selection

| Task profile | Preferred model |
|---|---|
| Short-window promoter/enhancer classification | DNABERT-2 / Nucleotide Transformer |
| Distal regulatory context (100kb+) | HyenaDNA |
| Prokaryotic sequence generation / scoring | Evo |
| Variant effect on expression tracks | Enformer / AlphaGenome |
| Splicing-focused interpretation | SpliceAI |

Use lightweight baselines first; escalate to foundation models when context or task complexity demands it.

## Key Points

- Tokenization choice changes both context length and biological granularity
- k-mer embeddings are a strong sanity-check baseline
- Long-range rules motivate long-context architectures like HyenaDNA
- Match model family to task: embedding, regulatory prediction, splicing, or generation

## References

- [Nucleotide Transformer repository](https://github.com/instadeepai/nucleotide-transformer)
- [HyenaDNA repository](https://github.com/HazyResearch/hyena-dna)
- [Evo repository](https://github.com/evo-design/evo)

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
