---
name: ai-science-genomic-llms
description: "**Tier 5 — Modern AI for Science | Module 05 · Notebook 1**"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/01_genomic_llms.ipynb"
---

# Genomic Foundation Models: Nucleotide Transformers, HyenaDNA, and Evo

*Source: Course notebook `Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/01_genomic_llms.ipynb`*

# Genomic Foundation Models: Nucleotide Transformers, HyenaDNA, and Evo

**Tier 5 — Modern AI for Science | Module 05 · Notebook 1**

*Prerequisites: Tier 3 Module 10 (Deep Learning), Tier 5 Module 01 (LLM Fine-tuning)*

---

**By the end of this notebook you will be able to:**
1. Explain tokenization choices for DNA language modeling
2. Build compact sequence embeddings from k-mer features
3. Evaluate embedding quality with a simple probe task
4. Decide when to use short-context transformers vs long-context models
5. Position NT, HyenaDNA, Evo, DNABERT-2 in practical workflows

## Why this notebook matters

The human genome is 3 billion base pairs of DNA — but the language of regulatory biology is written in much shorter motifs and their long-range combinations. Genomic foundation models learn this language by pretraining on DNA sequences from thousands of species, enabling transfer to diverse downstream tasks without task-specific feature engineering. Understanding how these models tokenize DNA, what context lengths they support, and what tasks each model family excels at is essential for applying them correctly.

## How to work through this notebook

1. Read the tokenization section (Section 1) carefully — the k-mer vs character-level tradeoff is a recurring decision point in genomic ML.
2. Build and run the k-mer embedding baseline (Section 2) before moving to foundation models. This gives you a concrete sense of what 64-dimensional 3-mer features capture.
3. The promoter probe task (Section 3) is a minimal example of the standard fine-tuning evaluation workflow.
4. Section 4 (long context) motivates HyenaDNA and similar models — run it to internalize why context length matters biologically.

## Common sticking points

- **k-mer vocabulary explosion**: for k=6, there are 4⁶ = 4096 possible k-mers. For k=8, there are 65,536. This is why most models use k≤6 for explicit k-mer tokenization, and character-level or BPE tokenization at larger scales.
- **Nucleotide Transformer tokenization**: NT uses 6-mer tokens with a stride of 1, giving a 4096-token vocabulary. A sequence of length L uses approximately L/6 tokens — much shorter than character-level, but loses single-nucleotide resolution.
- **HyenaDNA's long context**: HyenaDNA replaces self-attention with a state-space model (Hyena operator) to handle sequences up to 1 million nucleotides without the O(L²) cost of standard attention. This is essential for capturing enhancer-promoter interactions at 100kb+ distances.
- **Evo's scope**: Evo is trained on prokaryotic and viral genomes (not human). Use it for tasks involving microbial sequence design, not mammalian regulatory biology.

```python
import numpy as np
from collections import Counter

np.random.seed(7)
```

## 1. Tokenization in Genomic LMs

DNA models use different tokenization strategies:
- **Character-level** (`A,C,G,T,N`): highest resolution, long sequences
- **k-mer tokens** (e.g., k=6): compressed sequence representation
- **BPE/subword**: data-driven token units (used in some genomic LMs)

k-mers are often a strong baseline for lightweight analysis and sanity checks.

```python
def kmers(seq: str, k: int = 6):
    seq = seq.upper()
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]

sequence = "ATGCGTACGTTAGCGTATCGATCGGATCGA"
print("Sequence length:", len(sequence))
print("First 10 6-mers:", kmers(sequence, 6)[:10])
```

## 2. Build a Simple Embedding Baseline (k-mer frequency)

Before using large pretrained models, it is useful to benchmark against a compact embedding baseline.
Here we use normalized k-mer frequencies as sequence embeddings.

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

## 3. Probe Task: Promoter-like vs Background Sequences

We create synthetic sequences where promoter-like samples include a TATA motif.
Then we train a tiny nearest-centroid classifier using only k-mer embeddings.

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

## 4. Why Long Context Matters

Some regulatory logic depends on distal sequence elements. Short windows can miss interactions.
HyenaDNA and similar long-context models are useful when signals are distributed over large spans.

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

## 5. Practical Model Selection

| Task profile | Good first model |
|---|---|
| Short-window promoter/enhancer classification | DNABERT-2 / Nucleotide Transformer |
| Distal regulatory context (100kb+) | HyenaDNA |
| Prokaryotic sequence generation / scoring | Evo |
| Variant effect on expression tracks | Enformer / AlphaGenome |
| Splicing-focused interpretation | SpliceAI |

Use lightweight baselines first, then escalate to foundation models when context or task complexity demands it.

## Optional: Real Model Loading (commented for portability)

```python
# from transformers import AutoTokenizer, AutoModel
# model_name = 'InstaDeepAI/nucleotide-transformer-v2-500m-multi-species'
# tokenizer = AutoTokenizer.from_pretrained(model_name)
# model = AutoModel.from_pretrained(model_name).eval()
```

In this course version we keep runnable cells CPU-friendly and dependency-light.

## Summary

- Tokenization choice changes both context length and biological granularity.
- k-mer embeddings are a strong sanity-check baseline.
- Long-range rules motivate long-context architectures like HyenaDNA.
- Match model family to task: embedding, regulatory prediction, splicing, or generation.

## Source-backed Context

- Nucleotide Transformer is maintained as a genomics foundation-model hub by InstaDeep.
- HyenaDNA emphasizes long-context nucleotide modeling up to ~1M tokens.
- Evo reports long-context genome-scale modeling and design in prokaryotic settings.

## Validated Sources

Checked online during content expansion.

- [Nucleotide Transformer repository](https://github.com/instadeepai/nucleotide-transformer)
- [HyenaDNA repository](https://github.com/HazyResearch/hyena-dna)
- [Evo repository](https://github.com/evo-design/evo)
