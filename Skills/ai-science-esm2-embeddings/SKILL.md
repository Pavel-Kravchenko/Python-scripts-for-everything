---
name: ai-science-esm2-embeddings
description: ESM2 Embeddings and ESMFold with NumPy
tool_type: python
primary_tool: NumPy
---

# ESM2 Embeddings and ESMFold

## How to work through this notebook

1. Read the ESM2 architecture section (Section 1) to understand what information the embeddings encode.
2. The toy embedding baseline (Section 2) demonstrates that even amino acid composition gives useful signal — ESM2 embeddings are much richer.
3. The real ESM2 loading pattern (Section 3) is commented for portability; run it when you have a GPU available.
4. ESMFold (Section 4) shows the structure prediction API pattern.

## Common sticking points

- **Mean pooling vs per-residue embeddings**: mean pooling over all residues gives a single vector per protein (use for sequence-level tasks). Per-residue embeddings (one per amino acid) are used for residue-level tasks like variant scoring and annotation.
- **Which layer to extract from?** Later layers encode more abstract biological function; earlier layers encode more local sequence context. For most downstream tasks, the last or second-to-last layer works well. For structural prediction tasks, intermediate layers can be better.
- **ESMFold vs AF2**: ESMFold is ~10–60× faster (no MSA needed) but typically less accurate than AF2, especially for proteins with many homologs in databases. Use ESMFold for rapid screening and AF2 for final structural analysis.

```python
import numpy as np

np.random.seed(9)
```python

## ESM2 Architecture

ESM2 (Lin et al., 2023, Science) is a family of transformer-based protein language models trained on UniRef50 (~250M sequences). It uses the standard **masked language model (MLM)** pretraining objective — randomly masking ~15% of amino acids and training the model to predict them from context.

### Model variants

| Model | Layers | Embedding dim | Parameters |
|---|---|---|---|
| ESM2-8M | 6 | 320 | 8M |
| ESM2-35M | 12 | 480 | 35M |
| ESM2-150M | 30 | 640 | 150M |
| ESM2-650M | 33 | 1280 | 650M |
| ESM2-3B | 36 | 2560 | 3B |
| ESM2-15B | 48 | 5120 | 15B |

For most bioinformatics tasks, ESM2-650M is the sweet spot: large enough to capture complex structural signals, small enough to run on a single GPU.

### Tokenization

ESM2 uses a 33-token vocabulary:
- 20 standard amino acids
- 4 non-standard (B, U, Z, O)
- Special tokens: `<cls>`, `<eos>`, `<pad>`, `<mask>`, `<unk>`
- NO k-mer or subword tokenization — every amino acid is one token

This means a 500-residue protein uses 502 tokens (sequence + `<cls>` + `<eos>`), and the attention matrix is 502×502.

### What embeddings encode

Despite never seeing structural annotations during training, ESM2 representations:
- Predict secondary structure (helix/strand/coil) at ~80% accuracy from linear probes
- Predict contact maps (residue pairs in contact) from attention patterns
- Encode enzymatic function: sequences from the same EC number cluster together
- Capture evolutionary distance: sequence similarity correlates with embedding cosine similarity

### Real ESM2 usage pattern (requires GPU + `esm` package)

```python
# pip install fair-esm

import esm
import torch

# Load ESM2-650M
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()

data = [
    ("protein1", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"),
    ("protein2", "ACDEFGHIKLMNPQRSTVWY"),
]
batch_labels, batch_strs, batch_tokens = batch_converter(data)

# Extract representations from the last layer
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=False)

# Per-residue embeddings: shape (batch, seq_len, 1280)
token_representations = results["representations"][33]

# Mean-pool over sequence positions (exclude cls and eos)
sequence_representations = []
for i, (_, seq) in enumerate(data):
    seq_emb = token_representations[i, 1 : len(seq) + 1].mean(0)
    sequence_representations.append(seq_emb)
```python

### ESMFold: structure from sequence without MSA

ESMFold wraps ESM2-3B with a folding trunk to predict 3D structures. It achieves AlphaFold2-level accuracy on easy targets and is ~60× faster because it skips the MSA search step.

```python
# pip install fair-esm

import esm
model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()

sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
with torch.no_grad():
    output = model.infer_pdb(sequence)

# output is a PDB string
with open("predicted.pdb", "w") as f:
    f.write(output)
```python

The output PDB B-factor column contains pLDDT values (0–100), interpreted identically to AlphaFold2.

```python
AA = list('ACDEFGHIKLMNPQRSTVWY')
AA_INDEX = {a: i for i, a in enumerate(AA)}

def toy_embed(seq: str) -> np.ndarray:
    # 20-aa composition + simple length feature
    vec = np.zeros(21, dtype=float)
    for a in seq:
        if a in AA_INDEX:
            vec[AA_INDEX[a]] += 1
    if len(seq) > 0:
        vec[:20] /= len(seq)
    vec[20] = min(len(seq) / 500.0, 1.0)
    return vec

seq = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ'
emb = toy_embed(seq)
print('Embedding shape:', emb.shape)
print('Top composition entries:', np.argsort(emb[:20])[-5:][::-1])
```python

## Amino Acid Composition Baseline

Before using ESM2, it is always worth benchmarking against simple baselines. The `toy_embed` function below computes normalized amino acid composition — a 20-dimensional vector that captures global sequence statistics but ignores position and context.

ESM2 embeddings (1280 dimensions, context-aware) should substantially outperform this baseline on any real task. If they do not, consider whether the task is actually sequence-context-dependent or just composition-dependent.

```python
def random_protein(n: int) -> str:
    return ''.join(np.random.choice(AA, size=n))

def inject_motif(seq: str, motif: str, pos: int) -> str:
    return seq[:pos] + motif + seq[pos+len(motif):]

N = 200
motif = 'HGH'
seqs, y = [], []
for _ in range(N):
    s = random_protein(80)
    if np.random.rand() < 0.5:
        s = inject_motif(s, motif, 25)
        y.append(1)
    else:
        y.append(0)
    seqs.append(s)

X = np.vstack([toy_embed(s) for s in seqs])
y = np.array(y)

train, test = np.arange(150), np.arange(150, N)
c0 = X[train][y[train] == 0].mean(axis=0)
c1 = X[train][y[train] == 1].mean(axis=0)
pred = (((X[test]-c1)**2).sum(axis=1) < ((X[test]-c0)**2).sum(axis=1)).astype(int)
print('Probe accuracy:', float((pred == y[test]).mean()))
```python

## ESMFold Confidence Interpretation

The pLDDT thresholds for ESMFold output are the same as for AlphaFold2. The confidence_bucket function below encodes the standard interpretation:

```python
# Optional real API usage (network required)
# import requests
# sequence  'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ'
# r  requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', datasequence)
# with open('esmfold_prediction.pdb', 'w') as f:
#     f.write(r.text)

def confidence_bucket(mean_plddt: float) -> str:
    if mean_plddt >= 90:
        return 'very_high'
    if mean_plddt >= 70:
        return 'usable'
    if mean_plddt >= 50:
        return 'low'
    return 'very_low'

for v in [96, 82, 63, 41]:
    print(v, confidence_bucket(v))
```python

## Practical Guidance

- Use embeddings for large-scale screening and annotation tasks.
- Use structure prediction when mechanism likely depends on 3D context.
- Keep confidence-aware filters to avoid overinterpretation.

## Summary

- Embeddings provide compact sequence representations for downstream ML.
- ESMFold offers fast structure hypotheses.
- Combining embedding and structure workflows is often more robust than either alone.

## Source-backed Context

- ESM repository documentation groups ESM-2, ESMFold, ESM-1v, and ESM-IF1 as a coherent protein LM ecosystem.
- ESMFold is used in practice as a fast MSA-free structural hypothesis generator.

## Validated Sources

Checked online during content expansion.

- [ESM official repository](https://github.com/facebookresearch/esm)
- [ESM Atlas](https://esmatlas.com)

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
