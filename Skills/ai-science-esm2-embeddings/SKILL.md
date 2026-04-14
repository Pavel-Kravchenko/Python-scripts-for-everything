---
name: ai-science-esm2-embeddings
description: ESM2 Embeddings and ESMFold with NumPy
tool_type: python
primary_tool: NumPy
---

# ESM2 Embeddings and ESMFold

## Usage Notes

- **Mean pooling**: pool over all residues → single vector per protein (sequence-level tasks)
- **Per-residue**: one vector per AA → residue-level tasks (variant scoring, annotation)
- **Layer selection**: last/second-to-last for most tasks; intermediate layers can improve structural tasks
- **ESMFold vs AF2**: ESMFold ~60× faster (no MSA), but less accurate for proteins with many homologs — use for rapid screening, AF2 for final analysis

```python
import numpy as np

np.random.seed(9)
```

## ESM2 Architecture

ESM2 (Lin et al., 2023, Science) — transformer protein LM trained on UniRef50 (~250M sequences) with MLM objective (~15% masking).

### Model Variants

| Model | Layers | Embedding dim | Parameters |
|---|---|---|---|
| ESM2-8M | 6 | 320 | 8M |
| ESM2-35M | 12 | 480 | 35M |
| ESM2-150M | 30 | 640 | 150M |
| ESM2-650M | 33 | 1280 | 650M |
| ESM2-3B | 36 | 2560 | 3B |
| ESM2-15B | 48 | 5120 | 15B |

ESM2-650M is the sweet spot for most tasks: strong signal, single-GPU feasible.

### Tokenization

- 33-token vocabulary: 20 standard AAs + 4 non-standard (B, U, Z, O) + special tokens (`<cls>`, `<eos>`, `<pad>`, `<mask>`, `<unk>`)
- No k-mer or subword tokenization — one token per amino acid
- 500-residue protein → 502 tokens, attention matrix 502×502

### What Embeddings Encode

Despite no structural annotations during training:
- Secondary structure ~80% accuracy from linear probes
- Contact map prediction from attention patterns
- Enzymatic function: same EC number → cluster together
- Evolutionary distance correlates with embedding cosine similarity

### Real ESM2 Usage (requires GPU + `esm`)

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
```

## ESMFold

Wraps ESM2-3B with a folding trunk. ~60× faster than AF2 (no MSA), comparable accuracy on easy targets.

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
```

Output PDB B-factor column = pLDDT (0–100, same scale as AF2).

## Amino Acid Composition Baseline

Toy 20-dim composition embedding — always benchmark against this before using ESM2. If ESM2 doesn't outperform it, the task may be composition-dependent, not context-dependent.

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
```

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
```

## pLDDT Confidence Buckets

```python
# Optional REST API (no GPU):
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
```

## References

- [ESM official repository](https://github.com/facebookresearch/esm)
- [ESM Atlas](https://esmatlas.com)

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
