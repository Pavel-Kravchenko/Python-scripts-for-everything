---
name: ai-science-geneformer-scgpt
description: Geneformer and scGPT for Single-Cell Modeling
tool_type: python
primary_tool: Python
---

# Geneformer and scGPT for Single-Cell Modeling

## Architecture Overview

### Geneformer (Theodoris et al., Nature 2023)

Pretrained on ~30M human single-cell transcriptomes.

**Rank-value tokenization:**
- Rank all ~20,000 genes by expression per cell (highest = rank 1)
- Keep top-K genes (K=2048) — eliminates dropout problem
- Each gene → learned embedding vector; ordered sequence encodes cell state
- Raw count values never seen — only gene ordering

**Pretraining:** Masked LM — predict masked gene tokens from remaining ranked genes.

**Fine-tuning tasks:** cell-type classification, gene dosage sensitivity, in-silico perturbation (silence gene token → observe CLS shift), drug response prediction.

**Why ranks over raw counts?** Raw counts vary across cells (library size) and batches (technology/protocol). Rank ordering is inherently normalized without explicit batch correction.

### scGPT (Cui et al., Nature Methods 2024)

Pretrained on ~33M cells from CellxGene.

**Expression binning tokenization:**
- Gene vocab: ~60,000 Entrez gene IDs → learned embeddings
- Expression discretized into bins (e.g., 51 bins)
- Each cell = sequence of (gene_token, expression_bin_token) pairs
- Conditioned on batch/condition token for multi-batch training

**Perturbation module:** learned perturbation embedding injected into sequence; predicts post-perturbation expression profiles for in-silico knockouts and drug response.

### Comparison

| Feature | Geneformer | scGPT |
|---|---|---|
| Input representation | Gene ranks only | Gene identity + binned expression |
| Pretraining data | ~30M cells | ~33M cells |
| Context length | 2048 genes | variable |
| Perturbation modeling | Via embedding shift | Via explicit perturbation tokens |
| Output | CLS embedding, masked gene prediction | Expression profiles, perturbation responses |

## Toy Expression Matrix

```python
genes = [f'G{i:02d}' for i in range(30)]
n_cells = 240

types = np.random.choice(['Tcell', 'Bcell', 'Myeloid'], size=n_cells, p=[0.35, 0.30, 0.35])
X = np.random.gamma(shape=1.5, scale=1.0, size=(n_cells, len(genes)))

markers = {
    'Tcell': [1, 2, 3],
    'Bcell': [8, 9, 10],
    'Myeloid': [15, 16, 17],
}

for i, t in enumerate(types):
    X[i, markers[t]] += np.random.uniform(2.5, 4.0)

expr = pd.DataFrame(X, columns=genes)
expr.head(3)
```

## Rank-Token Embedding

Geneformer-style: represent each cell by ranked genes; weight by 1/(rank+1) — higher expressed genes contribute more.

```python
def topk_rank_tokens(row: np.ndarray, k: int = 8):
    idx = np.argsort(row)[::-1][:k]
    return tuple(idx.tolist())

tokens = [topk_rank_tokens(expr.iloc[i].values, k=8) for i in range(n_cells)]

def toy_cell_embedding(token_tuple, n_genes=30):
    v = np.zeros(n_genes, dtype=float)
    for rank, g in enumerate(token_tuple):
        v[g] += 1.0 / (rank + 1)
    return v / (np.linalg.norm(v) + 1e-9)

E = np.vstack([toy_cell_embedding(t, n_genes=len(genes)) for t in tokens])
print('Embedding matrix:', E.shape)
```

## Cell-Type Annotation (Nearest Centroid)

Analogous to fine-tuning Geneformer's CLS token for classification.

```python
y = np.array(types)
idx = np.arange(n_cells)
np.random.shuffle(idx)
train = idx[:180]
test = idx[180:]

centroids = {c: E[train][y[train] == c].mean(axis=0) for c in np.unique(y)}

def predict_one(vec):
    d = {c: np.linalg.norm(vec - centroids[c]) for c in centroids}
    return min(d, key=d.get)

pred = np.array([predict_one(E[i]) for i in test])
acc = float((pred == y[test]).mean())
print('Annotation accuracy:', round(acc, 3))
```

## Perturbation Prediction Prototype

In real scGPT, learned perturbation token embeddings replace these rule-based modifications.

```python
def perturb_predict(cell_vec: np.ndarray, pert: str) -> np.ndarray:
    out = cell_vec.copy()
    if pert == 'KO_G02':
        out[2] *= 0.2
    elif pert == 'KO_G09':
        out[9] *= 0.2
    elif pert == 'CYTOKINE_X':
        out[[1, 3, 15]] *= 1.3
    return out

cell0 = expr.iloc[0].values
for p in ['KO_G02', 'KO_G09', 'CYTOKINE_X']:
    pred_expr = perturb_predict(cell0, p)
    print(p, 'delta_mean=', round(float((pred_expr - cell0).mean()), 4))
```

## Key Points

- Rank-based representations are effective for cell-state encoding
- Embedding-space centroid methods give strong annotation baselines
- Perturbation prediction adds causal/extrapolative capability beyond annotation
- Transfer learning requires fewer labeled cells than training from scratch — verify pretraining data overlaps with your tissue/cell type

## References

- [Geneformer paper (Nature 2023)](https://www.nature.com/articles/s41586-023-06139-9)
- [scGPT paper (Nature Methods 2024)](https://www.nature.com/articles/s41592-024-02201-0)
- [CellxGene Census documentation](https://chanzuckerberg.github.io/cellxgene-census/)

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
