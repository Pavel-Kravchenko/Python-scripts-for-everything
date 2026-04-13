---
name: ai-science-geneformer-scgpt
description: "**Tier 5 — Modern AI for Science | Module 07 · Notebook 1**"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/07_Foundation_Models_Single_Cell/01_geneformer_scgpt.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Geneformer and scGPT for Single-Cell Modeling

*Source: Course notebook `Tier_5_Modern_AI_for_Science/07_Foundation_Models_Single_Cell/01_geneformer_scgpt.ipynb`*

# Geneformer and scGPT for Single-Cell Modeling

**Tier 5 — Modern AI for Science | Module 07 · Notebook 1**

*Prerequisites: Single-cell preprocessing basics, transformer fundamentals*

---

**By the end of this notebook you will be able to:**
1. Explain rank-based tokenization for single-cell profiles
2. Build toy cell embeddings from expression matrices
3. Perform simple cell-type annotation from embedding space
4. Prototype perturbation-response prediction workflow

## Why this notebook matters

Single-cell RNA sequencing generates expression profiles for thousands of individual cells, but the data is extremely sparse (most genes are not detected per cell), and each cell is measured by a different set of captured molecules. Foundation models for single-cell biology learn representations that are robust to this sparsity and that generalize across tissues, datasets, and perturbations. They are now used for cell-type annotation, perturbation response prediction, and cross-dataset integration at atlas scale.

## How to work through this notebook

1. Read the Geneformer and scGPT architecture section (Section 1) before the code — the rank-based tokenization idea is non-obvious and critical to understand.
2. The toy expression matrix (Section 2) uses only 30 genes; real datasets have 20,000+. The logic is identical.
3. The rank-token embedding (Section 3) is a faithful simplification of Geneformer's approach.
4. Section 5 (perturbation prediction) is conceptual — in real scGPT you would use learned perturbation token embeddings, not rule-based modifications.

## Common sticking points

- **Why rank instead of raw counts?** Raw count distributions are highly variable across cells (library size differences) and across batches (technology, protocol). Rank ordering is inherently normalized and comparable across cells without explicit batch correction.
- **Geneformer's rank-value encoding**: Geneformer ranks all ~20,000 genes by expression level in each cell, then encodes the top-K (K≈2048) as an ordered sequence of gene tokens. The position in the sequence encodes the rank. The model never sees raw count values — only the ordering.
- **scGPT's gene vocabulary**: scGPT uses a discrete token for each gene ID (from NCBI Entrez IDs). Expression levels are binned into discrete bins. The model then processes a sequence of (gene_token, expression_bin_token) pairs.
- **Transfer learning**: both models are pretrained on millions of cells and fine-tuned on small task-specific datasets. Fine-tuning requires fewer labeled cells than training from scratch — but always check that the pretraining data overlaps with your tissue/cell type of interest.

## 1. Geneformer and scGPT Architecture

### 1a. Geneformer

Geneformer (Theodoris et al., Nature 2023) was pretrained on ~30 million human single-cell transcriptomes. Its key innovations:

**Rank-value tokenization:**
1. For each cell, rank all ~20,000 genes by expression level (highest expression = rank 1)
2. Keep only the top-K genes (K=2048 by default) — this eliminates the dropout problem
3. Each gene is represented by a learned embedding vector (gene token)
4. The ordered sequence of gene tokens encodes the entire cell state

This is analogous to treating a cell as a "sentence" where the "words" are gene names and their "order" encodes relative expression.

**Pretraining objective:** Masked language modeling — randomly mask a fraction of gene tokens and train to predict them from the remaining ranked genes.

**What Geneformer learns:** Despite only seeing ranks, the model learns to predict masked genes correctly, which requires understanding gene regulatory relationships. The CLS token embedding provides a cell-state representation that encodes biology not captured by any individual gene.

**Fine-tuning tasks:**
- Cell-type classification (supervised, few-shot)
- Gene dosage sensitivity prediction
- In-silico perturbation (silence a gene token and see how CLS embedding shifts)
- Drug response prediction

### 1b. scGPT

scGPT (Cui et al., Nature Methods 2024) was pretrained on ~33 million cells from CellxGene. Key differences from Geneformer:

**Expression binning tokenization:**
- Gene vocabulary: ~60,000 Entrez gene IDs, each mapped to a learned embedding
- Expression is discretized into bins (e.g., 51 bins from 0 to max expression)
- Each cell is represented as a sequence of (gene_token, expression_bin_token) pairs
- The model processes both gene identity AND relative expression level

**Generative pretraining:**
- scGPT uses a masked-and-then-predict objective over gene expression
- The model is conditioned on a batch/condition token to handle multi-batch training

**scGPT perturbation module:**
- A learned perturbation embedding is injected into the sequence
- The model predicts post-perturbation expression profiles
- This enables in-silico gene knockouts and drug response simulation

### 1c. Comparison

| Feature | Geneformer | scGPT |
|---|---|---|
| Input representation | Gene ranks only | Gene identity + binned expression |
| Pretraining data | ~30M cells | ~33M cells |
| Context length | 2048 genes | variable |
| Perturbation modeling | Via embedding shift | Via explicit perturbation tokens |
| Output | CLS embedding, masked gene prediction | Expression profiles, perturbation responses |

## 2. Toy expression matrix

We create a synthetic single-cell expression matrix (cells × genes) with three cell types and marker genes. This illustrates the rank-tokenization workflow from Section 1.

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

## 3. Rank-token style embedding

Geneformer-style: represent each cell by ranked genes rather than raw count vectors. The `toy_cell_embedding` function below builds a weighted sum where earlier ranks (higher expressed genes) contribute more, using a 1/(rank+1) weight — a simplified version of Geneformer's positional encoding approach.

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

## 4. Cell-type annotation via nearest centroid

Using the rank-token embeddings, we train a nearest-centroid classifier to assign cell types. This is analogous to fine-tuning Geneformer's CLS token for classification — but simpler. The high accuracy (typically >95% here) reflects that our synthetic marker genes are strong discriminators.

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

## 5. Perturbation prediction prototype

This section illustrates the scGPT perturbation workflow conceptually. In real scGPT, a perturbation embedding (learned during pretraining on Perturb-seq datasets) is injected into the transformer, and the model predicts post-perturbation expression profiles. The rule-based modifications below are stand-ins for that learned prediction.

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

## 6. Scaling to atlas data

For real projects, use CellxGene Census and the official Geneformer/scGPT model APIs for million-cell scale workflows. The patterns demonstrated in this notebook — rank tokenization, centroid annotation, perturbation prediction — map directly to the real model APIs. Keep local prototypes like this for logic validation and QC before running at scale.

## Summary

- Rank-based representations are effective for cell-state encoding.
- Embedding-space centroid methods give strong annotation baselines.
- Perturbation prediction adds causal/extrapolative capability beyond annotation.

## Source-backed Context

- Geneformer and scGPT are both established foundation-model references for single-cell transcriptomics tasks.
- CellxGene Census is a practical large-scale resource for real-world atlas-scale integration workflows.

## Validated Sources

Checked online during content expansion.

- [Geneformer paper (Nature 2023)](https://www.nature.com/articles/s41586-023-06139-9)
- [scGPT paper (Nature Methods 2024)](https://www.nature.com/articles/s41592-024-02201-0)
- [CellxGene Census documentation](https://chanzuckerberg.github.io/cellxgene-census/)
