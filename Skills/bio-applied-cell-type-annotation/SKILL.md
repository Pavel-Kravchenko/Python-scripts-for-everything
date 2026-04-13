---
name: bio-applied-cell-type-annotation
description: "**Tier 3 — Applied Bioinformatics | Module 30 · Notebook 3**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/30_Single_Cell_RNA_seq/03_cell_type_annotation.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: anndata 0.10+, matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# scRNA-seq: Cell Type Annotation

*Source: Course notebook `Tier_3_Applied_Bioinformatics/30_Single_Cell_RNA_seq/03_cell_type_annotation.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 30 · Notebook 3**

*Prerequisites: Notebook 2 (Dimensionality Reduction & Clustering)*

---

**By the end of this notebook you will be able to:**
1. Manually annotate clusters using canonical marker genes from literature
2. Use automated tools (SingleR, CellTypist) for reference-based annotation
3. Assess annotation confidence scores and resolve ambiguous clusters
4. Integrate CellMarker and PanglaoDB marker databases
5. Visualize annotated cell types and produce publication-ready UMAP figures



**Key resources:**
- [SingleR documentation](https://bioconductor.org/packages/release/bioc/html/SingleR.html)
- [CellTypist](https://www.celltypist.org/)
- [CellMarker database](http://bio-bigdata.hrbmu.edu.cn/CellMarker/)
- [PanglaoDB](https://panglaodb.se/)

## Why cell type annotation matters

Clustering algorithms find groups of similar cells, but they have no concept of biology — they output numbers (Cluster 0, 1, 2...). Annotation maps these cluster numbers to biological cell identities. The quality of annotation determines whether your biological conclusions are valid. A mislabeled cluster contaminates every downstream analysis: differential abundance tests, trajectory inference, and ligand-receptor communication analysis all depend on correct cell type labels.

```python
# Setup: recreate clustered AnnData from previous notebooks
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

try:
    adata = ad.read_h5ad('/tmp/scrna_preprocessed.h5ad')
    print(f"Loaded: {adata.shape}")
except FileNotFoundError:
    rng = np.random.default_rng(42)
    n_cells, n_genes = 500, 1500
    cell_types = ["T_cell", "B_cell", "Monocyte", "NK_cell", "Dendritic"]
    cells_per = n_cells // len(cell_types)
    labels = np.repeat(cell_types, cells_per)
    counts = rng.negative_binomial(2, 0.9, (n_cells, n_genes)).astype(np.float32)
    type_idx = {ct: np.where(labels == ct)[0] for ct in cell_types}
    offsets = [50, 200, 400, 700, 1100]
    for ct, off in zip(cell_types, offsets):
        counts[np.ix_(type_idx[ct], [off, off+1, off+2])] += rng.poisson(20, (cells_per, 3))
    cell_totals = counts.sum(axis=1, keepdims=True)
    X_norm = np.log1p(counts / cell_totals * 1e4)
    adata = ad.AnnData(
        X=X_norm,
        obs=pd.DataFrame({'cell_type': labels, 'leiden': np.repeat(['0','1','2','3','4'], cells_per)},
                         index=[f'CELL_{i:04d}' for i in range(n_cells)]),
        var=pd.DataFrame(index=[f'GENE{i:04d}' for i in range(n_genes)])
    )
    from sklearn.decomposition import PCA
    from umap import UMAP
    X_scaled = (X_norm - X_norm.mean(0)) / (X_norm.std(0) + 1e-8)
    X_scaled = np.clip(X_scaled, -10, 10)
    X_pca = PCA(n_components=15, random_state=42).fit_transform(X_scaled)
    X_umap = UMAP(n_neighbors=15, min_dist=0.3, random_state=42).fit_transform(X_pca)
    adata.obsm['X_pca'] = X_pca
    adata.obsm['X_umap'] = X_umap
    print(f"Created synthetic dataset: {adata.shape}")

coords = adata.obsm['X_umap']
print(f"Clusters: {adata.obs['leiden'].unique().tolist()}")
```python

## 1. Canonical Marker-Based Manual Annotation

Manual annotation is the gold standard for well-characterized cell types with known markers. It requires domain knowledge but is transparent, reproducible, and does not depend on a reference dataset.

### PBMC canonical markers (widely agreed upon)

| Cell Type | Marker Genes | Notes |
|-----------|-------------|-------|
| T cells | CD3D, CD3E, CD3G | Pan-T marker; present in all T subsets |
| CD4+ T cells | CD4, IL7R | Helper T cells |
| CD8+ T cells | CD8A, CD8B | Cytotoxic T cells |
| B cells | CD79A, CD79B, MS4A1 (CD20) | MS4A1 is target of rituximab |
| NK cells | GNLY, NKG7, KLRD1 | No CD3, distinguishes from T cells |
| Monocytes (classical) | LYZ, S100A8, S100A9, CST3 | High LYZ = phagocytic activity |
| Monocytes (non-classical) | FCGR3A (CD16), MS4A7 | Patrolling monocytes |
| Dendritic cells | FCER1A, CLEC10A | Low count, high HLA-DR |
| Platelets | PPBP, PF4 | Often contaminants in PBMC prep |

### Annotation workflow
1. Plot known markers as feature plots on UMAP
2. Check violin plots per cluster to see distribution of expression
3. Assign cell type label if ≥2 concordant markers strongly expressed
4. Ambiguous clusters: check for intermediate expression (transitional states) or doublet markers (two cell type programs simultaneously)

```python
# Manual annotation: define marker dictionary and score each cluster
import numpy as np
import matplotlib.pyplot as plt

# Canonical PBMC marker dictionary
# In real data, these would be actual gene names; our synthetic genes use index-based proxies
# We'll use the clusters from the previous notebook and assign based on majority composition
marker_dict = {
    'T_cell':   ['GENE0050', 'GENE0051', 'GENE0052'],   # CD3D, CD3E, CD3G
    'B_cell':   ['GENE0200', 'GENE0201'],                # CD79A, MS4A1  
    'Monocyte': ['GENE0400', 'GENE0401'],                # LYZ, S100A8
    'NK_cell':  ['GENE0700', 'GENE0701'],                # GNLY, NKG7
    'Dendritic':['GENE1100', 'GENE1101'],                # FCER1A, CLEC10A
}

X = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
gene_names = adata.var_names.tolist()

# Compute per-cluster mean expression of each marker set
clusters = adata.obs['leiden'].values
cluster_ids = sorted(adata.obs['leiden'].unique())

# Score each cluster for each cell type
scores = {}
for ct, markers in marker_dict.items():
    valid_markers = [m for m in markers if m in gene_names]
    if valid_markers:
        marker_idx = [gene_names.index(m) for m in valid_markers]
        scores[ct] = {cid: X[clusters == cid][:, marker_idx].mean() for cid in cluster_ids}

# Show scores as table
score_df = pd.DataFrame(scores, index=cluster_ids)
print("Per-cluster marker gene scores (mean log-normalized expression):")
print(score_df.round(3))

# Assign annotation: cluster gets the cell type with highest score
cluster_annotation = score_df.idxmax(axis=1).to_dict()
print(f"\nCluster annotations:")
for cid, annot in cluster_annotation.items():
    truth = adata.obs[adata.obs['leiden'] == cid]['cell_type'].value_counts().index[0]
    print(f"  Cluster {cid} -> {annot} (true majority: {truth})")

# Assign annotations back to adata
adata.obs['cell_type_annotated'] = adata.obs['leiden'].map(cluster_annotation)

# Visualize annotated UMAP
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
colors_map = {'T_cell': '#e41a1c', 'B_cell': '#377eb8', 'Monocyte': '#4daf4a',
              'NK_cell': '#984ea3', 'Dendritic': '#ff7f00'}

for ct in sorted(adata.obs['cell_type_annotated'].unique()):
    mask = adata.obs['cell_type_annotated'] == ct
    axes[0].scatter(coords[mask, 0], coords[mask, 1],
                    c=colors_map.get(ct, 'gray'), label=ct, s=10, alpha=0.8)
axes[0].set_title('Annotated cell types')
axes[0].legend(bbox_to_anchor=(1.02, 1))

for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    axes[1].scatter(coords[mask, 0], coords[mask, 1],
                    c=colors_map.get(ct, 'gray'), label=ct, s=10, alpha=0.8)
axes[1].set_title('Ground truth cell types')
axes[1].legend(bbox_to_anchor=(1.02, 1))

for ax in axes:
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
plt.tight_layout()
plt.savefig('annotation_comparison.png', dpi=100, bbox_inches='tight')
plt.show()

# Accuracy
correct = (adata.obs['cell_type_annotated'] == adata.obs['cell_type']).mean()
print(f"\nAnnotation accuracy: {correct:.1%}")
```python

## 2. SingleR Automated Reference Annotation

SingleR (Aran et al. 2019, *Nature Methods*) performs reference-based annotation entirely in R (Bioconductor). It works by:

1. Computing Spearman correlations between each test cell's expression profile and each reference cell type's profile
2. Iteratively pruning low-correlating reference cell types until a single best match is found (fine-tuning step)
3. Assigning confidence scores; cells with ambiguous scores are labeled "pruned" (low confidence)

### Available reference datasets
- `HumanPrimaryCellAtlasData()` — 37 cell types from diverse human tissues
- `MonacoImmuneData()` — 29 immune cell subtypes (excellent for PBMC)
- `BlueprintEncodeData()` — hematopoietic and epithelial lineages
- Custom references: any labeled single-cell dataset

### Running SingleR in R
```r
library(SingleR)
library(celldex)

# Load reference
ref <- MonacoImmuneData()

# Predict (test is a SingleCellExperiment or matrix)
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

# Scores heatmap
plotScoreHeatmap(pred)

# Diagnostic: how confident is each annotation?
table(pred$pruned.labels)
```python

### Key output: `pred$pruned.labels`
Cells with low confidence across all reference types are labeled `NA` in `pruned.labels` (but have a label in `labels`). Inspect these cells carefully — they may be:
- Novel cell types not in the reference
- Transitional/intermediate states
- Poor-quality cells that survived QC

## 3. CellTypist: Python-Based Automated Annotation

CellTypist (Domínguez Conde et al. 2022, *Science*) is a logistic regression classifier trained on a curated atlas of 20+ million human cells. Unlike SingleR, it runs in Python and integrates directly with scanpy.

### Installation and usage
```bash
pip install celltypist
```python

```python
import celltypist
from celltypist import models

# Download pre-trained models (run once)
models.download_models(force_update=False)

# Available models:
# - Immune_All_Low.pkl     — 36 immune cell subtypes (low hierarchy)
# - Immune_All_High.pkl    — 9 broad immune categories
# - Pan_Fetal_Human.pkl    — fetal cell types across organs

model = models.Model.load(model='Immune_All_Low.pkl')

# IMPORTANT: CellTypist expects log1p-normalized data with target_sum=1e4
# This is exactly what sc.pp.normalize_total + sc.pp.log1p produces
predictions = celltypist.annotate(adata, model=model, majority_voting=True)

# majority_voting=True: each Leiden cluster gets the plurality label
# majority_voting=False: each cell gets an independent prediction
adata = predictions.to_adata()
```python

### Majority voting mode
With `majority_voting=True`, individual cell predictions within each cluster are pooled, and the cluster gets the most common prediction. This reduces noise from individual cell misclassifications. Confidence is measured by the fraction of cells in the cluster assigned to the winning label.

### When CellTypist fails
- For non-immune cells (neurons, hepatocytes), use tissue-specific models
- For non-human species, retrain on species-specific data
- If all clusters predict the same cell type, your normalization may be incorrect (check that `target_sum=1e4` was used)

```python
# CellTypist demo: simulate the annotation probability structure
# (actual celltypist requires downloading models; this shows the output format)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Simulate CellTypist-style probability output
cell_types_ref = ['CD4+ T cells', 'CD8+ T cells', 'B cells', 'NK cells',
                   'Classical monocytes', 'Dendritic cells']
n_cells = adata.n_obs

# Generate probability matrix (rows=cells, cols=cell types)
# Based on true cell type labels
rng = np.random.default_rng(42)
probs = rng.dirichlet(np.ones(6) * 0.1, size=n_cells)  # base noise

# Add true signal
true_to_ref = {
    'T_cell':    0, 'B_cell': 2, 'Monocyte': 4,
    'NK_cell':   3, 'Dendritic': 5
}
for i, true_ct in enumerate(adata.obs['cell_type']):
    if true_ct in true_to_ref:
        ref_idx = true_to_ref[true_ct]
        probs[i, ref_idx] += rng.uniform(0.5, 0.9)

# Normalize to sum to 1
probs = probs / probs.sum(axis=1, keepdims=True)
predicted_labels = [cell_types_ref[i] for i in probs.argmax(axis=1)]
max_prob = probs.max(axis=1)

adata.obs['celltypist_pred'] = predicted_labels
adata.obs['celltypist_conf'] = max_prob

print("CellTypist prediction summary:")
print(pd.Series(predicted_labels).value_counts())
print(f"\nMedian confidence: {max_prob.median():.2f}")
print(f"Cells with confidence < 0.5: {(max_prob < 0.5).sum()} ({(max_prob < 0.5).mean():.1%})")

# Probability heatmap (cells sorted by predicted type)
sort_idx = np.argsort(probs.argmax(axis=1))
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

im = axes[0].imshow(probs[sort_idx, :], aspect='auto', cmap='YlOrRd',
                    interpolation='nearest')
axes[0].set_xticks(range(len(cell_types_ref)))
axes[0].set_xticklabels(cell_types_ref, rotation=45, ha='right', fontsize=8)
axes[0].set_ylabel('Cells (sorted by prediction)')
axes[0].set_title('CellTypist prediction probabilities')
plt.colorbar(im, ax=axes[0], label='Probability')

# UMAP colored by CellTypist prediction
colors_ct = {'CD4+ T cells':'#e41a1c', 'CD8+ T cells':'#fc8d62',
             'B cells':'#377eb8', 'NK cells':'#984ea3',
             'Classical monocytes':'#4daf4a', 'Dendritic cells':'#ff7f00'}
for ct in set(predicted_labels):
    mask = np.array(predicted_labels) == ct
    axes[1].scatter(coords[mask, 0], coords[mask, 1],
                    c=colors_ct.get(ct, 'gray'), label=ct, s=10, alpha=0.8)
axes[1].set_title('CellTypist annotations on UMAP')
axes[1].set_xlabel('UMAP1'); axes[1].set_ylabel('UMAP2')
axes[1].legend(bbox_to_anchor=(1.02, 1), fontsize=8)
plt.tight_layout()
plt.savefig('celltypist_results.png', dpi=100, bbox_inches='tight')
plt.show()
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
