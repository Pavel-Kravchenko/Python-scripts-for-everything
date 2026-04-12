# scrna-seq-analysis

Single-cell RNA-seq analysis from count matrix to annotated cell types and trajectories.

## Quick Reference

| Step | Tool | Key Parameter |
|------|------|---------------|
| Load 10x data | `sc.read_10x_mtx()` | `var_names='gene_symbols'` |
| QC filter | `sc.pp.filter_cells/genes` | `min_genes=200`, `pct_counts_mt < 20` |
| Normalize | `sc.pp.normalize_total` | `target_sum=1e4` |
| Log transform | `sc.pp.log1p` | — |
| HVG selection | `sc.pp.highly_variable_genes` | `n_top_genes=2000` |
| PCA | `sc.tl.pca` | `n_comps=50` |
| Neighbors | `sc.pp.neighbors` | `n_neighbors=10, n_pcs=40` |
| UMAP | `sc.tl.umap` | `min_dist=0.3` |
| Cluster | `sc.tl.leiden` | `resolution=0.5` |
| Markers | `sc.tl.rank_genes_groups` | `method='wilcoxon'` |

## Standard Workflow

```python
import scanpy as sc
import numpy as np

# Load
adata = sc.read_10x_mtx('data/', var_names='gene_symbols')
adata.var_names_make_unique()

# Annotate mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter
adata = adata[adata.obs.n_genes_by_counts > 200]
adata = adata[adata.obs.pct_counts_mt < 20]
sc.pp.filter_genes(adata, min_cells=3)

# Normalize and log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # store raw for DE

# HVG → Scale → PCA → Neighbors → UMAP → Leiden
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# Plot
sc.pl.umap(adata, color=['leiden', 'CD3D', 'CD79A'])

# Marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)
```

## CellTypist Automated Annotation

```python
import celltypist
from celltypist import models

models.download_models(force_update=True)
model = models.Model.load(model='Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model=model, majority_voting=True)
adata = predictions.to_adata()
sc.pl.umap(adata, color='majority_voting')
```

## RNA Velocity (scVelo)

```python
import scvelo as scv

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap')
```

## Common Pitfalls
- **Forgetting `adata.raw = adata`** before filtering HVGs — stores full gene set for DE
- **Over-clustering**: start with `resolution=0.3–0.5`, increase if subtypes expected
- **Batch effects**: always color UMAP by batch/sample before interpreting clusters
- **MT% threshold**: 5–10% for neurons, 20–25% for heart/muscle; adjust by tissue
- **Doublets**: run DoubletFinder/scDblFinder before QC filtering

## Key Formats
- **AnnData (.h5ad)**: `adata.obs` (cell metadata), `adata.var` (gene metadata), `adata.X` (count matrix)
- **10x output**: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`

## Module
Tier 3 · Module 30 (Single-Cell RNA-seq)
