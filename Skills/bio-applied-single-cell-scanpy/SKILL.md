---
name: bio-applied-single-cell-scanpy
description: "Full scanpy scRNA-seq workflow: QC, normalization, HVG, PCA, UMAP, Leiden clustering, marker gene detection"
tool_type: python
primary_tool: scanpy
---

# Single-Cell RNA-seq Analysis with Scanpy

## Standard Workflow

| Step | Function | Key Parameters |
|------|----------|----------------|
| Load data | `sc.read_*()` | File format |
| QC metrics | `sc.pp.calculate_qc_metrics()` | `qc_vars=['mt']` |
| Filter cells | `sc.pp.filter_cells()` | `min_genes`, `max_genes` |
| Filter genes | `sc.pp.filter_genes()` | `min_cells` |
| Normalize | `sc.pp.normalize_total()` | `target_sum=1e4` |
| Log transform | `sc.pp.log1p()` | — |
| Find HVGs | `sc.pp.highly_variable_genes()` | `min_mean`, `min_disp` |
| Scale | `sc.pp.scale()` | `max_value=10` |
| PCA | `sc.tl.pca()` | `n_comps=50` |
| Neighbor graph | `sc.pp.neighbors()` | `n_neighbors`, `n_pcs` |
| UMAP | `sc.tl.umap()` | — |
| Cluster | `sc.tl.leiden()` | `resolution` |
| Find markers | `sc.tl.rank_genes_groups()` | `method='wilcoxon'` |

## QC Filtering

| Problem | Symptom | Filter |
|---------|---------|--------|
| Empty droplet | Very few genes | `n_genes_by_counts < 200` |
| Doublet | Unusually high gene count | `n_genes_by_counts > 5000` |
| Dying cell | High MT% | `pct_counts_mt > 20%` |

```python
import scanpy as sc

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 20, :]
```

## Normalization and HVG Selection

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
```

## Dimensionality Reduction and Clustering

```python
adata.raw = adata  # save for DE (needs unscaled data)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
```

Leiden `resolution`: 0.2-0.4 = few large clusters; 1.0-2.0 = many subtypes.

## Marker Gene Detection

```python
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)
```

### Canonical PBMC markers

```python
marker_genes = {
    'T cells (CD4+)': ['CD3D', 'CD4', 'IL7R'],
    'T cells (CD8+)': ['CD3D', 'CD8A', 'CD8B'],
    'B cells': ['CD79A', 'MS4A1', 'CD19'],
    'NK cells': ['GNLY', 'NKG7', 'NCAM1'],
    'Monocytes (CD14+)': ['CD14', 'LYZ', 'S100A8'],
    'Monocytes (FCGR3A+)': ['FCGR3A', 'MS4A7'],
    'Dendritic cells': ['FCER1A', 'CST3'],
    'Platelets': ['PPBP', 'PF4']
}
sc.pl.dotplot(adata, [g for gs in marker_genes.values() for g in gs],
              groupby='leiden', standard_scale='var')
```

## Pitfalls

- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **HVG on raw counts**: HVG selection must use normalized (not scaled) data
- **Scaling before HVG**: inflates dispersion estimates — normalize first, select HVGs, then scale
