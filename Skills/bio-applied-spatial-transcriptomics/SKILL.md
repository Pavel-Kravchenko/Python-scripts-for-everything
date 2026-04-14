---
name: bio-applied-spatial-transcriptomics
description: "Spatial transcriptomics: Visium data loading, spatial QC, spatial neighbor graphs, Moran's I for spatially variable genes, tissue visualization with Squidpy"
tool_type: python
primary_tool: NumPy
---

# Spatial Transcriptomics

## Platform Comparison

| Platform | Resolution | Type | Key feature |
|---|---|---|---|
| 10x Visium | 55 um (~10-20 cells/spot) | Capture-based | H&E co-registration |
| 10x Xenium | ~10 um (single-cell) | In situ sequencing | Targeted ~400 genes |
| MERFISH | Sub-cellular | In situ imaging | Error-robust barcoding |
| Slide-seq v2 | ~10 um | Capture-based | High-res bead array |
| Stereo-seq | 500 nm | Capture-based | Ultra-high resolution |

## AnnData Spatial Slots

- `adata.obsm["spatial"]` — (n_spots, 2) pixel coordinates
- `adata.uns["spatial"]` — tissue image and scale factors
- `adata.obsp["spatial_connectivities"]` — spatial neighbor graph

## Standard Workflow

```python
import squidpy as sq
import scanpy as sc

adata = sq.datasets.visium_hne_adata()

# QC + normalization (same as scRNA-seq)
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.filter_cells(adata, min_counts=200)
sc.pp.filter_genes(adata, min_cells=5)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Clustering
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5)

# Spatial neighbor graph + autocorrelation
sq.gr.spatial_neighbors(adata, coord_type="visium", n_rings=1)
sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names[:500], n_perms=100)
moran = adata.uns["moranI"].sort_values("I", ascending=False)

# Visualize
sq.pl.spatial_scatter(adata, color="leiden", size=1.4)
sq.pl.spatial_scatter(adata, color=moran.head(6).index.tolist(), ncols=3, img_alpha=0.4)
```

## Moran's I (Spatial Autocorrelation)

- **I -> +1**: strong positive autocorrelation (spatial clusters)
- **I -> 0**: random spatial pattern
- **I -> -1**: checkerboard pattern

Significance via permutation testing (`n_perms=100` approximate; `n_perms=1000` for publication).

## Cell-type Deconvolution (Visium)

Each spot captures ~10-20 cells. Deconvolution estimates cell-type composition per spot.

| Tool | Approach |
|---|---|
| RCTD | Poisson GLM |
| cell2location | Negative binomial + hierarchical |
| Stereoscope | Probabilistic |

## Pipeline Summary

| Step | Tool | Output |
|---|---|---|
| Load | `sq.datasets.visium_hne_adata()` | AnnData with spatial coords |
| QC | `sc.pp.filter_cells/genes` | Clean matrix |
| Normalize | `sc.pp.normalize_total + log1p` | Log-normalized expression |
| Cluster | `sc.tl.leiden` | Spot clusters |
| Spatial graph | `sq.gr.spatial_neighbors` | `obsp["spatial_connectivities"]` |
| SVG detection | `sq.gr.spatial_autocorr(mode="moran")` | Moran's I per gene |
| Visualization | `sq.pl.spatial_scatter` | Expression on tissue image |

## Pitfalls

- **Batch effects**: check for batch confounding before interpreting spatial patterns
- **Permutation count**: n_perms=100 gives approximate p-values only; use 1000 for publication
- **Deconvolution requires scRNA-seq reference**: results depend heavily on reference quality
