---
name: spatial-transcriptomics
description: Spatial transcriptomics analysis — AnnData with spatial coordinates, Visium/Xenium QC, normalization, spatial neighbor graphs, SVG detection, deconvolution concepts, Squidpy/Scanpy
---

## When to Use

Use this skill when:
- Working with spatially-resolved gene expression data (10x Visium, Xenium, Slide-seq)
- Building spatial neighbor graphs
- Detecting spatially variable genes (SVGs)
- Performing cell-type deconvolution on bulk-like spatial spots
- Visualizing expression patterns overlaid on tissue coordinates

## Quick Reference

| Task | Tool | Key Method |
|---|---|---|
| Load Visium data | `squidpy.read.visium()` | Returns AnnData with `obsm["spatial"]` |
| Normalize | `sc.pp.normalize_total()` + `sc.pp.log1p()` | standard scRNA-seq normalization |
| Spatial neighbor graph | `sq.gr.spatial_neighbors()` | `coord_type="visium"`, `n_rings=1` |
| SVG detection | `sq.gr.spatial_autocorr()` | Moran's I; `mode="moran"` |
| Visualization | `sq.pl.spatial_scatter()` | color by gene or cluster |
| Clustering | `sc.pp.neighbors()` + `sc.tl.leiden()` | same as scRNA-seq |
| Deconvolution | RCTD, cell2location, Stereoscope | concept; external tools |

## Key Patterns

**Pattern 1: Load and inspect spatial AnnData**
```python
import squidpy as sq
import scanpy as sc

adata = sq.datasets.visium_hne_adata()  # public mouse brain dataset
print(adata)
print(adata.obsm["spatial"][:5])  # x, y pixel coordinates
print(adata.uns["spatial"])        # image and scale metadata
```

**Pattern 2: Standard preprocessing**
```python
sc.pp.filter_cells(adata, min_counts=100)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)
```

**Pattern 3: Spatial neighbor graph**
```python
sq.gr.spatial_neighbors(adata, coord_type="visium", n_rings=1)
# adata.obsp["spatial_connectivities"] is now populated
print(adata.obsp["spatial_connectivities"].shape)
```

**Pattern 4: Moran's I (spatial autocorrelation)**
```python
sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names[:100])
svgs = adata.uns["moranI"].sort_values("I", ascending=False).head(20)
print(svgs)
```

**Pattern 5: Spatial scatter plot**
```python
sq.pl.spatial_scatter(adata, color=["Hpca", "leiden"], size=1.5, img_alpha=0.5)
```

## Code Templates

**Template 1: SVG detection and visualization**
```python
import squidpy as sq, scanpy as sc, matplotlib.pyplot as plt

adata = sq.datasets.visium_hne_adata()
sc.pp.normalize_total(adata); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sq.gr.spatial_neighbors(adata, coord_type="visium")
sq.gr.spatial_autocorr(adata, mode="moran")
top_svgs = adata.uns["moranI"].head(6).index.tolist()
sq.pl.spatial_scatter(adata, color=top_svgs, ncols=3, size=1.5)
plt.suptitle("Top spatially variable genes (Moran's I)")
plt.tight_layout()
```

**Pattern 6: Spatially variable gene ranking**
```python
# Run Moran's I on all genes after neighbor graph construction
sq.gr.spatial_autocorr(adata, mode="moran", n_perms=100, n_jobs=4)
svg_df = adata.uns["moranI"].sort_values("I", ascending=False)
# Top SVGs: high Moran's I (close to 1), low p-value
top_svgs = svg_df[svg_df["pval_sim_fdr_bh"] < 0.05].head(20)
print(top_svgs[["I", "pval_sim_fdr_bh"]])
```

**Pattern 7: Clustering + spatial overlay**
```python
# Standard scRNA-seq clustering on spatial data
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5, key_added="spatial_cluster")
sc.tl.umap(adata)

# Side-by-side: UMAP and spatial
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sc.pl.umap(adata, color="spatial_cluster", ax=axes[0], show=False)
sq.pl.spatial_scatter(adata, color="spatial_cluster", ax=axes[1], size=1.5)
plt.tight_layout()
```

**Pattern 8: NMF-based deconvolution (lightweight)**
```python
from sklearn.decomposition import NMF

# NMF on spot × gene matrix to find spatially coherent expression programs
X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
nmf = NMF(n_components=10, init="nndsvda", random_state=0)
W = nmf.fit_transform(X)   # spots × programs
H = nmf.components_        # programs × genes

# Add program scores to adata
for k in range(W.shape[1]):
    adata.obs[f"NMF_{k}"] = W[:, k]

sq.pl.spatial_scatter(adata, color="NMF_0", size=1.5,
                       cmap="magma", title="NMF program 0")
```

**Pattern 9: Ligand–receptor co-expression (spatial context)**
```python
# Compute neighborhood enrichment (co-occurrence of clusters)
sq.gr.nhood_enrichment(adata, cluster_key="spatial_cluster")
sq.pl.nhood_enrichment(adata, cluster_key="spatial_cluster",
                        title="Neighborhood enrichment score")

# For full L-R analysis: squidpy's ligrec or CellChat (R) / LIANA (Python)
sq.gr.ligrec(adata, n_perms=1000, cluster_key="spatial_cluster",
             interactions_params={"resources": "CellChatDB"})
```

## Code Templates

**Template 1: SVG detection and visualization**
```python
import squidpy as sq, scanpy as sc, matplotlib.pyplot as plt

adata = sq.datasets.visium_hne_adata()
sc.pp.normalize_total(adata); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sq.gr.spatial_neighbors(adata, coord_type="visium")
sq.gr.spatial_autocorr(adata, mode="moran")
top_svgs = adata.uns["moranI"].head(6).index.tolist()
sq.pl.spatial_scatter(adata, color=top_svgs, ncols=3, size=1.5)
plt.suptitle("Top spatially variable genes (Moran's I)")
plt.tight_layout()
```

**Template 2: Full QC + clustering + spatial layout**
```python
import squidpy as sq, scanpy as sc

adata = sq.read.visium("path/to/visium/")  # or sq.datasets.*

# QC: mitochondrial fraction
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Filter low-quality spots
adata = adata[adata.obs["pct_counts_mt"] < 20]
adata = adata[adata.obs["total_counts"] > 500]

# Normalize + cluster
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_pcs=20)
sc.tl.leiden(adata, resolution=0.5)
sq.gr.spatial_neighbors(adata, coord_type="visium")

# Visualize
sq.pl.spatial_scatter(adata, color=["leiden", "total_counts"], ncols=2, size=1.5)
```

**Template 3: Cell-type deconvolution concept (RCTD)**
```python
# RCTD requires R (spacexr package) — Python interface via rpy2 or run in R
# Conceptual Python workflow using pre-computed proportions

import pandas as pd

# Assume rctd_results.csv from R RCTD run:
# columns: spot_barcode, celltype, proportion
rctd = pd.read_csv("rctd_results.csv", index_col=0)
for ct in rctd.columns:
    adata.obs[f"RCTD_{ct}"] = rctd.reindex(adata.obs_names)[ct].fillna(0)

# Visualize dominant cell type per spot
adata.obs["dominant_ct"] = rctd.reindex(adata.obs_names).idxmax(axis=1)
sq.pl.spatial_scatter(adata, color="dominant_ct", size=1.5)
```

## Common Pitfalls

- **Coordinate units:** Squidpy expects pixel coordinates; scale factor in `adata.uns["spatial"]` converts to microns
- **Spot size:** `size` parameter in `sq.pl.spatial_scatter` is relative; tune for your tissue section
- **SVG inflation:** Moran's I is sensitive to normalization; always normalize before running
- **Deconvolution assumptions:** RCTD assumes each spot is a mixture of pure cell types from a reference scRNA-seq dataset; choose a matched reference
- **Xenium vs Visium:** Xenium/MERFISH are subcellular (single-cell resolution); use `sq.read.vizgen()` or parse CSV directly — no spot-level assumptions
- **n_rings vs n_neighs:** `coord_type="visium"` uses hexagonal grid rings; `coord_type="generic"` uses k-nearest neighbors (set `n_neighs=6`)
- **Moran's I permutation test:** always use `n_perms ≥ 100`; default without perms gives only analytical p-value which is less reliable for spatial data

## Related Skills

- `ngs-variant-calling` — upstream alignment and BAM processing
- `ml-deep-learning-bio` — cell-type annotation with deep learning embeddings
- `gwas-population-genetics` — spatial eQTL concept (spatial expression vs genotype)
- `motif-discovery` — regulatory motif analysis complementary to spatial gene programs
