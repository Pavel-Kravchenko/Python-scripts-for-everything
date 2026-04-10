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

## Common Pitfalls

- **Coordinate units:** Squidpy expects pixel coordinates; scale factor in `adata.uns["spatial"]` converts to microns
- **Spot size:** `size` parameter in `sq.pl.spatial_scatter` is relative; tune for your tissue section
- **SVG inflation:** Moran's I is sensitive to normalization; always normalize before running
- **Deconvolution assumptions:** RCTD assumes each spot is a mixture of pure cell types from a reference scRNA-seq dataset; choose a matched reference

## Related Skills

- `ngs-variant-calling` — upstream alignment and BAM processing
- `ml-deep-learning-bio` — cell-type annotation with deep learning embeddings
- `gwas-population-genetics` — spatial eQTL concept (spatial expression vs genotype)
