---
name: bio-applied-spatial-transcriptomics
description: "**By the end you will be able to:** - Load spatial AnnData with coordinates - Perform QC and normalization for spatial data - Build spatial neighbor graphs - Detect spatially variable genes (Moran's I"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/20_Spatial_Transcriptomics/20_spatial_transcriptomics.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scanpy 1.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 20: Spatial Transcriptomics

*Source: Course notebook `Tier_3_Applied_Bioinformatics/20_Spatial_Transcriptomics/20_spatial_transcriptomics.ipynb`*

**Tier 3 — Applied Bioinformatics | Module 20**
Prerequisites: Module 03 (RNA-seq), Module 07 (ML/Clustering), Module 12 (Single-cell Scanpy).

---

**By the end you will be able to:**
- Load spatial AnnData with coordinates
- Perform QC and normalization for spatial data
- Build spatial neighbor graphs
- Detect spatially variable genes (Moran's I)
- Visualize expression overlaid on tissue space

**Attribution:** *Patterns inspired by NGSchool 2023 spatial transcriptomics practical. Uses public 10x Visium mouse brain dataset via Squidpy.*

## Background: Spatial Transcriptomics

Unlike bulk or single-cell RNA-seq, spatial methods preserve the physical location of each measurement on a tissue section.

**Key platforms:**

| Platform | Resolution | Throughput | Type | Key feature |
|---|---|---|---|---|
| 10x Visium | 55 µm spot (~10–20 cells) | ~5,000 spots | Capture-based | H&E image co-registration |
| 10x Xenium | Single-cell (~10 µm) | ~200,000–1M+ cells | In situ sequencing | Targeted gene panel |
| MERFISH | Sub-cellular | ~100,000+ cells | In situ imaging | Error-robust barcoding of 100–10,000 genes |
| seqFISH+ | Sub-cellular | ~10,000 cells | In situ imaging | Sequential hybridization, ~10,000 genes |
| Slide-seq v2 | ~10 µm | ~50,000 beads | Capture-based | High spatial resolution bead array |
| Stereo-seq | 500 nm (DNA nanoball) | Whole tissue | Capture-based | Ultra-high resolution |

**Platform choice:**
- **Visium**: best for whole-tissue overview; well-supported by Squidpy/Seurat; H&E co-staining available
- **MERFISH/seqFISH**: subcellular resolution, in situ detection without capture bias; probe design limits gene panel size
- **Xenium**: targeted but up to ~400 genes; single-cell resolution without dissociation

**Moran's I — spatial autocorrelation:**
Moran's I measures whether nearby spots tend to have similar expression:
$$I = \frac{n}{\sum_{i}\sum_{j} w_{ij}} \cdot \frac{\sum_{i}\sum_{j} w_{ij}(x_i - \bar{x})(x_j - \bar{x})}{\sum_{i}(x_i - \bar{x})^2}$$
where $w_{ij}$ is the spatial weight between spots $i$ and $j$ (1 if neighbors, 0 otherwise) and $x_i$ is expression at spot $i$.
- **I → +1**: strong positive spatial autocorrelation (clusters of high/low expression)
- **I → 0**: random spatial pattern
- **I → −1**: perfect negative autocorrelation (checkerboard)
Spatially variable genes (SVGs) have high |I| with a significant permutation-based p-value.

**AnnData for spatial data:**
- `adata.obs` — spot/cell metadata (tissue_section, in_tissue flag)
- `adata.obsm["spatial"]` — (n_spots, 2) pixel coordinates
- `adata.uns["spatial"]` — tissue image and scale factor metadata
- `adata.X` — count matrix (n_spots × n_genes)

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

try:
    import squidpy as sq
    import scanpy as sc
    print(f"squidpy {sq.__version__}, scanpy {sc.__version__}")
except ImportError:
    print("Install: pip install squidpy scanpy")

plt.rcParams.update({"figure.dpi":110})
sc.settings.verbosity = 0
```python

```python
# Load public mouse brain H&E Visium dataset (shipped with Squidpy)
adata = sq.datasets.visium_hne_adata()
print(adata)
print(f"\nSpot coordinates shape: {adata.obsm['spatial'].shape}")
print(f"Genes: {adata.n_vars}")
print(f"\nobs columns: {adata.obs.columns.tolist()}")

# Show the tissue image
fig, ax = plt.subplots(figsize=(5,5))
sq.pl.spatial_scatter(adata, ax=ax, size=0.8, title="Raw spots on tissue")
plt.tight_layout(); plt.show()
```python

```python
# QC metrics
sc.pp.calculate_qc_metrics(adata, inplace=True)
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

fig, axes = plt.subplots(1, 3, figsize=(14,4))
axes[0].hist(adata.obs["total_counts"], bins=40, color="steelblue")
axes[0].set_xlabel("Total counts"); axes[0].set_title("Counts per spot")
axes[1].hist(adata.obs["n_genes_by_counts"], bins=40, color="steelblue")
axes[1].set_xlabel("Genes detected"); axes[1].set_title("Genes per spot")
axes[2].hist(adata.obs["pct_counts_mt"], bins=40, color="salmon")
axes[2].set_xlabel("% MT counts"); axes[2].set_title("Mitochondrial fraction")
plt.tight_layout(); plt.show()

# Filter
sc.pp.filter_cells(adata, min_counts=200)
sc.pp.filter_genes(adata, min_cells=5)
print(f"After QC: {adata.n_obs} spots, {adata.n_vars} genes")
```python

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)
print(f"Highly variable genes: {adata.var.highly_variable.sum()}")

sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0.5")
sc.tl.umap(adata)

fig, axes = plt.subplots(1, 2, figsize=(12,5))
sc.pl.umap(adata, color="leiden_0.5", ax=axes[0], show=False, title="UMAP clusters")
sq.pl.spatial_scatter(adata, color="leiden_0.5", ax=axes[1], size=1.4, title="Spatial clusters")
plt.tight_layout(); plt.show()
```python

```python
sq.gr.spatial_neighbors(adata, coord_type="visium", n_rings=1)
print("Spatial graph built.")
print(f"Connectivity matrix shape: {adata.obsp['spatial_connectivities'].shape}")
print(f"Mean neighbors per spot: {adata.obsp['spatial_connectivities'].sum(1).mean():.1f}")
```python

```python
sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names[:500], n_perms=100)
moran = adata.uns["moranI"].sort_values("I", ascending=False)
print("Top 10 spatially variable genes (Moran's I):")
print(moran.head(10))

fig, ax = plt.subplots(figsize=(6,4))
ax.scatter(range(len(moran)), moran["I"], s=5, alpha=0.5, color="steelblue")
ax.set_xlabel("Gene rank"); ax.set_ylabel("Moran's I")
ax.set_title("Spatial autocorrelation (Moran's I) — all genes")
plt.tight_layout(); plt.show()
```python

```python
top_svgs = moran.head(6).index.tolist()
print("Plotting top SVGs:", top_svgs)
fig = sq.pl.spatial_scatter(adata, color=top_svgs, ncols=3, size=1.4,
                            img_alpha=0.4, return_fig=True)
if fig: fig.suptitle("Top spatially variable genes", y=1.01)
plt.tight_layout(); plt.show()
```python

## Cell-type Deconvolution

Each Visium spot captures RNA from ~10–20 cells. To infer cell-type composition of each spot, **deconvolution** methods use a scRNA-seq reference.

**Key tools:**

| Tool | Approach | Reference needed |
|---|---|---|
| RCTD (robust cell type decomposition) | Poisson GLM | Yes (scRNA-seq) |
| cell2location | Negative binomial + hierarchical | Yes |
| Stereoscope | Probabilistic | Yes |
| NNLS | Non-negative least squares | Yes (marker genes) |

**Conceptual pattern:**
```python
# RCTD (spacexr package — R, or Python wrapper)
# 1. Build reference: scRNA-seq adata with cell_type labels
# 2. Run deconvolution on spatial adata
# 3. adata.obsm["cell_type_proportions"] shape: (n_spots, n_types)
# 4. sq.pl.spatial_scatter(adata, color=adata.obsm["cell_type_proportions"].columns)
```python

## Summary

| Step | Tool | Key Output |
|---|---|---|
| Load Visium | `sq.datasets.visium_hne_adata()` | AnnData with spatial coords |
| QC | `sc.pp.filter_cells/genes` | Clean count matrix |
| Normalize | `sc.pp.normalize_total + log1p` | Log-normalized expression |
| Cluster | `sc.tl.leiden` | Spot clusters |
| Spatial graph | `sq.gr.spatial_neighbors` | `obsp["spatial_connectivities"]` |
| Spatial autocorrelation | `sq.gr.spatial_autocorr(mode="moran")` | Moran's I per gene (−1 to +1) |
| SVG visualization | `sq.pl.spatial_scatter` | Gene expression on tissue image |
| Deconvolution | RCTD / cell2location / Stereoscope | Cell-type proportions per spot |

**Spatial statistics note:** Moran's I significance is assessed by permutation testing (shuffling spot labels) rather than a parametric test, because spatial data violates independence assumptions. `n_perms=100` gives approximate p-values; use `n_perms=1000` for publication-quality results.

**Related skill:** `spatial-transcriptomics.md`

```python
# Exercise 20
# 1. Identify the top 3 marker genes for each leiden cluster (sc.tl.rank_genes_groups).
#    Visualize them with sq.pl.spatial_scatter. Do clusters correspond to tissue regions?
# 2. Change the leiden resolution from 0.5 to 1.0 and 0.3. How does the number of
#    clusters change? Which resolution best matches the tissue anatomy?
# 3. Compute Moran's I for the top 5 marker genes from exercise 1.
#    Are highly spatially variable genes the same as cluster markers?
# 4. Challenge: implement a simple NNLS deconvolution using two "pure" reference
#    gene signatures (you can define them manually). Visualize the deconvolution
#    result spatially.
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
