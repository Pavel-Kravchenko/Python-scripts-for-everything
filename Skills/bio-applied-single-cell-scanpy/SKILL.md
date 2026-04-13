---
name: bio-applied-single-cell-scanpy
description: "Single-cell analysis with Scanpy: AnnData objects, preprocessing, clustering, UMAP visualization, and marker gene detection. Use when analyzing scRNA-seq with Scanpy."
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/12_Modern_Workflows/01_single_cell_scanpy.ipynb"
primary_tool: scanpy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scanpy 1.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Single-Cell RNA-seq Analysis with Scanpy

*Source: Course notebook `Tier_3_Applied_Bioinformatics/12_Modern_Workflows/01_single_cell_scanpy.ipynb`*


**Tier 3 -- Applied Bioinformatics**

---

## Learning Objectives

By the end of this notebook you will be able to:

1. Explain **why** single-cell RNA-seq reveals biology that bulk RNA-seq cannot
2. Navigate the **AnnData** data structure and understand its components
3. Apply standard **quality control** metrics to filter low-quality cells
4. Perform **normalization**, log transformation, and highly variable gene selection
5. Reduce dimensionality with **PCA** and visualize with **UMAP**
6. Cluster cells with the **Leiden algorithm** and identify **marker genes**
7. Annotate cell types using known marker gene expression patterns

---

## 1. Why Single-Cell RNA-seq?

### The Biological Motivation

Traditional **bulk RNA-seq** measures the average gene expression across millions of cells. This averaging hides critical biological information:

| Bulk RNA-seq limitation | scRNA-seq reveals |
|-------------------------|-------------------|
| Cannot identify rare cell populations | Even 1% rare cells are detected |
| Masks cell-to-cell variability | Captures full expression distribution |
| Cannot track differentiation paths | Enables pseudotime trajectory analysis |
| Averages across tumor heterogeneity | Maps clonal diversity |

### Applications

- **Developmental biology**: Mapping lineage trees from stem cells to mature cell types
- **Cancer research**: Identifying drug-resistant subclones
- **Immunology**: Characterizing immune cell states in response to infection
- **Neuroscience**: Building cell-type atlases of the brain

```python
# Standard imports for single-cell analysis
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Configure scanpy settings
sc.settings.verbosity = 3  # Show informative output
sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(6, 5))

print(f"Scanpy version: {sc.__version__}")
```python

---

## 2. The AnnData Structure

Scanpy uses **AnnData** (Annotated Data), a purpose-built data structure for single-cell data. Think of it as a container that keeps expression data, cell metadata, and gene metadata together.

```python
AnnData object
    ├── X          : Expression matrix (cells × genes) -- sparse or dense
    ├── obs        : Cell metadata (DataFrame) -- cell barcodes, QC metrics, clusters
    ├── var        : Gene metadata (DataFrame) -- gene symbols, highly variable flags
    ├── obsm       : Cell embeddings (dict) -- PCA, UMAP, tSNE coordinates
    ├── varm       : Gene embeddings (dict) -- PCA loadings
    ├── obsp       : Cell-cell graphs (dict) -- neighbor connectivities
    ├── uns        : Unstructured data (dict) -- colors, clustering parameters
    └── raw        : Copy of original data before filtering
```python

### Why This Structure?

Single-cell datasets are **large** (often 10,000--100,000 cells × 20,000 genes) and **sparse** (most genes have zero counts in most cells). AnnData handles this efficiently with sparse matrix storage and keeps all analysis results organized.

```python
# Load the classic PBMC 3k dataset (Peripheral Blood Mononuclear Cells)
# This dataset from 10X Genomics contains ~2,700 cells and ~32,000 genes
adata = sc.datasets.pbmc3k()

print(f"Dataset shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
print(f"\nExpression matrix type: {type(adata.X).__name__}")
print(f"Sparsity: {1 - (adata.X.nnz / np.prod(adata.X.shape)):.1%} zeros")

# Display the AnnData object
adata
```python

```python
# Explore the structure
print("Cell metadata (adata.obs):")
print(adata.obs.head())

print("\nGene metadata (adata.var):")
print(adata.var.head())

print("\nFirst cell's expression (top 10 nonzero genes):")
cell_expr = pd.Series(adata.X[0].toarray().flatten(), index=adata.var_names)
print(cell_expr[cell_expr > 0].sort_values(ascending=False).head(10))
```python

---

## 3. Quality Control

### Why QC Matters

Not all captured cells are high-quality. Common problems include:

| Problem | Symptom | Filter criterion |
|---------|---------|------------------|
| Empty droplet | Very few genes detected | `n_genes_by_counts < 200` |
| Doublet (two cells) | Unusually high gene count | `n_genes_by_counts > 5000` |
| Dying/stressed cell | High mitochondrial % | `pct_counts_mt > 20%` |
| Low sequencing depth | Few total counts | `total_counts < 500` |

**Mitochondrial genes** are especially important: dying cells lose cytoplasmic mRNA, but mitochondrial mRNA (inside the organelle) is retained, causing artificially high MT percentages.

```python
# Mark mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
)

# Visualize QC metrics
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], show=False)
sc.pl.violin(adata, 'total_counts', ax=axes[1], show=False)
sc.pl.violin(adata, 'pct_counts_mt', ax=axes[2], show=False)
plt.tight_layout()
plt.show()
```python

```python
# Filter cells
print(f"Before filtering: {adata.n_obs} cells")

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 20, :]

print(f"After filtering: {adata.n_obs} cells")
```python

---

## 4. Normalization and Highly Variable Genes

### Why Normalize?

Different cells have different **sequencing depths** (total UMI counts), which can mask biological variation. We normalize to make cells comparable:

1. **Library-size normalization**: Scale each cell to a target sum (e.g., 10,000 counts)
2. **Log transformation**: Compress the dynamic range -- some genes have 0 counts, others have thousands

### Why Highly Variable Genes?

Most genes are either:
- **Not expressed** (zeros everywhere) -- no information
- **Housekeeping genes** (constant across all cells) -- no discriminatory power

We select **highly variable genes (HVGs)** -- genes that vary more than expected given their mean expression. These capture biological differences between cell types.

```python
# Normalize to 10,000 counts per cell
# This makes gene expression comparable across cells with different sequencing depth
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform to stabilize variance and compress dynamic range
# log1p = log(x + 1) to handle zeros
sc.pp.log1p(adata)

# Find highly variable genes
# These genes will drive our clustering -- they distinguish cell types
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

n_hvg = sum(adata.var.highly_variable)
print(f"Selected {n_hvg} highly variable genes ({n_hvg/adata.n_vars:.1%} of total)")

# Visualize: x = mean expression, y = dispersion (variance/mean)
# Black dots = highly variable genes
sc.pl.highly_variable_genes(adata)
```python

---

## 5. Dimensionality Reduction

### The Curse of Dimensionality

With ~2,000 genes, each cell is a point in 2,000-dimensional space. This causes problems:
- **Visualization**: Cannot plot 2,000 dimensions
- **Distance metrics**: Euclidean distance becomes meaningless in high dimensions
- **Clustering**: Algorithms struggle with sparse, high-dimensional data

### The Solution: Two-Step Reduction

1. **PCA** (Principal Component Analysis): Linear reduction to 30--50 dimensions
   - Captures most variance in few components
   - Denoises data by discarding low-variance directions

2. **UMAP** (Uniform Manifold Approximation and Projection): Nonlinear reduction to 2D
   - Preserves local neighborhood structure
   - Great for visualization (clusters visible)

```python
# Save raw counts for differential expression (needs unscaled data)
adata.raw = adata

# Subset to highly variable genes only
adata = adata[:, adata.var.highly_variable]
print(f"Reduced to {adata.n_vars} HVG features")

# Scale to zero mean, unit variance (required for PCA)
sc.pp.scale(adata, max_value=10)  # Clip outliers to max_value

# PCA: reduce from ~2000 genes to 50 principal components
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# How many PCs to keep? Look at variance explained
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
```python

```python
# Build k-nearest neighbor graph in PCA space
# This graph is used for both clustering and UMAP
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

# Compute UMAP embedding from the neighbor graph
sc.tl.umap(adata)

# Visualize: color by QC metrics to check for batch effects
sc.pl.umap(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
           ncols=3, wspace=0.4)
```python

---

## 6. Clustering

### The Leiden Algorithm

Clustering assigns each cell to a discrete group. The **Leiden algorithm** (an improvement over Louvain) optimizes **modularity** -- finding groups of cells that are more connected to each other than to cells outside the group.

The `resolution` parameter controls cluster granularity:
- **Low resolution** (0.2--0.4): Few large clusters (major cell types)
- **High resolution** (1.0--2.0): Many small clusters (subtypes)

```python
# Leiden clustering on the neighbor graph
sc.tl.leiden(adata, resolution=0.5)

print(f"Found {adata.obs['leiden'].nunique()} clusters")
print(f"\nCells per cluster:\n{adata.obs['leiden'].value_counts()}")

# Visualize clusters on UMAP
sc.pl.umap(adata, color='leiden', legend_loc='on data', 
           title='Leiden Clustering (resolution=0.5)')
```python

---

## 7. Marker Gene Detection

### What Are Marker Genes?

Marker genes are differentially expressed between clusters. They help identify:
- **What cell type** each cluster represents
- **Biological differences** between clusters

We use the **Wilcoxon rank-sum test** (non-parametric, robust to outliers) to find genes significantly higher in each cluster compared to all others.

```python
# Find marker genes for each cluster vs all others
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Plot top 5 markers per cluster
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, fontsize=10)
```python

```python
# Use known markers to identify cell types
# These are canonical markers for human PBMC populations
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

# Flatten to list and filter to available genes
all_markers = [g for markers in marker_genes.values() for g in markers]
available = [g for g in all_markers if g in adata.raw.var_names]
print(f"Available markers: {len(available)}/{len(all_markers)}")

# Dot plot: size = fraction expressing, color = mean expression
sc.pl.dotplot(adata, available, groupby='leiden', standard_scale='var',
              figsize=(12, 4))
```python

---

## Exercises

### Exercise 1: Annotate Cell Types

Based on the marker gene dot plot above, assign cell type labels to each cluster. Create a new column `adata.obs['cell_type']` with your annotations.

```python
# Your code here
# Example:
# cluster_annotations = {
#     '0': 'CD4+ T cells',
#     '1': 'CD14+ Monocytes',
#     '2': 'B cells',
#     # ... continue for all clusters
# }
# adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)
# sc.pl.umap(adata, color='cell_type')
```python

### Exercise 2: Resolution Sensitivity

Re-cluster the data at resolution 0.2 and 1.0. How does cluster number change? What is the trade-off between over-splitting and under-clustering?

```python
# Your code here
# sc.tl.leiden(adata, resolution=0.2, key_added='leiden_low')
# sc.tl.leiden(adata, resolution=1.0, key_added='leiden_high')
# sc.pl.umap(adata, color=['leiden_low', 'leiden', 'leiden_high'], ncols=3)
```python

### Exercise 3: Compare Wilcoxon and t-test

Re-run marker detection with `method='t-test'`. Do the top markers change significantly? When might one method be preferred?

```python
# Your code here
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', key_added='ttest_markers')
# sc.pl.rank_genes_groups(adata, key='ttest_markers', n_genes=5)
```python

---

## Summary

| Step | Function | Key Parameters |
|------|----------|----------------|
| **Load data** | `sc.read_*()` | File format |
| **QC metrics** | `sc.pp.calculate_qc_metrics()` | `qc_vars=['mt']` |
| **Filter cells** | `sc.pp.filter_cells()` | `min_genes`, `max_genes` |
| **Filter genes** | `sc.pp.filter_genes()` | `min_cells` |
| **Normalize** | `sc.pp.normalize_total()` | `target_sum=1e4` |
| **Log transform** | `sc.pp.log1p()` | -- |
| **Find HVGs** | `sc.pp.highly_variable_genes()` | `min_mean`, `min_disp` |
| **Scale** | `sc.pp.scale()` | `max_value=10` |
| **PCA** | `sc.tl.pca()` | `n_comps=50` |
| **Neighbor graph** | `sc.pp.neighbors()` | `n_neighbors`, `n_pcs` |
| **UMAP** | `sc.tl.umap()` | -- |
| **Cluster** | `sc.tl.leiden()` | `resolution` |
| **Find markers** | `sc.tl.rank_genes_groups()` | `method='wilcoxon'` |

---

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
