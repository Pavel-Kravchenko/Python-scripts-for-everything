---
name: bio-applied-scrna-preprocessing
description: "**Tier 3 — Applied Bioinformatics | Module 30 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/30_Single_Cell_RNA_seq/01_scrna_preprocessing.ipynb"
primary_tool: scanpy
---

## Version Compatibility

Reference examples tested with: anndata 0.10+, matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# scRNA-seq: Quality Control and Preprocessing

*Source: Course notebook `Tier_3_Applied_Bioinformatics/30_Single_Cell_RNA_seq/01_scrna_preprocessing.ipynb`*

# scRNA-seq: Quality Control and Preprocessing

**Tier 3 — Applied Bioinformatics | Module 30 · Notebook 1**

*Prerequisites: Module 03 (RNA-seq Analysis), Module 07 (Machine Learning for Biology)*

---

**By the end of this notebook you will be able to:**
1. Explain the scRNA-seq experimental workflow (droplet-based 10x Genomics vs plate-based Smart-seq)
2. Generate a count matrix from raw FASTQ files using Cell Ranger or STARsolo
3. Perform quality control: filter cells by UMI count, gene count, and % mitochondrial reads
4. Normalize (library-size, scran) and log-transform count data
5. Identify highly variable genes (HVGs) for downstream analysis



**Key resources:**
- [Scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html)
- [Seurat tutorials](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
- [Harvard HBC scRNA-seq training](https://hbctraining.github.io/scRNA-seq_online/)
- [Best practices for scRNA-seq (Luecken & Theis 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746)

## Why single-cell preprocessing matters

Bulk RNA-seq averages gene expression over thousands of cells, masking rare populations and cellular heterogeneity. A tumor biopsy with 5% cancer-stem-like cells will show those cells' signatures diluted 20-fold. scRNA-seq resolves this by profiling each cell independently — but the raw data is exceptionally noisy: dropouts (zero counts for expressed genes), ambient RNA contamination, and doublets (two cells captured as one) all corrupt the signal if not addressed. This notebook walks through every decision in the preprocessing pipeline, explaining *why* each step is necessary and what goes wrong if you skip it.

## Environment setup and synthetic data

This notebook uses synthetic count matrix data that mimics a real 10x Genomics PBMC dataset. All code runs without downloading external files.

## Common preprocessing pitfalls

- **Mitochondrial threshold too strict**: setting pct_counts_mt < 5% in metabolically active cells (cardiomyocytes, hepatocytes) will discard real cells. Always look at the joint distribution with total counts.
- **Skipping `var_names_make_unique()`**: duplicate gene symbols (e.g., two ENSG IDs mapped to the same symbol) cause cryptic downstream errors in scanpy operations.
- **Normalizing before QC filtering**: including low-quality cells biases the normalization target. Filter first, then normalize.
- **HVG selection on raw counts**: highly variable gene selection must be done on normalized (but not scaled) data. Scaling before HVG selection inflates dispersion estimates.

## Quick demo: Build a synthetic scRNA-seq AnnData object

We create a 500-cell × 2000-gene count matrix that mimics a 10x Genomics PBMC experiment. Five cell populations are embedded with distinct marker gene programs.

```python
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

rng = np.random.default_rng(42)

n_cells = 500
n_genes = 2000
n_mt_genes = 13  # human mtDNA encodes 13 proteins

# Gene names: first 13 are mitochondrial (MT-), rest are nuclear
gene_names = [f"MT-{g}" for g in ["ND1","ND2","COX1","COX2","ATP8","ATP6","COX3",
                                    "ND3","ND4L","ND4","ND5","ND6","CYB"]]
gene_names += [f"GENE{i:04d}" for i in range(n_genes - n_mt_genes)]

# Five cell types with different baseline expression levels
cell_types = ["T_cell", "B_cell", "Monocyte", "NK_cell", "Dendritic"]
cells_per_type = n_cells // len(cell_types)
cell_type_labels = np.repeat(cell_types, cells_per_type)

# Base count matrix: negative binomial draws
counts = rng.negative_binomial(n=2, p=0.9, size=(n_cells, n_genes)).astype(np.float32)

# Add cell-type-specific marker programs
type_idx = {ct: np.where(cell_type_labels == ct)[0] for ct in cell_types}
# T cells: high CD3D (gene 50), CD3E (gene 51)
counts[np.ix_(type_idx["T_cell"], [50, 51, 52])] += rng.poisson(15, (cells_per_type, 3))
# B cells: high CD79A (gene 100), MS4A1 (gene 101)
counts[np.ix_(type_idx["B_cell"], [100, 101])] += rng.poisson(20, (cells_per_type, 2))
# Monocytes: high LYZ (gene 150), CST3 (gene 151)
counts[np.ix_(type_idx["Monocyte"], [150, 151])] += rng.poisson(25, (cells_per_type, 2))

# Inject 20 low-quality cells (high MT%, low counts)
low_q_idx = rng.choice(n_cells, 20, replace=False)
counts[low_q_idx, :n_mt_genes] += rng.poisson(50, (20, n_mt_genes))
counts[low_q_idx, n_mt_genes:] = counts[low_q_idx, n_mt_genes:] * 0.1

# Inject 5 doublets (very high total counts)
doublet_idx = rng.choice(n_cells, 5, replace=False)
counts[doublet_idx, :] *= 2.5

# Build AnnData object (the core data structure in scanpy)
obs = pd.DataFrame({'cell_type_truth': cell_type_labels}, 
                   index=[f"CELL_{i:04d}" for i in range(n_cells)])
var = pd.DataFrame({'gene_name': gene_names,
                    'is_mt': [g.startswith('MT-') for g in gene_names]},
                   index=gene_names)

adata = ad.AnnData(X=counts, obs=obs, var=var)
print(f"AnnData object: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"Mitochondrial genes: {adata.var['is_mt'].sum()}")
print(f"\nData type: {type(adata.X)}")
print(f"Memory: {adata.X.nbytes / 1e6:.1f} MB (dense); real 10x data uses sparse scipy.csr_matrix)")
```

## 2. Generating Count Matrices

### From FASTQ to count matrix
Raw reads from a 10x experiment are processed by **Cell Ranger** (commercial) or **STARsolo** (open-source, recommended):

1. **Barcode correction**: Cell barcodes matched against a whitelist (~6 million valid 10x barcodes), allowing 1 Hamming distance
2. **Alignment**: Reads aligned to the genome using STAR with splice-aware alignment
3. **UMI counting**: Reads sharing cell barcode + gene + UMI collapsed to 1 count (UMI deduplication)
4. **Cell calling**: EmptyDrops (Cell Ranger 3+) models the ambient RNA distribution and uses a statistical test to call cells vs empty droplets

### The three output files (MEX / Market Exchange format)
- `barcodes.tsv.gz` — one cell barcode per line
- `features.tsv.gz` — Ensembl ID, gene symbol, feature type per row
- `matrix.mtx.gz` — sparse matrix (row=gene, col=cell, val=UMI count)

Loading in scanpy is one line, but `var_names_make_unique()` is essential — multiple Ensembl IDs can map to the same gene symbol (e.g., MATR3 appears twice in hg38), which causes silent errors later.

**Read depth recommendations**: 20,000–50,000 reads/cell for robust clustering; 500,000+ reads/cell for rare isoform detection. Saturation curves in Cell Ranger QC report show whether sequencing depth is sufficient.

```python
# Demonstrate what AnnData looks like after loading real 10x data
# (using our synthetic object from above)

# Key AnnData slots:
# adata.X       -- count matrix (n_cells x n_genes), usually scipy sparse CSR
# adata.obs     -- per-cell metadata DataFrame  
# adata.var     -- per-gene metadata DataFrame
# adata.obsm    -- per-cell embeddings (PCA, UMAP stored here as 'X_pca', 'X_umap')
# adata.uns     -- unstructured metadata (color palettes, clustering parameters)
# adata.layers  -- additional count matrices (e.g. 'raw_counts' before normalization)

# Save raw counts before any modification (important for downstream DE testing)
import scipy.sparse as sp
adata.layers['raw_counts'] = sp.csr_matrix(adata.X.copy())

print("AnnData structure:")
print(f"  .X shape: {adata.X.shape}")
print(f"  .obs columns: {list(adata.obs.columns)}")
print(f"  .var columns: {list(adata.var.columns)}")
print(f"\nFirst 5 cell barcodes:")
print(adata.obs_names[:5].tolist())
print(f"\nFirst 5 gene names:")
print(adata.var_names[:5].tolist())
```

## 3. Quality Control Metrics

Three metrics capture the most common data quality problems:

| Metric | What it measures | Low = problem | High = problem |
|--------|-----------------|---------------|----------------|
| `n_genes_by_counts` | Number of distinct genes detected | Empty droplet or dead cell | Doublet (two cells) |
| `total_counts` | Total UMI count per cell (library size) | Dead/lysed cell | Doublet |
| `pct_counts_mt` | % UMIs from mitochondrial genes | — | Dying cell (cytoplasmic RNA leaked, MT RNA retained) |

### Why high MT% indicates dying cells
When a cell ruptures, cytoplasmic mRNA leaks out of the droplet before being captured, but mitochondria (membrane-enclosed) remain intact and their transcripts are captured. A cell with 40% MT reads has almost certainly lost its cytoplasmic transcriptome — it is not a biological MT-rich cell type.

### MAD-based filtering (preferred over fixed thresholds)
The median absolute deviation (MAD) approach adapts thresholds to each dataset's distribution:
- Cells with `total_counts` < median − 3×MAD are filtered
- Cells with `pct_counts_mt` > median + 3×MAD are filtered
This is more principled than arbitrary cutoffs like "< 200 genes" or "> 20% MT".

```python
# Compute QC metrics manually to understand what scanpy's sc.pp.calculate_qc_metrics does
import numpy as np
import scipy.sparse as sp

X = adata.X if not sp.issparse(adata.X) else adata.X.toarray()
mt_mask = adata.var['is_mt'].values

adata.obs['total_counts'] = X.sum(axis=1)
adata.obs['n_genes_by_counts'] = (X > 0).sum(axis=1)
mt_counts = X[:, mt_mask].sum(axis=1)
adata.obs['pct_counts_mt'] = (mt_counts / adata.obs['total_counts'] * 100)

print("QC metric summary:")
print(adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe().round(2))

# MAD-based thresholds (robust, adapts to each dataset)
def mad_outlier(series, n_mad=3, direction="high"):
    median = series.median()
    mad = (series - median).abs().median()
    if direction == "high":
        return series > median + n_mad * mad
    elif direction == "low":
        return series < median - n_mad * mad

fail_total = mad_outlier(adata.obs["total_counts"], direction="low")
fail_genes = mad_outlier(adata.obs["n_genes_by_counts"], direction="low")
fail_mt    = mad_outlier(adata.obs["pct_counts_mt"], direction="high")

quality_mask = ~(fail_total | fail_genes | fail_mt)
print(f"\nCells before QC: {adata.n_obs}")
print(f"  Failed total_counts (low): {fail_total.sum()}")
print(f"  Failed n_genes (low):      {fail_genes.sum()}")
print(f"  Failed MT% (high):         {fail_mt.sum()}")
print(f"Cells after QC: {quality_mask.sum()}")

# Joint scatter plot: key diagnostic
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

colors = ["red" if not q else "steelblue" for q in quality_mask]
axes[0].scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"],
                c=colors, alpha=0.4, s=8)
axes[0].set_xlabel("Total UMI counts")
axes[0].set_ylabel("Genes detected")
axes[0].set_title("QC: genes vs counts (red = filtered)")

axes[1].scatter(adata.obs["total_counts"], adata.obs["pct_counts_mt"],
                c=colors, alpha=0.4, s=8)
axes[1].set_xlabel("Total UMI counts")
axes[1].set_ylabel("% Mitochondrial counts")
axes[1].set_title("QC: MT% vs counts")

plt.tight_layout()
plt.savefig("qc_scatter.png", dpi=100, bbox_inches="tight")
plt.show()

# Apply filter (subset AnnData)
adata = adata[quality_mask].copy()
print(f"\nFiltered to {adata.n_obs} cells x {adata.n_vars} genes")
```

## 4. Normalization and Log-Transformation

### Why normalize?
Each cell has a different total UMI count (library size) due to technical variation — more cDNA molecules were captured from some cells than others. Without normalization, a gene appearing at 100 counts/cell in a 10,000-count cell has the same normalized expression as a gene at 200 counts in a 20,000-count cell.

### Library-size normalization (CPM-style)
Divide each cell's counts by its total, then multiply by a scale factor (10,000 by convention, making units "counts per 10k" or CP10K):

```
normalized_count = raw_count / cell_total * 10,000
```

This is `sc.pp.normalize_total(adata, target_sum=1e4)` in scanpy.

### Log1p transformation
After normalization, counts span several orders of magnitude (0.001 to 1000). Log1p compression `log(x+1)`:
- Brings the distribution closer to normal (required for PCA)
- The +1 prevents log(0) for zero counts
- Makes fold-changes between groups symmetric: doubling and halving are equal distances

### scran pooling-based normalization (alternative)
scran (Lun et al. 2016) pools cells to estimate size factors, then normalizes each cell by its estimated factor. More statistically rigorous than library-size normalization, particularly for datasets with strong composition bias (e.g. one cell type dominates). Implemented in `sc.external.pp.scran_normalization()`.

### SCTransform / analytic Pearson residuals (Seurat approach)
Fits a regularized negative binomial regression model to remove the mean-variance trend inherent to count data. Equivalent to variance-stabilizing normalization. For scanpy users: `sc.experimental.pp.normalize_pearson_residuals()`.
