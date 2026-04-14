---
name: bio-applied-scrna-preprocessing
description: "scRNA-seq QC and preprocessing: AnnData construction, QC metrics, MAD filtering, normalization strategies"
tool_type: python
primary_tool: scanpy
---

# scRNA-seq: Quality Control and Preprocessing

## Pitfalls

- **Mitochondrial threshold too strict**: pct_counts_mt < 5% discards real cells in metabolically active types (cardiomyocytes, hepatocytes). Always check the joint distribution with total counts.
- **Skipping `var_names_make_unique()`**: duplicate gene symbols cause cryptic downstream errors.
- **Normalizing before QC filtering**: low-quality cells bias normalization. Filter first.
- **HVG selection on raw counts**: must use normalized (not scaled) data. Scaling before HVG inflates dispersion estimates.

## Count Matrix Generation (FASTQ to counts)

**Cell Ranger** or **STARsolo** pipeline:
1. Barcode correction (1 Hamming distance against whitelist)
2. Splice-aware alignment (STAR)
3. UMI deduplication (cell barcode + gene + UMI collapsed)
4. Cell calling (EmptyDrops statistical test)

### MEX output files
- `barcodes.tsv.gz` — one cell barcode per line
- `features.tsv.gz` — Ensembl ID, gene symbol, feature type
- `matrix.mtx.gz` — sparse matrix (row=gene, col=cell, val=UMI count)

Read depth: 20k-50k reads/cell for clustering; 500k+ for rare isoforms.

## AnnData Structure

| Slot | Content |
|------|---------|
| `adata.X` | count matrix (n_cells x n_genes), usually sparse CSR |
| `adata.obs` | per-cell metadata |
| `adata.var` | per-gene metadata |
| `adata.obsm` | embeddings (X_pca, X_umap) |
| `adata.uns` | unstructured metadata |
| `adata.layers` | additional matrices (e.g. raw_counts) |

```python
# Save raw counts before modification (needed for DE testing)
adata.layers['raw_counts'] = scipy.sparse.csr_matrix(adata.X.copy())
```

## QC Metrics

| Metric | Low = problem | High = problem |
|--------|---------------|----------------|
| `n_genes_by_counts` | Empty droplet / dead cell | Doublet |
| `total_counts` | Dead/lysed cell | Doublet |
| `pct_counts_mt` | — | Dying cell (cytoplasmic RNA leaked) |

### MAD-based filtering (preferred over fixed thresholds)

```python
import scanpy as sc
import numpy as np

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

def mad_outlier(series, n_mad=3, direction="high"):
    median = series.median()
    mad = (series - median).abs().median()
    if direction == "high":
        return series > median + n_mad * mad
    return series < median - n_mad * mad

fail_total = mad_outlier(adata.obs["total_counts"], direction="low")
fail_genes = mad_outlier(adata.obs["n_genes_by_counts"], direction="low")
fail_mt    = mad_outlier(adata.obs["pct_counts_mt"], direction="high")

adata = adata[~(fail_total | fail_genes | fail_mt)].copy()
```

## Normalization

### Library-size (CP10K)
```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)  # log(x+1) — compresses dynamic range, symmetrizes fold-changes
```

### Alternatives
- **scran pooling**: `sc.external.pp.scran_normalization()` — better for strong composition bias
- **Pearson residuals** (SCTransform equivalent): `sc.experimental.pp.normalize_pearson_residuals()` — variance-stabilizing
