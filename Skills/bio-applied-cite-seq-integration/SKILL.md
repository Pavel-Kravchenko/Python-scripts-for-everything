---
name: bio-applied-cite-seq-integration
description: "CITE-seq and Multiome integration — ADT normalization (CLR/DSB), WNN graph construction, and paired RNA+ATAC analysis with muon"
tool_type: python
primary_tool: NumPy
---

# CITE-seq and Multiome Data Integration

- [Seurat WNN tutorial](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis)
- [muon Python package](https://muon.readthedocs.io/)

## ADT Normalization

ADT counts differ fundamentally from RNA: high background (~100-500 counts even in negative cells), no true zeros, bimodal distribution (negative + positive populations).

### CLR (Centered Log-Ratio)

```python
def clr_normalize(X):
    """X: cells x proteins matrix (raw counts)."""
    X_ps = X + 0.5
    geo_mean = np.exp(np.log(X_ps).mean(axis=1, keepdims=True))
    return np.log(X_ps / geo_mean)
```

In Seurat: `NormalizeData(assay='ADT', normalization.method='CLR', margin=2)`.

**Caveat**: CLR assumes balanced panel. If 80% of antibodies target one lineage, geometric mean is skewed.

### DSB (Denoised and Scaled by Background)

```python
def dsb_normalize(X, background):
    """Subtract per-protein background, scale by noise. Requires empty droplet data."""
    bg_mean = np.log1p(background).mean(axis=0)
    bg_std  = np.log1p(background).std(axis=0) + 1e-6
    return (np.log1p(X) - bg_mean) / bg_std
```

DSB values > 3-4 = clearly positive cell. Requires `raw_feature_bc_matrix.h5` (empty droplets).

| Method | Needs empty droplets | Best for |
|--------|---------------------|---------|
| CLR | No | Quick analysis, <50 proteins |
| DSB | Yes | Rigorous analysis, high background proteins |

## WNN Integration (Hao et al. 2021)

Per-cell, per-modality weights learned automatically: RNA gets higher weight where RNA is informative, protein gets higher weight where protein is clearer (e.g., CD4 protein >> CD4 mRNA for T cell subtyping).

### muon Implementation

```python
import muon as mu

mdata = mu.MuData({'rna': adata_rna, 'prot': adata_prot})

# RNA preprocessing
sc.pp.normalize_total(mdata['rna'], target_sum=1e4)
sc.pp.log1p(mdata['rna'])
sc.pp.pca(mdata['rna'], n_comps=30)
sc.pp.neighbors(mdata['rna'])

# Protein: CLR + PCA
mu.prot.pp.clr(mdata['prot'])
sc.pp.pca(mdata['prot'], n_comps=15)
sc.pp.neighbors(mdata['prot'])

# WNN integration
mu.pp.neighbors(mdata, key_added='wnn')
sc.tl.umap(mdata, neighbors_key='wnn')
```

## 10x Multiome (RNA + ATAC)

| Aspect | CITE-seq | 10x Multiome |
|--------|----------|-------------|
| Modalities | RNA + protein | RNA + ATAC |
| Cell vs nucleus | Cells | Nuclei only |
| Regulatory inference | Limited | Direct (peak -> gene) |
| Read requirement | 20k/cell RNA + 1k/cell ADT | 20k/cell RNA + 10k/cell ATAC |

```python
mdata = mu.read_10x_h5('filtered_feature_bc_matrix.h5')

# RNA: standard preprocessing
sc.pp.normalize_total(mdata['rna']); sc.pp.log1p(mdata['rna']); sc.pp.pca(mdata['rna'])

# ATAC: TF-IDF + LSI (drop first component — depth-correlated)
mu.atac.pp.tfidf(mdata['atac'], scale_factor=1e4)
sc.tl.pca(mdata['atac'])

mu.pp.neighbors(mdata)
sc.tl.umap(mdata)
```

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
