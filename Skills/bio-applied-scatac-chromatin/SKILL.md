---
name: bio-applied-scatac-chromatin
description: "Single-Cell ATAC-seq: Chromatin Accessibility with NumPy"
tool_type: python
primary_tool: NumPy
---

# Single-Cell ATAC-seq: Chromatin Accessibility

**References:** [Signac](https://stuartlab.org/signac/) | [ArchR](https://www.archrproject.com/) | [chromVAR](https://www.nature.com/articles/nmeth.4401)

scATAC-seq measures chromatin accessibility at single-cell resolution — direct readout of where TFs bind and regulatory decisions are made, upstream of expression changes.

## Technology

**Tn5 tagmentation:** Hyperactive Tn5 inserts into nucleosome-free (accessible) regions, fragmenting DNA and ligating sequencing adapters. For 10x Genomics: isolate nuclei first (removes mtDNA), then add Tn5 inside droplets.

**Fragment file format** (TSV.gz, Tabix-indexed):
```
chr1    10000    10200    ACGTCAGTACGT-1    1
```
Columns: chrom, start, end, cell barcode, read count. Compact and cell-indexed — faster than BAM for per-cell operations.

**QC metrics:**
1. **TSS enrichment score** — ratio of coverage at TSS vs. flanking background. Score > 4 acceptable; > 8 high quality.
2. **Nucleosome banding** — fragment length histogram shows mono- (147 bp), di- (294 bp), tri-nucleosome peaks + sub-nucleosomal (<147 bp) free DNA.
3. **Fragments per cell** — 10,000–50,000 unique fragments typical for 10x.

| Feature | Bulk ATAC-seq | scATAC-seq |
|---------|--------------|------------|
| Input | 50,000–500,000 cells | Single nuclei |
| Sparsity | Dense (~100% at peaks) | Sparse (binary: 0 or 1) |
| Peak calling | Directly on BAM | Requires pseudo-bulk |

## Peak Calling and TF-IDF Normalization

**Why not call peaks per cell:** Each cell has only 0–1 reads per position — insufficient signal. Solution: pseudo-bulk peak calling (cluster cells → merge fragments per cluster → run MACS3 → create master peak set).

**TF-IDF** handles sparse binary ATAC data where library-size normalization is meaningless:
- TF = fragments in peak / total fragments in cell (depth normalization)
- IDF = log(1 + n_cells / cells_with_peak_accessible) (downweights ubiquitous peaks)
- Final = TF × IDF — emphasizes cell-type-specific peaks

SnapATAC2 uses raw fragment count with TF-IDF; ArchR uses binary (0/1) by default. Both work; binary is more robust to PCR duplicates.

## LSI Dimensionality Reduction

**Why PCA fails:** scATAC matrices are binary/sparse; PC1 correlates with read depth (technical artifact), not biology.

**LSI** (SVD on TF-IDF matrix):
1. Compute TF-IDF
2. Apply SVD
3. **Discard component 1** — almost always correlated with total fragment count
4. Use components 2–30 for UMAP + clustering

### SnapATAC2 (Python)
```python
import snapatac2 as snap

data = snap.read('fragments.tsv.gz', backed='r')
snap.pp.make_tile_matrix(data, bin_size=500)
snap.pp.select_features(data, n_features=50000)
snap.tl.spectral(data, n_comps=30)   # LSI via spectral embedding
snap.tl.umap(data)
snap.tl.leiden(data)
snap.pl.umap(data, color='leiden')
```

### Signac (R)
```r
library(Signac); library(Seurat)

chrom_assay <- CreateChromatinAssay(
  counts = peak_counts, fragments = 'fragments.tsv.gz',
  genome = 'hg38', min.cells = 10
)
seurat_obj <- CreateSeuratObject(counts = chrom_assay)
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q75')
seurat_obj <- RunSVD(seurat_obj)
DepthCor(seurat_obj)  # LSI1 should show high depth correlation -> drop it
```

### Python TF-IDF + LSI demo
```python
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from scipy.stats import pearsonr

rng = np.random.default_rng(42)
n_cells, n_peaks = 400, 5000
cell_types = ['CD4_T', 'CD8_T', 'B_cell', 'Monocyte']
cells_per = n_cells // len(cell_types)
cell_labels = np.repeat(cell_types, cells_per)

counts = np.zeros((n_cells, n_peaks), dtype=np.float32)
peak_offsets = {'CD4_T': 0, 'CD8_T': 500, 'B_cell': 1500, 'Monocyte': 3000}
for ct, offset in peak_offsets.items():
    ct_mask = cell_labels == ct
    specific_peaks = np.arange(offset, offset + 500)
    counts[np.ix_(np.where(ct_mask)[0], specific_peaks)] = 1
    for i in np.where(ct_mask)[0]:
        random_peaks = rng.choice(n_peaks, 500, replace=False)
        counts[i, random_peaks] += 1
counts = np.clip(counts, 0, 1)

depth_factors = rng.lognormal(0, 0.5, n_cells)
counts_noisy = (counts * depth_factors[:, None] > rng.uniform(0, 1, (n_cells, n_peaks))).astype(np.float32)

# TF-IDF
cell_totals = counts_noisy.sum(axis=1, keepdims=True) + 1e-6
TF = counts_noisy / cell_totals
cells_with_peak = (counts_noisy > 0).sum(axis=0) + 1
IDF = np.log1p(n_cells / cells_with_peak)
TFIDF = TF * IDF[np.newaxis, :]

# LSI
svd = TruncatedSVD(n_components=30, random_state=42)
X_lsi = svd.fit_transform(TFIDF)

depth_arr = counts_noisy.sum(axis=1)
corr1, _ = pearsonr(X_lsi[:, 0], depth_arr)
corr2, _ = pearsonr(X_lsi[:, 1], depth_arr)
print(f"LSI1 correlation with depth: {corr1:.3f}  -> DISCARD")
print(f"LSI2 correlation with depth: {corr2:.3f}")

X_lsi_clean = X_lsi[:, 1:]   # drop component 1
```

## Co-accessibility and Peak-Gene Links

**Cicero** (Pliner et al. 2018): trains graphical LASSO on pseudo-cells to compute co-accessibility scores (0–1) within 500 kb. Score > 0.25 suggests regulatory connection.

### Peak-gene linking (Signac)
```r
seurat_obj <- LinkPeaks(
  object = seurat_obj,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  distance = 5e5
)
```

**Validation:** Hi-C/HiChIP (3D proximity), eQTL overlap, ChromHMM enhancer states, ENCODE cCREs.

**Limitations:** Co-accessibility requires >500 cells; correlation ≠ causation; cell-type mixing creates spurious co-accessibility.

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing causes off-by-one errors.
- **Batch effects**: Always check for batch confounding before interpreting biological signal.
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously.
