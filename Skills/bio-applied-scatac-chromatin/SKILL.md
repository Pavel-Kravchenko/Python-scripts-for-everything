---
name: bio-applied-scatac-chromatin
description: "Single-Cell ATAC-seq: Chromatin Accessibility with NumPy"
tool_type: python
primary_tool: NumPy
---

# Single-Cell ATAC-seq: Chromatin Accessibility

- [Signac documentation](https://stuartlab.org/signac/)
- [ArchR documentation](https://www.archrproject.com/)
- [chromVAR (Schep et al. 2017)](https://www.nature.com/articles/nmeth.4401)

## Why scATAC-seq matters

scRNA-seq measures gene expression, but expression is downstream of chromatin accessibility. Accessible chromatin (open regions) is where transcription factors bind and where regulatory decisions are made. scATAC-seq provides direct measurement of chromatin accessibility at single-cell resolution, revealing enhancer activity, TF binding landscapes, and regulatory states that precede changes in gene expression — making it essential for understanding gene regulation and cell fate decisions.

## scATAC-seq Technology

### Tn5 tagmentation
ATAC-seq (Assay for Transposase-Accessible Chromatin with sequencing) uses a hyperactive Tn5 transposase loaded with sequencing adapters. Tn5 preferentially inserts into accessible (nucleosome-free) regions of chromatin, simultaneously fragmenting the DNA and ligating sequencing adapters. Closed chromatin (wrapped around nucleosomes) is protected from Tn5 access.

For single-cell ATAC-seq (10x Genomics Chromium), nuclei are isolated first (critical: cytoplasm is removed to prevent mtDNA contamination), then Tn5 is added to the transposition reaction inside droplets.

### Fragment file format
The primary output is a fragments file (TSV.gz, Tabix-indexed):
```python
chr1    10000    10200    ACGTCAGTACGT-1    1
chr1    10050    10300    TGCAGTACGACC-1    2
```python
Columns: chromosome, start, end, cell barcode, read count for this fragment.

**Why fragments, not BAM?** Fragment files are compact and cell-indexed, enabling fast per-cell operations without loading entire BAM files.

### Quality metrics specific to scATAC-seq
1. **TSS enrichment score**: ratio of read coverage at TSS to flanking background regions. TSS are strongly enriched in ATAC signal because promoters are consistently accessible. Score > 4 is acceptable; > 8 is high quality.
2. **Nucleosome banding pattern**: histogram of fragment lengths shows mono- (147 bp), di- (294 bp), tri-nucleosome peaks plus sub-nucleosomal (< 147 bp) free DNA peaks.
3. **Fragments per cell**: 10,000–50,000 unique fragments per cell is typical for 10x scATAC-seq.

### Comparison with bulk ATAC-seq
| Feature | Bulk ATAC-seq | scATAC-seq |
|---------|--------------|------------|
| Input | 50,000–500,000 cells | Single nuclei |
| Resolution | Population average | Per-cell |
| Sparsity | Dense (~100% coverage at peaks) | Sparse (binary: 0 or 1) |
| Sensitivity | High | Lower per cell |
| Peak quality | Sharp, clear | Requires pseudo-bulk for calling |

## Peak Calling and Count Matrix

### Why you can't just call peaks per cell
Single-cell ATAC data is extremely sparse — each cell has only 0 or 1 reads at any given position. Peak callers like MACS3 require sufficient read depth to distinguish signal from noise, which a single cell cannot provide.

**Solution: pseudo-bulk peak calling**
1. Assign cells to clusters (based on LSI + Leiden)
2. Merge all fragments from cells in the same cluster → pseudo-bulk ATAC-seq sample
3. Run MACS3 on each pseudo-bulk sample
4. Merge peaks across all clusters using a master peak set

This gives high-quality peaks while preserving cell-type specificity.

### TF-IDF normalization for sparse binary ATAC data
The count matrix (cells × peaks) is mostly binary (0 or 1 read per peak per cell). Standard normalization fails because:
- Library-size normalization is meaningless for binary data
- Highly accessible peaks are accessible in nearly all cells → low information content

**TF-IDF (Term Frequency-Inverse Document Frequency)** addresses both:
- TF = fragment count in peak / total fragments in cell (normalizes for depth)
- IDF = log(1 + n_cells / cells_with_peak_accessible) (downweights ubiquitous peaks)
- Final score = TF × IDF

This emphasizes peaks that are accessible in a subset of cells (cell-type specific) over peaks accessible in all cells (constitutive).

### Binary vs count-based approaches
SnapATAC2 (Python) and ArchR (R) differ here:
- SnapATAC2 uses the raw fragment count (0, 1, 2...) with TF-IDF
- ArchR uses a binary accessibility score (0/1) as default
Both approaches work; binary is more robust to PCR duplicate issues.

## LSI Dimensionality Reduction

### Why PCA fails for scATAC data
PCA assumes variables follow a multivariate normal distribution and looks for linear combinations that maximize variance. scATAC count matrices are:
- Binary (0/1): not normally distributed
- Highly sparse: > 99% zeros
- The first PC is usually correlated with read depth (sequencing depth artifact), not biology

### LSI (Latent Semantic Indexing)
LSI is a dimensionality reduction method originally from information retrieval (finding documents similar in topic). It applies SVD (Singular Value Decomposition) to the TF-IDF transformed matrix:

1. Compute TF-IDF matrix (cells × peaks)
2. Apply SVD to extract latent dimensions
3. **Discard the first LSI component**: LSI component 1 is almost always correlated with total fragment count (sequencing depth), similar to how PC1 in bulk RNA-seq correlates with read depth. Dropping it removes this technical artifact.
4. Use LSI components 2–30 for UMAP + clustering

### SnapATAC2 Python implementation
```python
import snapatac2 as snap

# Load fragments file
data = snap.read('fragments.tsv.gz', backed='r')
snap.pp.make_tile_matrix(data, bin_size=500)  # 500bp bins
snap.pp.select_features(data, n_features=50000)  # top variable bins
snap.tl.spectral(data, n_comps=30)  # LSI via spectral embedding
snap.tl.umap(data)
snap.tl.leiden(data)
snap.pl.umap(data, color='leiden')
```python

### Signac R implementation
```r
library(Signac)
library(Seurat)

# Create ChromatinAssay from fragment file
chrom_assay <- CreateChromatinAssay(
  counts = peak_counts, fragments = 'fragments.tsv.gz',
  genome = 'hg38', min.cells = 10
)
seurat_obj <- CreateSeuratObject(counts = chrom_assay)

# Normalize with TF-IDF
seurat_obj <- RunTFIDF(seurat_obj)
# Feature selection: select top variable peaks
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q75')
# LSI
seurat_obj <- RunSVD(seurat_obj)
# Plot: check that LSI1 is not correlated with depth
DepthCor(seurat_obj)  # should show high correlation for LSI1, low for others
```python

```python
# Demonstrate TF-IDF normalization and LSI for scATAC-seq
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
import warnings
warnings.filterwarnings('ignore')

rng = np.random.default_rng(42)

# Simulate a sparse binary scATAC-seq count matrix
n_cells = 400
n_peaks = 5000
n_peak_per_cell = 1000  # ~20% accessible

cell_types = ['CD4_T', 'CD8_T', 'B_cell', 'Monocyte']
cells_per = n_cells // len(cell_types)
cell_labels = np.repeat(cell_types, cells_per)

# Each cell type has different accessible peaks
counts = np.zeros((n_cells, n_peaks), dtype=np.float32)
peak_offsets = {'CD4_T': 0, 'CD8_T': 500, 'B_cell': 1500, 'Monocyte': 3000}
for ct, offset in peak_offsets.items():
    ct_mask = cell_labels == ct
    # Type-specific peaks (constitutive for that type)
    specific_peaks = np.arange(offset, offset + 500)
    counts[np.ix_(np.where(ct_mask)[0], specific_peaks)] = 1
    # Random accessible peaks (variable)
    for i in np.where(ct_mask)[0]:
        random_peaks = rng.choice(n_peaks, n_peak_per_cell - 500, replace=False)
        counts[i, random_peaks] += 1

counts = np.clip(counts, 0, 1)  # binary

# Add varying depth (technical noise)
depth_factors = rng.lognormal(0, 0.5, n_cells)
counts_noisy = (counts * depth_factors[:, None] > rng.uniform(0, 1, (n_cells, n_peaks))).astype(np.float32)

print(f"scATAC-seq matrix: {n_cells} cells x {n_peaks} peaks")
print(f"Sparsity: {1 - counts_noisy.mean():.3f} ({(counts_noisy == 0).mean()*100:.1f}% zeros)")
print(f"Mean accessible peaks per cell: {counts_noisy.sum(axis=1).mean():.0f}")

# TF-IDF normalization ----
# TF: normalize each cell by its total fragment count (depth normalization)
cell_totals = counts_noisy.sum(axis=1, keepdims=True) + 1e-6  # avoid division by zero
TF = counts_noisy / cell_totals

# IDF: log(1 + n_cells / cells_with_peak_open)
cells_with_peak = (counts_noisy > 0).sum(axis=0) + 1  # +1 to avoid log(0)
IDF = np.log1p(n_cells / cells_with_peak)

TFIDF = TF * IDF[np.newaxis, :]
print(f"\nTF-IDF matrix range: [{TFIDF.min():.4f}, {TFIDF.max():.4f}]")

# LSI via SVD ----
n_components = 30
svd = TruncatedSVD(n_components=n_components, random_state=42)
X_lsi = svd.fit_transform(TFIDF)

# Check: LSI component 1 correlation with sequencing depth
from scipy.stats import pearsonr
depth_arr = counts_noisy.sum(axis=1)
corr_comp1, _ = pearsonr(X_lsi[:, 0], depth_arr)
corr_comp2, _ = pearsonr(X_lsi[:, 1], depth_arr)
print(f"\nLSI component 1 correlation with depth: {corr_comp1:.3f}")
print(f"LSI component 2 correlation with depth: {corr_comp2:.3f}")
print(f"-> Component 1 is dominated by depth, DISCARD it (use LSI 2-30)")

# Use components 2+ for downstream analysis
X_lsi_clean = X_lsi[:, 1:]  # drop first component

# UMAP on LSI (components 2-16)
try:
    from umap import UMAP
    X_umap = UMAP(n_neighbors=15, min_dist=0.3, random_state=42).fit_transform(X_lsi_clean[:, :15])
except ImportError:
    X_umap = X_lsi_clean[:, :2]

adata_atac = ad.AnnData(
    X=sp.csr_matrix(counts_noisy),
    obs=pd.DataFrame({'cell_type': cell_labels}, index=[f'CELL_{i:04d}' for i in range(n_cells)]),
    var=pd.DataFrame(index=[f'peak_{i:05d}' for i in range(n_peaks)])
)
adata_atac.obsm['X_lsi'] = X_lsi_clean
adata_atac.obsm['X_umap'] = X_umap

# Visualize
fig, axes = plt.subplots(1, 3, figsize=(16, 4))
colors_map = {'CD4_T':'#e41a1c','CD8_T':'#ff7f00','B_cell':'#377eb8','Monocyte':'#4daf4a'}

for ct in cell_types:
    mask = cell_labels == ct
    axes[0].scatter(X_lsi[mask, 0], X_lsi[mask, 1], c=colors_map[ct], s=8, alpha=0.7, label=ct)
axes[0].set_title('LSI (all components)\nComponent 1 = depth artifact')
axes[0].set_xlabel('LSI1'); axes[0].set_ylabel('LSI2')
axes[0].legend(fontsize=7)

for ct in cell_types:
    mask = cell_labels == ct
    axes[1].scatter(X_lsi_clean[mask, 0], X_lsi_clean[mask, 1], c=colors_map[ct], s=8, alpha=0.7, label=ct)
axes[1].set_title('LSI (components 2+)\nAfter depth correction')
axes[1].set_xlabel('LSI2'); axes[1].set_ylabel('LSI3')
axes[1].legend(fontsize=7)

for ct in cell_types:
    mask = cell_labels == ct
    axes[2].scatter(X_umap[mask, 0], X_umap[mask, 1], c=colors_map[ct], s=8, alpha=0.7, label=ct)
axes[2].set_title('UMAP from LSI 2-16')
axes[2].set_xlabel('UMAP1'); axes[2].set_ylabel('UMAP2')
axes[2].legend(fontsize=7)

plt.tight_layout()
plt.savefig('scatac_lsi.png', dpi=100, bbox_inches='tight')
plt.show()
```python

## Co-accessibility and Peak-Gene Links

### Cicero co-accessibility
Cicero (Pliner et al. 2018) identifies co-accessible peak pairs — peaks that tend to be accessible in the same cells. This reveals regulatory interactions: a distal enhancer that co-activates with a gene's promoter peak is likely a functional regulatory element for that gene.

**How it works**: Cicero trains a graphical LASSO model on pseudo-cells (groups of similar cells) to compute co-accessibility scores (0 to 1) between all pairs of peaks within 500 kb. A score > 0.25 suggests a regulatory connection.

### Peak-gene linking in Signac
```r
# Link peaks to genes using correlation within Seurat/Signac
seurat_obj <- LinkPeaks(
  object = seurat_obj,
  peak.assay = "ATAC",
  expression.assay = "RNA",  # requires matched RNA
  distance = 5e5  # 500 kb window
)
# Result: each peak gets a list of correlated genes
```python

### Validation approaches
- **Hi-C / HiChIP**: 3D genome contact maps confirm spatial proximity
- **eQTL overlap**: variants in co-accessible peaks that affect gene expression
- **ChromHMM enhancer states**: overlap with H3K27ac ChIP-seq enhancer marks
- **ENCODE cCREs**: catalog of candidate cis-regulatory elements

### Limitations
- Co-accessibility requires sufficient cells (> 500) to reliably estimate correlations
- Correlation ≠ causation: co-accessible peaks may both respond to a third regulatory factor
- Cell type mixing can create spurious co-accessibility if not properly controlled

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
