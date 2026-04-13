---
name: bio-applied-cite-seq-integration
description: "**Tier 3 — Applied Bioinformatics | Module 31 · Notebook 2**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/31_Single_Cell_Multi_Omics/02_cite_seq_integration.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, muon 0.1+, numpy 1.26+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# CITE-seq and Multiome Data Integration

*Source: Course notebook `Tier_3_Applied_Bioinformatics/31_Single_Cell_Multi_Omics/02_cite_seq_integration.ipynb`*

# CITE-seq and Multiome Data Integration

**Tier 3 — Applied Bioinformatics | Module 31 · Notebook 2**

*Prerequisites: Notebook 1 (scATAC-seq)*

---

**By the end of this notebook you will be able to:**
1. Describe CITE-seq (RNA + protein), 10x Multiome (RNA + ATAC), and Spatial Transcriptomics co-modality designs
2. Process ADT (antibody-derived tags) counts and perform DSB normalization
3. Perform weighted nearest neighbor (WNN) joint embedding of RNA and protein
4. Integrate matched RNA and ATAC profiles from 10x Multiome data
5. Identify cell-type-specific regulatory elements linking chromatin to expression



**Key resources:**
- [Seurat WNN tutorial](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis)
- [CITE-seq paper (Stoeckius et al. 2017)](https://www.nature.com/articles/nmeth.4380)
- [muon Python package](https://muon.readthedocs.io/)

## 1. CITE-seq Experimental Design

### What is CITE-seq?
CITE-seq (Cellular Indexing of Transcriptomes and Epitopes by sequencing, Stoeckius et al. 2017) simultaneously measures **RNA** and **surface protein** expression in the same cell. This is achieved by:

1. Incubating cells with a cocktail of antibodies each conjugated to a unique **DNA oligonucleotide barcode** (Feature Barcode Technology)
2. Performing standard 10x Genomics droplet capture — the antibody-DNA conjugates co-encapsulate with mRNA from the same cell
3. During library preparation, RNA and ADT (Antibody-Derived Tag) libraries are separated by size and sequenced separately

### Antibody-Derived Tags (ADT) vs Hashtag Oligos (HTO)
- **ADT**: antibodies targeting cell-surface proteins (CD3, CD4, CD8, CD19...). Measures protein expression.
- **HTO**: antibodies conjugated to sample-specific barcodes. Used to multiplex multiple samples into a single 10x run, then "demultiplex" computationally. Reduces cost; 4–16 samples per lane.

### Why measure protein AND RNA?
- **Protein is more stable**: surface protein levels reflect the "display" state of the cell, while RNA reflects production. They are correlated but capture different biology.
- **Classic markers are protein-based**: FACS-defined T cell subsets (CD4/CD8) are far cleaner by protein than RNA due to low CD4/CD8 mRNA expression.
- **Sub-population resolution**: B cell subsets (naive, memory, plasmablast) show overlapping RNA signatures but distinct protein profiles. 
- **Cross-validate clusters**: an RNA-based cluster that lacks the expected protein marker is suspect.

### Panel design considerations
- Start with known lineage markers (CD3, CD19, CD14, CD56 for PBMC)
- Include "anchoring proteins" that are highly enriched in specific cell types to validate clustering
- Include both specific and broad markers — CD45 (pan-immune) serves as a positive control
- Maximum ~200 antibodies per 10x run; prioritize functionally informative proteins
- Avoid antibodies with high non-specific background (test empirically in bulk first)

### The ADT count matrix structure
ADTs produce a separate count matrix: rows = protein targets, columns = cell barcodes. This is loaded alongside the RNA matrix in the Cell Ranger `filtered_feature_bc_matrix.h5` file, accessible as feature_type == "Antibody Capture".

## 2. ADT Processing and Normalization

### The protein background problem
ADT counts have a fundamentally different distribution from RNA:
- **High background**: even cells that don't express a protein will show hundreds of counts due to non-specific antibody binding. This creates a bimodal distribution: negative population (background ~100–500) + positive population (signal ~1000–10,000).
- **No true zeros**: unlike RNA where unexpressed genes genuinely have 0 reads, every protein shows at least background counts in every cell.
- **Log-normal, not Poisson**: the negative population follows a log-normal distribution, not the Poisson/NB model that fits RNA.

### CLR (Centered Log-Ratio) normalization
The simplest and most widely used approach. For each cell, compute the geometric mean of all ADT counts, then divide each protein count by the geometric mean and take the log:

```
CLR(x_i) = log(x_i / geometric_mean(x))
```

This centers the distribution around 0, making comparison across cells and proteins meaningful. CLR is implemented in Seurat as `NormalizeData(assay='ADT', normalization.method='CLR', margin=2)` (margin=2 = normalize across features within each cell).

**Caveat**: CLR assumes the panel is roughly balanced — if 80% of your antibodies target T cell markers and you run on B cells, the geometric mean is dominated by negatives and CLR will over-correct positives.

### DSB (Denoised and Scaled by Background) normalization
DSB (Mulè et al. 2022) is more principled:
1. Estimate per-protein background level from **empty droplets** (droplets with no cell but still containing ambient antibody)
2. Subtract the expected background from each protein's count
3. Divide by the technical noise component (estimated from the data)

**Result**: DSB-normalized values are interpretable as signal-to-noise ratios above background. Values > 3–4 indicate a clearly positive cell.

**Requirement**: you need empty droplet data (from Cell Ranger `raw_feature_bc_matrix.h5`, not just the filtered version).

### When to use which normalization
| Method | Needs empty droplets | Speed | Best for |
|--------|---------------------|-------|---------|
| CLR | No | Fast | Quick analysis, <50 proteins |
| DSB | Yes | Moderate | Rigorous analysis, high background proteins |
| Raw log | No | Fastest | Initial exploration only |

```python
# ADT normalization: CLR and DSB-style approaches with visualization
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# --- CLR normalization ---
def clr_normalize(X):
    """Centered log-ratio normalization for ADT counts.
    X: cells x proteins matrix (raw counts, no zeros assumed but we add pseudocount)
    """
    # Add pseudocount to avoid log(0)
    X_ps = X + 0.5
    # Geometric mean per cell (axis=1 = across proteins within a cell)
    geo_mean = np.exp(np.log(X_ps).mean(axis=1, keepdims=True))
    return np.log(X_ps / geo_mean)

adt_clr = clr_normalize(adt_counts)

# --- DSB-style normalization (simplified without empty droplets) ---
# Simulate empty droplet background: lower count range
n_empty = 200
empty_adt = rng.poisson(80, (n_empty, n_proteins)).astype(np.float32)  # ambient only

def dsb_normalize(X, background):
    """Simplified DSB normalization.
    Subtracts per-protein background mean, scales by background std.
    """
    bg_mean = np.log1p(background).mean(axis=0)
    bg_std  = np.log1p(background).std(axis=0) + 1e-6
    return (np.log1p(X) - bg_mean) / bg_std

adt_dsb = dsb_normalize(adt_counts, empty_adt)

print("CLR normalized ADT:")
print(f"  Range: [{adt_clr.min():.2f}, {adt_clr.max():.2f}]")
print(f"  Mean: {adt_clr.mean():.3f}")

print("\nDSB normalized ADT:")
print(f"  Range: [{adt_dsb.min():.2f}, {adt_dsb.max():.2f}]")
print(f"  DSB > 3 (positive cells): {(adt_dsb > 3).mean():.3f}")

# Visualize: bimodal distribution for CD3 protein
fig, axes = plt.subplots(2, 3, figsize=(14, 8))

# Raw counts for CD3 and CD4
for col, prot_idx, prot_name in [(0, 0, 'CD3'), (1, 1, 'CD4'), (2, 3, 'CD19')]:
    raw_vals = np.log1p(adt_counts[:, prot_idx])
    clr_vals = adt_clr[:, prot_idx]
    dsb_vals = adt_dsb[:, prot_idx]
    
    axes[0, col].hist(raw_vals, bins=40, edgecolor='k', alpha=0.7, color='steelblue')
    axes[0, col].set_title(f'{prot_name}: log1p raw counts')
    axes[0, col].set_xlabel('log1p(ADT count)')
    axes[0, col].set_ylabel('Cells')
    
    axes[1, col].hist(dsb_vals, bins=40, edgecolor='k', alpha=0.7, color='tomato')
    axes[1, col].axvline(3, color='black', linestyle='--', label='DSB=3 (positive cutoff)')
    axes[1, col].set_title(f'{prot_name}: DSB normalized')
    axes[1, col].set_xlabel('DSB score')
    axes[1, col].legend(fontsize=8)

plt.suptitle('ADT distribution: raw counts vs DSB normalization\n'
             'Bimodal distribution indicates positive/negative cell populations', y=1.02)
plt.tight_layout()
plt.savefig('adt_normalization.png', dpi=100, bbox_inches='tight')
plt.show()

# Show that CLR-normalized values separate cell types for CD3
print("\nCD3 CLR values by cell type:")
for ct in ['CD4_T', 'CD8_T', 'B_cell', 'Monocyte', 'NK_cell']:
    mask = np.array(all_labels) == ct
    print(f"  {ct}: mean={adt_clr[mask, 0].mean():.2f}  (positive should be T/NK cells)")
```

## 3. Weighted Nearest Neighbor (WNN) Integration

### The core problem: combining two noisy views
RNA and protein measure partially overlapping cell identity signals. Naively concatenating normalized RNA and ADT matrices and running PCA would be dominated by RNA (2000+ features vs 20 proteins). A simple average of RNA and protein distances would give equal weight regardless of which modality is more informative for a given cell.

**WNN (Hao et al. 2021, Cell)** solves this by learning per-cell, per-modality weights:
- In cells where RNA gives a confident, well-clustered position → RNA gets higher weight
- In cells where RNA is ambiguous but protein is clear → protein gets higher weight

The result is a **within-cell adaptive weighting** that extracts the most information from each cell across modalities.

### WNN algorithm steps
1. Compute a **within-modality k-NN graph** (RNA k-NN and protein k-NN independently)
2. For each cell, estimate the **predictive power** of each modality by asking: how well can RNA neighbors predict protein expression, and vice versa?
3. Combine k-NN graphs using per-cell weights: 
   `WNN_affinity(i,j) = w_RNA(i) * RNA_affinity(i,j) + w_protein(i) * protein_affinity(i,j)`
4. UMAP is run on the WNN graph (not on individual modality embeddings)

### Why WNN is superior for CITE-seq
- **T cell sub-typing**: CD4 mRNA is low-expressed and noisy; CD4 protein is clear. WNN gives protein high weight for T cells.
- **Plasma cells**: very high RNA output (immunoglobulins), distinct from memory B cells in RNA space, but both are CD19+. RNA dominates appropriately here.
- **Monocyte/DC boundary**: CD11c protein distinguishes plasmacytoid DC from myeloid DC; RNA doesn't always separate them cleanly.

### muon Python implementation
```python
import muon as mu

# MuData object holds multiple AnnData modalities
mdata = mu.MuData({'rna': adata_rna, 'prot': adata_prot})

# RNA: standard scanpy preprocessing
mu.pp.filter_var(mdata['rna'], 'n_cells_by_counts', lb=10)
sc.pp.normalize_total(mdata['rna'], target_sum=1e4)
sc.pp.log1p(mdata['rna'])
sc.pp.pca(mdata['rna'], n_comps=30)
sc.pp.neighbors(mdata['rna'])

# Protein: CLR normalization + PCA
mu.prot.pp.clr(mdata['prot'])
sc.pp.pca(mdata['prot'], n_comps=15)
sc.pp.neighbors(mdata['prot'])

# WNN integration
mu.pp.neighbors(mdata, key_added='wnn')
sc.tl.umap(mdata, neighbors_key='wnn')
sc.pl.umap(mdata, color='leiden_wnn')
```

## 4. 10x Multiome: Paired RNA + ATAC from the Same Cell

### What is 10x Multiome?
The 10x Chromium Multiome ATAC + Gene Expression kit captures both **nuclear RNA** and **chromatin accessibility** from the same single nucleus in one experiment. This is fundamentally different from CITE-seq (RNA + protein) — here you get two views of the regulatory genome: what's accessible AND what's being transcribed.

### Key differences from CITE-seq
| Aspect | CITE-seq | 10x Multiome |
|--------|----------|-------------|
| Modalities | RNA + protein | RNA + ATAC |
| Cell vs nucleus | Cells (cytoplasm+nucleus) | Nuclei only |
| Protein panels | Up to ~200 targets | N/A |
| Regulatory inference | Limited | Direct (peak → gene) |
| Read requirement | 20k/cell RNA + 1k/cell ADT | 20k/cell RNA + 10k/cell ATAC |

### Why paired RNA+ATAC is powerful
Standard scATAC-seq gives you open chromatin peaks, but assigning these peaks to regulatory target genes requires assumptions. With paired RNA+ATAC from the same cells:
- Directly correlate peak accessibility with target gene expression (no assumptions)
- Identify **cell-type-specific regulatory elements**: peaks open in cell type A that correlate with a gene highly expressed in cell type A
- **TF footprinting validation**: find TF motif in open peak AND see TF expression in same cell
- **Gene regulation directionality**: peak opens before gene activates → regulatory (not reactive)

### muon for Multiome analysis
```python
import muon as mu

# Load 10x Multiome h5 file (contains both RNA and ATAC)
mdata = mu.read_10x_h5('filtered_feature_bc_matrix.h5')
# mdata['rna']  — AnnData with RNA counts
# mdata['atac'] — AnnData with peak counts

# Standard RNA preprocessing
mu.pp.filter_var(mdata['rna'], 'n_cells_by_counts', lb=10)
sc.pp.normalize_total(mdata['rna'])
sc.pp.log1p(mdata['rna'])
sc.pp.pca(mdata['rna'])

# ATAC: TF-IDF + LSI (drop first component)
mu.atac.pp.tfidf(mdata['atac'], scale_factor=1e4)
sc.tl.pca(mdata['atac'])  # LSI component 1 will be depth-correlated

# Joint UMAP (WNN concept applied to RNA + ATAC)
mu.pp.neighbors(mdata)
sc.tl.umap(mdata)
```

### Gene activity scores (bridging ATAC to RNA)
In the absence of paired RNA, scATAC clusters can be annotated using **gene activity scores**: sum the ATAC signal in the gene body + promoter (2kb upstream) for each gene, creating a pseudo-expression matrix. This proxy expression is much noisier than actual RNA but enables:
- Coarse cell type annotation using well-known marker genes
- Reference mapping to RNA atlases for label transfer
