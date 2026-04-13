---
name: bio-applied-sc-integration
description: "**Tier 3 — Applied Bioinformatics | Module 31 · Notebook 3**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/31_Single_Cell_Multi_Omics/03_sc_integration.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Single-Cell Batch Correction and Dataset Integration

*Source: Course notebook `Tier_3_Applied_Bioinformatics/31_Single_Cell_Multi_Omics/03_sc_integration.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 31 · Notebook 3**

*Prerequisites: Module 30 (scRNA-seq), Notebook 2 (CITE-seq)*

---

**By the end of this notebook you will be able to:**
1. Diagnose batch effects in scRNA-seq data using mixing metrics and kBET
2. Apply Harmony, scVI, and Seurat CCA integration to correct batch effects
3. Evaluate integration quality: cell type mixing vs biological signal preservation
4. Perform label transfer between reference atlas and query datasets
5. Map query cells onto a reference atlas (Azimuth / CELLxGENE Census)



**Key resources:**
- [scib benchmarking (Luecken et al. 2022)](https://www.nature.com/articles/s41592-021-01336-8)
- [Harmony documentation](https://portals.broadinstitute.org/harmony/)
- [scVI-tools documentation](https://docs.scvi-tools.org/)
- [Azimuth reference mapping](https://azimuth.hubmapconsortium.org/)

## 1. Diagnosing Batch Effects

### What is a batch effect?
A batch effect is any systematic technical variation that is correlated with a non-biological variable: different days of library preparation, different operators, different 10x kits, different sequencing runs, or different laboratories. Batch effects manifest as cells from different batches separating in UMAP space even when they should be the same cell type.

### Why batch effects are severe in scRNA-seq
- Minute differences in cell viability at capture time affect thousands of stress-response genes simultaneously
- Different ambient RNA concentrations between batches contaminate different cell fractions
- Protocol differences (enzyme lots, temperature) affect cDNA amplification non-uniformly
- Even samples processed on the same day but on different chips can show separation

### Visual diagnosis
The first diagnostic: color your UMAP by batch. If clusters separate by batch (a distinct "batch shape" overlapping true biology), you have a batch effect requiring correction.

**Red flags:**
- Cluster that contains cells from only one batch
- Same cell type in two batches forms two separate UMAP blobs
- Leiden clustering gives clusters that are 90%+ from one batch

### Quantitative metrics
Simply looking at UMAP is subjective. Quantitative metrics give objective measurements:

**LISI (Local Inverse Simpson's Index, Korsunsky et al. 2019)**:
- For each cell, look at its k nearest neighbors
- Compute the inverse Simpson's diversity index of batch labels among those neighbors
- High LISI (close to number of batches) = well-mixed = good integration
- Low LISI (close to 1) = neighbors are all from same batch = segregated batches

**ASW (Average Silhouette Width)**:
- For batch integration: lower ASW-batch = more mixed (desired)
- For cell type preservation: higher ASW-celltype = cleaner separation (desired)
- These two objectives conflict — good integration must balance both

**kBET (k-nearest neighbor Batch Effect Test)**:
- For each cell, test whether batch composition in its k-NN matches the global batch composition (chi-square test)
- Reports rejection rate: high = significant batch effect; low = well-integrated

## 2. Harmony Integration

### How Harmony works (Korsunsky et al. 2019, *Nature Methods*)
Harmony operates in PCA space and iteratively adjusts the PCA embedding to remove batch effects while preserving biological structure:

1. **Initialize**: start with the standard PCA embedding
2. **Cluster** cells into soft clusters (fuzzy k-means) based on current embedding
3. **Compute correction**: for each cluster, compute the batch centroid offsets — how much cells from Batch B are shifted relative to Batch A in that cluster
4. **Apply correction**: shift cells to remove the batch offset, merging cells from different batches that belong to the same cluster
5. **Iterate** until convergence (typically 10 iterations)

The result is a corrected PCA embedding (`adata.obsm['X_pca_harmony']`) that retains cell type structure while removing batch separation.

### Key properties of Harmony
- **In-memory**: operates on the PCA matrix, not the full count matrix → very fast (seconds to minutes)
- **Preserves sparse cell types**: cells from rare populations in one batch won't be overcorrected because the soft clusters ensure they influence only their own cluster's correction
- **No data imputation**: Harmony does NOT modify the count matrix — it only adjusts the embedding. Differential expression must still be done with the raw counts.
- **Python package**: `harmonypy` (install: `pip install harmonypy`)

### When Harmony may fail
- **Compositional imbalance**: if Batch A has only T cells and Batch B has only B cells, Harmony cannot find corresponding cell types to align — the two batches simply cannot be integrated
- **Very strong batch effects**: technical variation larger than biological variation may require deeper models (scVI)
- **Dataset-specific effects**: when batches have different cell type frequencies AND strong batch effects, over-correction may merge truly different populations

### Harmony pitfall: over-correction
Setting `theta` too high (the diversity penalty) can force batch mixing even when batches genuinely contain different cell types (e.g., disease vs healthy with truly different immune compositions). Recommended `theta` range: 1–3 (default = 2).

```python
# Harmony batch correction: implemented from scratch to show the core idea
# (Production: use harmonypy — pip install harmonypy)
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

def harmony_correct(X_pca, batch_labels, n_clusters=20, max_iter=10, theta=2.0, sigma=0.1):
    """
    Simplified Harmony: iterative soft-clustering batch correction in PCA space.
    X_pca: (n_cells, n_pcs) — initial PCA embedding
    batch_labels: array of batch identifiers per cell
    """
    X = X_pca.copy()
    unique_batches = np.unique(batch_labels)
    n_batches = len(unique_batches)
    batch_matrix = np.column_stack([(batch_labels == b).astype(float) for b in unique_batches])
    
    for iteration in range(max_iter):
        # --- Step 1: Soft cluster assignment (fuzzy k-means) ---
        kmeans = KMeans(n_clusters=n_clusters, random_state=iteration, n_init=1)
        hard_assign = kmeans.fit_predict(X)
        centroids = kmeans.cluster_centers_
        
        # Soft assignment via inverse distance
        dists = np.array([np.linalg.norm(X - c, axis=1) for c in centroids]).T  # (n_cells, K)
        soft_assign = np.exp(-dists**2 / (2 * sigma**2))
        soft_assign /= soft_assign.sum(axis=1, keepdims=True) + 1e-10
        
        # --- Step 2: Compute batch centroids per cluster ---
        # For each cluster k and batch b: mean position of cells from batch b in cluster k
        correction = np.zeros_like(X)
        for k in range(n_clusters):
            weights = soft_assign[:, k]  # per-cell weight for cluster k
            cluster_center = (X * weights[:, None]).sum(0) / (weights.sum() + 1e-10)
            
            for b_idx, b in enumerate(unique_batches):
                b_mask = batch_labels == b
                b_weights = weights * b_mask
                if b_weights.sum() < 1e-6:
                    continue
                b_center = (X * b_weights[:, None]).sum(0) / (b_weights.sum() + 1e-10)
                # Correction: move cells toward cluster center, away from batch-specific center
                offset = cluster_center - b_center
                correction[b_mask] += (weights[b_mask, None] * offset[None, :] * theta / n_batches)
        
        X = X + correction * 0.3  # learning rate
    
    return X

# Apply simplified Harmony
print("Running simplified Harmony correction...")
X_pca_harmony = harmony_correct(
    adata.obsm['X_pca'][:, :20], 
    batches_arr, 
    n_clusters=20, max_iter=8, theta=2.0
)
adata.obsm['X_pca_harmony'] = X_pca_harmony
print("Done.")

# UMAP after correction
try:
    from umap import UMAP
    X_umap_post = UMAP(n_neighbors=15, min_dist=0.3, random_state=42).fit_transform(X_pca_harmony)
    adata.obsm['X_umap_harmony'] = X_umap_post

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    batch_colors = {'Batch_A':'#e41a1c','Batch_B':'#377eb8','Batch_C':'#4daf4a'}
    ct_colors = {'T_cell':'#e41a1c','B_cell':'#377eb8','Monocyte':'#4daf4a','NK_cell':'#984ea3'}
    
    X_umap_pre = adata.obsm['X_umap_uncorrected']
    
    for batch in batch_effects:
        mask = batches_arr == batch
        axes[0,0].scatter(X_umap_pre[mask,0], X_umap_pre[mask,1],
                          c=batch_colors[batch], label=batch, s=6, alpha=0.6)
        axes[1,0].scatter(X_umap_post[mask,0], X_umap_post[mask,1],
                          c=batch_colors[batch], label=batch, s=6, alpha=0.6)
    
    for ct in cell_types:
        mask = labels_arr == ct
        axes[0,1].scatter(X_umap_pre[mask,0], X_umap_pre[mask,1],
                          c=ct_colors[ct], label=ct, s=6, alpha=0.6)
        axes[1,1].scatter(X_umap_post[mask,0], X_umap_post[mask,1],
                          c=ct_colors[ct], label=ct, s=6, alpha=0.6)
    
    axes[0,0].set_title('BEFORE correction: by batch'); axes[0,0].legend(fontsize=8)
    axes[0,1].set_title('BEFORE correction: by cell type'); axes[0,1].legend(fontsize=8)
    axes[1,0].set_title('AFTER Harmony: by batch (should be mixed)'); axes[1,0].legend(fontsize=8)
    axes[1,1].set_title('AFTER Harmony: by cell type (should be preserved)'); axes[1,1].legend(fontsize=8)
    for ax in axes.flatten():
        ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    plt.tight_layout()
    plt.savefig('harmony_correction.png', dpi=100, bbox_inches='tight')
    plt.show()
except ImportError:
    print("umap-learn not installed")

# Quantify improvement
lisi_batch_post = compute_lisi(X_pca_harmony, batches_arr, n_neighbors=30)
lisi_ct_post    = compute_lisi(X_pca_harmony, labels_arr, n_neighbors=30)
print(f"\nBefore Harmony: LISI-batch={lisi_batch.mean():.3f}, LISI-celltype={lisi_ct.mean():.3f}")
print(f"After  Harmony: LISI-batch={lisi_batch_post.mean():.3f}, LISI-celltype={lisi_ct_post.mean():.3f}")
print(f"  -> LISI-batch increase = better batch mixing")
print(f"  -> LISI-celltype should remain similar = cell types preserved")
```python

## 3. scVI: Deep Generative Integration

### What is scVI? (Lopez et al. 2018, *Nature Methods*)
scVI (Single-Cell Variational Inference) is a **variational autoencoder (VAE)** trained on raw count data. Unlike Harmony (which works in PCA space), scVI learns a probabilistic model directly from counts, using a negative binomial likelihood to account for overdispersion.

### Architecture
- **Encoder**: 2-layer MLP (gene expression → 128-dim latent) with batch and library size as conditional inputs
- **Latent space**: 10-dimensional Gaussian distribution (mean + variance)
- **Decoder**: reconstructs expected counts per gene, with batch-specific dispersion

The batch label is fed as a covariate to the encoder and decoder. During training, the model is penalized when the latent space retains batch information — it must encode biology, not batch.

### scVI advantages over Harmony
- **Probabilistic**: gives uncertainty estimates for each cell's position
- **Handles raw counts**: no need to normalize/log-transform before scVI
- **Native differential expression**: scVI can compute DE on the latent space using normalized expression estimates, properly accounting for batch
- **Scalable**: supports GPU acceleration; handles millions of cells with data loaders

### scVI installation and usage
```bash
pip install scvi-tools
```python

```python
import scvi

# Prepare data: scVI works on raw integer counts
scvi.model.SCVI.setup_anndata(adata, batch_key='batch', layer='raw_counts')

# Build and train the model
model = scvi.model.SCVI(adata, n_layers=2, n_latent=10, gene_likelihood='nb')
model.train(max_epochs=400, early_stopping=True)

# Get batch-corrected latent representation
adata.obsm['X_scVI'] = model.get_latent_representation()

# Differential expression (handles batch properly)
de_df = model.differential_expression(
    adata,
    groupby='cell_type',
    group1='T_cell', group2='B_cell',
    delta=0.25  # minimum effect size (log2FC)
)
```python

### scANVI: labeled integration
scANVI (scVI + semi-supervision) extends scVI by using available cell type labels as training signal. Labeled cells help anchor the latent space, making integration more biologically meaningful. Used for label transfer between datasets.

### BBKNN: graph-based integration
BBKNN (Batch Balanced k-NN, Polański et al. 2020) works differently:
1. For each cell, find k nearest neighbors from EACH batch separately
2. Merge these batch-balanced neighborhoods into a single graph
3. Run UMAP on this graph

This ensures each cell has neighbors from all batches, forcing cross-batch connections. Very fast (faster than Harmony), works directly in PCA space.

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
