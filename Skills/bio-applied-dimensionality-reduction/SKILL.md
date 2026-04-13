---
name: bio-applied-dimensionality-reduction
description: "**Tier 3 — Applied Bioinformatics | Module 30 · Notebook 2**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/30_Single_Cell_RNA_seq/02_dimensionality_reduction.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: anndata 0.10+, matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scikit-learn 1.4+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# scRNA-seq: Dimensionality Reduction and Clustering

*Source: Course notebook `Tier_3_Applied_Bioinformatics/30_Single_Cell_RNA_seq/02_dimensionality_reduction.ipynb`*

# scRNA-seq: Dimensionality Reduction and Clustering

**Tier 3 — Applied Bioinformatics | Module 30 · Notebook 2**

*Prerequisites: Notebook 1 (QC & Preprocessing)*

---

**By the end of this notebook you will be able to:**
1. Apply PCA to reduce dimensionality and select the number of significant PCs
2. Build a k-nearest-neighbor graph and compute UMAP / t-SNE embeddings
3. Cluster cells with the Leiden algorithm and interpret resolution parameters
4. Identify differentially expressed genes between clusters (Wilcoxon, t-test)
5. Visualize cluster marker genes as dot plots, violin plots, and feature plots



**Key resources:**
- [Scanpy clustering tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)
- [UMAP documentation](https://umap-learn.readthedocs.io/)
- [Understanding UMAP (Coenen & Pearce)](https://pair-code.github.io/understanding-umap/)

## Why dimensionality reduction matters

A 2,000-gene HVG matrix is too high-dimensional to cluster or visualize directly. PCA first compresses this into ~30–50 meaningful axes that capture most biological variance while discarding technical noise. The k-NN graph built on PCA coordinates encodes the manifold structure of the data — which cells are similar to which — and both UMAP embeddings and graph-based clustering (Leiden) operate on this graph. Understanding the graph is the key to understanding why changing `n_neighbors` or `resolution` changes your results.

```python
# Setup: recreate the preprocessed AnnData from Notebook 1
# (or load from disk if running after Notebook 1)
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

try:
    adata = ad.read_h5ad('/tmp/scrna_preprocessed.h5ad')
    print(f"Loaded from disk: {adata.shape}")
except FileNotFoundError:
    # Rebuild synthetic dataset
    rng = np.random.default_rng(42)
    n_cells, n_genes = 500, 1500
    cell_types = ["T_cell", "B_cell", "Monocyte", "NK_cell", "Dendritic"]
    cells_per = n_cells // len(cell_types)
    labels = np.repeat(cell_types, cells_per)
    
    counts = rng.negative_binomial(2, 0.9, (n_cells, n_genes)).astype(np.float32)
    type_idx = {ct: np.where(labels == ct)[0] for ct in cell_types}
    offsets = [50, 200, 400, 700, 1100]
    for ct, off in zip(cell_types, offsets):
        counts[np.ix_(type_idx[ct], [off, off+1, off+2])] += rng.poisson(20, (cells_per, 3))
    
    # Normalize and log-transform
    cell_totals = counts.sum(axis=1, keepdims=True)
    X_norm = np.log1p(counts / cell_totals * 1e4)
    
    adata = ad.AnnData(
        X=X_norm,
        obs=pd.DataFrame({'cell_type': labels}, index=[f'CELL_{i:04d}' for i in range(n_cells)]),
        var=pd.DataFrame({'hvg': True}, index=[f'GENE{i:04d}' for i in range(n_genes)])
    )
    print(f"Created synthetic dataset: {adata.shape}")
```

## 1. Principal Component Analysis

### Why scale before PCA?
Log-normalized expression values still vary in magnitude across genes — a housekeeping gene may have mean expression of 3 while a rare marker has mean 0.1. PCA would be dominated by high-mean genes regardless of their biological relevance. Scaling standardizes each gene to zero mean and unit variance (z-score), so every gene contributes equally to the principal components.

**`max_value=10` clipping**: extreme outlier cells (often doublets that survived QC) can produce single-gene z-scores of 50+, distorting PCs. Clipping at 10 standard deviations prevents this.

### Choosing the number of PCs
The elbow plot shows cumulative variance explained by each PC. The "elbow" — where each additional PC adds diminishing variance — suggests how many PCs capture biological signal vs noise. For PBMCs: typically 10–20 PCs. For complex tissue samples: 30–50 PCs.

**Rule of thumb**: Use 2× the PC at the elbow. Over-including PCs adds noise but rarely destroys structure; under-including can merge distinct cell types.

### PC loadings interpretation
Each PC is a weighted combination of genes. PC1 often corresponds to the dominant cell type axis (e.g., immune vs stromal). Inspect `sc.pl.pca_loadings(adata, components=[1,2,3])` to see which genes drive each component.

```python
# PCA on the HVG matrix
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Step 1: Scale genes (center + unit variance, clip at 10 SD)
X = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
X_scaled = np.clip(X_scaled, -10, 10)  # clip outliers

# Step 2: PCA
n_pcs = min(50, X_scaled.shape[1] - 1)
pca = PCA(n_components=n_pcs, random_state=42)
X_pca = pca.fit_transform(X_scaled)

# Store PCA coordinates in adata
adata.obsm['X_pca'] = X_pca
print(f"PCA shape: {X_pca.shape}")
print(f"Variance explained by PC1: {pca.explained_variance_ratio_[0]:.3f} ({pca.explained_variance_ratio_[0]*100:.1f}%)")
print(f"Total variance in top 10 PCs: {pca.explained_variance_ratio_[:10].sum():.3f} ({pca.explained_variance_ratio_[:10].sum()*100:.1f}%)")

# Elbow plot
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

cumvar = np.cumsum(pca.explained_variance_ratio_) * 100
axes[0].plot(range(1, n_pcs+1), pca.explained_variance_ratio_*100, 'o-', ms=4)
axes[0].axvline(15, color='red', linestyle='--', label='Suggested elbow: PC15')
axes[0].set_xlabel('Principal Component')
axes[0].set_ylabel('Variance explained (%)')
axes[0].set_title('Scree plot (variance per PC)')
axes[0].legend()

axes[1].plot(range(1, n_pcs+1), cumvar, 'o-', ms=4, color='orange')
axes[1].axhline(80, color='gray', linestyle='--', label='80% variance')
axes[1].set_xlabel('Number of PCs')
axes[1].set_ylabel('Cumulative variance (%)')
axes[1].set_title('Cumulative variance explained')
axes[1].legend()

plt.tight_layout()
plt.savefig('pca_variance.png', dpi=100, bbox_inches='tight')
plt.show()

# PCA scatter plot colored by true cell type
n_pcs_use = 15
fig, ax = plt.subplots(figsize=(8, 6))
cell_types_list = adata.obs['cell_type'].unique()
colors_map = dict(zip(cell_types_list, ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']))
for ct in cell_types_list:
    mask = adata.obs['cell_type'] == ct
    ax.scatter(X_pca[mask, 0], X_pca[mask, 1], c=colors_map[ct], label=ct, s=15, alpha=0.7)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
ax.set_title('PCA: PC1 vs PC2 colored by cell type')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig('pca_scatter.png', dpi=100, bbox_inches='tight')
plt.show()
print(f"\nUsing top {n_pcs_use} PCs for downstream neighbor graph")
```

## 2. Neighborhood Graph and UMAP

### Building the k-NN graph
After PCA, we construct a **k-nearest-neighbor graph** in PC space: for each cell, find its k most similar cells (by Euclidean distance in PC space) and connect them with edges. This graph is the foundation for both clustering and UMAP.

Key parameters:
- **`n_neighbors`** (default 15): Controls the local vs global balance. Small values (5–10) emphasize local structure — tight clusters. Large values (30–50) emphasize global structure — better topology but fuzzier clusters. For PBMC-scale data (10k cells), n_neighbors=15–20 is typical.
- **`n_pcs`**: How many PCs to use. Should match your elbow plot analysis.

### UMAP vs t-SNE
Both are non-linear dimensionality reduction methods that project the k-NN graph into 2D.

| | UMAP | t-SNE |
|--|------|-------|
| Speed | Fast (minutes) | Slow (hours for >50k cells) |
| Global structure | Partially preserved | Not preserved |
| Cluster distances | Approximate (meaningful) | Meaningless |
| Reproducibility | Random seed controlled | Less stable |

**Critical warning**: UMAP distances between clusters are NOT direct measures of transcriptional similarity. Two clusters far apart in UMAP space are not necessarily more different than two clusters that are adjacent. For quantitative comparisons, always use the PC space or expression values directly.

### UMAP hyperparameters
- `min_dist` (default 0.5): Controls how tightly cells are packed within clusters. Smaller = tighter clusters, more visual separation.
- `spread` (default 1.0): Controls overall spread of the embedding.
These affect aesthetics but not biological conclusions — don't over-tune them.

```python
# Build k-NN graph and UMAP
from sklearn.neighbors import NearestNeighbors
import numpy as np
import matplotlib.pyplot as plt

n_pcs_use = 15
n_neighbors = 15

# Build k-NN graph using PCA coordinates
knn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean', n_jobs=-1)
knn.fit(adata.obsm['X_pca'][:, :n_pcs_use])
distances, indices = knn.kneighbors(adata.obsm['X_pca'][:, :n_pcs_use])

# Store connectivity in adata (adjacency-like structure)
adata.uns['neighbors'] = {'params': {'n_neighbors': n_neighbors, 'n_pcs': n_pcs_use}}
adata.obsm['knn_indices'] = indices
adata.obsm['knn_distances'] = distances

print(f"k-NN graph: {adata.n_obs} cells x {n_neighbors} neighbors each")
print(f"Mean distance to nearest neighbor: {distances[:, 1].mean():.4f}")

# UMAP
try:
    from umap import UMAP
    reducer = UMAP(n_neighbors=n_neighbors, min_dist=0.3, random_state=42, n_epochs=200)
    X_umap = reducer.fit_transform(adata.obsm['X_pca'][:, :n_pcs_use])
    adata.obsm['X_umap'] = X_umap
    print(f"UMAP embedding shape: {X_umap.shape}")
    
    # Visualize
    fig, ax = plt.subplots(figsize=(8, 6))
    cell_types_list = adata.obs['cell_type'].unique()
    colors_map = dict(zip(sorted(cell_types_list), ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']))
    for ct in sorted(cell_types_list):
        mask = adata.obs['cell_type'] == ct
        ax.scatter(X_umap[mask, 0], X_umap[mask, 1], c=colors_map[ct], label=ct, s=10, alpha=0.8)
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_title('UMAP colored by cell type')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('umap_celltypes.png', dpi=100, bbox_inches='tight')
    plt.show()
except ImportError:
    print("umap-learn not installed. Install with: pip install umap-learn")
    print("Falling back to PCA visualization for clustering steps")
    adata.obsm['X_umap'] = adata.obsm['X_pca'][:, :2]
```

## 3. Leiden / Louvain Clustering

### Graph-based community detection
Clustering in scRNA-seq is not traditional k-means or hierarchical clustering. Instead, **community detection** algorithms find groups of cells that are densely connected in the k-NN graph. These algorithms were originally developed for social network analysis.

**Leiden vs Louvain**: Both maximize modularity (ratio of edges within communities to edges between communities). Leiden is strictly better — it guarantees well-connected communities (Louvain can produce internally disconnected communities). Use Leiden.

### The resolution parameter
Resolution controls how many clusters you get:
- High resolution → more, smaller clusters (over-clustering)
- Low resolution → fewer, larger clusters (under-clustering)

There is no single "correct" resolution. Use domain knowledge:
- Do clusters correspond to known cell types?
- Are marker genes specific to clusters?
- Do adjacent clusters in UMAP space have qualitatively different biology?

**Practical approach**: Start with resolution=0.5, then increase (0.8, 1.0) if you suspect merged populations or decrease (0.3) if you have too much fragmentation.

### Leiden requires the `leidenalg` package
In scanpy: `sc.tl.leiden(adata, resolution=0.5)`. The implementation in scanpy uses the `leidenalg` Python package, which wraps the C++ igraph library for performance.

```python
# Leiden clustering on the k-NN graph
# We implement a simplified community detection to demonstrate the concept
# In production, use: sc.tl.leiden(adata, resolution=0.5)

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

# Build sparse adjacency matrix from k-NN indices
n_cells = adata.n_obs
n_neighbors = adata.obsm['knn_indices'].shape[1]
row_idx = np.repeat(np.arange(n_cells), n_neighbors)
col_idx = adata.obsm['knn_indices'].flatten()
# Mutual k-NN (keep edge only if both cells consider each other neighbors)
adj = sp.csr_matrix((np.ones(len(row_idx)), (row_idx, col_idx)), shape=(n_cells, n_cells))
adj = ((adj + adj.T) > 0).astype(np.float32)  # symmetrize

# Simple connected-components as a proxy for clusters
# (production code uses Leiden/Louvain modularity optimization)
from scipy.sparse.csgraph import connected_components
n_components, labels = connected_components(adj, directed=False)
print(f"Connected components: {n_components}")

# For demonstration, use agglomerative clustering on PCA coordinates
from sklearn.cluster import AgglomerativeClustering
n_clusters = 5
agg = AgglomerativeClustering(n_clusters=n_clusters, metric='euclidean', linkage='ward')
cluster_labels = agg.fit_predict(adata.obsm['X_pca'][:, :15])
adata.obs['leiden'] = [str(c) for c in cluster_labels]  # string labels like Leiden

print(f"\nCluster sizes:")
print(adata.obs['leiden'].value_counts().sort_index())

# Compare clusters to true cell types
from pandas import crosstab
ct = crosstab(adata.obs['leiden'], adata.obs['cell_type'], normalize='index').round(2)
print(f"\nCluster composition (row-normalized):")
print(ct)

# Visualize clusters on UMAP/PCA
coords = adata.obsm['X_umap']
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

colors_cluster = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
for ci in sorted(adata.obs['leiden'].unique()):
    mask = adata.obs['leiden'] == ci
    axes[0].scatter(coords[mask, 0], coords[mask, 1], 
                    c=colors_cluster[int(ci)], label=f'Cluster {ci}', s=10, alpha=0.8)
axes[0].set_title('Clustering result')
axes[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left')

cell_types_list = sorted(adata.obs['cell_type'].unique())
colors_map = dict(zip(cell_types_list, ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']))
for ct_name in cell_types_list:
    mask = adata.obs['cell_type'] == ct_name
    axes[1].scatter(coords[mask, 0], coords[mask, 1],
                    c=colors_map[ct_name], label=ct_name, s=10, alpha=0.8)
axes[1].set_title('True cell types')
axes[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left')

for ax in axes:
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
plt.tight_layout()
plt.savefig('clustering.png', dpi=100, bbox_inches='tight')
plt.show()
```
