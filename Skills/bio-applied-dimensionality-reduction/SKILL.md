---
name: bio-applied-dimensionality-reduction
description: "scRNA-seq dimensionality reduction and clustering: PCA, k-NN graph, UMAP, Leiden. Parameter selection guide, implementation patterns, and pitfalls."
tool_type: python
primary_tool: NumPy
---

# scRNA-seq: Dimensionality Reduction and Clustering

## Pipeline Overview

```
HVG matrix (cells × ~2000 genes)
  → Scale (z-score, clip ±10)
  → PCA (top 10–50 PCs)
  → k-NN graph (n_neighbors=15–20)
  → UMAP (visualization only)
  → Leiden clustering (on the graph)
```

## PCA

```python
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

X = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
X_scaled = StandardScaler().fit_transform(X)
X_scaled = np.clip(X_scaled, -10, 10)  # prevent outlier/doublet dominance

pca = PCA(n_components=50, random_state=42)
X_pca = pca.fit_transform(X_scaled)
adata.obsm['X_pca'] = X_pca
```

**Choosing n_pcs:** Plot cumulative variance; use 2× the "elbow" PC. Under-including merges cell types; over-including adds noise but rarely destroys structure. Typical: 10–20 PCs for PBMC, 30–50 for complex tissue.

## k-NN Graph

```python
from sklearn.neighbors import NearestNeighbors

knn = NearestNeighbors(n_neighbors=15, metric='euclidean', n_jobs=-1)
knn.fit(X_pca[:, :n_pcs])
distances, indices = knn.kneighbors(X_pca[:, :n_pcs])
```

**`n_neighbors` effect:**
- Small (5–10): emphasizes local structure, tight clusters, may fragment continuous populations
- Large (30–50): emphasizes global structure, smoother but may merge distinct types
- Default 15: good starting point for 5k–50k cells

## UMAP

```python
from umap import UMAP

X_umap = UMAP(n_neighbors=15, min_dist=0.3, random_state=42).fit_transform(X_pca[:, :n_pcs])
adata.obsm['X_umap'] = X_umap
```

**UMAP vs t-SNE:**

| | UMAP | t-SNE |
|--|------|-------|
| Speed | Fast (minutes) | Slow (hours for >50k cells) |
| Global structure | Partially preserved | Not preserved |
| Between-cluster distances | Approximate (directionally meaningful) | Meaningless |
| Reproducibility | Seed-controlled | Less stable |

## Leiden Clustering (scanpy)

```python
import scanpy as sc

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
sc.tl.leiden(adata, resolution=0.5)
# result: adata.obs['leiden'] — string cluster labels
```

**Resolution tuning:**
- Start at 0.5; increase to 0.8–1.0 if known cell types are merged; decrease to 0.3 if over-fragmented
- Validate: do clusters have specific marker genes? Do adjacent UMAP clusters have different biology?

## Marker Gene Identification

```python
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)
```

## Parameter Decision Table

| Parameter | Default | Increase when | Decrease when |
|-----------|---------|--------------|---------------|
| `n_pcs` | 15–20 | Complex tissue, many cell types | Simple dataset, elbow is early |
| `n_neighbors` | 15 | Need global structure | Need fine local detail |
| `min_dist` | 0.5 | Clusters too compressed | Want tighter visual clusters |
| Leiden `resolution` | 0.5 | Known subtypes being merged | Too many spurious clusters |

## Pitfalls

- **UMAP distances are not transcriptional distances**: two clusters far apart in UMAP are not necessarily more different. Use PCA space or expression values for quantitative comparisons.
- **Scaling required before PCA**: high-mean housekeeping genes dominate PCs without scaling, burying rare-marker variation.
- **Clip at ±10 SD**: doublets that survived QC can have extreme z-scores (50+) that distort PCs. Clipping at 10 is standard.
- **Leiden requires `leidenalg`**: `pip install leidenalg`. Louvain can produce internally disconnected communities — use Leiden.
- **Resolution is not portable**: resolution=0.5 on one dataset ≠ same number of clusters on another; always re-tune per dataset.
- **Batch effects appear as clusters**: run batch correction (Harmony, scVI) before clustering if samples from different batches.
- **Multiple testing**: apply FDR correction (Benjamini-Hochberg) for marker gene testing
