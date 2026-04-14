---
name: bio-applied-cancer-transcriptomics
description: Cancer transcriptomics — melanoma subtype classification (Tirosh/Harbst), preprocessing pipeline, PCA/t-SNE, hierarchical clustering, random forest, and Kaplan-Meier survival analysis
tool_type: python
primary_tool: Matplotlib
---

# Cancer Transcriptomics: Subtype Classification

## Key References
- Tirosh et al. (2016) *Science* — single-cell dissection of melanoma intratumoral heterogeneity
- Harbst et al. (2016) *Clin Cancer Res* — four molecular subtypes from bulk RNA-seq (TCGA-SKCM)
- Data: [cBioPortal TCGA-SKCM](https://www.cbioportal.org/study/summary?id=skcm_tcga)

## Melanoma Subtype Marker Genes

| Subtype | Key Markers | Classification |
|---------|-------------|----------------|
| Pigmentation (MITF-high) | MITF, DCT, TYRP1, MLANA | Tirosh + Harbst |
| Keratin | KRT5, KRT14, KRT6A, EGFR | Tirosh |
| Immune | CD3D, CD8A, GZMB, PRF1 | Tirosh + Harbst |
| Proliferative | MKI67, TOP2A, CDK1 | Harbst |
| Normal-like | VIM, CDH2, FN1 | Harbst |

## Preprocessing Pipeline

```python
import numpy as np
import pandas as pd

# 1. Load: genes as rows, samples as columns (cBioPortal TSV format)
expr_df = pd.read_csv("data_mrna_seq_v2_rsem.txt", sep="\t", index_col=0)

# 2. log1p transform — stabilizes variance, handles zeros
expr_log = np.log1p(expr_df)

# 3. Z-score per gene across samples
expr_z = (expr_log
    .subtract(expr_log.mean(axis=1), axis=0)
    .divide(expr_log.std(axis=1) + 1e-8, axis=0))

# 4. Variance filtering — retain top N most variable genes
top_genes = expr_z.var(axis=1).nlargest(1500).index
expr_top = expr_z.loc[top_genes]  # shape: (1500, n_samples)
```

## Exploratory Analysis

```python
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import seaborn as sns

X = expr_top.values.T  # (n_samples, n_genes)

# PCA
pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X)
# pca.explained_variance_ratio_ — report cumulative variance

# t-SNE (use PCA output as input for large datasets)
tsne = TSNE(n_components=2, perplexity=30, random_state=42)
X_tsne = tsne.fit_transform(X_pca[:, :50] if X.shape[0] > 500 else X)

# Clustermap of top 50 genes
subtype_palette = {'MITF-low': '#4C72B0', 'Keratin': '#DD8452', 'Immune': '#55A868'}
col_colors = pd.Series(subtype_labels, index=sample_names).map(subtype_palette)
sns.clustermap(expr_z.loc[top_genes[:50]], col_colors=col_colors,
               cmap='RdBu_r', center=0, vmin=-3, vmax=3,
               xticklabels=False, figsize=(12, 8))
```

## Tirosh 3-class: Random Forest Classification

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

X_train, X_test, y_train, y_test = train_test_split(
    X, subtype_labels, test_size=0.30, random_state=42, stratify=subtype_labels
)
rf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
rf.fit(X_train, y_train)
print(classification_report(y_test, rf.predict(X_test)))

# Feature importance — top marker genes
importance = pd.Series(rf.feature_importances_, index=top_genes).nlargest(20)
```

## Harbst 4-class: Hierarchical Clustering

```python
from sklearn.cluster import AgglomerativeClustering

harbst_markers = ['MITF','DCT','TYRP1','MLANA',
                  'MKI67','TOP2A','CDK1',
                  'VIM','CDH2','FN1',
                  'CD3D','CD8A','GZMB','PRF1']

X_harbst = expr_z.loc[harbst_markers].values.T  # (n_samples, 14)
cluster_ids = AgglomerativeClustering(n_clusters=4, linkage='ward').fit_predict(X_harbst)

# Assign labels by cluster centroid expression
centroids = pd.DataFrame(X_harbst, columns=harbst_markers)
centroids['cluster'] = cluster_ids
ctrs = centroids.groupby('cluster').mean()

subtype_map = {
    ctrs[['MITF','DCT','TYRP1','MLANA']].mean(axis=1).idxmax(): 'Pigmentation',
    ctrs[['MKI67','TOP2A','CDK1']].mean(axis=1).idxmax(): 'Proliferative',
    ctrs[['CD3D','CD8A','GZMB','PRF1']].mean(axis=1).idxmax(): 'High-immune',
}
remaining = [c for c in range(4) if c not in subtype_map]
if remaining:
    subtype_map[remaining[0]] = 'Normal-like'
harbst_labels = [subtype_map.get(c, f'Cluster_{c}') for c in cluster_ids]
```

## Kaplan-Meier Survival Analysis

```python
# With lifelines (preferred)
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

fig, ax = plt.subplots(figsize=(9, 5))
for subtype in ['MITF-low', 'Keratin', 'Immune']:
    mask = subtype_labels == subtype
    kmf = KaplanMeierFitter()
    kmf.fit(survival_times[mask], events[mask], label=subtype)
    kmf.plot_survival_function(ax=ax, ci_show=True)

# Log-rank test between two groups
res = logrank_test(t1, t2, e1, e2)
print(f"p = {res.p_value:.4f}")
```

## Pitfalls

- **log1p before z-score, not after**: log1p compresses dynamic range first; z-scoring after ensures comparable gene scales; reversing the order gives distorted variance estimates
- **Variance filter on z-scored data is circular**: compute gene variance on log1p (not z-scored) data for filtering, then z-score the filtered matrix
- **t-SNE is not deterministic**: always set `random_state`; t-SNE axes have no biological meaning — distances between clusters are not interpretable, only local neighborhood structure is
- **Ward linkage requires Euclidean distance**: `AgglomerativeClustering(linkage='ward')` assumes Euclidean; for correlation-based clustering use `linkage='average'` with precomputed distance matrix
- **Kaplan-Meier requires censoring indicator**: `events` must be boolean (True = event observed); patients lost to follow-up are censored (False), not dead — mislabeling censored as events inflates mortality estimates
- **Batch effects in TCGA**: TCGA-SKCM has plate/batch effects; always check PC1/PC2 coloring by batch variable before biological interpretation; use ComBat or limma's `removeBatchEffect` if needed
