---
name: cancer-transcriptomics
description: Cancer transcriptomics subtype classification — variance filtering, hierarchical clustering, semi-supervised Random Forest, Kaplan-Meier survival analysis, TCGA data patterns
---

## When to Use

Use this skill when:
- Classifying tumors into molecular subtypes from bulk RNA-seq expression matrices
- Applying the Tirosh semi-supervised framework (labeled reference set → Random Forest → predict unlabeled)
- Running the Harbst hierarchical clustering approach (marker genes + Ward linkage → 4 subtypes)
- Performing Kaplan–Meier survival analysis to compare subtype prognoses
- Assessing concordance between two classification schemes with NMI and cross-tabulation
- Loading and preprocessing TCGA/cBioPortal expression data (log1p + z-score normalization)

## Quick Reference

| Method | Use Case | Key Parameter |
|---|---|---|
| Variance filtering | Reduce noise before clustering/PCA | top 1,500–2,000 genes by var(axis=1) |
| seaborn clustermap | Visualize gene–sample expression patterns | `cmap='RdBu_r'`, `col_colors` for subtypes |
| PCA | Linear structure, explained variance | `n_components=2`, `explained_variance_ratio_` |
| t-SNE | Nonlinear local structure | `perplexity=30`, `random_state=42` |
| Random Forest (Tirosh) | Semi-supervised 3-class assignment | `n_estimators=200`, `stratify=y` in split |
| AgglomerativeClustering (Harbst) | 4-class from marker genes | `linkage='ward'`, `n_clusters=4` |
| KaplanMeierFitter | Per-subtype survival curves | `ci_show=True` |
| normalized_mutual_info_score | Concordance between two schemes | returns 0–1 |

| Tool | Purpose |
|---|---|
| `lifelines` | Kaplan–Meier, log-rank test, Cox PH |
| `sklearn.ensemble.RandomForestClassifier` | Semi-supervised subtype prediction |
| `sklearn.cluster.AgglomerativeClustering` | Ward linkage hierarchical clustering |
| `sklearn.metrics.normalized_mutual_info_score` | Compare two label assignments |
| `seaborn.clustermap` | Heatmap with row/column dendrograms |
| `sklearn.decomposition.PCA` | Principal component analysis |
| `sklearn.manifold.TSNE` | t-distributed stochastic neighbor embedding |

## Key Patterns

**Pattern 1: Loading cBioPortal expression data**
```python
import pandas as pd
import numpy as np

# cBioPortal tab-delimited format: Hugo_Symbol + Entrez_Gene_Id columns, then samples
df = pd.read_csv('data_mrna_seq_v2_rsem.txt', sep='\t', index_col=0)
df = df.drop(columns=['Entrez_Gene_Id'], errors='ignore')
df.index.name = 'gene'

# Remove duplicate gene symbols (keep highest mean)
df = df.loc[~df.index.duplicated(keep='first')]
print(f'Loaded: {df.shape[0]} genes × {df.shape[1]} samples')
```

**Pattern 2: Variance filtering**
```python
def filter_high_variance_genes(expr_df: pd.DataFrame, n: int = 1500) -> pd.DataFrame:
    """Keep top-n genes by variance across samples (axis=1)."""
    gene_vars = expr_df.var(axis=1)
    top_genes = gene_vars.nlargest(n).index
    return expr_df.loc[top_genes]

# Log1p + z-score normalization
def log_zscore_normalize(expr_df: pd.DataFrame) -> pd.DataFrame:
    log_df = np.log1p(expr_df)
    return log_df.subtract(log_df.mean(axis=1), axis=0).divide(
        log_df.std(axis=1) + 1e-8, axis=0
    )

expr_norm = log_zscore_normalize(expr_df)
expr_top  = filter_high_variance_genes(expr_norm, n=1500)
```

**Pattern 3: Clustermap with subtype color bars**
```python
import seaborn as sns
import matplotlib.pyplot as plt

palette = {'MITF-low': '#4C72B0', 'Keratin': '#DD8452', 'Immune': '#55A868'}
col_colors = pd.Series(subtype_labels, index=sample_names).map(palette)

top50 = expr_norm.var(axis=1).nlargest(50).index
cg = sns.clustermap(
    expr_norm.loc[top50],
    col_colors=col_colors,
    cmap='RdBu_r', center=0, vmin=-3, vmax=3,
    yticklabels=True, xticklabels=False,
    figsize=(12, 8)
)
cg.fig.suptitle('Top 50 high-variance genes', y=1.02)
plt.show()
```

**Pattern 4: PCA and t-SNE embedding**
```python
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

X = expr_top.values.T  # shape: (n_samples, n_genes)

pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X)
print(f'PC1+PC2 variance: {pca.explained_variance_ratio_.sum():.1%}')

tsne = TSNE(n_components=2, perplexity=30, random_state=42)
X_tsne = tsne.fit_transform(X)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for label, color in palette.items():
    mask = np.array(subtype_labels) == label
    axes[0].scatter(X_pca[mask, 0], X_pca[mask, 1], c=color, label=label, alpha=0.7, s=30)
    axes[1].scatter(X_tsne[mask, 0], X_tsne[mask, 1], c=color, label=label, alpha=0.7, s=30)
for ax, title in zip(axes, ['PCA', 't-SNE']):
    ax.set_title(title); ax.legend()
plt.tight_layout(); plt.show()
```

**Pattern 5: Kaplan–Meier survival curves**
```python
try:
    from lifelines import KaplanMeierFitter
    fig, ax = plt.subplots(figsize=(9, 5))
    for subtype, color in palette.items():
        mask = np.array(subtype_labels) == subtype
        kmf = KaplanMeierFitter()
        kmf.fit(survival_times[mask], events[mask], label=subtype)
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color)
    ax.set_xlabel('Time (months)'); ax.set_ylabel('Survival probability')
    ax.set_ylim(0, 1.05); plt.tight_layout(); plt.show()
except ImportError:
    # Manual KM estimator
    def km_estimate(times, evts):
        order = np.argsort(times)
        t_s, e_s = times[order], evts[order]
        unique_t = np.unique(t_s[e_s])
        S, t_out = [1.0], [0.0]
        for t in unique_t:
            d = np.sum((t_s == t) & e_s)
            at_risk = np.sum(t_s >= t)
            S.append(S[-1] * (1 - d / at_risk))
            t_out.append(t)
        return np.array(t_out), np.array(S)
```

**Pattern 6: Comparing two classification schemes**
```python
from sklearn.metrics import normalized_mutual_info_score
import pandas as pd

# Cross-tabulation
crosstab = pd.crosstab(
    pd.Series(tirosh_labels, name='Tirosh'),
    pd.Series(harbst_labels, name='Harbst')
)
print(crosstab)

# Normalized Mutual Information
nmi = normalized_mutual_info_score(tirosh_labels, harbst_labels)
print(f'NMI: {nmi:.4f}')  # 0 = independent, 1 = perfect agreement

# Heatmap of row-normalized crosstab
import seaborn as sns
crosstab_norm = crosstab.div(crosstab.sum(axis=1), axis=0)
sns.heatmap(crosstab_norm, annot=True, fmt='.2f', cmap='Blues')
plt.title(f'Subtype concordance (NMI = {nmi:.3f})')
plt.show()
```

## Common Pitfalls

- **Axis confusion in z-score**: use `axis=1` to normalize each *gene* across samples; `axis=0` would normalize each *sample* across genes
- **Duplicate gene symbols**: cBioPortal files often have multiple rows per gene — deduplicate before analysis
- **Ward linkage requires Euclidean distance**: `AgglomerativeClustering(linkage='ward')` only works with `affinity='euclidean'` (default); other metrics will raise an error
- **Tn5 offset vs. RNA-seq**: log1p is appropriate for read counts; do not apply it twice
- **Survival data censoring**: KM estimator ignores censored observations when computing S(t) but counts them toward the at-risk set — never exclude censored samples entirely
- **NMI is symmetric**: `nmi(A, B) == nmi(B, A)`, so label ordering does not matter for concordance
- **Marker-based centroid labeling**: after AgglomerativeClustering, cluster IDs are arbitrary integers — always inspect centroids to map cluster numbers to biological subtypes
- **t-SNE is stochastic**: always set `random_state` for reproducibility; results may differ across runs with different perplexity values
