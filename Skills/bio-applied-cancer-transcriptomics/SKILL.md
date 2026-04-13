---
name: bio-applied-cancer-transcriptomics
description: "*Prerequisites: Modules 01–23. Module 03 (RNA-seq), Module 07 (Machine Learning), Module 06 (Statistics) strongly recommended.*"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/24_Cancer_Transcriptomics/01_cancer_transcriptomics.ipynb"
---

# Module 24: Cancer Transcriptomics: Subtype Classification

*Source: Course notebook `Tier_3_Applied_Bioinformatics/24_Cancer_Transcriptomics/01_cancer_transcriptomics.ipynb`*

# Module 24: Cancer Transcriptomics: Subtype Classification

**Tier 3 — Applied Bioinformatics | Module 24**

*Prerequisites: Modules 01–23. Module 03 (RNA-seq), Module 07 (Machine Learning), Module 06 (Statistics) strongly recommended.*

---

Cancer transcriptomics uses bulk or single-cell RNA-seq to stratify tumors into molecular subtypes with distinct prognoses and therapeutic vulnerabilities. This module focuses on melanoma (TCGA-SKCM) using synthetic data mirroring the TCGA portal, covering variance filtering, hierarchical clustering, semi-supervised Random Forest classification, Kaplan–Meier survival analysis, and subtype comparison.

**Key references:**
- Tirosh et al. (2016) *Science* — single-cell dissection of intratumoral heterogeneity in melanoma
- Harbst et al. (2016) *Clin Cancer Res* — four molecular subtypes of metastatic melanoma from bulk RNA-seq (TCGA-SKCM)
- TCGA-SKCM data available via [cBioPortal](https://www.cbioportal.org/study/summary?id=skcm_tcga)

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import classification_report, normalized_mutual_info_score

%matplotlib inline
plt.rcParams['figure.dpi'] = 100
sns.set_style('whitegrid')
```

## 1. Introduction: Cancer Transcriptomics & Melanoma Subtypes

### What is Cancer Transcriptomics?

Cancer transcriptomics applies bulk or single-cell RNA sequencing to tumors to:
- **Classify** tumors into biologically meaningful subtypes
- **Identify** subtype-specific vulnerabilities and therapy targets
- **Predict** patient outcomes (prognosis) from gene-expression profiles

### TCGA-SKCM: Melanoma as a Model

The Cancer Genome Atlas (TCGA) Skin Cutaneous Melanoma (SKCM) cohort contains RNA-seq profiles for ~470 metastatic tumors, accessible via [cBioPortal](https://www.cbioportal.org/study/summary?id=skcm_tcga). Two landmark studies defined the major transcriptional subtypes:

| Study | Subtypes | Method |
|---|---|---|
| Tirosh et al. 2016 *Science* | MITF-low / Keratin / Immune | Semi-supervised, scRNA-seq |
| Harbst et al. 2016 *Clin Cancer Res* | Pigmentation / Proliferative / Normal-like / High-immune | Hierarchical clustering, bulk RNA-seq |

### Marker Genes Used in This Module

| Subtype | Key Markers |
|---|---|
| Pigmentation (MITF-high) | MITF, DCT, TYRP1, MLANA |
| Keratin | KRT5, KRT14, KRT6A, EGFR |
| Immune | CD3D, CD8A, GZMB, PRF1 |
| Proliferative | MKI67, TOP2A, CDK1 |
| Normal-like | VIM, CDH2, FN1 |

## 2. Data Loading & Preprocessing

Real cBioPortal data arrives as tab-delimited files with genes as rows and samples as columns. We generate a synthetic matrix with the same structure and realistic subtype-specific expression patterns.

**Preprocessing steps:**
1. `log1p` transform to stabilize variance (handles zero counts)
2. Z-score normalization per gene (mean=0, std=1 across samples)

```python
np.random.seed(42)
N_SAMPLES = 200
N_GENES = 5000

# Ground-truth subtype labels (Tirosh 3-class)
subtype_labels = np.array(['MITF-low'] * 70 + ['Keratin'] * 65 + ['Immune'] * 65)
np.random.shuffle(subtype_labels)

# Base expression matrix (lognormal)
expr = np.random.lognormal(mean=2, sigma=1.5, size=(N_GENES, N_SAMPLES))

gene_names = [f'GENE_{i:04d}' for i in range(N_GENES)]
marker_genes = {
    'MITF': 0, 'DCT': 1, 'TYRP1': 2, 'MLANA': 3,
    'KRT5': 4, 'KRT14': 5, 'KRT6A': 6, 'EGFR': 7,
    'CD3D': 8, 'CD8A': 9, 'GZMB': 10, 'PRF1': 11,
    'MKI67': 12, 'TOP2A': 13, 'CDK1': 14,
    'VIM': 15, 'CDH2': 16, 'FN1': 17,
}
for gene, idx in marker_genes.items():
    gene_names[idx] = gene

# Amplify marker gene signals per subtype
mitf_mask = subtype_labels == 'MITF-low'
krt_mask  = subtype_labels == 'Keratin'
imm_mask  = subtype_labels == 'Immune'

for gene in ['MITF', 'DCT', 'TYRP1', 'MLANA']:
    idx = gene_names.index(gene)
    expr[idx, mitf_mask] *= 4.0
for gene in ['KRT5', 'KRT14', 'KRT6A', 'EGFR']:
    idx = gene_names.index(gene)
    expr[idx, krt_mask] *= 4.0
for gene in ['CD3D', 'CD8A', 'GZMB', 'PRF1']:
    idx = gene_names.index(gene)
    expr[idx, imm_mask] *= 4.0

sample_names = [f'SAMPLE_{i:03d}' for i in range(N_SAMPLES)]
expr_df = pd.DataFrame(expr, index=gene_names, columns=sample_names)
print(f'Expression matrix: {expr_df.shape[0]} genes × {expr_df.shape[1]} samples')
expr_df.iloc[:5, :5]
```

```python
# log1p transform
expr_log = np.log1p(expr_df)

# Z-score normalize per gene (axis=1: across samples)
expr_z = expr_log.subtract(expr_log.mean(axis=1), axis=0).divide(
    expr_log.std(axis=1) + 1e-8, axis=0
)

print('After log1p+zscore:')
print(f'  Mean per gene (first 3): {expr_z.iloc[:3].mean(axis=1).values.round(3)}')
print(f'  Std  per gene (first 3): {expr_z.iloc[:3].std(axis=1).values.round(3)}')
```

## 3. Exploratory Analysis

Before classification, we explore global structure:
- **Variance filtering** — retain the 1,500 most variable genes to reduce noise
- **Clustermap** — hierarchically clustered heatmap of the top-50 genes
- **PCA** — linear projection to visualize dominant variance axes
- **t-SNE** — nonlinear embedding preserving local neighborhood structure

```python
# ── Variance filtering ──────────────────────────────
gene_vars = expr_z.var(axis=1)
top_genes = gene_vars.nlargest(1500).index
expr_top  = expr_z.loc[top_genes]
print(f'Retained {len(top_genes)} high-variance genes')

# ── Clustermap (top 50 genes for readability) ────────
top50 = gene_vars.nlargest(50).index
subtype_palette = {'MITF-low': '#4C72B0', 'Keratin': '#DD8452', 'Immune': '#55A868'}
col_colors = pd.Series(subtype_labels, index=sample_names).map(subtype_palette)

cg = sns.clustermap(
    expr_z.loc[top50],
    col_colors=col_colors,
    cmap='RdBu_r', center=0, vmin=-3, vmax=3,
    yticklabels=True, xticklabels=False,
    figsize=(12, 8)
)
cg.fig.suptitle('Clustermap: Top 50 high-variance genes (TCGA-SKCM synthetic)', y=1.02)
plt.show()
```

```python
# ── PCA ─────────────────────────────────────────────
X = expr_top.values.T  # samples × genes
pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (label, color) in zip(
    [axes[0]] * 3,
    [(s, subtype_palette[s]) for s in ['MITF-low', 'Keratin', 'Immune']]
):
    mask = subtype_labels == label
    axes[0].scatter(X_pca[mask, 0], X_pca[mask, 1], c=color, label=label, alpha=0.7, s=30)
axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
axes[0].set_title('PCA — Tirosh subtypes')
axes[0].legend()

# ── t-SNE ────────────────────────────────────────────
tsne = TSNE(n_components=2, random_state=42, perplexity=30)
X_tsne = tsne.fit_transform(X)

for label, color in [(s, subtype_palette[s]) for s in ['MITF-low', 'Keratin', 'Immune']]:
    mask = subtype_labels == label
    axes[1].scatter(X_tsne[mask, 0], X_tsne[mask, 1], c=color, label=label, alpha=0.7, s=30)
axes[1].set_xlabel('t-SNE 1'); axes[1].set_ylabel('t-SNE 2')
axes[1].set_title('t-SNE — Tirosh subtypes')
axes[1].legend()

plt.tight_layout()
plt.show()
print(f'PCA cumulative variance (2 PCs): {pca.explained_variance_ratio_.sum():.1%}')
```

## 4. Semi-supervised Classification (Tirosh 3-class)

In the Tirosh framework, a **labeled training set** of cells with known phenotypes (MITF-high pigmented, Keratin, Immune-enriched) is used to train a classifier that is then applied to unlabeled samples.

We simulate this by holding out 30% of samples as 'unlabeled' and predicting their subtypes using a **Random Forest** trained on the labeled 70%.

```python
from sklearn.model_selection import train_test_split

X_clf = expr_top.values.T  # 200 samples × 1500 genes
y = subtype_labels

X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
    X_clf, y, np.arange(N_SAMPLES), test_size=0.30, random_state=42, stratify=y
)

rf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
rf.fit(X_train, y_train)
y_pred = rf.predict(X_test)

print('=== Tirosh 3-class Random Forest ===')
print(classification_report(y_test, y_pred))

# Store predictions for all samples (train uses ground truth, test uses RF)
tirosh_labels = y.copy()
tirosh_labels[idx_test] = y_pred  # replace test labels with predictions
```

## 5. Hierarchical Clustering Subtypes (Harbst 4-class)

Harbst et al. identified four subtypes using Ward-linkage hierarchical clustering on a signature of melanoma-relevant marker genes:

| Subtype | Primary Signature |
|---|---|
| Pigmentation | MITF, DCT, TYRP1, MLANA |
| Proliferative | MKI67, TOP2A, CDK1 |
| Normal-like | VIM, CDH2, FN1 |
| High-immune | CD3D, CD8A, GZMB, PRF1 |

We extract expression of these marker genes, then apply `AgglomerativeClustering` with Ward linkage.

```python
harbst_markers = ['MITF', 'DCT', 'TYRP1', 'MLANA',
                  'MKI67', 'TOP2A', 'CDK1',
                  'VIM', 'CDH2', 'FN1',
                  'CD3D', 'CD8A', 'GZMB', 'PRF1']

# Extract marker gene expression (z-scored)
X_harbst = expr_z.loc[harbst_markers].values.T  # 200 × 14

# Ward-linkage hierarchical clustering into 4 clusters
hc = AgglomerativeClustering(n_clusters=4, linkage='ward')
cluster_ids = hc.fit_predict(X_harbst)

# Assign biologically meaningful labels by inspecting cluster centroids
cluster_df = pd.DataFrame(X_harbst, columns=harbst_markers)
cluster_df['cluster'] = cluster_ids
centroids = cluster_df.groupby('cluster').mean()

subtype_map = {
    centroids[['MITF', 'DCT', 'TYRP1', 'MLANA']].mean(axis=1).idxmax(): 'Pigmentation',
    centroids[['MKI67', 'TOP2A', 'CDK1']].mean(axis=1).idxmax(): 'Proliferative',
    centroids[['CD3D', 'CD8A', 'GZMB', 'PRF1']].mean(axis=1).idxmax(): 'High-immune',
}
remaining = [c for c in range(4) if c not in subtype_map]
if remaining:
    subtype_map[remaining[0]] = 'Normal-like'

harbst_labels = np.array([subtype_map.get(c, f'Cluster_{c}') for c in cluster_ids])

# Heatmap of marker gene expression by Harbst subtype
plot_df = pd.DataFrame(X_harbst, columns=harbst_markers)
plot_df['Subtype'] = harbst_labels
plot_df_sorted = plot_df.sort_values('Subtype')

harbst_palette = {'Pigmentation': '#E07B54', 'Proliferative': '#9B59B6',
                  'Normal-like': '#2ECC71', 'High-immune': '#3498DB'}
row_colors = plot_df_sorted['Subtype'].map(harbst_palette)

g = sns.clustermap(
    plot_df_sorted.drop(columns='Subtype').T,
    col_colors=row_colors.values,
    cmap='RdBu_r', center=0, vmin=-3, vmax=3,
    row_cluster=False, col_cluster=False,
    yticklabels=True, xticklabels=False, figsize=(12, 6)
)
g.fig.suptitle('Harbst 4-class subtypes — marker gene expression', y=1.02)
plt.show()

print('Harbst subtype counts:')
print(pd.Series(harbst_labels).value_counts())
```

## 6. Survival Analysis: Kaplan–Meier Curves

We assess whether the Tirosh subtypes differ in survival using **Kaplan–Meier estimator** — a non-parametric method that accounts for censored observations (patients lost to follow-up).

**Lifelines** is the standard Python library for survival analysis. A manual fallback is provided.

```python
# Synthetic survival data correlated with Tirosh subtypes
np.random.seed(42)
survival_times = np.random.exponential(scale=24, size=N_SAMPLES)  # months
survival_times[subtype_labels == 'Immune'] *= 1.8  # Immune subtype: better prognosis
events = np.random.binomial(1, 0.6, size=N_SAMPLES).astype(bool)

fig, ax = plt.subplots(figsize=(9, 5))
colors = {'MITF-low': '#4C72B0', 'Keratin': '#DD8452', 'Immune': '#55A868'}

try:
    from lifelines import KaplanMeierFitter
    for subtype, color in colors.items():
        mask = subtype_labels == subtype
        kmf = KaplanMeierFitter()
        kmf.fit(survival_times[mask], events[mask], label=subtype)
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color)
    ax.set_title('Kaplan–Meier survival curves (lifelines)')
except ImportError:
    # Manual KM implementation
    def km_estimate(times, evts):
        order = np.argsort(times)
        t_sorted, e_sorted = times[order], evts[order]
        unique_times = np.unique(t_sorted[e_sorted])
        n = len(times)
        S, t_out = [1.0], [0.0]
        for t in unique_times:
            d = np.sum((t_sorted == t) & e_sorted)
            at_risk = np.sum(t_sorted >= t)
            S.append(S[-1] * (1 - d / at_risk))
            t_out.append(t)
        return np.array(t_out), np.array(S)

    for subtype, color in colors.items():
        mask = subtype_labels == subtype
        t_km, s_km = km_estimate(survival_times[mask], events[mask])
        ax.step(t_km, s_km, where='post', label=subtype, color=color, lw=2)
    ax.set_title('Kaplan–Meier survival curves (manual)')

ax.set_xlabel('Time (months)')
ax.set_ylabel('Survival probability')
ax.legend()
ax.set_ylim(0, 1.05)
plt.tight_layout()
plt.show()
```
