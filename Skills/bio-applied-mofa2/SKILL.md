---
name: bio-applied-mofa2
description: "**Tier 3 — Applied Bioinformatics | Module 27 · Notebook 2**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/27_Multi_Omics_Integration/02_mofa2.ipynb"
---

# MOFA2: Multi-Omics Factor Analysis

*Source: Course notebook `Tier_3_Applied_Bioinformatics/27_Multi_Omics_Integration/02_mofa2.ipynb`*

# MOFA2: Multi-Omics Factor Analysis

**Tier 3 — Applied Bioinformatics | Module 27 · Notebook 2**

*Prerequisites: Notebook 1 (Data Harmonization), Module 07 (Machine Learning)*

---

**By the end of this notebook you will be able to:**
1. Explain MOFA2's probabilistic model and latent factor concept
2. Train a MOFA2 model on a multi-omics dataset
3. Interpret variance explained per factor per omics layer
4. Identify top features (weights) driving each latent factor
5. Correlate latent factors with sample metadata (age, treatment, subtype)



**Key resources:**
- [MOFA2 documentation and tutorials](https://biofam.github.io/MOFA2/)
- [MOFA2 paper (Argelaguet et al., 2020)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1)
- [mofapy2 Python package](https://github.com/bioFAM/mofapy2)

## 1. MOFA2 Model Overview

**MOFA2** (Multi-Omics Factor Analysis v2) is a probabilistic framework that learns a low-dimensional representation of multi-omics data.

### The generative model

```
                          Z (shared factors)
                          |    |    |
                         W_1  W_2  W_3
                          |    |    |
                         X_1  X_2  X_3
                    (RNA)  (Prot) (Meth)
```

- **Z** (N × K): Latent factor matrix (N samples, K factors) — shared across views
- **W_m** (D_m × K): Feature weight matrix for view m — view-specific
- **X_m = Z @ W_m.T + noise** for each omics view m

### Key advantages over single-omics PCA
- Identifies factors that are shared across multiple omics layers
- Identifies factors unique to a single layer (view-specific)
- Handles **missing views** natively (samples lacking one omics layer)
- Produces interpretable, sparse weights for feature importance

### MOFA2 installation and usage
```python
# Install
pip install mofapy2

# Python API
from mofapy2.run.entry_point import entry_point
ent = entry_point()
ent.set_data_options(scale_groups=False, scale_views=True)
ent.set_data_df(data_df)  # long-format dataframe
ent.set_model_options(factors=10)
ent.set_train_options(iter=1000, convergence_mode="fast")
ent.build()
ent.run()
ent.save("mofa_model.hdf5")
```

Or via muon (recommended):
```python
import muon as mu
mdata = mu.MuData({'rna': rna_adata, 'prot': prot_adata, 'meth': meth_adata})
mu.tl.mofa(mdata, n_factors=10, outfile="model.hdf5")
```

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

np.random.seed(42)

# Generate multi-omics data with known latent structure
n_samples = 80
n_rna = 500
n_prot = 150
n_meth = 300

sample_ids = [f'S{i:03d}' for i in range(n_samples)]
sample_type = (['TypeA'] * 20 + ['TypeB'] * 20 + ['TypeC'] * 20 + ['TypeD'] * 20)
group = ['Group1'] * 40 + ['Group2'] * 40

# ---- True latent factors ----
# Factor 1: discriminates TypeA/B from C/D (biological signal)
# Factor 2: continuous gradient (tumor purity, patient age)
# Factor 3: batch/technical factor
# Factor 4: TypeA vs B within group
n_factors_true = 5
Z_true = np.zeros((n_samples, n_factors_true))
# Factor 1: type-driven
type_map = {'TypeA': -1.5, 'TypeB': -0.5, 'TypeC': 0.5, 'TypeD': 1.5}
Z_true[:, 0] = [type_map[t] + np.random.normal(0, 0.2) for t in sample_type]
# Factor 2: continuous gradient
Z_true[:, 1] = np.linspace(-2, 2, n_samples) + np.random.normal(0, 0.3, n_samples)
# Factor 3: group effect
Z_true[:, 2] = [1.5 if g == 'Group1' else -1.5 for g in group]
Z_true[:, 2] += np.random.normal(0, 0.2, n_samples)
# Factor 4: noise
Z_true[:, 3] = np.random.normal(0, 1, n_samples)
# Factor 5: noise
Z_true[:, 4] = np.random.normal(0, 1, n_samples)

# ---- Weights per omics layer ----
# Factor 1 is captured by all 3 layers
# Factor 2 only by RNA and methylation
# Factor 3 only by proteomics
W_rna = np.zeros((n_rna, n_factors_true))
W_rna[:100, 0] = np.random.randn(100) * 0.5   # Factor 1 captured by first 100 RNA features
W_rna[100:200, 1] = np.random.randn(100) * 0.4  # Factor 2
W_rna[200:250, 3] = np.random.randn(50) * 0.3

W_prot = np.zeros((n_prot, n_factors_true))
W_prot[:50, 0] = np.random.randn(50) * 0.5
W_prot[50:80, 2] = np.random.randn(30) * 0.6   # Factor 3 (proteomics only)

W_meth = np.zeros((n_meth, n_factors_true))
W_meth[:80, 0] = np.random.randn(80) * 0.4
W_meth[80:140, 1] = np.random.randn(60) * 0.4

# Generate data: X = Z @ W.T + noise
rna_data = Z_true @ W_rna.T + np.random.normal(0, 0.5, (n_samples, n_rna))
prot_data = Z_true @ W_prot.T + np.random.normal(0, 0.5, (n_samples, n_prot))
meth_data = Z_true @ W_meth.T + np.random.normal(0, 0.5, (n_samples, n_meth))

# Standardize each view
rna_scaled = StandardScaler().fit_transform(rna_data)
prot_scaled = StandardScaler().fit_transform(prot_data)
meth_scaled = StandardScaler().fit_transform(meth_data)

print("Multi-omics dataset (MOFA2 input format):")
print(f"  RNA-seq:     {rna_scaled.shape} (samples x features)")
print(f"  Proteomics:  {prot_scaled.shape}")
print(f"  Methylation: {meth_scaled.shape}")
print(f"\nTrue latent factors: {n_factors_true}")
print(f"Sample types: {dict(pd.Series(sample_type).value_counts())}")
```

## 2. Training the MOFA2 Model

### Data input format
MOFA2 expects a **long-format DataFrame** with columns:
- `sample`: sample ID
- `group`: group label (optional, for multi-group MOFA)
- `feature`: feature name
- `view`: omics layer name
- `value`: numeric value (pre-normalized)

### Model parameters

| Parameter | Default | Notes |
|---|---|---|
| `n_factors` | 10 | Number of factors to infer |
| `convergence_mode` | "fast" | fast/medium/slow — controls ELBO tolerance |
| `spikeslab_factors` | False | Sparse factors (binary inclusion) |
| `spikeslab_weights` | True | Sparse weights (ARD prior) |
| `scale_views` | True | Normalize each view to unit variance |

### Convergence
MOFA2 optimizes the **Evidence Lower Bound (ELBO)** via variational Bayes. Monitor convergence:
```python
model.plot_ELBO()  # ELBO should plateau
```
Typically converges in 100-500 iterations.

```python
# ----- MOFA2-style factor analysis (implemented with sklearn/numpy) -----
# Note: in practice you would use mofapy2 or muon.tl.mofa
# Here we implement a similar matrix factorization to demonstrate the concept

from sklearn.decomposition import TruncatedSVD
from scipy.stats import pearsonr

def mofa_like(views, n_factors=5, max_iter=50):
    """
    Simple multi-view matrix factorization mimicking MOFA2.
    Views: list of (samples x features) matrices.
    Returns: Z (samples x factors), W_list (features x factors per view)
    """
    n_samples = views[0].shape[0]
    
    # Initialize Z by PCA on concatenated data
    concat = np.hstack(views)
    svd = TruncatedSVD(n_components=n_factors)
    Z = svd.fit_transform(concat)
    Z /= np.std(Z, axis=0, keepdims=True)
    
    W_list = []
    var_explained = []
    
    for v_idx, X in enumerate(views):
        # Regress X on Z to get weights
        W = np.linalg.lstsq(Z, X, rcond=None)[0].T  # features x factors
        W_list.append(W)
        
        # Variance explained per factor per view
        X_recon = Z @ W.T
        ss_total = np.sum((X - X.mean(axis=0))**2)
        r2_per_factor = []
        for k in range(n_factors):
            X_k = Z[:, k:k+1] @ W[:, k:k+1].T
            r2 = 1 - np.sum((X - X_k)**2) / ss_total
            r2_per_factor.append(max(0, r2))
        var_explained.append(r2_per_factor)
    
    return Z, W_list, np.array(var_explained)

views = [rna_scaled, prot_scaled, meth_scaled]
view_names = ['RNA-seq', 'Proteomics', 'Methylation']
n_factors = 5

Z_est, W_est, var_explained = mofa_like(views, n_factors=n_factors)

print(f"MOFA2-style model trained with {n_factors} factors")
print(f"\nVariance explained per factor per view (R2):")
df_var = pd.DataFrame(var_explained * 100,
                       index=view_names,
                       columns=[f'Factor{k+1}' for k in range(n_factors)])
print(df_var.round(2).to_string())
```

## 3. Variance Decomposition

The central output of MOFA2 is the **R² (variance explained)** per factor per view.

### Interpretation guide

| Pattern | Meaning |
|---|---|
| High R² in all views | Shared biological signal (master regulator, global effect) |
| High R² in RNA only | Transcription-specific factor (e.g., splicing) |
| High R² in methylation only | Epigenetic reprogramming |
| High R² in proteomics only | Post-translational regulatory mechanism |

### Determining optimal number of factors
- **Elbow** in cumulative variance explained curve
- **Robustness**: factors that appear consistently across runs with different random seeds
- Rule of thumb: start with 15-20, prune inactive factors (var explained < 2%)

```python
# ----- Variance explained visualization -----

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('MOFA2 Model Results: Variance Decomposition', fontsize=13, fontweight='bold')

# Panel 1: Heatmap of variance explained per factor per view
sns.heatmap(df_var, ax=axes[0], cmap='YlOrRd', annot=True, fmt='.1f',
            vmin=0, cbar_kws={'label': 'Variance explained (%)'})
axes[0].set_title('R² per Factor per View (%)')
axes[0].set_xlabel('Factor')

# Panel 2: Total variance per factor (stacked bar)
factor_labels = [f'Factor{k+1}' for k in range(n_factors)]
bottom_arr = np.zeros(n_factors)
colors = ['steelblue', 'orange', 'green']
for i, (vname, color) in enumerate(zip(view_names, colors)):
    axes[1].bar(factor_labels, var_explained[i] * 100, bottom=bottom_arr,
                color=color, alpha=0.8, label=vname)
    bottom_arr += var_explained[i] * 100
axes[1].set_ylabel('Total variance explained (%)')
axes[1].set_title('Variance per Factor (across all views)')
axes[1].legend()

# Panel 3: Per-view total variance
total_per_view = var_explained.sum(axis=1) * 100
axes[2].bar(view_names, total_per_view, color=colors, alpha=0.8)
axes[2].set_ylabel('Total R² (%)')
axes[2].set_title('Total Variance Explained per View')
for bar, val in zip(axes[2].patches, total_per_view):
    axes[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                 f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')

plt.tight_layout()
plt.show()

print("\nInterpretation:")
print("- Factors with high R² in multiple views = shared biological signal")
print("- Factors with high R² in only one view = view-specific signal")
print("- Factor 1 captures most shared variance (present in all views)")
```

## 4. Factor Interpretation

Each factor is characterized by:
1. **Sample scores** (Z values): how strongly each sample expresses the factor
2. **Feature weights** (W values): which features drive the factor

### Biological interpretation workflow
1. **Correlate factor scores** with metadata (sample type, clinical variables)
2. **Top weighted features** → candidate genes/proteins/CpGs driving the factor
3. **Gene set enrichment** on top RNA weights → pathway interpretation
4. **TF motif enrichment** on top ATAC weights → regulatory interpretation

### MOFA2 vs PCA
| Aspect | PCA | MOFA2 |
|---|---|---|
| Input | Single matrix | Multiple matrices |
| Shared signal | Implicit | Explicit |
| View-specific signal | Mixed | Separated |
| Missing views | Cannot handle | Handles natively |
| Sparse weights | No | Yes (ARD prior) |

```python
# ----- Factor interpretation -----
import matplotlib.patches as mpatches

type_colors = {'TypeA': '#E74C3C', 'TypeB': '#3498DB', 'TypeC': '#2ECC71', 'TypeD': '#F39C12'}
group_colors = {'Group1': '#8E44AD', 'Group2': '#1ABC9C'}

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('MOFA2 Factor Interpretation', fontsize=13, fontweight='bold')

# Row 1: Factor scatter plots
factor_pairs = [(0, 1), (0, 2), (1, 2)]
for col, (f1, f2) in enumerate(factor_pairs):
    for s in range(n_samples):
        axes[0, col].scatter(Z_est[s, f1], Z_est[s, f2],
                              c=type_colors[sample_type[s]], s=25, alpha=0.8)
    axes[0, col].set_xlabel(f'Factor {f1+1}')
    axes[0, col].set_ylabel(f'Factor {f2+1}')
    axes[0, col].set_title(f'Factor {f1+1} vs Factor {f2+1}')

type_handles = [mpatches.Patch(color=v, label=k) for k, v in type_colors.items()]
axes[0, -1].legend(handles=type_handles, loc='lower right', fontsize=8)

# Row 2: Feature weights for Factor 1 (top features per view)
for col, (vname, W) in enumerate(zip(view_names, W_est)):
    factor_weights = W[:, 0]  # Factor 1 weights
    top_idx = np.argsort(np.abs(factor_weights))[-20:]
    weights_top = factor_weights[top_idx]
    
    pos_colors = ['red' if w > 0 else 'blue' for w in weights_top]
    axes[1, col].barh(range(20), weights_top, color=pos_colors, alpha=0.7)
    axes[1, col].axvline(0, color='black', linewidth=0.8)
    axes[1, col].set_yticks(range(20))
    axes[1, col].set_yticklabels([f'Feature_{idx}' for idx in top_idx], fontsize=6)
    axes[1, col].set_xlabel('Weight')
    axes[1, col].set_title(f'Factor 1 Top Weights\n({vname})')

plt.tight_layout()
plt.show()

# Correlation of estimated factors with true factors
print("\nCorrelation of estimated factors with true latent factors:")
header = ''.join([f'  TrueF{k+1}' for k in range(5)])
print(f'{"":10s}' + header)
for i in range(n_factors):
    row = f'EstF{i+1}:    '
    for j in range(n_factors_true):
        r, _ = pearsonr(Z_est[:, i], Z_true[:, j])
        row += f'   {abs(r):.2f} '
    print(row)
```
