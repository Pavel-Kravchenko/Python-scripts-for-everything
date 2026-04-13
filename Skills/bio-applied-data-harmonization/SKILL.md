---
name: bio-applied-data-harmonization
description: "**Tier 3 — Applied Bioinformatics | Module 27 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/27_Multi_Omics_Integration/01_data_harmonization.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scikit-learn 1.4+, scipy 1.12+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Data Harmonization for Multi-Omics

*Source: Course notebook `Tier_3_Applied_Bioinformatics/27_Multi_Omics_Integration/01_data_harmonization.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 27 · Notebook 1**

*Prerequisites: Module 03 (RNA-seq), Module 07 (Machine Learning), Module 18 (Proteomics)*

---

**By the end of this notebook you will be able to:**
1. Describe the challenges of integrating heterogeneous omics data
2. Apply sample-centric and feature-centric scaling strategies
3. Handle missing data across omics layers (imputation strategies)
4. Detect and visualize batch effects between omics platforms
5. Prepare a harmonized multi-omics matrix for downstream integration



**Key resources:**
- [MOFA2 documentation](https://biofam.github.io/MOFA2/)
- [mixOmics documentation](http://mixomics.org/)
- [Galaxy Training — Multi-Omics](https://training.galaxyproject.org/)

## 1. The Multi-Omics Data Landscape

Multi-omics integration combines measurements from different molecular layers to gain a holistic view of biological systems.

### Data modalities and their properties

| Omics layer | Measures | Typical features | Data type |
|---|---|---|---|
| **Genomics** (WGS/WES) | DNA variants, CNVs | ~4M SNPs | Binary/integer |
| **Transcriptomics** (RNA-seq) | mRNA abundance | ~20,000 genes | Continuous (counts) |
| **Proteomics** (LC-MS) | Protein abundance | 3,000-10,000 proteins | Continuous, many missing |
| **Metabolomics** | Metabolite levels | 500-5,000 metabolites | Continuous |
| **Epigenomics** (ATAC/ChIP) | Chromatin accessibility | ~100k-1M peaks | Binary/continuous |
| **Methylomics** (RRBS/WGBS) | DNA methylation | ~500k-28M CpGs | Beta (0-1) |

### Integration strategies

**Concatenation (early integration)**: merge all matrices horizontally → dimension reduction
**Meta-analysis (late integration)**: analyze each omics separately → combine results
**Factor analysis (intermediate)**: MOFA2, iCluster — learn shared latent factors
**Supervised (DIABLO)**: use sample labels to guide integration

### Key challenges
- Vastly different feature dimensionalities (p >> n for each layer)
- Missing data in proteomics (~20-40% MNAR)
- Batch effects from different processing dates/platforms
- Scale differences: log-counts vs beta values vs raw intensities

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.preprocessing import StandardScaler, MinMaxScaler, QuantileTransformer
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer

np.random.seed(42)

# ----- Generate synthetic multi-omics dataset -----
# 60 samples: 20 normal, 20 tumor type A, 20 tumor type B
n_samples = 60
n_rna = 2000   # RNA-seq features (genes)
n_prot = 300   # proteomics features
n_meth = 500   # methylation features (CpG sites)

sample_ids = [f'S{i:03d}' for i in range(n_samples)]
sample_type = ['Normal']*20 + ['TumorA']*20 + ['TumorB']*20

# Two batches (different sequencing dates)
batch = ['Batch1']*30 + ['Batch2']*30

# Simulate RNA-seq (log-normalized counts)
rna_matrix = np.random.randn(n_samples, n_rna)
# Add true biological signal
for i in range(20, 40):   # TumorA signal
    rna_matrix[i, :200] += 2.0
for i in range(40, 60):   # TumorB signal
    rna_matrix[i, 200:400] += 2.0
# Add batch effect
for i in range(30):  # Batch 2
    rna_matrix[i] += np.random.normal(2.0, 0.5, n_rna)

# Simulate proteomics (MS intensities, log2-transformed)
prot_matrix = np.random.randn(n_samples, n_prot)
for i in range(20, 40):
    prot_matrix[i, :50] += 1.5
for i in range(40, 60):
    prot_matrix[i, 50:100] += 1.5
prot_matrix[:30] += np.random.normal(1.5, 0.3)  # batch effect

# Simulate methylation (beta values 0-1)
meth_matrix = np.random.beta(2, 5, (n_samples, n_meth))
for i in range(20, 40):
    meth_matrix[i, :80] = np.clip(meth_matrix[i, :80] + 0.3, 0, 1)
for i in range(40, 60):
    meth_matrix[i, 80:160] = np.clip(meth_matrix[i, 80:160] + 0.3, 0, 1)

# Metadata
metadata = pd.DataFrame({
    'SampleID': sample_ids,
    'Type': sample_type,
    'Batch': batch,
    'Age': np.random.randint(30, 75, n_samples),
    'Sex': np.random.choice(['M', 'F'], n_samples)
})

print("Multi-omics dataset created:")
print(f"  RNA-seq: {rna_matrix.shape[0]} samples x {rna_matrix.shape[1]} features")
print(f"  Proteomics: {prot_matrix.shape[0]} samples x {prot_matrix.shape[1]} features")
print(f"  Methylation: {meth_matrix.shape[0]} samples x {meth_matrix.shape[1]} features")
print(f"\nMetadata:")
print(metadata.groupby(['Type', 'Batch']).size().to_string())
```python

## 2. Normalization and Scaling

Different omics layers require different normalization strategies before integration.

### Normalization by data type

| Data type | Recommended normalization | Rationale |
|---|---|---|
| RNA-seq counts | TMM/DESeq2 VST → z-score | Library size + biological variation |
| Proteomics intensities | Median centering + z-score | Inter-sample variation |
| Methylation beta | logit transform (M-value) → z-score | Beta values are bounded 0-1 |
| Metabolomics | Probabilistic quotient + log | Dilution effects |

### Z-score vs quantile normalization

**Z-score**: preserves relative differences between features; assumes approximately normal distribution.

**Quantile normalization**: forces identical distributions across samples; removes global sample-level effects but can distort signal.

**Key principle**: normalize *within* each omics layer before integration. Applying the same normalization across different omics would not be appropriate.

```python
# ----- Normalization and scaling strategies -----

def compare_normalizations(data, title, sample_labels, n_show=20):
    """Compare z-score, min-max, and quantile normalization."""
    scalers = {
        'Raw': data,
        'Z-score': StandardScaler().fit_transform(data),
        'Min-Max': MinMaxScaler().fit_transform(data),
        'Quantile\n(normal)': QuantileTransformer(output_distribution='normal',
                                                    random_state=42).fit_transform(data),
    }
    
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    fig.suptitle(f'Normalization Comparison: {title}', fontsize=12, fontweight='bold')
    
    type_colors = {'Normal': 'green', 'TumorA': 'orange', 'TumorB': 'red'}
    colors = [type_colors[t] for t in sample_type]
    
    for j, (norm_name, norm_data) in enumerate(scalers.items()):
        # Distribution of first 5 features
        axes[0, j].boxplot([norm_data[:, k] for k in range(5)],
                           patch_artist=True,
                           boxprops=dict(facecolor='steelblue', alpha=0.5))
        axes[0, j].set_title(f'{norm_name}\n(feature distributions)')
        axes[0, j].set_xlabel('Feature index')
        axes[0, j].set_xticks(range(1, 6))
        axes[0, j].set_xticklabels([f'F{k}' for k in range(5)])
        
        # PCA plot
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(norm_data)
        for s, (x, y) in enumerate(pcs):
            axes[1, j].scatter(x, y, c=type_colors[sample_type[s]], s=20, alpha=0.7)
        axes[1, j].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
        axes[1, j].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
        axes[1, j].set_title(f'{norm_name}\nPCA')
    
    handles = [mpatches.Patch(color=v, label=k) for k, v in type_colors.items()]
    axes[1, -1].legend(handles=handles, loc='lower right', fontsize=8)
    plt.tight_layout()
    plt.show()

compare_normalizations(rna_matrix, 'RNA-seq', sample_type)
print("Key observation: Quantile normalization equalizes feature distributions")
print("but can over-correct if batch effect is confounded with biology.")
```python

## 3. Missing Data in Multi-Omics

Missing data mechanisms affect imputation strategy:

### Missing data mechanisms
- **MCAR** (Missing Completely at Random): random technical failures → KNN or mean imputation
- **MAR** (Missing at Random): depends on observed data → KNN imputation
- **MNAR** (Missing Not at Random): missing because value is too low to detect → min/2 imputation

### Proteomics missing data
In LC-MS proteomics, missing values predominantly follow **MNAR** (low-abundance proteins below detection limit). Best practices:
- Filter proteins missing in >50% of samples
- For remaining: half-minimum imputation for MNAR; KNN for MAR
- Flag imputed values for downstream sensitivity analysis

### RNA-seq: no traditional missing values
Zeros in RNA-seq are true zeros (not detected), not missing data. Handle with count-based models (DESeq2, edgeR) that account for zero-inflation.

```python
# ----- Missing data imputation for proteomics -----

# Introduce missing values (MS data has lots of missing - ~30%)
prot_missing = prot_matrix.copy()
missing_mask = np.random.random(prot_missing.shape) < 0.25  # 25% missing
prot_missing[missing_mask] = np.nan

print(f"Missing values introduced: {missing_mask.sum()} / {missing_mask.size}")
print(f"Missing rate: {missing_mask.mean()*100:.1f}%")

# Missing per feature
miss_per_feat = np.isnan(prot_missing).mean(axis=0)
miss_per_sample = np.isnan(prot_missing).mean(axis=1)

fig, axes = plt.subplots(2, 3, figsize=(14, 8))
fig.suptitle('Missing Data Analysis and Imputation Strategies', fontsize=12, fontweight='bold')

# 1. Missing data heatmap (subset)
sns.heatmap(np.isnan(prot_missing[:, :50]).T, ax=axes[0,0],
            cmap='RdYlGn_r', xticklabels=False, yticklabels=False,
            cbar_kws={'label': 'Missing (1=yes)'})
axes[0,0].set_xlabel('Samples'); axes[0,0].set_ylabel('Features (first 50)')
axes[0,0].set_title('Missing Data Pattern')

# 2. Missing rate distribution per feature
axes[0,1].hist(miss_per_feat * 100, bins=25, color='steelblue', alpha=0.8)
axes[0,1].axvline(50, color='red', linestyle='--', label='50% threshold')
axes[0,1].set_xlabel('Missing rate per feature (%)')
axes[0,1].set_ylabel('Count'); axes[0,1].set_title('Feature Missing Rate Distribution')
axes[0,1].legend()

# 3. Missing rate per sample
axes[0,2].bar(range(n_samples), miss_per_sample * 100, color='salmon', alpha=0.8)
axes[0,2].set_xlabel('Sample index'); axes[0,2].set_ylabel('Missing rate (%)')
axes[0,2].set_title('Sample Missing Rate')
axes[0,2].axhline(40, color='red', linestyle='--', label='40% threshold')
axes[0,2].legend()

# 4-6. Imputation methods comparison
methods = {'KNN\n(k=5)': KNNImputer(n_neighbors=5).fit_transform(prot_missing),
           'Mean': np.where(np.isnan(prot_missing),
                            np.nanmean(prot_missing, axis=0), prot_missing),
           'Min/2\n(MNAR)': np.where(np.isnan(prot_missing),
                                      np.nanmin(prot_missing, axis=0) / 2, prot_missing)}

for ax, (method_name, imputed) in zip(axes[1, :], methods.items()):
    # Compare original vs imputed at previously-observed positions
    non_missing_idx = ~missing_mask
    original_vals = prot_matrix[non_missing_idx][:500]
    imputed_vals = imputed[non_missing_idx][:500]
    ax.scatter(original_vals, imputed_vals, alpha=0.3, s=8)
    lo, hi = original_vals.min(), original_vals.max()
    ax.plot([lo, hi], [lo, hi], 'r--', linewidth=1)
    from scipy.stats import pearsonr
    r, _ = pearsonr(original_vals, imputed_vals)
    ax.set_title(f'{method_name}\nR = {r:.3f}')
    ax.set_xlabel('Original value'); ax.set_ylabel('Imputed value')

plt.tight_layout()
plt.show()
print("\nKNN imputation generally performs best for MAR data.")
print("Min/2 (half-minimum) is preferred for MNAR (missing not at random) in proteomics.")
```python

## 4. Batch Effect Detection and Correction

Batch effects arise from technical variation: different sequencing runs, laboratory days, operators, reagent lots.

### Detection methods
- **PCA**: samples cluster by batch, not biology → clear batch effect
- **PVCA** (Principal Variance Component Analysis): quantifies proportion of variance from batch vs. biology
- **RLE plot** (Relative Log Expression): systematic shifts indicate batch effects

### Correction methods

| Tool | Method | Use case |
|---|---|---|
| **ComBat** | Empirical Bayes | Single covariate batch correction |
| **ComBat-seq** | Negative binomial | RNA-seq count data |
| **Harmony** | Iterative correction | scRNA-seq integration |
| **limma removeBatchEffect** | Linear model | After normalization |
| **MNN** | Mutual nearest neighbors | Single-cell multi-batch |

### Important caveat
**Never correct batch when it is confounded with biology** (e.g., all tumor samples in Batch1, all normals in Batch2). Always check confounding with `table(batch, type)` first.

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
