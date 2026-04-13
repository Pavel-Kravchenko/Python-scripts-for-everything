---
name: bio-applied-mixomics
description: "**Tier 3 — Applied Bioinformatics | Module 27 · Notebook 3**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/27_Multi_Omics_Integration/03_mixomics.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scikit-learn 1.4+, scipy 1.12+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# mixOmics: Supervised Multi-Omics Integration

*Source: Course notebook `Tier_3_Applied_Bioinformatics/27_Multi_Omics_Integration/03_mixomics.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 27 · Notebook 3**

*Prerequisites: Notebook 2 (MOFA2), Module 07 (Machine Learning)*

---

**By the end of this notebook you will be able to:**
1. Apply PLS-DA for supervised classification with a single omics layer
2. Use sparse PLS-DA (sPLS-DA) for simultaneous classification and feature selection
3. Build a DIABLO model for supervised multi-block integration
4. Evaluate model performance with cross-validation
5. Interpret loadings plots and correlation circles



**Key resources:**
- [mixOmics documentation](http://mixomics.org/)
- [mixOmics paper (Rohart et al., 2017)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005752)
- [DIABLO tutorial](http://mixomics.org/mixdiablo/)

## 1. mixOmics Framework Overview

**mixOmics** is an R package (with Python concepts applicable) for supervised multi-omics integration. Its core methods all derive from **Partial Least Squares (PLS)**.

### The PLS framework

PLS finds **latent components** that maximize the covariance between X (features) and Y (outcome). This makes it ideal for high-dimensional data where p >> n.

```python
X (samples × features)  →  T (samples × components)  →  Y (samples × outcome)
         "scores"                "predictors"
         "loadings"
```python

### mixOmics method family

| Method | Data | Task |
|---|---|---|
| **PLS-DA** | Single omics | Supervised classification |
| **sPLS-DA** | Single omics | Classification + feature selection |
| **MINT** | Multi-study | Cross-study integration |
| **DIABLO** | Multi-omics | Supervised multi-block classification |
| **PLS** | Single omics | Regression |
| **sPLS** | Single omics | Regression + feature selection |

### Key concepts
- **Components**: latent variables summarizing variation (like PCA components)
- **Loadings**: contribution of each feature to each component
- **Scores**: sample projection onto components
- **Balanced Error Rate (BER)**: classification error accounting for class imbalance

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import balanced_accuracy_score
from scipy.linalg import svd

np.random.seed(42)

# ---- Generate multi-omics classification dataset ----
# Task: classify 3 cancer subtypes from RNA + proteomics + methylation
n_per_class = 25
n_classes = 3
n_samples = n_per_class * n_classes

n_rna = 300
n_prot = 100
n_meth = 200

class_labels = ['Luminal_A', 'Luminal_B', 'TNBC'] * n_per_class
np.random.shuffle(class_labels)

def make_data_with_signal(n_samples, n_features, class_labels, class_signals, noise=0.8):
    X = np.random.randn(n_samples, n_features) * noise
    for i, label in enumerate(class_labels):
        if label in class_signals:
            feat_idx, direction, strength = class_signals[label]
            X[i, feat_idx] += direction * strength
    return X

# Define which features define each class
rna_signals = {
    'Luminal_A':  (slice(0, 60),   1,  2.0),
    'Luminal_B':  (slice(60, 120), 1,  2.0),
    'TNBC':       (slice(120, 180), 1, 2.5),
}
prot_signals = {
    'Luminal_A':  (slice(0, 20),  1, 2.0),
    'Luminal_B':  (slice(20, 40), 1, 2.0),
    'TNBC':       (slice(40, 60), 1, 2.5),
}
meth_signals = {
    'Luminal_A':  (slice(0, 40),  1, 1.5),
    'Luminal_B':  (slice(40, 80), 1, 1.5),
    'TNBC':       (slice(80, 120),1, 2.0),
}

rna_data  = make_data_with_signal(n_samples, n_rna,  class_labels, rna_signals)
prot_data = make_data_with_signal(n_samples, n_prot, class_labels, prot_signals)
meth_data = make_data_with_signal(n_samples, n_meth, class_labels, meth_signals)

rna_scaled  = StandardScaler().fit_transform(rna_data)
prot_scaled = StandardScaler().fit_transform(prot_data)
meth_scaled = StandardScaler().fit_transform(meth_data)

le = LabelEncoder()
y = le.fit_transform(class_labels)
y_onehot = pd.get_dummies(class_labels).values

print(f"Dataset: {n_samples} samples, {n_classes} classes")
print(f"  RNA-seq:     {rna_scaled.shape}")
print(f"  Proteomics:  {prot_scaled.shape}")
print(f"  Methylation: {meth_scaled.shape}")
print(f"  Classes: {list(le.classes_)}")
print(f"  Balanced: {dict(pd.Series(class_labels).value_counts())}")
```python

## 2. PLS-DA: Single Omics Classification

**PLS-DA** (Partial Least Squares Discriminant Analysis) is the workhorse for supervised omics analysis.

### How PLS-DA works
1. Encode class labels as dummy matrix Y (one-hot encoding)
2. Find components that maximize covariance between X and Y
3. Project samples onto first K components for visualization
4. Classify new samples by soft-max on predicted Y scores

### Number of components
- Typically 2-5 components (each extracts progressively less signal)
- Determined by cross-validation: minimize BER
- Rule: first 2 components for visualization; add more for prediction

### Assumptions and caveats
- Sensitive to overfitting when p >> n (use sPLS-DA instead)
- Requires standardized features (mean=0, var=1)
- Multi-class: one-vs-all approach in component extraction

```python
# ----- PLS-DA: Partial Least Squares Discriminant Analysis -----

class PLSDA:
    """PLS-DA implementation: PLS regression with one-hot Y."""
    def __init__(self, n_components=2):
        self.n_comp = n_components
        self.pls = PLSRegression(n_components=n_components, scale=True)
    
    def fit(self, X, Y):
        self.pls.fit(X, Y)
        return self
    
    def transform(self, X):
        return self.pls.transform(X)[0]
    
    def predict_class(self, X):
        Y_pred = self.pls.predict(X)
        return np.argmax(Y_pred, axis=1)

# Apply PLS-DA to RNA-seq data
plsda_rna = PLSDA(n_components=3)
plsda_rna.fit(rna_scaled, y_onehot)
scores_rna = plsda_rna.transform(rna_scaled)

plsda_prot = PLSDA(n_components=3)
plsda_prot.fit(prot_scaled, y_onehot)
scores_prot = plsda_prot.transform(prot_scaled)

class_colors = {'Luminal_A': '#3498DB', 'Luminal_B': '#F39C12', 'TNBC': '#E74C3C'}
point_colors = [class_colors[c] for c in class_labels]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('PLS-DA per Omics Layer', fontsize=13, fontweight='bold')

for ax, scores, title in zip(axes[:2], [scores_rna, scores_prot], ['RNA-seq', 'Proteomics']):
    ax.scatter(scores[:, 0], scores[:, 1], c=point_colors, s=40, alpha=0.8,
               edgecolors='black', linewidths=0.3)
    ax.set_xlabel('Component 1'); ax.set_ylabel('Component 2')
    ax.set_title(f'PLS-DA: {title}')

# Cross-validation for each view
def cv_plsda(X, y_oh, y, n_comp=3, cv=5):
    skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=42)
    accs = []
    for train, test in skf.split(X, y):
        model = PLSDA(n_components=n_comp)
        model.fit(X[train], y_oh[train])
        y_pred = model.predict_class(X[test])
        accs.append(balanced_accuracy_score(y[test], y_pred))
    return np.mean(accs), np.std(accs)

results = {}
for name, X in [('RNA-seq', rna_scaled), ('Proteomics', prot_scaled), ('Methylation', meth_scaled)]:
    mean_acc, std_acc = cv_plsda(X, y_onehot, y)
    results[name] = (mean_acc, std_acc)

bars = axes[2].bar(results.keys(), [v[0]*100 for v in results.values()],
                   yerr=[v[1]*100 for v in results.values()],
                   color=['steelblue', 'orange', 'green'], alpha=0.8, capsize=5)
axes[2].axhline(100/n_classes, color='red', linestyle='--', label=f'Random ({100/n_classes:.0f}%)')
axes[2].set_ylabel('Balanced Accuracy (%)')
axes[2].set_title('5-fold CV Accuracy per Omics Layer')
axes[2].legend()

handles = [mpatches.Patch(color=v, label=k) for k, v in class_colors.items()]
axes[1].legend(handles=handles, fontsize=8, loc='lower right')

plt.tight_layout()
plt.show()
for name, (m, s) in results.items():
    print(f"  {name}: {m*100:.1f} +/- {s*100:.1f}%")
```python

## 3. sPLS-DA: Sparse Feature Selection

**sPLS-DA** (sparse PLS-DA) adds an L1-like penalty to select informative features. This is critical for omics data where most features are noise.

### Sparsity in mixOmics
The `keepX` parameter controls how many features to retain per component:
- Tuned by cross-validation: choose `keepX` that minimizes BER
- Different `keepX` values per component allowed
- Features selected consistently across bootstrap iterations = "stable"

### Feature stability
Not all selected features are equally reliable. Bootstrap stability measures:
- Selection frequency > 50%: reliably selected
- Selection frequency > 80%: very stable (high confidence biomarker candidates)

### Tuning keepX in R/Python
```r
# R mixOmics
test.keepX <- c(5, 10, 20, 50, 100)
tune.splsda <- tune.splsda(X, Y, ncomp=3, validation='Mfold', 
                            folds=5, dist='max.dist', 
                            test.keepX=test.keepX, nrepeat=50)
optimal.keepX <- tune.splsda$choice.keepX
```python

```python
# ----- sPLS-DA: Sparse PLS-DA for feature selection -----

class SPLSDA:
    """Sparse PLS-DA: penalized regression coefficients."""
    def __init__(self, n_components=3, keepX=50):
        self.n_comp = n_components
        self.keepX = keepX  # features to keep per component
        self.models = []
        self.selected_features = {}
    
    def fit(self, X, Y):
        self.models = []
        self.selected_features = {}
        X_deflated = X.copy()
        
        for k in range(self.n_comp):
            # PLS on current deflated X
            pls = PLSRegression(n_components=1, scale=False)
            pls.fit(X_deflated, Y)
            
            # Get weights (loadings)
            weights = np.abs(pls.x_weights_[:, 0])
            
            # Keep only top keepX features (soft thresholding)
            threshold_idx = np.argsort(weights)[-self.keepX:]
            sparse_weights = np.zeros(X.shape[1])
            sparse_weights[threshold_idx] = weights[threshold_idx]
            
            self.selected_features[f'comp_{k+1}'] = threshold_idx
            self.models.append((pls, sparse_weights))
            
            # Deflate X
            t = X_deflated @ sparse_weights
            t = t / (t @ t + 1e-10)
            X_deflated = X_deflated - np.outer(t, X_deflated.T @ t)
        
        return self
    
    def get_selected_features(self):
        return self.selected_features

# Apply sPLS-DA to RNA data
splsda = SPLSDA(n_components=3, keepX=40)
splsda.fit(rna_scaled, y_onehot)

# Feature selection stability via bootstrap
n_boot = 30
all_selected = {f'comp_{k+1}': [] for k in range(3)}
for b in range(n_boot):
    idx = np.random.choice(n_samples, n_samples, replace=True)
    m = SPLSDA(n_components=3, keepX=40)
    m.fit(rna_scaled[idx], y_onehot[idx])
    for k in range(3):
        all_selected[f'comp_{k+1}'].extend(m.get_selected_features()[f'comp_{k+1}'])

# Stability scores
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('sPLS-DA Feature Selection', fontsize=13, fontweight='bold')

for col, (comp_name, selected) in enumerate(all_selected.items()):
    # Count selection frequency
    counts = np.bincount(selected, minlength=n_rna) / n_boot
    top20_idx = np.argsort(counts)[-20:]
    
    colors_bar = ['red' if counts[i] > 0.5 else 'steelblue' for i in top20_idx]
    axes[col].barh(range(20), counts[top20_idx], color=colors_bar, alpha=0.8)
    axes[col].axvline(0.5, color='red', linestyle='--', label='50% stability')
    axes[col].set_yticks(range(20))
    axes[col].set_yticklabels([f'Gene_{i}' for i in top20_idx], fontsize=7)
    axes[col].set_xlabel('Selection frequency')
    axes[col].set_title(f'{comp_name} Stable Features')
    axes[col].legend(fontsize=8)

plt.tight_layout()
plt.show()

stable_genes = {k: np.where(np.bincount(v, minlength=n_rna)/n_boot > 0.5)[0]
                for k, v in all_selected.items()}
print("Stably selected features (>50% bootstrap frequency):")
for comp, genes in stable_genes.items():
    print(f"  {comp}: {len(genes)} features")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
