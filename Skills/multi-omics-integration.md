---
name: multi-omics-integration
description: Multi-omics integration — data harmonization, MOFA2 latent factor analysis, mixOmics PLS-DA and DIABLO, batch correction, missing data imputation
---

# Multi-Omics Integration

## When to Use

Use this skill when:
- Integrating RNA-seq with proteomics, metabolomics, or methylation data
- Performing unsupervised multi-view dimensionality reduction (MOFA2)
- Building supervised multi-omics classifiers (DIABLO/mixOmics)
- Harmonizing data across platforms (batch correction, scaling)
- Identifying cross-omics biomarkers

## Quick Reference

| Task | Tool | Key Function |
|------|------|-------------|
| Scaling per feature | scikit-learn | `StandardScaler().fit_transform()` |
| Missing value imputation | scikit-learn | `KNNImputer(n_neighbors=5)` |
| Batch correction | pyComBat / sva (R) | `pycombat(data, batch)` |
| Multi-omics factor analysis | mofapy2 / MOFA2 | `entry_point()` + `ent.run()` |
| Load MOFA2 model | muon | `mu.read_h5mu('model.hdf5')` |
| PLS-DA (single view) | mixOmics (R) | `plsda(X, Y, ncomp=3)` |
| sPLS-DA feature selection | mixOmics (R) | `splsda(X, Y, keepX=c(50,30))` |
| Multi-block DIABLO | mixOmics (R) | `block.splsda(X=list(...), Y=y)` |
| Cross-validation | mixOmics (R) | `perf(model, validation='Mfold')` |

## Key Patterns

**Pattern 1: Data harmonization**
```python
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer

# Per-feature z-score scaling
scaler = StandardScaler()
rna_scaled = pd.DataFrame(
    scaler.fit_transform(rna.T).T,
    index=rna.index, columns=rna.columns
)

# KNN imputation for proteomics missing values
imputer = KNNImputer(n_neighbors=5)
prot_imputed = pd.DataFrame(
    imputer.fit_transform(prot.T).T,
    index=prot.index, columns=prot.columns
)
```

**Pattern 2: MOFA2 model training (Python)**
```python
from mofapy2.run.entry_point import entry_point

ent = entry_point()
ent.set_data_options(scale_groups=False, scale_views=True)
ent.set_model_options(factors=10, spikeslab_weights=True, ard_factors=True)
ent.set_train_options(iter=1000, convergence_mode='fast', seed=42)

# Data: list of views (samples × features matrices)
ent.set_data_matrix([[rna_scaled.values, prot_imputed.values, meth_scaled.values]],
                    views_names=['RNA', 'Proteomics', 'Methylation'])
ent.build()
ent.run()
ent.save('mofa2_model.hdf5')
```

**Pattern 3: MOFA2 analysis with muon**
```python
import muon as mu

model = mu.read_h5mu('mofa2_model.hdf5')
# Variance explained heatmap
mu.pl.mofa(model)
# Factor scores
factor_scores = model.obsm['X_mofa']  # samples × n_factors
```

**Pattern 4: DIABLO (R)**
```r
library(mixOmics)

# Design matrix: 0 = independent, 1 = fully correlated
design <- matrix(0.1, nrow=3, ncol=3,
                 dimnames=list(c('RNA','Prot','Meth'), c('RNA','Prot','Meth')))
diag(design) <- 0

diablo <- block.splsda(
    X = list(RNA=rna_scaled, Prot=prot_imputed, Meth=meth_scaled),
    Y = subtype,
    ncomp = 2,
    keepX = list(RNA=c(25,25), Prot=c(25,25), Meth=c(25,25)),
    design = design
)

# Evaluate
set.seed(42)
perf_res <- perf(diablo, validation='Mfold', folds=10, nrepeat=10, auc=TRUE)
plot(perf_res)

# Visualize
plotIndiv(diablo)
plotVar(diablo, overlap=FALSE)
cimDiablo(diablo)
```

## Method Selection Guide

| Scenario | Recommended Method |
|----------|-------------------|
| Exploratory, no outcome labels | MOFA2 |
| Classification with feature selection | mixOmics sPLS-DA or DIABLO |
| Survival prediction | MOFA2 factors + Cox model |
| 2 views, correlation structure | sPLS (mixOmics) |
| Many views (≥4), sparse data | MOFA2 |
| Multi-study integration | mixOmics MINT |

## Common Pitfalls

- **Scale before integration** — omics layers have very different value ranges; always scale per feature
- **Sample matching** — all views must share the same samples; handle missing samples carefully in MOFA2 (supports NAs)
- **Feature selection** — reduce to high-variance features (top 5000 per layer) before MOFA2 to reduce noise
- **DIABLO keepX tuning** — use `tune.block.splsda()` to select the optimal number of features per component
- **Biological interpretation** — MOFA2 factors and PLS loadings need GO enrichment to interpret; top weights alone are insufficient

## Code Templates

### Batch Correction with pyComBat
```python
import pandas as pd
from inmoose.pycombat import pycombat_norm

# data: features × samples DataFrame
# batch: list of batch labels per sample
data_corrected = pycombat_norm(data, batch)
# Returns features × samples DataFrame with batch effect removed
```

### Top-Variance Feature Selection
```python
import pandas as pd
import numpy as np

def select_top_variance(df, n_features=5000):
    """Select top N most variable features (features in rows)."""
    variances = df.var(axis=1)
    top_idx = variances.nlargest(n_features).index
    return df.loc[top_idx]

rna_hv = select_top_variance(rna_matrix, n_features=5000)
prot_hv = select_top_variance(prot_matrix, n_features=3000)
```

### Load MOFA2 Factor Scores
```python
import h5py
import numpy as np
import pandas as pd

def load_mofa2_factors(h5_path):
    """Extract factor scores from MOFA2 HDF5 output."""
    with h5py.File(h5_path, 'r') as f:
        # Factor scores: groups × samples × factors
        factors = f['expectations/Z/group1'][:]  # samples × factors
        factor_names = [f'Factor{i+1}' for i in range(factors.shape[1])]
        samples = [s.decode() for s in f['samples/group1'][:]]
    return pd.DataFrame(factors, index=samples, columns=factor_names)

factors_df = load_mofa2_factors('mofa2_model.hdf5')
print(factors_df.shape)  # (n_samples, n_factors)
```

### Variance Explained per Factor/View
```python
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def mofa2_variance_explained(h5_path):
    with h5py.File(h5_path, 'r') as f:
        r2 = f['variance_explained/r2_per_factor/group1'][:]
        views = [v.decode() for v in f['views/views'][:]]
        n_factors = r2.shape[1]
    df = pd.DataFrame(r2.T, columns=views,
                      index=[f'Factor{i+1}' for i in range(n_factors)])
    df.plot(kind='bar', figsize=(10, 4), title='Variance Explained by Factor/View')
    plt.tight_layout()
    return df
```

## Related Skills
- `rnaseq` — RNA-seq differential expression (DESeq2, edgeR)
- `ml-deep-learning-bio` — dimensionality reduction (PCA, UMAP, autoencoders)
- `numpy-pandas-wrangling` — DataFrame merging, scaling, handling missing data
- `cancer-transcriptomics` — tumor subtype classification, multi-omics cancer studies
