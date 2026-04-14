---
name: bio-applied-mixomics
description: "mixOmics PLS-DA and DIABLO for supervised multi-omics integration and feature selection. Use when classifying samples or selecting biomarkers from multi-omics data."
tool_type: python
primary_tool: Matplotlib
---

# mixOmics: Supervised Multi-Omics Integration

- [mixOmics documentation](http://mixomics.org/)
- [DIABLO tutorial](http://mixomics.org/mixdiablo/)

## Method Family

| Method | Data | Task |
|---|---|---|
| **PLS-DA** | Single omics | Supervised classification |
| **sPLS-DA** | Single omics | Classification + feature selection |
| **MINT** | Multi-study | Cross-study integration |
| **DIABLO** | Multi-omics | Supervised multi-block classification |

## PLS-DA

PLS finds latent components maximizing covariance between X (features) and Y (outcome). Ideal for p >> n.

1. One-hot encode class labels as Y
2. Find components maximizing cov(X, Y)
3. Classify by soft-max on predicted Y scores

Key parameters: 2-5 components (CV to minimize BER). Requires standardized features.

```python
from sklearn.cross_decomposition import PLSRegression

class PLSDA:
    def __init__(self, n_components=2):
        self.pls = PLSRegression(n_components=n_components, scale=True)
    def fit(self, X, Y_onehot):
        self.pls.fit(X, Y_onehot)
    def transform(self, X):
        return self.pls.transform(X)[0]
    def predict_class(self, X):
        return np.argmax(self.pls.predict(X), axis=1)
```

## sPLS-DA: Sparse Feature Selection

Adds L1-like penalty via `keepX` parameter (features retained per component).

- Tune keepX by CV: minimize BER
- Bootstrap stability: selection frequency >50% = reliable, >80% = high confidence biomarker

```r
# R mixOmics
tune.splsda <- tune.splsda(X, Y, ncomp=3, validation='Mfold',
                            folds=5, test.keepX=c(5,10,20,50,100), nrepeat=50)
```

## Pitfalls

- **Overfitting**: PLS-DA overfits when p >> n; use sPLS-DA with CV-tuned keepX
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Balanced Error Rate**: Use BER (not accuracy) for imbalanced classes
