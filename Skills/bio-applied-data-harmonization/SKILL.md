---
name: bio-applied-data-harmonization
description: "Multi-omics data harmonization — normalization strategies, missing data imputation, batch correction, and integration approaches (MOFA2, DIABLO)"
tool_type: python
primary_tool: Matplotlib
---

# Data Harmonization for Multi-Omics

- [MOFA2 documentation](https://biofam.github.io/MOFA2/)
- [mixOmics documentation](http://mixomics.org/)

## Integration Strategies

- **Concatenation (early)**: merge all matrices -> dimension reduction
- **Meta-analysis (late)**: analyze each omics separately -> combine results
- **Factor analysis (intermediate)**: MOFA2, iCluster — learn shared latent factors
- **Supervised (DIABLO)**: use sample labels to guide integration

## Normalization by Data Type

| Data type | Recommended normalization | Rationale |
|---|---|---|
| RNA-seq counts | TMM/DESeq2 VST -> z-score | Library size + biological variation |
| Proteomics intensities | Median centering + z-score | Inter-sample variation |
| Methylation beta | logit (M-value) -> z-score | Beta bounded 0-1 |
| Metabolomics | Probabilistic quotient + log | Dilution effects |

**Key principle**: normalize *within* each omics layer before integration.

## Missing Data Imputation

| Mechanism | Description | Strategy |
|---|---|---|
| MCAR | Random technical failures | KNN or mean imputation |
| MAR | Depends on observed data | KNN imputation |
| MNAR | Below detection limit | min/2 imputation |

Proteomics missing values (~20-40%) are predominantly MNAR. Best practice:
1. Filter proteins missing in >50% of samples
2. Half-minimum imputation for MNAR; KNN for MAR
3. Flag imputed values for sensitivity analysis

```python
from sklearn.impute import KNNImputer

# KNN imputation (MAR)
imputed_knn = KNNImputer(n_neighbors=5).fit_transform(prot_missing)

# Min/2 imputation (MNAR)
imputed_mnar = np.where(np.isnan(prot_missing),
                        np.nanmin(prot_missing, axis=0) / 2, prot_missing)
```

## Batch Effect Correction

**Detection**: PCA (samples cluster by batch, not biology), PVCA (variance decomposition), RLE plot.

| Tool | Method | Use case |
|---|---|---|
| ComBat | Empirical Bayes | Single covariate batch correction |
| ComBat-seq | Negative binomial | RNA-seq count data |
| Harmony | Iterative correction | scRNA-seq integration |
| limma removeBatchEffect | Linear model | After normalization |
| MNN | Mutual nearest neighbors | Single-cell multi-batch |

**Never correct batch when confounded with biology** (all tumor in Batch1, all normal in Batch2). Always check with `table(batch, type)` first.

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
