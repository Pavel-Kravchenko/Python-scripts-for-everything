---
name: bio-applied-mofa2
description: "MOFA2 unsupervised multi-omics factor analysis: variance decomposition, factor interpretation, and shared/view-specific signal separation. Use when integrating multiple omics layers."
tool_type: python
primary_tool: Matplotlib
---

# MOFA2: Multi-Omics Factor Analysis

- [MOFA2 documentation](https://biofam.github.io/MOFA2/)
- [mofapy2 Python package](https://github.com/bioFAM/mofapy2)

## Generative Model

X_m = Z @ W_m.T + noise, for each omics view m.
- **Z** (N x K): shared latent factors across views
- **W_m** (D_m x K): view-specific feature weights

## Usage

```python
# Via mofapy2
from mofapy2.run.entry_point import entry_point
ent = entry_point()
ent.set_data_options(scale_groups=False, scale_views=True)
ent.set_data_df(data_df)  # long-format: sample, group, feature, view, value
ent.set_model_options(factors=10)
ent.set_train_options(iter=1000, convergence_mode="fast")
ent.build(); ent.run(); ent.save("mofa_model.hdf5")

# Via muon (recommended)
import muon as mu
mdata = mu.MuData({'rna': rna_adata, 'prot': prot_adata})
mu.tl.mofa(mdata, n_factors=10, outfile="model.hdf5")
```

## Model Parameters

| Parameter | Default | Notes |
|---|---|---|
| `n_factors` | 10 | Number of latent factors |
| `convergence_mode` | "fast" | fast/medium/slow |
| `spikeslab_weights` | True | Sparse weights (ARD prior) |
| `scale_views` | True | Normalize each view to unit variance |

## Variance Decomposition (R^2 per Factor per View)

| Pattern | Interpretation |
|---|---|
| High R^2 in all views | Shared biological signal |
| High R^2 in one view only | View-specific signal (e.g., epigenetic reprogramming) |

Optimal factors: start with 15-20, prune inactive (R^2 < 2%). Robustness: factors consistent across random seeds.

## Factor Interpretation Workflow

1. Correlate factor scores with metadata (sample type, clinical variables)
2. Top weighted features -> candidate genes/proteins/CpGs
3. Gene set enrichment on top RNA weights -> pathway interpretation

## MOFA2 vs PCA

| Aspect | PCA | MOFA2 |
|---|---|---|
| Input | Single matrix | Multiple matrices |
| Shared vs specific | Mixed | Explicit separation |
| Missing views | Cannot handle | Native support |
| Sparse weights | No | Yes (ARD prior) |

## Pitfalls

- **Factor 1 depth artifact**: Check whether Factor 1 correlates with library size rather than biology
- **Batch effects**: If batch confounds biology, MOFA will capture batch as a factor -- check factor-metadata correlations
- **View scaling**: Always set `scale_views=True` so high-dimensional views don't dominate
