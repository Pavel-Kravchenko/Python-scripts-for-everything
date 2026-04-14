---
name: bio-applied-gene-regulatory-networks
description: "GRN inference methods: correlation, mutual information (ARACNE), and random forest (GENIE3). Decision table for method selection, evaluation patterns, and key pitfalls."
tool_type: python
primary_tool: Matplotlib
---

# Gene Regulatory Network Inference

## GRN vs PPI Networks

| Aspect | PPI Network | GRN |
|--------|-------------|-----|
| Edges | Physical protein binding | Transcriptional regulation |
| Direction | Undirected | Directed (TF → target) |
| Data source | Y2H, co-IP | Expression + ChIP-seq |
| Databases | STRING, BioGRID | TRRUST, ChEA, DoRothEA |

## Method Comparison

| Method | Captures | Pros | Cons |
|--------|----------|------|------|
| Pearson/Spearman correlation | Linear co-variation | Fast, interpretable | Symmetric; indirect edges |
| Partial correlation (GeneNet) | Direct linear links | Fewer false positives | Requires n >> p |
| ARACNE (mutual information + DPI) | Non-linear relationships | Prunes indirect edges | Symmetric |
| GENIE3 (random forest) | Complex dependencies | Directional importance scores | Slow, no sign (activation/repression) |
| SCENIC+ | Cell-type-specific regulation | Integrates chromatin + expression | Complex setup |

**Choosing a method:** Start with correlation for speed. Use GENIE3 when you need directionality scores. Use ARACNE when non-linearity is suspected. Use SCENIC+ for single-cell multiome data.

## Correlation-Based Inference

```python
from scipy.stats import pearsonr
import numpy as np

def infer_grn_correlation(tf_expr, target_expr, tf_names, target_names,
                           r_thresh=0.25, p_thresh=0.01):
    """Return DataFrame of TF-target edges above thresholds."""
    edges = []
    for i, tf in enumerate(tf_names):
        for j, tgt in enumerate(target_names):
            r, p = pearsonr(tf_expr[:, i], target_expr[:, j])
            if abs(r) > r_thresh and p < p_thresh:
                edges.append({'TF': tf, 'Target': tgt, 'r': r, 'p': p})
    return edges
```

## Mutual Information / ARACNE

```python
from sklearn.feature_selection import mutual_info_regression

def compute_mi_matrix(tf_expr, target_expr):
    """MI between each TF and all targets."""
    n_tfs = tf_expr.shape[1]
    mi = np.zeros((n_tfs, target_expr.shape[1]))
    for i in range(n_tfs):
        mi[i] = mutual_info_regression(
            tf_expr[:, i:i+1], target_expr,
            discrete_features=False, random_state=42
        )
    return mi

# ARACNE DPI: for triplet (A, B, C), remove weakest MI edge
# if MI(A,C) < min(MI(A,B), MI(B,C)): remove edge A-C
```

## Partial Correlation (Graphical Lasso)

```python
from sklearn.covariance import GraphicalLassoCV
import numpy as np

glasso = GraphicalLassoCV(cv=5).fit(tf_expr)
prec = glasso.precision_
d = np.sqrt(np.diag(prec))
partial_corr = -prec / np.outer(d, d)
np.fill_diagonal(partial_corr, 1.0)
```

## Evaluation Pattern

```python
def evaluate_grn(predicted_binary, true_binary):
    tp = (predicted_binary & true_binary).sum()
    fp = (predicted_binary & ~true_binary).sum()
    fn = (~predicted_binary & true_binary).sum()
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return {'precision': precision, 'recall': recall, 'f1': f1}
```

## Key Concepts

- **Regulon**: set of genes regulated by one TF
- **Feed-forward loop (FFL)**: A→B→C and A→C — most common 3-node network motif
- **VIPER**: infers TF *activity* (not expression) from regulon weights — more informative than TF mRNA level
- **PAGA score**: (observed − expected inter-cluster edges) / expected; used to assess trajectory connectivity

## Pitfalls

- **Correlation ≠ regulation**: two genes co-expressed due to shared upstream regulator will appear as direct edges
- **Symmetric methods (correlation, MI) cannot infer direction**: need additional evidence (ChIP-seq, perturbation) for A→B vs B→A
- **High false positive rate**: GRN inference on bulk RNA-seq is noisy; always validate top edges with ChIP-seq or ENCODE data
- **Sample size**: reliable GRN inference requires n >> number of TFs; partial correlation fails when n < p
- **Batch effects**: always check for batch confounding before interpreting co-expression signal
- **Multiple testing**: apply FDR correction when testing thousands of TF-target pairs
