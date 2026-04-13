---
name: bio-applied-gene-regulatory-networks
description: "**Tier 3 — Applied Bioinformatics | Module 28 · Notebook 3**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/28_Network_Biology/03_gene_regulatory_networks.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, networkx 3.2+, numpy 1.26+, pandas 2.1+, scikit-learn 1.4+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Gene Regulatory Network Inference

*Source: Course notebook `Tier_3_Applied_Bioinformatics/28_Network_Biology/03_gene_regulatory_networks.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 28 · Notebook 3**

*Prerequisites: Notebook 2 (Network Modules), Module 03 (RNA-seq Analysis)*

---

**By the end of this notebook you will be able to:**
1. Distinguish GRN inference from PPI network construction
2. Apply GENIE3 (Random Forest) to infer regulatory edges from expression data
3. Identify network motifs (feed-forward loops, autoregulatory motifs)
4. Validate inferred regulations against known TF binding data (JASPAR/ENCODE)
5. Visualize a subnetwork around a master regulator TF



**Key resources:**
- [GENIE3 documentation](https://bioconductor.org/packages/release/bioc/html/GENIE3.html)
- [ARACNE web server](http://califano.c2b2.columbia.edu/aracne)
- [TRRUST — curated TF-target database](https://www.grnpedia.org/trrust/)

## 1. Gene Regulatory Networks vs. PPI Networks

**Gene Regulatory Networks (GRNs)** differ from PPI networks in fundamental ways:

| Aspect | PPI Network | GRN |
|---|---|---|
| **Edges** | Physical protein-protein binding | Transcriptional regulation |
| **Direction** | Undirected (usually) | Directed (TF → target) |
| **Weight** | Interaction confidence | Regulation strength/direction |
| **Data source** | Experimental (Y2H, co-IP) | Expression data + ChIP-seq |
| **Database** | STRING, BioGRID | TRRUST, ChEA, DoRothEA |

### GRN inference approaches

| Method | Algorithm | Captures |
|---|---|---|
| **Correlation** | Pearson/Spearman r | Linear co-variation |
| **Mutual information** (ARACNE) | MI + DPI pruning | Non-linear relationships |
| **Random forest** (GENIE3) | Feature importance | Complex dependencies |
| **Bayesian network** (BANJO) | Structure learning | Causal structure |
| **scGRN** (SCENIC+) | chromatin + expression | Cell-type-specific regulation |

### Key concepts
- **Transcription factor (TF)**: protein that binds DNA to regulate transcription
- **Regulon**: set of genes regulated by a TF
- **Network motif**: recurring sub-graph pattern with functional significance
- **Feed-forward loop (FFL)**: A → B → C and A → C (most common motif)

```python
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from collections import defaultdict

np.random.seed(42)

# ----- Simulate expression data for GRN inference -----
n_samples = 100   # cells or samples
n_tfs = 15        # transcription factors
n_targets = 50    # target genes

# TF names
tf_names = ['TP53', 'MYC', 'E2F1', 'FOXM1', 'SP1', 'NF1', 'YAP1',
            'TWIST1', 'SNAI1', 'ZEB1', 'HIF1A', 'STAT3', 'NFkB1', 'AR', 'ESR1']
target_names = [f'Gene_{i:03d}' for i in range(1, n_targets + 1)]
all_genes = tf_names + target_names
n_all = len(all_genes)

# True GRN: sparse TF-target adjacency matrix
# Each TF regulates 3-8 targets (activation or repression)
true_grn = np.zeros((n_tfs, n_targets))
for i in range(n_tfs):
    n_reg = np.random.randint(2, 8)
    targets = np.random.choice(n_targets, n_reg, replace=False)
    directions = np.random.choice([-1, 1], n_reg)  # repression or activation
    true_grn[i, targets] = directions

# Simulate expression: targets driven by TFs + noise
tf_expression = np.random.randn(n_samples, n_tfs)
# Add co-regulation structure
tf_expression[:50, :5] += 1.5   # group1 high TF1-5
tf_expression[50:, 5:10] += 1.5  # group2 high TF6-10

target_expression = tf_expression @ true_grn.T
target_expression += np.random.randn(n_samples, n_targets) * 0.8

# Full expression matrix (samples x genes)
expr_matrix = np.hstack([tf_expression, target_expression])
df_expr = pd.DataFrame(expr_matrix, columns=all_genes)

print(f"Expression matrix: {df_expr.shape} (samples x genes)")
print(f"Transcription factors: {n_tfs}")
print(f"Target genes: {n_targets}")
print(f"True regulatory edges: {int(np.abs(true_grn) > 0).sum()}")
print(f"Sparsity: {1 - (np.abs(true_grn)>0).mean():.2%}")
```python

## 2. Correlation-Based GRN Inference

The simplest approach: compute **Pearson or Spearman correlation** between TF and target gene expression.

### Advantages
- Fast and interpretable
- Easily scaled to genome-wide analysis
- No hyperparameters (just threshold on r)

### Limitations
- Cannot distinguish direct from indirect regulation
- Symmetric: does not infer direction (A→B or B→A)
- Sensitive to outliers (Spearman more robust)
- High false positive rate in dense expression data

### Partial correlation
Partial correlation between TF and target **conditioned on all other TFs** removes indirect associations. More computationally expensive but fewer false positives. Used in **GeneNet**, **PCIT** packages.

```python
from sklearn.covariance import GraphicalLassoCV
glasso = GraphicalLassoCV(cv=5).fit(tf_expr)
partial_corr = -glasso.precision_ / np.outer(np.sqrt(np.diag(glasso.precision_)),
                                               np.sqrt(np.diag(glasso.precision_)))
```python

```python
# ----- Correlation-based GRN inference -----
from scipy.stats import pearsonr, spearmanr

# Compute TF-target correlation matrix
n_tfs_loc = len(tf_names)
tf_expr = df_expr[tf_names].values
tgt_expr = df_expr[target_names].values

# Pearson correlation for each TF-target pair
corr_matrix = np.zeros((n_tfs_loc, n_targets))
pval_matrix = np.zeros((n_tfs_loc, n_targets))

for i in range(n_tfs_loc):
    for j in range(n_targets):
        r, p = pearsonr(tf_expr[:, i], tgt_expr[:, j])
        corr_matrix[i, j] = r
        pval_matrix[i, j] = p

# Threshold: |r| > 0.25 and p < 0.01
r_threshold = 0.25
p_threshold = 0.01
significant = (np.abs(corr_matrix) > r_threshold) & (pval_matrix < p_threshold)

print(f"Correlation-based edges (|r|>0.25, p<0.01): {significant.sum()}")

# Evaluate against true GRN
true_binary = np.abs(true_grn) > 0
tp = (significant & true_binary).sum()
fp = (significant & ~true_binary).sum()
fn = (~significant & true_binary).sum()
tn = (~significant & ~true_binary).sum()

precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0
f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

print(f"Evaluation vs. true GRN:")
print(f"  Precision: {precision:.3f}")
print(f"  Recall:    {recall:.3f}")
print(f"  F1:        {f1:.3f}")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Correlation-Based GRN Inference', fontsize=12, fontweight='bold')

# Correlation heatmap
im = axes[0].imshow(corr_matrix, cmap='RdBu_r', vmin=-0.7, vmax=0.7, aspect='auto')
axes[0].set_title('TF-Target Correlation Matrix (Pearson r)')
axes[0].set_xlabel('Target genes'); axes[0].set_ylabel('TFs')
axes[0].set_yticks(range(n_tfs_loc)); axes[0].set_yticklabels(tf_names, fontsize=7)
plt.colorbar(im, ax=axes[0], label='Pearson r')

# True GRN
im2 = axes[1].imshow(true_grn, cmap='RdBu_r', vmin=-1.5, vmax=1.5, aspect='auto')
axes[1].set_title('True GRN (simulated)')
axes[1].set_xlabel('Target genes'); axes[1].set_ylabel('TFs')
axes[1].set_yticks(range(n_tfs_loc)); axes[1].set_yticklabels(tf_names, fontsize=7)
plt.colorbar(im2, ax=axes[1], label='Regulation direction')

# Precision-recall at varying thresholds
thresholds = np.linspace(0, 0.8, 40)
precisions, recalls = [], []
for t in thresholds:
    sig = np.abs(corr_matrix) > t
    tp_t = (sig & true_binary).sum()
    fp_t = (sig & ~true_binary).sum()
    fn_t = (~sig & true_binary).sum()
    prec = tp_t / (tp_t + fp_t) if (tp_t + fp_t) > 0 else 1
    rec = tp_t / (tp_t + fn_t) if (tp_t + fn_t) > 0 else 0
    precisions.append(prec); recalls.append(rec)

axes[2].plot(recalls, precisions, 'b-', linewidth=2)
axes[2].fill_between(recalls, precisions, alpha=0.2)
axes[2].axhline(true_binary.mean(), color='red', linestyle='--', label='Random (baseline)')
axes[2].scatter([recall], [precision], color='orange', s=100, zorder=5, label=f'r>0.25 threshold')
axes[2].set_xlabel('Recall'); axes[2].set_ylabel('Precision')
axes[2].set_title('Precision-Recall Curve\n(threshold = |r|)')
axes[2].legend(fontsize=8)

plt.tight_layout()
plt.show()
```python

## 3. Mutual Information: ARACNE Algorithm

**ARACNE** (Algorithm for the Reconstruction of Accurate Cellular Networks) uses **mutual information** to infer edges and then prunes indirect edges using the **Data Processing Inequality (DPI)**.

### Mutual Information vs. Correlation
- MI captures **non-linear** relationships (correlation only captures linear)
- MI = 0 means statistical independence
- MI is symmetric (like correlation)

### Data Processing Inequality (DPI)
For three variables A, B, C forming a chain A → B → C:
```python
MI(A, B) >= MI(A, C)  and  MI(B, C) >= MI(A, C)
```python
ARACNE removes the weakest edge in any triangle of 3 highly connected genes, keeping only direct regulatory relationships.

### VIPER: Protein activity inference from GRN
VIPER uses the inferred regulon (TF → target weights from ARACNE) to infer TF **activity** (not just expression) from gene expression data:
```python
# In R: viper(expr_matrix, regulon=aracne_output, method="ttest")
```python

```python
# ----- Mutual Information-based GRN (ARACNE concept) -----
from sklearn.feature_selection import mutual_info_regression

# ARACNE: Algorithm for the Reconstruction of Accurate Cellular Networks
# Step 1: Compute mutual information between all TF-target pairs
# Step 2: Data Processing Inequality (DPI): remove indirect edges
#         if MI(A,B) > MI(A,C) and MI(B,C) > 0: remove the weakest edge A-C

# Compute MI for each TF-target pair
mi_matrix = np.zeros((n_tfs, n_targets))
for i in range(n_tfs):
    # mutual_info_regression: MI between TF_i and all targets
    mi_vals = mutual_info_regression(tf_expr[:, i:i+1], tgt_expr,
                                      discrete_features=False, random_state=42)
    mi_matrix[i] = mi_vals

# Threshold at top 20% MI values per TF
edges_mi = []
for i in range(n_tfs):
    threshold_mi = np.percentile(mi_matrix[i], 80)
    for j in range(n_targets):
        if mi_matrix[i, j] >= threshold_mi:
            edges_mi.append({'TF': tf_names[i], 'Target': target_names[j],
                              'MI': mi_matrix[i, j]})

df_mi_edges = pd.DataFrame(edges_mi)

# Evaluate
mi_binary = np.zeros_like(true_binary)
for _, row in df_mi_edges.iterrows():
    i = tf_names.index(row['TF'])
    j = target_names.index(row['Target'])
    mi_binary[i, j] = 1

tp_mi = (mi_binary & true_binary).sum()
fp_mi = (mi_binary & ~true_binary).sum()
fn_mi = (~mi_binary & true_binary).sum()
prec_mi = tp_mi / (tp_mi + fp_mi) if (tp_mi + fp_mi) > 0 else 0
rec_mi = tp_mi / (tp_mi + fn_mi) if (tp_mi + fn_mi) > 0 else 0

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Mutual Information (ARACNE-style) GRN Inference', fontsize=12, fontweight='bold')

# MI matrix heatmap
im = axes[0].imshow(mi_matrix, cmap='YlOrRd', aspect='auto')
axes[0].set_title('MI Matrix (TF vs Targets)')
axes[0].set_xlabel('Target genes'); axes[0].set_ylabel('TFs')
axes[0].set_yticks(range(n_tfs)); axes[0].set_yticklabels(tf_names, fontsize=7)
plt.colorbar(im, ax=axes[0], label='Mutual Information')

# Compare precision-recall: Pearson vs MI
# Method comparison bar chart
methods_data = {
    'Pearson\nCorrelation': (precision, recall, f1),
    'Mutual\nInformation': (prec_mi, rec_mi,
                             2*prec_mi*rec_mi/(prec_mi+rec_mi) if (prec_mi+rec_mi)>0 else 0),
}

x = np.arange(len(methods_data))
width = 0.25
bar_prec = [v[0] for v in methods_data.values()]
bar_rec = [v[1] for v in methods_data.values()]
bar_f1 = [v[2] for v in methods_data.values()]

axes[1].bar(x - width, bar_prec, width, label='Precision', color='steelblue', alpha=0.8)
axes[1].bar(x, bar_rec, width, label='Recall', color='orange', alpha=0.8)
axes[1].bar(x + width, bar_f1, width, label='F1', color='green', alpha=0.8)
axes[1].set_xticks(x); axes[1].set_xticklabels(methods_data.keys())
axes[1].set_ylabel('Score')
axes[1].set_title('GRN Inference Method Comparison')
axes[1].legend()
axes[1].set_ylim(0, 1)

plt.tight_layout()
plt.show()

print(f"MI-based: Precision={prec_mi:.3f}, Recall={rec_mi:.3f}")
print(f"Pearson:  Precision={precision:.3f}, Recall={recall:.3f}")
print(f"\nARACNE applies Data Processing Inequality (DPI) to prune indirect edges.")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
