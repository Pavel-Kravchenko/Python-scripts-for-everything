---
name: bio-applied-statistics-for-bioinformatics
description: "Statistical testing for bioinformatics: distributions, hypothesis testing, multiple testing correction (Bonferroni, BH/FDR), parametric vs non-parametric test selection"
tool_type: python
primary_tool: NumPy
---

# Statistics for Bioinformatics

## Distribution Selection Guide

| Distribution | Biology use case | Example |
|---|---|---|
| Normal | Continuous measurements | Log-transformed expression |
| Poisson | Count data (rare events) | Mutations per gene |
| Negative binomial | Overdispersed counts | RNA-seq reads (DESeq2) |
| Binomial | Success/failure in n trials | Methylated CpGs |
| Exponential | Waiting times | Distance between mutations |
| Uniform | p-values under null | Expected p-value distribution when H0 true |

## Test Selection Guide

| Situation | Parametric | Non-parametric |
|---|---|---|
| Two independent groups | Independent t-test | Mann-Whitney U |
| Two paired groups | Paired t-test | Wilcoxon signed-rank |
| 3+ independent groups | One-way ANOVA | Kruskal-Wallis |
| Categorical association | — | Chi-squared / Fisher's exact |

**Use parametric** when: ~normal distribution (or n>30), roughly equal variances, continuous data.
**Use non-parametric** when: small samples, skewed/outlier-heavy, ordinal data.

```python
from scipy import stats

# Two-group comparison
t_stat, t_p = stats.ttest_ind(treated, control)           # assumes equal var
tw_stat, tw_p = stats.ttest_ind(treated, control, equal_var=False)  # Welch's
u_stat, u_p = stats.mannwhitneyu(treated, control, alternative='two-sided')
```

## Multiple Testing Correction

Testing 20,000 genes at alpha=0.05 yields ~1,000 false positives.

### Bonferroni
- alpha_adj = alpha / m — controls family-wise error rate (FWER)
- Very conservative; misses many true positives

### Benjamini-Hochberg (FDR)
- Controls false discovery rate (expected proportion of FP among discoveries)
- Sort p-values, find largest k where p_(k) <= k/m * q
- Standard choice for genomics (q=0.05)

```python
from statsmodels.stats.multitest import multipletests

reject_bonf, pvals_bonf, _, _ = multipletests(p_values, alpha=0.05, method='bonferroni')
reject_bh, pvals_bh, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
```

## Pitfalls

- **Multiple testing**: Always apply FDR correction when testing thousands of features
- **Skewed data + t-test**: outliers distort t-test results; Mann-Whitney is robust
- **P-value histogram diagnostic**: uniform = all null; spike near 0 = true signal present
