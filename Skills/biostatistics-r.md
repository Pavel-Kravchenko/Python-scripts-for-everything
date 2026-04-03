---
name: biostatistics-r
description: R programming fundamentals, statistical distributions, hypothesis testing, and multiple testing correction for biological data
---

# Biostatistics & R for Bioinformatics

## When to Use
- Differential expression analysis (RNA-seq, microarray)
- Variant/mutation association testing
- Comparing gene expression across sample groups
- Planning experiments (power analysis, sample size)
- Reporting results with effect sizes + corrected p-values

---

## Quick Reference

### R vs Python Gotchas
| Feature | Python | R |
|---|---|---|
| Index start | 0 | **1** |
| `x[-1]` | last element | **all except first** |
| NA/None | `None` | `NA` |
| Boolean | `True/False` | `TRUE/FALSE` |
| Assignment | `=` | `<-` (preferred) |

### R Distribution Functions (prefix convention)
| Prefix | Purpose | Example |
|---|---|---|
| `d` | PDF value | `dnorm(x, mean, sd)` |
| `p` | CDF / cumulative prob | `pnorm(q, mean, sd)` |
| `q` | Quantile (inverse CDF) | `qnorm(p, mean, sd)` |
| `r` | Random sample | `rnorm(n, mean, sd)` |
Distributions: `norm`, `t`, `binom`, `pois`, `nbinom`, `chisq`, `f`, `unif`

### Test Selection Decision Tree
```
Data type?
├── Continuous
│   ├── 2 groups
│   │   ├── Normal + equal variance  → t-test (ttest_ind)
│   │   ├── Normal + unequal var     → Welch's t-test (equal_var=False)
│   │   ├── Paired samples           → paired t-test (ttest_rel)
│   │   └── Non-normal or small n    → Mann-Whitney U (mannwhitneyu)
│   └── 3+ groups
│       ├── Normal                   → ANOVA (f_oneway)
│       └── Non-normal               → Kruskal-Wallis (kruskal)
├── Categorical (contingency table)
│   ├── Small expected counts (<5)   → Fisher's exact (fisher_exact)
│   └── Large samples                → Chi-squared (chi2_contingency)
└── Correlation
    ├── Linear, normal               → Pearson (pearsonr)
    └── Monotonic or non-normal      → Spearman (spearmanr)
```

### Multiple Testing Correction Guide
| Context | Method | Threshold |
|---|---|---|
| RNA-seq DE (exploratory) | BH/FDR | FDR < 0.10 |
| RNA-seq DE (standard) | BH/FDR | FDR < 0.05 |
| GWAS | Bonferroni | p < 5×10⁻⁸ |
| Drug target (few tests, high stakes) | Bonferroni | α/n |
| Gene set enrichment | BH/FDR | FDR < 0.05 |

### Effect Size Interpretation (Cohen's d)
| d | Interpretation |
|---|---|
| 0.2 | Small |
| 0.5 | Medium |
| 0.8 | Large |

---

## Key Patterns

### Distributions Matched to Biology
- **Normal**: log-transformed expression values, continuous measurements
- **Poisson**: read counts (simple model), mutations per gene (variance = mean)
- **Negative Binomial**: RNA-seq counts with overdispersion — variance > mean (DESeq2/edgeR use this)
- **Binomial**: variant calling, methylated CpGs out of n trials
- **Exponential**: waiting times, distances between features

### Why NB over Poisson for RNA-seq
Poisson assumes variance = mean. Biological replicate variability causes overdispersion (variance > mean). NB adds a dispersion parameter: `rnbinom(n, size=dispersion, mu=mean_count)`.

### Multiple Testing Problem
Testing 20,000 genes at α=0.05 → ~1,000 false positives from null genes alone.
- **Bonferroni** controls FWER: zero tolerance for any false positive. Very conservative.
- **BH/FDR** controls expected proportion of false positives among discoveries. Preferred for genomics.

### CLT in Practice
t-tests and ANOVA test the *mean*. By CLT, sample means are approximately normal for n ≥ 30, even when individual observations are non-normal. With small n and skewed data, use non-parametric alternatives.

---

## Code Templates

### R: Vectors, Data Frames, Filtering
```r
# Vectors
gene_expr <- c(2.5, 3.1, 8.2, 1.8, 5.6)
log2(gene_expr + 1)                        # pseudocount for log(0) safety

# Named vector
gc <- c(BRCA1=0.42, TP53=0.38, EGFR=0.55)
gc[gc > 0.5]                               # logical subsetting

# Data frame
df <- data.frame(
    gene    = c("BRCA1", "TP53", "EGFR"),
    log2fc  = c(1.2, -0.5, 3.8),
    pvalue  = c(0.01, 0.15, 0.001),
    stringsAsFactors = FALSE
)
df[df$log2fc > 1 & df$pvalue < 0.05, ]    # filter DE genes

# Apply corrections
df$padj_bonf <- p.adjust(df$pvalue, method = "bonferroni")
df$fdr       <- p.adjust(df$pvalue, method = "BH")
df[order(df$pvalue), ]                     # sort by p-value
```

### R: Matrix Operations (expression data)
```r
expr <- matrix(counts, nrow=n_genes, ncol=n_samples,
               dimnames=list(gene_names, sample_names))
apply(expr, 1, mean)   # row-wise (per gene)
apply(expr, 2, mean)   # col-wise (per sample)
```

### R: Distributions
```r
set.seed(42)
rnorm(100, mean=10, sd=2)          # normal expression values
rpois(100, lambda=20)              # Poisson read counts
rnbinom(100, size=5, mu=20)        # NB counts (overdispersed)
rbinom(100, size=30, prob=0.5)     # variant allele counts

pnorm(8, mean=10, sd=2)            # P(X < 8)
qnorm(0.95, mean=10, sd=2)         # 95th percentile
```

### Python: Hypothesis Tests
```python
from scipy import stats

# Two-group comparisons
t_stat, p = stats.ttest_ind(tumor, normal)                    # equal variance
t_stat, p = stats.ttest_ind(tumor, normal, equal_var=False)   # Welch's
t_stat, p = stats.ttest_rel(after, before)                    # paired
u_stat, p = stats.mannwhitneyu(a, b, alternative='two-sided') # non-parametric

# Multi-group
f_stat, p = stats.f_oneway(group1, group2, group3)            # ANOVA
h_stat, p = stats.kruskal(group1, group2, group3)             # Kruskal-Wallis

# Categorical
odds_r, p = stats.fisher_exact(contingency_2x2)
chi2, p, dof, expected = stats.chi2_contingency(table)

# Correlation
r, p = stats.pearsonr(x, y)
rho, p = stats.spearmanr(x, y)

# Normality check (use before choosing parametric/non-parametric)
stat, p = stats.shapiro(data)   # p > 0.05 → normal
```

### Python: Multiple Testing Correction
```python
from statsmodels.stats.multitest import multipletests

reject, pvals_adj, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
reject, pvals_adj, _, _ = multipletests(p_values, alpha=0.05, method='bonferroni')

# Volcano plot threshold: FDR < 0.05 AND |log2FC| > 1
is_sig = (pvals_adj < 0.05) & (np.abs(log2fc) > 1)
```

### Python: Effect Size and Power Analysis
```python
import numpy as np
from statsmodels.stats.power import TTestIndPower

def cohens_d(g1, g2):
    n1, n2 = len(g1), len(g2)
    sp = np.sqrt(((n1-1)*g1.std(ddof=1)**2 + (n2-1)*g2.std(ddof=1)**2) / (n1+n2-2))
    return (g1.mean() - g2.mean()) / sp

# Required n per group
pa = TTestIndPower()
n = pa.solve_power(effect_size=0.5, alpha=0.05, power=0.80)  # → ~64 per group

# Achieved power given n
power = pa.power(effect_size=0.5, nobs1=30, alpha=0.05, ratio=1.0)
```

### Complete DE Analysis Workflow (Python)
```python
# 1. t-test per gene
pvalues, log2fcs = [], []
for gene_row in expression_matrix:
    ctrl = gene_row[:n_ctrl]
    treat = gene_row[n_ctrl:]
    _, p = stats.ttest_ind(treat, ctrl)
    pvalues.append(p)
    log2fcs.append(treat.mean() - ctrl.mean())

# 2. Correct for multiple testing
_, fdr, _, _ = multipletests(pvalues, method='fdr_bh')

# 3. Significant genes: FDR < 0.05 AND |log2FC| > 1
is_sig = (np.array(fdr) < 0.05) & (np.abs(np.array(log2fcs)) > 1)
```

---

## Common Pitfalls

- **Forgetting pseudocount**: `log2(count)` fails on zeros; use `log2(count + 1)`.
- **Unpaired test on paired data**: loses statistical power; always use `ttest_rel` for matched samples.
- **Parametric test on clearly non-normal data**: use Shapiro-Wilk first; switch to Mann-Whitney/Kruskal if violated.
- **Raw p-values in genomics**: always apply BH correction; raw p < 0.05 across 20k genes yields ~1k false positives.
- **Statistical vs biological significance**: a tiny fold change (e.g., 1.01×) can be p < 0.001 with large n; always report log2FC alongside FDR.
- **Poisson for RNA-seq**: real data is overdispersed; use negative binomial (DESeq2, edgeR).
- **R negative indexing**: `x[-1]` removes first element (not the last as in Python).
- **R 1-based indexing**: `x[1]` is the first element; `x[0]` returns `numeric(0)`.

---

## Related Skills
- `rnaseq-metagenomics` — DESeq2/edgeR differential expression using NB model
- `ngs-variant-calling` — genome-wide association, variant filtering, Manhattan plots
- `ml-deep-learning-bio` — when to use ML vs statistical tests for classification/prediction
