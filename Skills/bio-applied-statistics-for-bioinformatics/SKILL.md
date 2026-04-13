---
name: bio-applied-statistics-for-bioinformatics
description: "Bioinformatics analyses generate thousands to millions of measurements. Proper statistical methods are essential to distinguish genuine biological signals from noise. This notebook covers the statisti"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/06_Statistics_for_Bioinformatics/01_statistics_for_bioinformatics.ipynb"
---

# Statistics for Bioinformatics

*Source: Course notebook `Tier_3_Applied_Bioinformatics/06_Statistics_for_Bioinformatics/01_statistics_for_bioinformatics.ipynb`*

# Statistics for Bioinformatics

**Tier 3 -- Applied Bioinformatics**

Bioinformatics analyses generate thousands to millions of measurements. Proper statistical methods are essential to distinguish genuine biological signals from noise. This notebook covers the statistical toolkit every bioinformatician needs: from hypothesis testing and multiple testing correction to survival analysis.

**Prerequisites:** Tier 2 (basic Python, NumPy, pandas)  
**Libraries:** `numpy`, `pandas`, `scipy.stats`, `statsmodels`, `matplotlib`

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

%matplotlib inline
plt.rcParams['figure.figsize'] = (12, 5)
plt.rcParams['font.size'] = 12
np.random.seed(42)
```

---
## 1. Distributions in Biology

Choosing the right statistical test starts with understanding your data's distribution.

| Distribution | Use in biology | Example |
|-------------|---------------|--------|
| **Normal** | Continuous measurements | Height, weight, log-transformed expression |
| **Poisson** | Count data (rare events) | Mutations per gene, read counts |
| **Negative binomial** | Overdispersed counts | RNA-seq read counts (DESeq2 uses this) |
| **Binomial** | Success/failure in n trials | SNPs in a region, methylated CpGs |
| **Exponential** | Waiting times | Distance between mutations |
| **Uniform** | p-values under the null | Expected distribution of p-values when H0 is true |

```python
fig, axes = plt.subplots(2, 3, figsize=(16, 9))

# Normal -- gene expression (log2)
x_norm = np.random.normal(8, 2, 5000)
axes[0, 0].hist(x_norm, bins=50, color='steelblue', edgecolor='white', density=True)
axes[0, 0].set_title('Normal\n(Log2 expression)')

# Poisson -- mutations per gene
x_pois = np.random.poisson(3, 5000)
axes[0, 1].hist(x_pois, bins=range(15), color='coral', edgecolor='white', density=True, align='left')
axes[0, 1].set_title('Poisson\n(Mutations per gene)')

# Negative binomial -- RNA-seq counts
x_nb = np.random.negative_binomial(5, 0.01, 5000)
axes[0, 2].hist(x_nb, bins=50, color='mediumseagreen', edgecolor='white', density=True)
axes[0, 2].set_title('Negative Binomial\n(RNA-seq read counts)')

# Binomial -- methylated CpGs out of 20
x_binom = np.random.binomial(20, 0.3, 5000)
axes[1, 0].hist(x_binom, bins=range(22), color='mediumpurple', edgecolor='white', density=True, align='left')
axes[1, 0].set_title('Binomial\n(Methylated CpGs / 20 sites)')

# Exponential -- inter-mutation distance
x_exp = np.random.exponential(1000, 5000)
axes[1, 1].hist(x_exp, bins=50, color='goldenrod', edgecolor='white', density=True)
axes[1, 1].set_title('Exponential\n(Distance between mutations, bp)')

# Uniform -- p-values under null
x_unif = np.random.uniform(0, 1, 5000)
axes[1, 2].hist(x_unif, bins=50, color='gray', edgecolor='white', density=True)
axes[1, 2].set_title('Uniform\n(p-values under H0)')

for ax in axes.flat:
    ax.set_ylabel('Density')

plt.suptitle('Common Distributions in Bioinformatics', fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.show()
```

---
## 2. Hypothesis Testing Review

### The framework

1. **Null hypothesis (H0)**: No effect / no difference (e.g., gene expression is the same in treated and control)
2. **Alternative hypothesis (H1)**: There is an effect
3. **Test statistic**: Summarizes the data (e.g., t-statistic, chi-squared)
4. **p-value**: Probability of observing the test statistic (or more extreme) if H0 is true
5. **Decision**: Reject H0 if p-value < significance level (typically 0.05)

### Type I and Type II errors

| | H0 is true | H0 is false |
|---|---|---|
| **Reject H0** | Type I error (false positive, rate = alpha) | Correct (power = 1 - beta) |
| **Fail to reject H0** | Correct | Type II error (false negative, rate = beta) |

```python
# Demonstration: t-test on gene expression data
np.random.seed(42)

# Simulated gene expression (log2 scale): treated vs control
control = np.random.normal(loc=8.0, scale=1.5, size=30)
treated = np.random.normal(loc=9.2, scale=1.5, size=30)  # 1.2-fold increase in log2

t_stat, p_val = stats.ttest_ind(treated, control)

fig, ax = plt.subplots(figsize=(8, 5))
ax.boxplot([control, treated], labels=['Control', 'Treated'], widths=0.5)
ax.scatter(np.ones(30) + np.random.normal(0, 0.03, 30), control, alpha=0.5, color='steelblue')
ax.scatter(2 * np.ones(30) + np.random.normal(0, 0.03, 30), treated, alpha=0.5, color='coral')
ax.set_ylabel('Log2 expression')
ax.set_title(f'Gene Expression: Control vs Treated\nt={t_stat:.2f}, p={p_val:.4f}')
plt.tight_layout()
plt.show()

print(f"Control: mean={control.mean():.2f}, std={control.std():.2f}")
print(f"Treated: mean={treated.mean():.2f}, std={treated.std():.2f}")
print(f"Welch's t-test: t={t_stat:.3f}, p={p_val:.4e}")
```

---
## 3. The Multiple Testing Problem

### Why it matters in genomics

In a typical RNA-seq experiment, you test ~20,000 genes for differential expression. If you use alpha = 0.05:

$$\text{Expected false positives} = 20{,}000 \times 0.05 = 1{,}000 \text{ genes!}$$

You would report 1,000 genes as differentially expressed purely by chance. This is why **multiple testing correction** is critical in genomics.

```python
# Demonstrate the multiple testing problem
np.random.seed(42)

n_genes = 20000
n_true_de = 500       # truly differentially expressed genes
n_null = n_genes - n_true_de

# Generate p-values
# Null genes: p-values are uniform
p_null = np.random.uniform(0, 1, n_null)
# True DE genes: p-values are small (drawn from beta distribution)
p_de = np.random.beta(0.3, 5, n_true_de)

p_values = np.concatenate([p_null, p_de])
is_true_de = np.array([False] * n_null + [True] * n_true_de)

# Without correction: how many false positives?
alpha = 0.05
sig_uncorrected = p_values < alpha
fp_uncorrected = sig_uncorrected & ~is_true_de
tp_uncorrected = sig_uncorrected & is_true_de

print(f"Total genes tested: {n_genes}")
print(f"Truly DE genes: {n_true_de}")
print(f"\nWithout correction (alpha = {alpha}):")
print(f"  Significant calls: {sig_uncorrected.sum()}")
print(f"  True positives: {tp_uncorrected.sum()}")
print(f"  FALSE POSITIVES: {fp_uncorrected.sum()}")
print(f"  False discovery rate: {fp_uncorrected.sum() / max(sig_uncorrected.sum(), 1):.3f}")
```

```python
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# P-value distribution
axes[0].hist(p_null, bins=50, alpha=0.7, color='gray', label=f'Null genes (n={n_null})', density=True)
axes[0].hist(p_de, bins=50, alpha=0.7, color='red', label=f'True DE genes (n={n_true_de})', density=True)
axes[0].axvline(0.05, color='black', linestyle='--', label='alpha = 0.05')
axes[0].set_xlabel('p-value')
axes[0].set_ylabel('Density')
axes[0].set_title('P-value Distributions')
axes[0].legend()

# Combined histogram (what you actually see)
axes[1].hist(p_values, bins=50, color='steelblue', edgecolor='white', density=True)
axes[1].axvline(0.05, color='red', linestyle='--', linewidth=2, label='alpha = 0.05')
axes[1].set_xlabel('p-value')
axes[1].set_ylabel('Density')
axes[1].set_title('Combined P-value Distribution (all 20,000 genes)\nNote the spike near 0 = true signal')
axes[1].legend()

plt.tight_layout()
plt.show()
```

---
## 4. Multiple Testing Corrections

### 4.1 Bonferroni correction

The simplest and most conservative approach: divide alpha by the number of tests.

$$\alpha_{\text{Bonferroni}} = \frac{\alpha}{m}$$

For 20,000 genes at alpha = 0.05: $\alpha_{\text{Bonf}} = 2.5 \times 10^{-6}$

**Pros:** Controls the family-wise error rate (FWER) -- the probability of even one false positive.  
**Cons:** Very conservative. Many true positives are missed.

### 4.2 Benjamini-Hochberg (FDR control)

Controls the **false discovery rate (FDR)** -- the expected proportion of false positives among all discoveries.

Algorithm:
1. Sort p-values: $p_{(1)} \le p_{(2)} \le \ldots \le p_{(m)}$
2. For each rank $k$, compute threshold: $\frac{k}{m} \times q$ (where q is the desired FDR)
3. Find the largest $k$ where $p_{(k)} \le \frac{k}{m} \times q$
4. Reject all hypotheses with rank $\le k$

```python
from statsmodels.stats.multitest import multipletests

# Apply corrections
reject_bonf, pvals_bonf, _, _ = multipletests(p_values, alpha=0.05, method='bonferroni')
reject_bh, pvals_bh, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

# Compare results
methods = {
    'No correction': sig_uncorrected,
    'Bonferroni': reject_bonf,
    'Benjamini-Hochberg (FDR)': reject_bh,
}

print(f"{'Method':<28} {'Significant':>11} {'TP':>6} {'FP':>6} {'FN':>6} {'FDR':>8} {'Sensitivity':>12}")
print("-" * 85)
for name, rejected in methods.items():
    tp = (rejected & is_true_de).sum()
    fp = (rejected & ~is_true_de).sum()
    fn = (~rejected & is_true_de).sum()
    total_sig = rejected.sum()
    fdr = fp / max(total_sig, 1)
    sensitivity = tp / n_true_de
    print(f"{name:<28} {total_sig:>11} {tp:>6} {fp:>6} {fn:>6} {fdr:>8.3f} {sensitivity:>12.3f}")
```

```python
# Visualize the BH procedure
sorted_indices = np.argsort(p_values)
sorted_p = p_values[sorted_indices]
ranks = np.arange(1, len(sorted_p) + 1)
bh_threshold = ranks / len(sorted_p) * 0.05

fig, ax = plt.subplots(figsize=(10, 6))

# Zoom in on the first 2000 genes for visibility
n_show = 2000
ax.scatter(ranks[:n_show], sorted_p[:n_show], s=2, alpha=0.5,
           c=['red' if is_true_de[sorted_indices[i]] else 'gray' for i in range(n_show)])
ax.plot(ranks[:n_show], bh_threshold[:n_show], 'b-', linewidth=2, label='BH threshold (FDR=0.05)')
ax.set_xlabel('Rank')
ax.set_ylabel('p-value')
ax.set_title('Benjamini-Hochberg Procedure\n(red = true DE genes, gray = null genes)')
ax.legend()
plt.tight_layout()
plt.show()
```

---
## 5. Parametric vs Non-parametric Tests

### When to use which?

| Situation | Parametric test | Non-parametric alternative |
|-----------|----------------|---------------------------|
| Two independent groups, continuous data | Independent t-test | Mann-Whitney U |
| Two paired groups | Paired t-test | Wilcoxon signed-rank |
| 3+ independent groups | One-way ANOVA | Kruskal-Wallis |
| Association between categorical variables | -- | Chi-squared / Fisher's exact |

**Use parametric tests when:**
- Data is approximately normally distributed (or n > 30 by CLT)
- Variances are roughly equal between groups
- Data is continuous and measured on an interval/ratio scale

**Use non-parametric tests when:**
- Small sample sizes
- Skewed distributions, outliers
- Ordinal data or ranks
- Normality assumption clearly violated

```python
# Demonstrate: when do parametric and non-parametric tests disagree?
np.random.seed(42)

# Case 1: Normal data (both tests should agree)
normal_a = np.random.normal(10, 2, 25)
normal_b = np.random.normal(12, 2, 25)

# Case 2: Skewed data with outliers (non-parametric more robust)
skewed_a = np.random.exponential(2, 25)
skewed_b = np.random.exponential(3, 25)
# Add outliers to group a
skewed_a[0] = 50
skewed_a[1] = 45

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (a, b), title in zip(axes,
                              [(normal_a, normal_b), (skewed_a, skewed_b)],
                              ['Normal data', 'Skewed data with outliers']):
    t_stat, t_p = stats.ttest_ind(a, b)
    u_stat, u_p = stats.mannwhitneyu(a, b, alternative='two-sided')
    
    bp = ax.boxplot([a, b], labels=['Group A', 'Group B'], widths=0.5)
    ax.set_title(f"{title}\nt-test p={t_p:.4f} | Mann-Whitney p={u_p:.4f}")
    ax.set_ylabel('Value')

plt.tight_layout()
plt.show()

print("With normal data, both tests give similar p-values.")
print("With skewed/outlier data, the t-test may be misleading while Mann-Whitney is robust.")
```

---
## 6. Common Statistical Tests in Detail

### 6.1 Two-group comparisons

```python
# Scenario: comparing gene expression between tumor and normal samples
np.random.seed(42)

normal_expr = np.random.normal(7.5, 1.8, 40)    # normal tissue
tumor_expr = np.random.normal(9.0, 2.0, 35)     # tumor tissue

# t-test (assumes normality, roughly equal variance)
t_stat, t_p = stats.ttest_ind(tumor_expr, normal_expr)

# Welch's t-test (does not assume equal variance)
tw_stat, tw_p = stats.ttest_ind(tumor_expr, normal_expr, equal_var=False)

# Mann-Whitney U (non-parametric)
u_stat, u_p = stats.mannwhitneyu(tumor_expr, normal_expr, alternative='two-sided')

print("=== Comparing Tumor vs Normal Expression ===")
print(f"Normal: n={len(normal_expr)}, mean={normal_expr.mean():.2f}, std={normal_expr.std():.2f}")
print(f"Tumor:  n={len(tumor_expr)}, mean={tumor_expr.mean():.2f}, std={tumor_expr.std():.2f}")
print(f"\nStudent's t-test: t={t_stat:.3f}, p={t_p:.4e}")
print(f"Welch's t-test:   t={tw_stat:.3f}, p={tw_p:.4e}")
print(f"Mann-Whitney U:   U={u_stat:.0f}, p={u_p:.4e}")
```
