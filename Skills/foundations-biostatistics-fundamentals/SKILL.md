---
name: foundations-biostatistics-fundamentals
description: "Biostatistics fundamentals: descriptive statistics, distributions, hypothesis testing, and confidence intervals for biological data. Use when analyzing experimental results."
tool_type: python
source_notebook: "Tier_0_Computational_Foundations/06_Biostatistics/01_biostatistics_fundamentals.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, scipy 1.12+, statsmodels 0.14+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Biostatistics Fundamentals

*Source: Course notebook `Tier_0_Computational_Foundations/06_Biostatistics/01_biostatistics_fundamentals.ipynb`*


**Tier 0 -- Computational Foundations | Module 6**

---

## Why Statistics Matters in Bioinformatics

Every major question in modern biology is ultimately a statistical question:

- **Differential expression**: Is this gene significantly up-regulated in tumor vs. normal tissue?
- **Variant calling**: Is this SNP a real mutation or a sequencing error?
- **Enrichment analysis**: Are DNA repair genes overrepresented in our hit list?
- **Clinical trials**: Does this drug improve patient survival?
- **Quality control**: Is this sample an outlier that should be removed?

Without solid statistical reasoning, you cannot distinguish signal from noise. This module provides the statistical foundations you will use throughout the course.

### Learning Objectives

By the end of this module you will be able to:
- Distinguish populations from samples and choose appropriate statistics
- Compute and interpret descriptive statistics for biological data
- Explain and visualize common probability distributions
- Apply the Central Limit Theorem to understand sampling distributions
- Formulate null and alternative hypotheses and interpret p-values correctly
- Choose the right statistical test (parametric vs. non-parametric)
- Apply multiple testing correction (Bonferroni, FDR) -- essential for genomics
- Perform power analysis to plan experiments

## How to use this notebook
1. Run imports first (cell below the section headers), then run each section top-to-bottom.
2. Focus on the visualizations — the histograms, boxplots, and power curves convey the concepts better than the formulas.
3. The final worked example (Section 11) ties everything together in a simulated RNA-seq analysis — read it carefully before the exercises.
4. For each statistical test, note which assumptions it requires and what happens when those assumptions are violated.

## Common stumbling points

- **P-value is not the probability that the null hypothesis is true**: It is the probability of observing data this extreme *if* the null hypothesis were true. These are completely different statements.
- **Statistical significance ≠ biological significance**: With n = 10,000 samples, even a 0.001-fold change will be "significant." Always look at effect sizes (log2 fold change, Cohen's d) alongside p-values.
- **Multiple testing explosion**: Testing 20,000 genes at α = 0.05 expects 1,000 false positives even with no real effects. Always apply BH/FDR correction in genomics analyses.
- **Paired vs. unpaired tests**: Using an unpaired test on paired data (before/after on same patients) throws away information and loses statistical power. Always match the test to the study design.
- **The Shapiro-Wilk test is not the right way to choose between t-test and Mann-Whitney**: With large n, Shapiro-Wilk rejects normality for trivially small deviations. With small n, it has low power to detect real non-normality. Look at the data visually and use domain knowledge.

```python
# Core imports for this module
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.power import TTestIndPower

# Reproducibility
np.random.seed(42)

# Plot styling
plt.rcParams['figure.figsize'] = (10, 4)
plt.rcParams['figure.dpi'] = 100
```python

---

## 1. Populations vs. Samples

| Concept | Population | Sample |
|---------|-----------|--------|
| Definition | All possible observations | The subset we actually measure |
| Size | Usually infinite or impractical to measure | Finite and known |
| Goal | Has true **parameters** | Provides **estimates** (statistics) |
| Notation (mean) | mu (population mean) | x-bar (sample mean) |
| Notation (SD) | sigma (population SD) | s (sample SD) |

### Biological examples

| Population | Sample |
|-----------|--------|
| All human genomes that ever existed or will exist | 1000 Genomes Project (~2,500 individuals) |
| All cancer cells in a patient's body | Cells from a needle biopsy |
| All possible measurements of a gene's expression | Your 3 biological replicates |

**Key insight**: We never know the true population parameters. Statistics lets us make rigorous inferences about the population from our sample.

---

## 2. Types of Data

The type of data determines which statistical tests are appropriate.

### Categorical Data
- **Nominal**: No natural ordering. Examples: blood type (A, B, AB, O), mutation type (missense, nonsense, frameshift), tissue type
- **Ordinal**: Has a natural order, but intervals are not equal. Examples: cancer stage (I, II, III, IV), quality score (low/medium/high)

### Numerical Data
- **Discrete**: Countable integers. Examples: read counts, number of mutations, number of exons
- **Continuous**: Any value in a range. Examples: gene expression (TPM), GC content (%), protein concentration

---

## 3. Descriptive Statistics

Before running any test, always look at your data. Descriptive statistics summarize the center and spread of a distribution.

```python
# Simulated gene expression data with outliers (realistic for biology)
np.random.seed(42)
expression = np.random.normal(loc=10, scale=2, size=100)

# Add a few outliers (e.g., batch effects or highly expressed genes)
expression = np.append(expression, [20, 22, 25])

# Measures of central tendency
print("=== Measures of Central Tendency ===")
print(f"Mean:   {np.mean(expression):.2f}  (sensitive to outliers)")
print(f"Median: {np.median(expression):.2f}  (robust to outliers)")
print(f"Mode:   {stats.mode(np.round(expression), keepdims=True).mode[0]:.1f}  (most frequent value)")

# Measures of spread
print("\n=== Measures of Spread ===")
print(f"Variance (s^2): {np.var(expression, ddof=1):.2f}")
print(f"Std Dev (s):    {np.std(expression, ddof=1):.2f}")
print(f"Range:          {np.ptp(expression):.2f}  (max - min)")

# IQR -- Interquartile Range (robust to outliers)
q25, q75 = np.percentile(expression, [25, 75])
iqr = q75 - q25
print(f"IQR:            {iqr:.2f}  (Q3 - Q1)")
print(f"Q1 (25th pctl): {q25:.2f}")
print(f"Q3 (75th pctl): {q75:.2f}")

# Standard error of the mean
sem = np.std(expression, ddof=1) / np.sqrt(len(expression))
print(f"\nStandard Error of Mean: {sem:.2f}  (SD / sqrt(n))")
print(f"  This tells us how precisely we know the mean.")
print(f"  n = {len(expression)}, so SEM << SD")
```python

```python
# Why median > mean when outliers are high
print("Effect of outliers:")
print(f"  Without outliers: mean = {np.mean(expression[:100]):.2f}, median = {np.median(expression[:100]):.2f}")
print(f"  With 3 outliers:  mean = {np.mean(expression):.2f}, median = {np.median(expression):.2f}")
print(f"  The mean shifted by {np.mean(expression) - np.mean(expression[:100]):.2f}, "
      f"but median only by {np.median(expression) - np.median(expression[:100]):.2f}")
```python

```python
# Visualizing distributions
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Histogram
axes[0].hist(expression, bins=25, edgecolor='black', alpha=0.7, color='steelblue')
axes[0].axvline(np.mean(expression), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(expression):.1f}')
axes[0].axvline(np.median(expression), color='green', linestyle='--', linewidth=2, label=f'Median: {np.median(expression):.1f}')
axes[0].set_xlabel('Expression Level')
axes[0].set_ylabel('Frequency')
axes[0].set_title('Distribution with Outliers')
axes[0].legend(fontsize=9)

# Boxplot
bp = axes[1].boxplot(expression, vert=True, patch_artist=True,
                     boxprops=dict(facecolor='lightblue'))
axes[1].set_ylabel('Expression Level')
axes[1].set_title('Boxplot (outliers shown as circles)')

# Q-Q plot (for normality assessment)
stats.probplot(expression, dist='norm', plot=axes[2])
axes[2].set_title('Q-Q Plot (deviation at tails = outliers)')

plt.tight_layout()
plt.show()
```python

### Understanding the Boxplot

```python
          o           <- Outliers (beyond 1.5 * IQR from box edge)
          |
    ------+------     <- Upper whisker (last data point within 1.5 * IQR)
    |            |
    |   Q3 (75%) |    <- Top of box
    |            |
    |----Median---|
    |            |
    |   Q1 (25%) |    <- Bottom of box
    |            |
    ------+------     <- Lower whisker
```python

The box contains 50% of the data (the IQR). Points beyond 1.5 * IQR from the box edges are flagged as potential outliers.

---

## 4. Probability Distributions

### 4.1 The Normal (Gaussian) Distribution

The most important distribution in statistics, defined by two parameters: mean (mu) and standard deviation (sigma).

**The empirical rule (68-95-99.7)**:
- ~68% of data falls within 1 SD of the mean
- ~95% within 2 SD
- ~99.7% within 3 SD

```python
# Visualize the normal distribution and the empirical rule
x = np.linspace(-4, 4, 1000)
y = stats.norm.pdf(x)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(x, y, 'k-', linewidth=2)

# Fill regions
ax.fill_between(x, y, where=(x >= -1) & (x <= 1), alpha=0.3, color='blue', label='68.3% (1 SD)')
ax.fill_between(x, y, where=((x >= -2) & (x < -1)) | ((x > 1) & (x <= 2)), alpha=0.3, color='green', label='+ 27.2% (2 SD)')
ax.fill_between(x, y, where=((x >= -3) & (x < -2)) | ((x > 2) & (x <= 3)), alpha=0.3, color='orange', label='+ 4.3% (3 SD)')

# Mark standard deviations
for sd_val in [-3, -2, -1, 0, 1, 2, 3]:
    label = f'{sd_val} SD' if sd_val != 0 else 'mean'
    ax.axvline(sd_val, color='gray', linestyle=':', alpha=0.5)
    ax.text(sd_val, -0.02, label, ha='center', fontsize=8)

ax.set_xlabel('Standard Deviations from Mean')
ax.set_ylabel('Probability Density')
ax.set_title('Normal Distribution and the Empirical Rule')
ax.legend()
ax.set_xlim(-4, 4)
plt.tight_layout()
plt.show()
```python

### 4.2 Distributions Used in Bioinformatics

```python
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
np.random.seed(42)

# 1. Normal -- log-transformed expression data
data_norm = np.random.normal(10, 2, 2000)
axes[0, 0].hist(data_norm, bins=40, density=True, alpha=0.7, color='steelblue')
axes[0, 0].set_title('Normal\n(log-expression values)')

# 2. Poisson -- read counts (simple model)
data_pois = np.random.poisson(lam=20, size=2000)
axes[0, 1].hist(data_pois, bins=range(0, 45), density=True, alpha=0.7, color='coral')
axes[0, 1].set_title('Poisson (lambda=20)\n(simple read count model)')

# 3. Negative Binomial -- RNA-seq counts (accounts for overdispersion)
data_nb = np.random.negative_binomial(n=5, p=5/(5+20), size=2000)
axes[0, 2].hist(data_nb, bins=40, density=True, alpha=0.7, color='mediumpurple')
axes[0, 2].set_title('Negative Binomial\n(RNA-seq counts, DESeq2 model)')

# 4. Binomial -- variant allele frequency
data_binom = np.random.binomial(n=100, p=0.5, size=2000)
axes[1, 0].hist(data_binom, bins=30, density=True, alpha=0.7, color='darkseagreen')
axes[1, 0].set_title('Binomial (n=100, p=0.5)\n(variant allele counts)')

# 5. Uniform -- p-values under the null
data_unif = np.random.uniform(0, 1, 2000)
axes[1, 1].hist(data_unif, bins=20, density=True, alpha=0.7, color='goldenrod')
axes[1, 1].set_title('Uniform(0,1)\n(p-values when H0 is true)')

# 6. Chi-squared -- goodness of fit
data_chi2 = np.random.chisquare(df=5, size=2000)
axes[1, 2].hist(data_chi2, bins=40, density=True, alpha=0.7, color='tomato')
axes[1, 2].set_title('Chi-squared (df=5)\n(goodness-of-fit test)')

for ax in axes.flat:
    ax.set_ylabel('Density')

plt.suptitle('Distributions in Bioinformatics', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.show()
```python

**Why Negative Binomial for RNA-seq?**

The Poisson distribution assumes variance = mean. In real RNA-seq data, the variance is typically *larger* than the mean (overdispersion) due to biological variability between replicates. The Negative Binomial distribution adds a dispersion parameter to model this extra variability. This is why DESeq2 and edgeR use NB, not Poisson.

---

## 5. The Central Limit Theorem (CLT)

One of the most important theorems in statistics:

> Regardless of the original distribution, the **distribution of sample means** approaches a normal distribution as the sample size increases.

This is why t-tests and ANOVA work even when individual observations are not perfectly normal -- they test the *mean*, and the CLT guarantees that sample means are approximately normal for moderate n.

The standard error of the mean (SEM) quantifies how precisely we know the sample mean:

**SEM = SD / sqrt(n)**

More samples = smaller SEM = more precise estimate of the true mean.

## Common Pitfalls

- **Assumption violations**: Check normality and homoscedasticity before parametric tests; use nonparametric alternatives when assumptions fail
- **Multiple comparisons**: Bonferroni is conservative; use Benjamini-Hochberg FDR for large-scale testing
- **Correlation ≠ causation**: Statistical association does not imply biological mechanism
