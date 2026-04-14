---
name: foundations-statistics-python
description: Statistics with Python
tool_type: python
primary_tool: Python
---

# Statistics with Python

## Pitfalls

- **`ttest_ind` assumes equal or unequal variance**: By default, scipy uses Welch's t-test (unequal variance), which is appropriate for most bioinformatics comparisons. If you have reason to believe variances are equal, pass `equal_var=True`.
- **`mannwhitneyu` requires the `alternative` argument**: Omitting it uses a two-sided test by default in newer scipy, but older versions behave differently. Always specify `alternative='two-sided'` or `'greater'` explicitly.
- **Multiple testing correction input is an array**: `multipletests` from statsmodels takes an array of p-values and returns corrected p-values and boolean reject arrays. The order of returned values is `(reject, pvals_corrected, alphacSidak, alphacBonf)`.
- **`pearsonr` returns (r, p) not just r**: A common mistake is to write `r = stats.pearsonr(x, y)` and then try to use `r` as a number. Unpack as `r, p = stats.pearsonr(x, y)`.


## The scipy.stats API

Before applying statistical tests, you need to know how to work with probability distributions in Python. `scipy.stats` provides every distribution through a consistent interface — once you learn it for one, you know it for all.

`scipy.stats` provides every common distribution through a consistent interface. Once you learn it for one distribution, you know it for all of them:

| Method | What it returns | When to use it |
|--------|----------------|----------------|
| `.pdf(x)` / `.pmf(x)` | Probability density (or mass) at x | Plotting the distribution shape |
| `.cdf(x)` | P(X ≤ x) | Computing tail probabilities |
| `.sf(x)` | P(X > x) = 1 - CDF | One-tailed p-values directly |
| `.ppf(q)` | Quantile: x such that P(X ≤ x) = q | Finding critical values |
| `.rvs(size)` | Random samples from the distribution | Simulation and bootstrap |
| `.fit(data)` | MLE parameter estimates | Fitting to real data |

### Normal Distribution: Gene Expression

Log-transformed gene expression values are approximately normally distributed. Here we model a housekeeping gene measured across 200 samples.

```python
# Normal distribution: log2 expression of a housekeeping gene
# Parameters chosen to match realistic RNA-seq log2(TPM+1) values
mu_expr = 8.5      # mean log2 expression
sigma_expr = 1.2   # standard deviation across samples

gene_dist = stats.norm(loc=mu_expr, scale=sigma_expr)

# Generate 200 synthetic expression measurements
expression_values = gene_dist.rvs(size=200)

# Key probabilities
prob_above_10 = gene_dist.sf(10)          # P(expression > 10)
percentile_95  = gene_dist.ppf(0.95)      # 95th percentile
cdf_at_7       = gene_dist.cdf(7)         # P(expression <= 7)

print(f"Mean expression:         {mu_expr}")
print(f"P(expression > 10):      {prob_above_10:.4f}")
print(f"95th percentile:         {percentile_95:.3f}")
print(f"P(expression <= 7):      {cdf_at_7:.4f}")
```

```python
# Plot PDF and CDF side by side
x = np.linspace(mu_expr - 4*sigma_expr, mu_expr + 4*sigma_expr, 300)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# PDF
ax1.plot(x, gene_dist.pdf(x), 'steelblue', lw=2, label='PDF')
ax1.fill_between(x[x > 10], gene_dist.pdf(x[x > 10]), alpha=0.3, color='tomato', label=f'P(X>10) = {prob_above_10:.3f}')
ax1.axvline(mu_expr, ls='--', color='gray', lw=1)
ax1.set_xlabel('log2 expression')
ax1.set_ylabel('Density')
ax1.set_title('Normal PDF — housekeeping gene')
ax1.legend()

# CDF
ax2.plot(x, gene_dist.cdf(x), 'steelblue', lw=2, label='CDF')
ax2.axhline(0.95, ls='--', color='tomato', lw=1, label=f'95th pct = {percentile_95:.2f}')
ax2.axvline(percentile_95, ls='--', color='tomato', lw=1)
ax2.set_xlabel('log2 expression')
ax2.set_ylabel('P(X ≤ x)')
ax2.set_title('Normal CDF — housekeeping gene')
ax2.legend()

plt.tight_layout()
plt.show()
```

### Binomial Distribution: Variant Allele Counts

When sequencing a diploid organism at a heterozygous SNP position, each read independently has a 50% chance of carrying the alternate allele. With 30x coverage, the number of alt-supporting reads follows Binomial(n=30, p=0.5). A somatic mutation in a tumor with 20% variant allele frequency follows Binomial(n=30, p=0.2).

```python
# Binomial: variant allele read counts at 30x coverage
coverage = 30
germline_vaf = 0.5   # expected heterozygous variant allele frequency
somatic_vaf  = 0.2   # somatic mutation in tumor (20% of cells)

germ_dist   = stats.binom(n=coverage, p=germline_vaf)
somatic_dist = stats.binom(n=coverage, p=somatic_vaf)

k_values = np.arange(0, coverage + 1)

fig, ax = plt.subplots(figsize=(10, 4))
ax.bar(k_values - 0.2, germ_dist.pmf(k_values), width=0.4,
       label=f'Germline (p={germline_vaf})', color='steelblue', alpha=0.8)
ax.bar(k_values + 0.2, somatic_dist.pmf(k_values), width=0.4,
       label=f'Somatic (p={somatic_vaf})', color='tomato', alpha=0.8)
ax.set_xlabel('Alt-supporting reads out of 30')
ax.set_ylabel('Probability')
ax.set_title('Binomial PMF: germline vs somatic variant allele counts')
ax.legend()
plt.tight_layout()
plt.show()

# Probability of observing exactly 6 alt reads for each model
print(f"P(alt=6 | germline model): {germ_dist.pmf(6):.5f}")
print(f"P(alt=6 | somatic model):  {somatic_dist.pmf(6):.5f}")
print(f"Likelihood ratio:          {somatic_dist.pmf(6)/germ_dist.pmf(6):.1f}x more likely under somatic model")
```

### Poisson Distribution: Mutations per Gene

The number of somatic mutations observed in a gene across cancer samples approximates a Poisson distribution when mutations are rare and independent. The expected count (lambda) depends on gene length and the background mutation rate.

```python
# Poisson: somatic mutation count per gene across 100 tumor samples
# Background rate: 1 mutation per 1 Mb, gene  3 kb  lambda  0.3
# A driver gene will have higher observed counts
lambda_background = 0.3
lambda_driver     = 3.5   # 10x enrichment — likely driver gene

bg_dist     = stats.poisson(mu=lambda_background)
driver_dist = stats.poisson(mu=lambda_driver)

k = np.arange(0, 12)

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(k, bg_dist.pmf(k), 'o-', color='steelblue', label=f'Background λ={lambda_background}')
ax.plot(k, driver_dist.pmf(k), 's-', color='tomato', label=f'Driver gene λ={lambda_driver}')
ax.set_xlabel('Mutation count')
ax.set_ylabel('Probability')
ax.set_title('Poisson PMF: somatic mutation counts per gene')
ax.legend()
plt.tight_layout()
plt.show()

# P-value: observing 5+ mutations under background model
pval_5_plus = bg_dist.sf(4)   # P(X > 4) = P(X >= 5)
print(f"P(mutations >= 5 | background): {pval_5_plus:.6f}")
print(f"This gene would be flagged as significantly mutated (p = {pval_5_plus:.2e})")
```

### Negative Binomial: RNA-seq Overdispersion

Raw RNA-seq read counts are *overdispersed* — their variance exceeds the mean. The Poisson model (variance = mean) underestimates this spread. The **Negative Binomial** adds a dispersion parameter that captures biological variability between replicates. DESeq2 and edgeR both use NB models.

scipy.stats uses the `nbinom(n, p)` parameterisation. To convert from the bioinformatics convention (mean μ, dispersion α where var = μ + αμ²):
- `p = 1 / (1 + α*μ)`
- `n = 1 / α`

```python
# Negative Binomial: RNA-seq count distribution
# Typical low-expressed gene: mean50, dispersion0.1 (variance  50 + 0.1*502  300)
mu_rnaseq   = 50
alpha_disp  = 0.1   # dispersion parameter

n_nb = 1 / alpha_disp
p_nb = 1 / (1 + alpha_disp * mu_rnaseq)

nb_dist     = stats.nbinom(n=n_nb, p=p_nb)
pois_dist   = stats.poisson(mu=mu_rnaseq)   # for comparison

print(f"NB mean:     {nb_dist.mean():.1f}")
print(f"NB variance: {nb_dist.var():.1f}")
print(f"Poisson var (=mean): {pois_dist.var():.1f}")
print(f"Overdispersion factor: {nb_dist.var() / nb_dist.mean():.1f}x")

k_rna = np.arange(0, 130)
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(k_rna, pois_dist.pmf(k_rna), color='steelblue', lw=2, label=f'Poisson (μ={mu_rnaseq})')
ax.plot(k_rna, nb_dist.pmf(k_rna), color='tomato', lw=2, label=f'NegBinom (μ={mu_rnaseq}, α={alpha_disp})')
ax.set_xlabel('Read count')
ax.set_ylabel('Probability')
ax.set_title('Poisson vs Negative Binomial: RNA-seq count overdispersion')
ax.legend()
plt.tight_layout()
plt.show()
```

### Visualising Multiple Distributions Together

A side-by-side comparison of the four distributions we use most often in bioinformatics, each with a biological interpretation.

```python
# Four-panel distribution gallery
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# 1. Normal -- gene expression
x_norm = np.linspace(4, 13, 300)
axes[0, 0].plot(x_norm, stats.norm(8.5, 1.2).pdf(x_norm), 'steelblue', lw=2)
axes[0, 0].set_title('Normal: log2 gene expression')
axes[0, 0].set_xlabel('log2(TPM+1)')

# 2. Binomial -- variant allele reads
k_bin = np.arange(0, 31)
axes[0, 1].bar(k_bin, stats.binom(30, 0.5).pmf(k_bin), color='steelblue', alpha=0.8)
axes[0, 1].set_title('Binomial: alt reads (n=30, p=0.5)')
axes[0, 1].set_xlabel('Alt-supporting reads')

# 3. Poisson -- somatic mutations
k_poi = np.arange(0, 12)
axes[1, 0].bar(k_poi, stats.poisson(2.0).pmf(k_poi), color='steelblue', alpha=0.8)
axes[1, 0].set_title('Poisson: mutation count (λ=2)')
axes[1, 0].set_xlabel('Mutations per gene')

# 4. Negative Binomial -- RNA-seq counts
k_nb = np.arange(0, 130)
axes[1, 1].plot(k_nb, stats.nbinom(10, 10/(10+50)).pmf(k_nb), color='steelblue', lw=2)
axes[1, 1].set_title('NegBinom: RNA-seq counts (μ=50, α=0.1)')
axes[1, 1].set_xlabel('Read count')

for ax in axes.flat:
    ax.set_ylabel('Probability')

plt.suptitle('Probability distributions in bioinformatics', fontsize=13, y=1.01)
plt.tight_layout()
plt.show()
```


## Descriptive Statistics and Visualization

Before running any test, explore your data visually. The goal is to understand the distribution shape, spot outliers, and decide whether parametric or non-parametric methods are appropriate.

We simulate a gene expression dataset: 60 samples across three tissue types (liver, kidney, brain), with 500 genes each.

```python
# Generate synthetic expression dataset: 60 samples x 500 genes
n_samples_per_tissue = 20
n_genes = 500
tissues = ['liver', 'kidney', 'brain']

# Each tissue has a slightly different mean expression profile
tissue_offsets = {'liver': 0.0, 'kidney': 0.3, 'brain': -0.2}

expr_data = {}
for tissue in tissues:
    offset = tissue_offsets[tissue]
    # Base expression drawn from N(8, 1.5), tissue-specific shift added
    expr_data[tissue] = np.random.normal(loc=8.0 + offset, scale=1.5,
                                          size=(n_samples_per_tissue, n_genes))

# Focus on one gene (index 0) for univariate analysis
gene_liver  = expr_data['liver'][:, 0]
gene_kidney = expr_data['kidney'][:, 0]
gene_brain  = expr_data['brain'][:, 0]

# Summary statistics with pandas
summary_df = pd.DataFrame({'liver': gene_liver, 'kidney': gene_kidney, 'brain': gene_brain})
print("Summary statistics for gene_0 across tissues:")
print(summary_df.describe().round(3))
```

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
