---
name: foundations-probability
description: Probability for Bioinformatics with NumPy
tool_type: python
primary_tool: NumPy
---

# Probability for Bioinformatics

## Pitfalls

- **PDF vs PMF vs CDF**: `pdf(x)` is a density, not a probability. `pmf(k)` IS a probability. CDF `P(X ≤ x)` is always a probability.
- **`scipy.stats` uses `loc`/`scale`, not `mean`/`std`**: For exponential, `scale` = 1/rate.
- **Negative Binomial parameterization varies**: NumPy, SciPy, and R all differ — check docs when switching.
- **Independence is a strong assumption**: Gene expression levels are correlated; treating them as independent leads to incorrect inference.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.power import TTestIndPower

np.random.seed(42)
plt.rcParams['figure.figsize'] = (10, 4)
plt.rcParams['figure.dpi'] = 100
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
```

## scipy.stats Unified API

| Method | Returns | Use when |
|--------|---------|----------|
| `.pdf(x)` / `.pmf(x)` | Probability density (or mass) at x | Plotting distribution shape |
| `.cdf(x)` | P(X ≤ x) | Computing tail probabilities |
| `.sf(x)` | P(X > x) = 1 - CDF | One-tailed p-values directly |
| `.ppf(q)` | x such that P(X ≤ x) = q | Finding critical values |
| `.rvs(size)` | Random samples | Simulation and bootstrap |
| `.fit(data)` | MLE parameter estimates | Fitting to real data |

## Normal Distribution: Gene Expression

```python
mu_expr, sigma_expr = 8.5, 1.2   # log2(TPM+1) housekeeping gene
gene_dist = stats.norm(loc=mu_expr, scale=sigma_expr)
expression_values = gene_dist.rvs(size=200)

prob_above_10 = gene_dist.sf(10)       # P(expression > 10)
percentile_95  = gene_dist.ppf(0.95)   # 95th percentile
cdf_at_7       = gene_dist.cdf(7)      # P(expression <= 7)

x = np.linspace(mu_expr - 4*sigma_expr, mu_expr + 4*sigma_expr, 300)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
ax1.plot(x, gene_dist.pdf(x), 'steelblue', lw=2)
ax1.fill_between(x[x > 10], gene_dist.pdf(x[x > 10]), alpha=0.3, color='tomato', label=f'P(X>10) = {prob_above_10:.3f}')
ax1.set_title('Normal PDF — housekeeping gene'); ax1.legend()
ax2.plot(x, gene_dist.cdf(x), 'steelblue', lw=2)
ax2.axhline(0.95, ls='--', color='tomato', lw=1, label=f'95th pct = {percentile_95:.2f}')
ax2.set_title('Normal CDF — housekeeping gene'); ax2.legend()
plt.tight_layout(); plt.show()
```

## Binomial Distribution: Variant Allele Counts

Diploid heterozygous SNP at 30x coverage: alt-reads ~ Binomial(n=30, p=0.5). Somatic mutation at 20% VAF: Binomial(n=30, p=0.2).

```python
coverage = 30
germ_dist   = stats.binom(n=coverage, p=0.5)
somatic_dist = stats.binom(n=coverage, p=0.2)

k_values = np.arange(0, coverage + 1)
fig, ax = plt.subplots(figsize=(10, 4))
ax.bar(k_values - 0.2, germ_dist.pmf(k_values), width=0.4, label='Germline (p=0.5)', color='steelblue', alpha=0.8)
ax.bar(k_values + 0.2, somatic_dist.pmf(k_values), width=0.4, label='Somatic (p=0.2)', color='tomato', alpha=0.8)
ax.set_title('Binomial PMF: germline vs somatic variant allele counts'); ax.legend()
plt.tight_layout(); plt.show()

print(f"P(alt=6 | germline): {germ_dist.pmf(6):.5f}")
print(f"P(alt=6 | somatic):  {somatic_dist.pmf(6):.5f}")
print(f"Likelihood ratio:    {somatic_dist.pmf(6)/germ_dist.pmf(6):.1f}x more likely under somatic model")
```

## Poisson Distribution: Mutations per Gene

Somatic mutation count per gene ~ Poisson(λ), where λ = gene_length × background_rate.

```python
lambda_background, lambda_driver = 0.3, 3.5
bg_dist, driver_dist = stats.poisson(mu=lambda_background), stats.poisson(mu=lambda_driver)

k = np.arange(0, 12)
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(k, bg_dist.pmf(k), 'o-', color='steelblue', label=f'Background λ={lambda_background}')
ax.plot(k, driver_dist.pmf(k), 's-', color='tomato', label=f'Driver gene λ={lambda_driver}')
ax.set_title('Poisson PMF: somatic mutation counts per gene'); ax.legend()
plt.tight_layout(); plt.show()

pval_5_plus = bg_dist.sf(4)   # P(X >= 5)
print(f"P(mutations >= 5 | background): {pval_5_plus:.2e}")
```

## Negative Binomial: RNA-seq Overdispersion

RNA-seq counts are overdispersed (variance > mean). NB adds a dispersion parameter α. DESeq2/edgeR use NB models.

SciPy `nbinom(n, p)` parameterization: given mean μ and dispersion α (var = μ + αμ²):
- `p = 1 / (1 + α*μ)`
- `n = 1 / α`

```python
mu_rnaseq, alpha_disp = 50, 0.1
n_nb = 1 / alpha_disp
p_nb = 1 / (1 + alpha_disp * mu_rnaseq)
nb_dist   = stats.nbinom(n=n_nb, p=p_nb)
pois_dist = stats.poisson(mu=mu_rnaseq)

print(f"NB variance: {nb_dist.var():.1f}  Poisson variance: {pois_dist.var():.1f}")
print(f"Overdispersion: {nb_dist.var() / nb_dist.mean():.1f}x")

k_rna = np.arange(0, 130)
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(k_rna, pois_dist.pmf(k_rna), color='steelblue', lw=2, label=f'Poisson (μ={mu_rnaseq})')
ax.plot(k_rna, nb_dist.pmf(k_rna), color='tomato', lw=2, label=f'NegBinom (μ={mu_rnaseq}, α={alpha_disp})')
ax.set_title('Poisson vs Negative Binomial: RNA-seq count overdispersion'); ax.legend()
plt.tight_layout(); plt.show()
```

## Descriptive Statistics

```python
n_samples_per_tissue, n_genes = 20, 500
tissues = ['liver', 'kidney', 'brain']
tissue_offsets = {'liver': 0.0, 'kidney': 0.3, 'brain': -0.2}

expr_data = {}
for tissue in tissues:
    offset = tissue_offsets[tissue]
    expr_data[tissue] = np.random.normal(loc=8.0 + offset, scale=1.5,
                                          size=(n_samples_per_tissue, n_genes))

gene_liver  = expr_data['liver'][:, 0]
gene_kidney = expr_data['kidney'][:, 0]
gene_brain  = expr_data['brain'][:, 0]

summary_df = pd.DataFrame({'liver': gene_liver, 'kidney': gene_kidney, 'brain': gene_brain})
print(summary_df.describe().round(3))
```

## Pitfalls (Genomic)

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors.
- **Batch effects**: Always check for batch confounding before interpreting biological signal.
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously.
