---
name: gwas-population-genetics
description: GWAS study design, genotype QC, PCA for population stratification, per-SNP association testing, Manhattan and QQ plots, LD clumping, fine-mapping concepts
---

## When to Use

Use this skill when:
- Designing a genome-wide association study
- Performing genotype quality control (SNP and sample filtering)
- Detecting population stratification via PCA
- Running per-SNP association tests (logistic or linear regression)
- Generating Manhattan and QQ plots
- Applying multiple testing correction (genome-wide significance threshold)
- Interpreting GWAS results: LD blocks, clumping, fine-mapping

## Quick Reference

| Step | Tool / Method | Key Parameters |
|---|---|---|
| Genotype QC | Filter by call rate, HWE, MAF | call_rate > 0.95, HWE p > 1e-6, MAF > 0.01 |
| Sample QC | Missingness, sex check, outliers | sample missingness < 0.05 |
| Population stratification | PCA on genotype matrix | top 10 PCs as covariates |
| Association test | Logistic regression (binary) / Linear (quantitative) | adjust for PCs + age + sex |
| Multiple testing | Genome-wide threshold | p < 5×10⁻⁸ (Bonferroni for ~1M SNPs) |
| Clumping | LD-based pruning of significant hits | r² < 0.1, window 250 kb |
| Fine-mapping | Credible sets (SuSiE, FINEMAP) | 95% credible set |
| GWAS catalog | `requests` to GWAS Catalog REST API | `https://www.ebi.ac.uk/gwas/rest/api/associations` |

## Key Patterns

**Pattern 1: Simulate genotype data**
```python
import numpy as np

def simulate_genotypes(n_samples=1000, n_snps=5000, mafs=None, seed=42):
    rng = np.random.default_rng(seed)
    if mafs is None:
        mafs = rng.uniform(0.05, 0.5, n_snps)
    # Each SNP: sample from {0,1,2} with HWE frequencies
    p = mafs
    freqs = np.column_stack([(1-p)**2, 2*p*(1-p), p**2])
    G = np.array([rng.choice([0,1,2], size=n_samples, p=freqs[j]) for j in range(n_snps)]).T
    return G, mafs
```

**Pattern 2: Genotype QC**
```python
def genotype_qc(G, min_maf=0.01, min_call_rate=0.95):
    """Filter SNPs by MAF and call rate. G: (n_samples, n_snps)."""
    maf = np.minimum(G.mean(0)/2, 1 - G.mean(0)/2)
    call_rate = (G >= 0).mean(0)
    keep = (maf >= min_maf) & (call_rate >= min_call_rate)
    return G[:, keep], keep
```

**Pattern 3: PCA for stratification**
```python
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def pca_stratification(G, n_components=10):
    G_std = StandardScaler().fit_transform(G)
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(G_std)
    return pcs, pca.explained_variance_ratio_
```

**Pattern 4: Per-SNP logistic regression**
```python
from sklearn.linear_model import LogisticRegression
from scipy.stats import chi2
import numpy as np

def gwas_logistic(G, phenotype, covariates=None):
    """Run logistic regression per SNP. Returns p-values."""
    n, p = G.shape
    pvals = np.ones(p)
    X_cov = covariates if covariates is not None else np.ones((n, 1))
    for j in range(p):
        X = np.column_stack([G[:, j], X_cov])
        try:
            lr = LogisticRegression(max_iter=200, solver="lbfgs")
            lr.fit(X, phenotype)
            # Wald test approximation via likelihood ratio
            lr0 = LogisticRegression(max_iter=200, solver="lbfgs")
            lr0.fit(X_cov, phenotype)
            ll1 = -lr.loss(X, phenotype) * n   # approx
            ll0 = -lr0.loss(X_cov, phenotype) * n
            lrt = 2 * (ll1 - ll0)
            pvals[j] = chi2.sf(max(lrt, 0), df=1)
        except Exception:
            pass
    return pvals
```

**Pattern 5: Manhattan plot**
```python
import matplotlib.pyplot as plt
import pandas as pd

def manhattan_plot(df, chrom_col="chrom", pos_col="pos", pval_col="pval"):
    df = df.copy()
    df["-log10p"] = -np.log10(df[pval_col].clip(1e-300))
    chroms = sorted(df[chrom_col].unique(), key=lambda x: int(x.replace("chr","")))
    offset = 0; xticks = []; colors = ["#1f77b4","#ff7f0e"]
    fig, ax = plt.subplots(figsize=(14,4))
    for i, c in enumerate(chroms):
        sub = df[df[chrom_col]==c]
        ax.scatter(sub[pos_col] + offset, sub["-log10p"],
                   c=colors[i%2], s=2, alpha=0.7)
        xticks.append(offset + sub[pos_col].mean())
        offset += sub[pos_col].max() + 5e6
    ax.axhline(-np.log10(5e-8), color="red", lw=1, ls="--", label="5×10⁻⁸")
    ax.set_xticks(xticks); ax.set_xticklabels([c.replace("chr","") for c in chroms], fontsize=7)
    ax.set_xlabel("Chromosome"); ax.set_ylabel("-log₁₀(p)")
    ax.set_title("Manhattan plot"); ax.legend(frameon=False)
    return fig, ax
```

## Code Templates

**Template 1: Complete GWAS pipeline**
```python
# 1. Simulate
G, mafs = simulate_genotypes(n_samples=1000, n_snps=5000)
# 2. Simulate phenotype with 3 causal SNPs
causal = [100, 500, 1500]
beta = [0.5, -0.4, 0.3]
logit = sum(b * G[:, c] for b, c in zip(beta, causal))
prob = 1 / (1 + np.exp(-logit))
phenotype = (np.random.default_rng(0).random(1000) < prob).astype(int)
# 3. QC
G_qc, keep = genotype_qc(G)
# 4. PCA
pcs, _ = pca_stratification(G_qc)
# 5. GWAS (demo on subset)
pvals = gwas_logistic(G_qc[:, :200], phenotype, covariates=pcs[:, :5])
```

## Common Pitfalls

- **Population stratification ignored:** always include PC covariates; use QQ plot λ (genomic inflation) to diagnose
- **Multiple testing threshold:** 5×10⁻⁸ is not always correct; depends on number of independent tests
- **HWE filtering in cases only:** HWE test should be run on controls only (causal SNPs deviate from HWE in cases)
- **MAF threshold:** rare variants need larger sample sizes; standard GWAS focuses on MAF > 0.01–0.05
- **LD clumping vs pruning:** clumping retains the most significant SNP; pruning retains a random one — always clump for GWAS results

## Related Skills

- `population-genetics-evolution` — drift, selection, Fst, dN/dS
- `ngs-variant-calling` — upstream VCF generation pipeline
- `bayesian-python` — Bayesian fine-mapping (SuSiE framework)
