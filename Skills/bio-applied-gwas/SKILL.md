---
name: bio-applied-gwas
description: Genome-Wide Association Studies (GWAS) with NumPy
tool_type: python
primary_tool: NumPy
---

# Genome-Wide Association Studies (GWAS)

## Study Design

GWAS tests ~1–10M SNPs for association with a phenotype by comparing allele frequencies (binary traits) or correlating allele dosage with trait values (quantitative traits).

| Factor | Recommendation |
|---|---|
| Sample size | n > 5,000 for common variants; n > 50,000 for small effects |
| Confounders | Age, sex, ancestry, batch — include as covariates |
| Multiple testing | ~1M independent tests → threshold p < 5×10⁻⁸ |

**Why p < 5×10⁻⁸:** ~1M approximately independent SNPs after LD pruning. Bonferroni: α = 0.05 / 10⁶ = 5×10⁻⁸.

**Genomic inflation factor λ:** λ = median(χ²_obs) / median(χ²_expected_null). λ > 1.1 indicates stratification, cryptic relatedness, or polygenicity. LDSC disentangles these.

**Association tests:** Binary → logistic regression; quantitative → linear regression. Both include top PCs as ancestry covariates.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings("ignore")
plt.rcParams.update({"figure.dpi":120,"axes.spines.top":False,"axes.spines.right":False})
```

```python
rng = np.random.default_rng(42)
N_SAMPLES, N_SNPS = 1000, 3000

mafs = rng.uniform(0.05, 0.45, N_SNPS)
p = mafs
freqs = np.column_stack([(1-p)**2, 2*p*(1-p), p**2])
G = np.vstack([rng.choice([0,1,2], size=N_SAMPLES, p=freqs[j]) for j in range(N_SNPS)]).T

# Simulate binary phenotype with 3 causal SNPs
causal_idx = [300, 1200, 2500]
betas      = [0.6, -0.5,  0.4]
logit = sum(b * G[:, c] for b, c in zip(betas, causal_idx))
logit += rng.normal(0, 0.5, N_SAMPLES)
prob = 1 / (1 + np.exp(-logit))
phenotype = (rng.random(N_SAMPLES) < prob).astype(int)

chrom_ids = np.repeat(["chr1","chr2","chr3"], [1000,1000,1000])
positions = np.concatenate([rng.integers(1e6, 250e6, 1000),
                            rng.integers(1e6, 240e6, 1000),
                            rng.integers(1e6, 198e6, 1000)])
snp_df = pd.DataFrame({"chrom": chrom_ids, "pos": positions, "maf": mafs})
```

```python
def qc_snps(G, snp_df, min_maf=0.01, min_call_rate=0.95):
    maf = np.minimum(G.mean(0)/2, 1 - G.mean(0)/2)
    call_rate = np.ones(G.shape[1])
    controls = phenotype == 0
    G_ctrl = G[controls]
    obs_hets = (G_ctrl == 1).mean(0)
    p_ctrl = G_ctrl.mean(0) / 2
    exp_hets = 2 * p_ctrl * (1 - p_ctrl)
    chi2_hwe = (obs_hets - exp_hets)**2 / (exp_hets + 1e-9) * controls.sum()
    hwe_p = stats.chi2.sf(chi2_hwe, df=1)
    keep = (maf >= min_maf) & (call_rate >= min_call_rate) & (hwe_p >= 1e-6)
    return G[:, keep], snp_df[keep].reset_index(drop=True), keep

G_qc, snp_qc, keep_mask = qc_snps(G, snp_df)
```

```python
G_std = StandardScaler().fit_transform(G_qc)
pca = PCA(n_components=10)
PCs = pca.fit_transform(G_std)
```

```python
# Per-SNP logistic regression (Wald test) with PC covariates
from scipy.stats import chi2 as chi2_dist

pvals = np.ones(G_qc.shape[1])
covariates = PCs[:, :5]

for j in range(G_qc.shape[1]):
    snp = G_qc[:, j].astype(float)
    X_full = np.column_stack([np.ones(len(phenotype)), snp, covariates])
    try:
        lr = LogisticRegression(solver='lbfgs', max_iter=300, C=1e9,
                                fit_intercept=False, random_state=0)
        lr.fit(X_full, phenotype)
        mu = lr.predict_proba(X_full)[:, 1]
        w = mu * (1 - mu)
        XtWX = (X_full.T * w) @ X_full
        cov_beta = np.linalg.inv(XtWX)
        se = np.sqrt(max(cov_beta[1, 1], 1e-15))
        z = lr.coef_[0][1] / se
        pvals[j] = chi2_dist.sf(z**2, df=1)
    except Exception:
        pvals[j] = 1.0

snp_qc["pval"] = pvals
snp_qc["-log10p"] = -np.log10(snp_qc["pval"].clip(1e-300))
```

```python
fig, axes = plt.subplots(2, 1, figsize=(14, 8))

# Manhattan plot
chroms = ["chr1","chr2","chr3"]
offset = 0; xticks = []; xticklabels = []
for i, c in enumerate(chroms):
    sub = snp_qc[snp_qc["chrom"]==c].copy()
    x = sub["pos"] + offset
    axes[0].scatter(x, sub["-log10p"], c=["#1f77b4","#ff7f0e","#2ca02c"][i], s=2, alpha=0.5)
    xticks.append(offset + sub["pos"].median()); xticklabels.append(c)
    offset += sub["pos"].max() + 10e6
axes[0].axhline(-np.log10(5e-8), color="red", lw=1.2, ls="--", label="5×10⁻⁸")
axes[0].axhline(-np.log10(1e-5), color="orange", lw=0.8, ls="--", label="1×10⁻⁵ (suggestive)")
axes[0].set_xticks(xticks); axes[0].set_xticklabels(xticklabels)
axes[0].set_ylabel("-log₁₀(p)"); axes[0].set_title("Manhattan plot"); axes[0].legend(frameon=False)

# QQ plot
observed = np.sort(pvals)
expected = (np.arange(1, len(observed)+1)) / (len(observed)+1)
lam = np.median(stats.chi2.ppf(1-observed, 1)) / stats.chi2.ppf(0.5, 1)
axes[1].scatter(-np.log10(expected), -np.log10(observed), s=3, alpha=0.5, color="steelblue")
axes[1].plot([0, -np.log10(expected[-1])], [0, -np.log10(expected[-1])], "r--", lw=1)
axes[1].set_title(f"QQ plot  (λ = {lam:.3f})")
plt.tight_layout(); plt.show()
```

## LD Clumping

Significant SNPs are not independent — nearby LD-linked SNPs share the same signal. Clumping retains the most significant SNP in each LD block:
1. Sort by p-value
2. For each lead SNP, remove all SNPs within ±250 kb with r² > 0.1

```python
chr1_G = G_qc[snp_qc["chrom"]=="chr1", :][:, :50].T
r2_matrix = np.corrcoef(chr1_G) ** 2
```

## Fine-mapping

Statistical fine-mapping (SuSiE, FINEMAP) models the locus as L causal signals, computing a 95% credible set. Use functional annotations (ENCODE, GTEx eQTLs) to prioritize credible set members.

```python
# GWAS Catalog lookup pattern
import urllib.request, json

def query_gwas_catalog(rsid):
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations?projection=associationBySnp"
    try:
        with urllib.request.urlopen(url, timeout=5) as r:
            data = json.load(r)
        assocs = data.get("_embedded", {}).get("associations", [])
        return [(a.get("pvalue",""), [t["trait"] for t in a.get("efoTraits",[])]) for a in assocs[:3]]
    except Exception:
        return [("(offline)", ["GWAS Catalog requires internet access"])]
```

## GWAS Checklist

| Step | Method | Key Threshold |
|---|---|---|
| SNP QC | MAF, call rate, HWE | MAF > 0.01, HWE p > 10⁻⁶ |
| Sample QC | Missingness, ancestry outliers | missingness < 5% |
| Stratification | PCA → include PCs as covariates | top 10 PCs standard |
| Association | Logistic / linear regression | adjusted for covariates |
| Multiple testing | Bonferroni / GW threshold | p < 5×10⁻⁸ |
| LD / Clumping | PLINK --clump | r² < 0.1, window 250 kb |
| Fine-mapping | SuSiE / FINEMAP | 95% credible set |
| Lookup | GWAS Catalog REST API | validate prior associations |

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing causes off-by-one errors.
- **Batch effects**: Always check for batch confounding before interpreting biological signal.
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously.
