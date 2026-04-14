---
name: bio-applied-gwas
description: Genome-Wide Association Studies (GWAS) with NumPy
tool_type: python
primary_tool: NumPy
---

# Genome-Wide Association Studies (GWAS)

Prerequisites: Modules 01–18. Module 15 (Population Genetics) strongly recommended.


**By the end you will be able to:**
- Apply SNP and sample QC filters
- Detect population stratification via PCA
- Run per-SNP association tests
- Generate Manhattan and QQ plots
- Interpret significance thresholds, LD, and clumping

**Attribution:** *Patterns inspired by NGSchool 2023 GWAS practical. Uses simulated genotype data.*

## Background: GWAS Study Design

GWAS tests ~1–10 million SNPs for association with a phenotype by comparing allele frequencies between cases and controls (binary traits) or correlating allele dosage with trait values (quantitative traits).

**Key design considerations:**

| Factor | Recommendation |
|---|---|
| Sample size | n > 5,000 for common variants; n > 50,000 for small effects |
| Phenotype definition | Binary (cases/controls) or quantitative; avoid mixing |
| Confounders | Age, sex, ancestry, batch; include as covariates |
| Population | Stratification is a major confounder; use ancestry-matched cohorts |
| Multiple testing | ~1M independent tests → threshold p < 5×10⁻⁸ |

**Why p < 5×10⁻⁸?**
The human genome has ~10 million common SNPs, of which ~1 million are approximately independent (due to linkage disequilibrium). Bonferroni correction: α = 0.05 / 10⁶ = 5×10⁻⁸. This is the genome-wide significance threshold established by Risch & Merikangas (1996) and confirmed empirically.

**Genomic inflation factor (λ):**
λ = median(χ²_observed) / median(χ²_expected_under_null)

λ is computed from the QQ plot. Under the null hypothesis (no associations, no confounding), λ = 1.0. λ > 1.1 indicates systematic inflation, which can be caused by:
- Population stratification (ancestry differences between cases and controls)
- Cryptic relatedness
- Model misspecification
- Actual polygenicity (many small true effects)

**LD Score Regression (LDSC):** If λ scales with the LD score (a measure of local LD), inflation is likely polygenic; if it does not scale, stratification is more likely. LDSC disentangles these.

**Association test:**
- Binary traits: logistic regression of case/control status on SNP dosage (0/1/2 copies of minor allele)
- Quantitative traits: linear regression of trait value on SNP dosage
- Both include principal components (PCs) of the genotype matrix as covariates to correct for ancestry

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
print("Libraries loaded.")
```

```python
rng = np.random.default_rng(42)
N_SAMPLES, N_SNPS = 1000, 3000

mafs = rng.uniform(0.05, 0.45, N_SNPS)
p = mafs
freqs = np.column_stack([(1-p)**2, 2*p*(1-p), p**2])
G = np.vstack([rng.choice([0,1,2], size=N_SAMPLES, p=freqs[j]) for j in range(N_SNPS)]).T
print(f"Genotype matrix: {G.shape} (samples × SNPs)")

# Simulate phenotype (binary) with 3 causal SNPs
causal_idx = [300, 1200, 2500]
betas      = [0.6, -0.5,  0.4]
logit = sum(b * G[:, c] for b, c in zip(betas, causal_idx))
logit += rng.normal(0, 0.5, N_SAMPLES)  # noise
prob = 1 / (1 + np.exp(-logit))
phenotype = (rng.random(N_SAMPLES) < prob).astype(int)
print(f"Cases: {phenotype.sum()}, Controls: {(phenotype==0).sum()}")

# SNP positions (fake chromosome coordinates for 3 chromosomes)
chrom_ids = np.repeat(["chr1","chr2","chr3"], [1000,1000,1000])
positions = np.concatenate([rng.integers(1e6, 250e6, 1000),
                            rng.integers(1e6, 240e6, 1000),
                            rng.integers(1e6, 198e6, 1000)])
snp_df = pd.DataFrame({"chrom": chrom_ids, "pos": positions, "maf": mafs})
print(snp_df.head())
```

```python
def qc_snps(G, snp_df, min_maf=0.01, min_call_rate=0.95):
    maf = np.minimum(G.mean(0)/2, 1 - G.mean(0)/2)
    call_rate = np.ones(G.shape[1])  # no missing in simulated data
    # HWE test (controls only)
    controls = phenotype == 0
    G_ctrl = G[controls]
    obs_hets = (G_ctrl == 1).mean(0)
    p_ctrl = G_ctrl.mean(0) / 2
    exp_hets = 2 * p_ctrl * (1 - p_ctrl)
    chi2_hwe = (obs_hets - exp_hets)**2 / (exp_hets + 1e-9) * controls.sum()
    hwe_p = stats.chi2.sf(chi2_hwe, df=1)
    keep = (maf >= min_maf) & (call_rate >= min_call_rate) & (hwe_p >= 1e-6)
    print(f"SNPs before QC: {G.shape[1]}")
    print(f"SNPs after QC:  {keep.sum()} (removed {(~keep).sum()})")
    return G[:, keep], snp_df[keep].reset_index(drop=True), keep

G_qc, snp_qc, keep_mask = qc_snps(G, snp_df)
print(f"\nFinal dataset: {G_qc.shape[0]} samples × {G_qc.shape[1]} SNPs")
```

```python
G_std = StandardScaler().fit_transform(G_qc)
pca = PCA(n_components=10)
PCs = pca.fit_transform(G_std)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].scatter(PCs[phenotype==0,0], PCs[phenotype==0,1], alpha=0.5, s=10, label="Control")
axes[0].scatter(PCs[phenotype==1,0], PCs[phenotype==1,1], alpha=0.5, s=10, label="Case")
axes[0].set_xlabel("PC1"); axes[0].set_ylabel("PC2")
axes[0].set_title("PC1 vs PC2 (population stratification check)")
axes[0].legend(markerscale=3)

axes[1].bar(range(1,11), pca.explained_variance_ratio_*100, color="steelblue")
axes[1].set_xlabel("PC"); axes[1].set_ylabel("% variance explained")
axes[1].set_title("PCA scree plot")
plt.tight_layout(); plt.show()
```

```python
# Per-SNP logistic regression with PC covariates (additive model)
# Standard GWAS tests SNP dosage (0/1/2) using logistic regression (binary trait)
# adjusted for principal components as ancestry covariates.
# chi2_contingency on binned genotypes ignores the additive model and cannot
# incorporate continuous covariates — replaced here with the Wald test.
from scipy.stats import chi2 as chi2_dist

pvals = np.ones(G_qc.shape[1])
covariates = PCs[:, :5]  # top 5 PCs as ancestry covariates

for j in range(G_qc.shape[1]):
    snp = G_qc[:, j].astype(float)
    X_full = np.column_stack([np.ones(len(phenotype)), snp, covariates])
    try:
        lr = LogisticRegression(solver='lbfgs', max_iter=300, C=1e9,
                                fit_intercept=False, random_state=0)
        lr.fit(X_full, phenotype)
        mu = lr.predict_proba(X_full)[:, 1]
        w = mu * (1 - mu)  # Fisher information weights
        XtWX = (X_full.T * w) @ X_full  # Information matrix
        try:
            cov_beta = np.linalg.inv(XtWX)
            se = np.sqrt(max(cov_beta[1, 1], 1e-15))
            z = lr.coef_[0][1] / se
            pvals[j] = chi2_dist.sf(z**2, df=1)  # Wald chi^2 with 1 df
        except np.linalg.LinAlgError:
            pvals[j] = 1.0
    except Exception:
        pvals[j] = 1.0

snp_qc["pval"] = pvals
snp_qc["-log10p"] = -np.log10(snp_qc["pval"].clip(1e-300))
print(f"SNPs with p < 5e-8: {(pvals < 5e-8).sum()}")
print(snp_qc.nsmallest(5, "pval")[["chrom","pos","-log10p","pval"]])
```

```python
fig, axes = plt.subplots(2, 1, figsize=(14, 8))

# Manhattan plot
chroms = ["chr1","chr2","chr3"]
offset = 0; xticks = []; xticklabels = []
colors = ["#1f77b4","#ff7f0e","#2ca02c"]
for i, c in enumerate(chroms):
    sub = snp_qc[snp_qc["chrom"]==c].copy()
    x = sub["pos"] + offset
    axes[0].scatter(x, sub["-log10p"], c=colors[i], s=2, alpha=0.5)
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
axes[1].set_xlabel("Expected -log₁₀(p)"); axes[1].set_ylabel("Observed -log₁₀(p)")
axes[1].set_title(f"QQ plot  (λ = {lam:.3f})")
plt.tight_layout(); plt.show()
print(f"Genomic inflation factor λ = {lam:.3f}  (ideal: 1.0; > 1.1 suggests confounding)")
```

## LD and Clumping

Significant SNPs are not independent — nearby SNPs in LD share the same signal. **Clumping** retains the most significant SNP in each LD block:

1. Sort hits by p-value (most significant first)
2. For each lead SNP, remove all SNPs within window (e.g. ±250 kb) with r² > 0.1
3. Remaining SNPs are independent signals


```python
# LD (r²) between adjacent SNPs on chr1 (illustrative)
chr1_G = G_qc[snp_qc["chrom"]=="chr1", :][:, :50].T  # 50 SNPs × samples
r2_matrix = np.corrcoef(chr1_G) ** 2

fig, ax = plt.subplots(figsize=(7,6))
im = ax.matshow(r2_matrix[:30,:30], cmap="Reds", vmin=0, vmax=1)
plt.colorbar(im, ax=ax, label="r²")
ax.set_title("LD (r²) matrix — chr1 SNPs 1–30")
ax.set_xlabel("SNP index"); ax.set_ylabel("SNP index")
plt.tight_layout(); plt.show()

# Identify hits above threshold
hits = snp_qc[snp_qc["pval"] < 5e-8].copy()
print(f"\nSignificant hits (p < 5e-8): {len(hits)}")
if len(hits) > 0:
    print(hits[["chrom","pos","pval"]].to_string(index=False))
```

## Fine-mapping Concept

After identifying a GWAS locus, **fine-mapping** narrows down which variant(s) actually drive the association. Because many SNPs in a locus are in LD with the lead SNP, it is often impossible to pinpoint the causal variant from association p-values alone.

**Statistical fine-mapping (SuSiE, FINEMAP):**
1. Model the locus as a sum of L causal signals (default L = 10)
2. Compute a **credible set**: the smallest set of SNPs containing the causal variant with 95% probability
3. Use functional annotations (ENCODE, GTEx eQTLs) to prioritize credible set members

**GWAS Catalog lookup:** submit your lead SNP rsIDs to [https://www.ebi.ac.uk/gwas/](https://www.ebi.ac.uk/gwas/) or use the REST API to find prior associations with the same locus.

```python
# Conceptual GWAS Catalog query (requires internet shown as pattern)
import urllib.request, json

def query_gwas_catalog(rsid):
    """Query GWAS Catalog for associations with a given rsID."""
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations?projection=associationBySnp"
    try:
        with urllib.request.urlopen(url, timeout=5) as r:
            data = json.load(r)
        assocs = data.get("_embedded", {}).get("associations", [])
        return [(a.get("pvalue",""), [t["trait"] for t in a.get("efoTraits",[])]) for a in assocs[:3]]
    except Exception as e:
        return [("(offline)", ["GWAS Catalog requires internet access"])]

# Example (rs429358  APOE ε4, known Alzheimer's variant)
results = query_gwas_catalog("rs429358")
for pval, traits in results:
    print(f"p={pval}: {', '.join(traits[:2])}")
```

## Summary: GWAS Checklist

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

**Related skill:** `gwas-population-genetics.md`


## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
