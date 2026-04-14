---
name: bio-applied-epigenetic-clocks
description: Epigenetic Clocks and Aging Analysis with Matplotlib
tool_type: python
primary_tool: Matplotlib
---

# Epigenetic Clocks and Aging Analysis

- [Horvath clock (2013, Genome Biology)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115)
- [GrimAge (Lu et al. 2019)](https://www.aging-us.com/article/101414/)
- [methylclock R package](https://bioconductor.org/packages/release/bioc/html/methylclock.html)

## Epigenetic Clock Theory

### The core observation

In 2013, Steve Horvath published a landmark paper showing that **353 CpG sites** scattered across the genome — when their methylation values are combined in an elastic-net regression model — can predict a person's chronological age with remarkable accuracy (median absolute error ~3.6 years across tissues). This was surprising because:

1. The same model works across diverse tissues (blood, brain, breast, kidney, etc.)
2. It uses only methylation data — no other clinical variables
3. The CpGs selected by the model are distributed genome-wide, not clustered

### Why does methylation track age?

Methylation changes with age through several mechanisms:
- **Epigenetic drift**: Random errors during DNA methylation maintenance after each cell division accumulate, causing gradual changes in bulk tissue methylation.
- **Programmatic changes**: Some CpGs show directional, reproducible methylation increases or decreases with age — these are what clocks exploit.
- **Developmental programs**: Some age-associated changes reflect the continued activity of developmental gene regulatory programs into adulthood.

The exact causal mechanism is debated, but the empirical predictive accuracy is robust.

### Generations of epigenetic clocks

| Generation | Examples | Trained to predict | Key advance |
|---|---|---|---|
| First | Horvath (2013), Hannum (2013) | Chronological age | Proof of concept; pan-tissue |
| Second | PhenoAge (2018), GrimAge (2019) | Biological/mortality age | Better predictor of health outcomes |
| Third | DunedinPACE (2022), PCClocks | Rate of aging | Longitudinal validation |

### Horvath clock construction

Horvath used 8,000 samples from 82 Illumina datasets. He trained an **elastic-net regression** (a mixture of L1/Lasso and L2/Ridge penalties) with a transformed age as the response:

$$\text{age\_transformed} = \begin{cases} \log(\text{age}+1) - \log(21) & \text{if age} < 20 \\ \text{age} & \text{if age} \geq 20 \end{cases}$$

This transformation handles the rapid methylation changes in childhood differently from slower adult changes. The elastic net selected 353 non-zero CpGs from ~21,000 candidates.

## CpG Site Selection and Elastic Net Regression

### Why elastic net for clock building?

Standard linear regression cannot handle the scenario where the number of predictors (850,000 CpGs) far exceeds the number of observations (n samples). Elastic net adds two regularization penalties to the loss function:

$$\text{Loss} = \underbrace{\sum_i (y_i - \hat{y}_i)^2}_{\text{MSE}} + \lambda \left[ \underbrace{\alpha \sum_j |\beta_j|}_{\text{L1 / Lasso}} + \underbrace{(1-\alpha) \sum_j \beta_j^2}_{\text{L2 / Ridge}} \right]$$

- **L1 (Lasso)**: drives many coefficients exactly to zero → automatic feature selection (picks a sparse set of CpGs)
- **L2 (Ridge)**: shrinks non-zero coefficients → stability when correlated CpGs compete
- **α**: mixing parameter (Horvath used α ≈ 0.5); **λ**: regularization strength (chosen by cross-validation)

### Key properties of Horvath clock CpGs

- 193 CpGs increase methylation with age (gain methylation)
- 160 CpGs decrease methylation with age (lose methylation)
- Enriched for: PRC2 polycomb targets, developmental TF binding sites, bivalent chromatin domains
- Not enriched for: cell-type-specific enhancers (which explains tissue robustness)

The PRC2/polycomb enrichment suggests that age-associated methylation gain occurs at genomic loci already poised for silencing — CpG island promoters with H3K27me3 marks in stem cells.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from sklearn.linear_model import ElasticNetCV, LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict

np.random.seed(42)

#
# Simulate a methylation dataset:
# 200 samples with ages 18-85 500 CpG features
# 50 clock CpGs have a real age effect 450 are noise
#
N_SAMPLES   = 200
N_CPGS      = 500
N_CLOCK     = 50    # CpGs with true age association

chronological_ages = np.random.uniform(18, 85, N_SAMPLES)

# Horvath age transformation
def transform_age(age):
    return np.where(age < 20, np.log(age + 1) - np.log(21), age)

def inv_transform_age(t):
    """Inverse of Horvath age transformation."""
    return np.where(t < 0, np.exp(t + np.log(21)) - 1, t)

transformed_ages = transform_age(chronological_ages)

# Build feature matrix
# Clock CpGs: linear relationship with transformed age + noise
clock_slopes = np.random.uniform(-0.005, 0.005, N_CLOCK)   # small per-unit-age effect
beta_clock = (
    np.outer(transformed_ages, clock_slopes)  # age effect
    + np.random.normal(0, 0.05, (N_SAMPLES, N_CLOCK))   # biological noise
    + np.random.uniform(0.1, 0.9, N_CLOCK)   # intercept per CpG
)
beta_clock = np.clip(beta_clock, 0.01, 0.99)

# Noise CpGs: no age relationship
beta_noise = np.random.beta(2, 3, (N_SAMPLES, N_CPGS - N_CLOCK))

X = np.hstack([beta_clock, beta_noise])

print(f"Feature matrix shape: {X.shape}  ({N_SAMPLES} samples × {N_CPGS} CpGs)")
print(f"  {N_CLOCK} CpGs have true age associations")
print(f"  {N_CPGS - N_CLOCK} CpGs are noise")
print(f"Age range: {chronological_ages.min():.0f}–{chronological_ages.max():.0f} years")
print(f"\nTransformed age range: {transformed_ages.min():.2f}–{transformed_ages.max():.2f}")
```

## Age Acceleration Analysis

**Epigenetic age acceleration** (EAA) is defined as the residual of DNAmAge regressed on chronological age:

$$\text{EAA} = \text{DNAmAge} - \text{chronological age}$$

A positive EAA means a person is biologically "older" than their chronological age; negative EAA means "younger." This residual is what epidemiologists associate with lifestyle and disease.

### Intrinsic vs extrinsic age acceleration

- **Intrinsic EAA (IEAA)**: residual after also adjusting for blood cell composition (CD4T, CD8T, NK, B cells, monocytes, granulocytes). Reflects aging within cell types.
- **Extrinsic EAA (EEAA)**: not adjusted for cell type; driven partly by immune aging (changes in blood cell proportions). Tracks immune system aging.

### Known associations with age acceleration

Studies have found EAA (typically using GrimAge or PhenoAge) to be associated with:
- **Smoking**: 3–5 years acceleration per pack-year in blood
- **BMI/obesity**: ~1–2 years acceleration per 5 BMI units
- **Physical activity**: deceleration in athletes / active individuals
- **Socioeconomic status**: lower SES → higher acceleration
- **Disease**: cancer, Alzheimer's, HIV, lupus all show accelerated EAA
- **Longevity**: centenarians and their offspring show decelerated EAA

### Horvath vs PhenoAge vs GrimAge

| Clock | Prediction target | Best mortality predictor? |
|---|---|---|
| Horvath | Chronological age | Weak |
| Hannum | Chronological age (blood) | Moderate |
| PhenoAge | Biological/phenotypic age (clinical markers) | Strong |
| GrimAge | Mortality risk (10-year survival) | Strongest |
| DunedinPACE | Pace of aging (longitudinal) | Best for intervention studies |

## GrimAge and Mortality Prediction

### From chronological to mortality-based clocks

GrimAge (Lu et al. 2019, *Aging*) was developed by training on time-to-death data rather than chronological age. It consists of:

1. **DNAm surrogate biomarkers**: 7 plasma proteins whose circulating levels are predictable from blood DNA methylation (including PAI-1, adrenomedullin, GDF15, leptin, tissue plasminogen activator, β-2 microglobulin, and TIMP-1).
2. **DNAmPAck-years**: a methylation-based predictor of smoking pack-years (even in never-smokers, reflecting secondhand smoke).
3. **Weighted combination**: GrimAge = baseline + weighted sum of 8 components.

GrimAge outperforms all other clocks for predicting:
- All-cause mortality
- Coronary heart disease
- Cancer incidence
- Physical disability
- Cognitive decline

### PhenoAge

PhenoAge (Levine et al. 2018) was trained on a "phenotypic age" composite score derived from 9 clinical biomarkers (albumin, creatinine, glucose, CRP, lymphocyte %, mean cell volume, red cell distribution width, alkaline phosphatase, white blood cell count) using the Klemera-Doubal biological age method.

Unlike Horvath/Hannum which were trained directly on chronological age, PhenoAge captures the biological heterogeneity at any given chronological age — some 60-year-olds are as biologically "old" as 70-year-olds.

### Key formula: GrimAge intercept adjustment

Because GrimAge components predict plasma protein levels (not age directly), the final age estimate involves:

$$\text{GrimAge} = \beta_0 + \sum_k \beta_k \cdot \hat{p}_k + \beta_{\text{smoke}} \cdot \text{DNAmPack-years}$$

where $\hat{p}_k$ = methylation-predicted plasma protein level for component $k$.

```python
np.random.seed(33)

# Simulate comparison of four clocks on the same cohort
# Clocks differ in how well they predict mortality (simulated)
n_cohort = 120
chron  = np.random.uniform(30, 80, n_cohort)
noise  = np.random.normal(0, 3, n_cohort)

# Simulate 4 clock estimates
horvath  = chron + 0.8 * noise + np.random.normal(0, 2, n_cohort)       # good chronological fit
hannum   = chron + 0.7 * noise + np.random.normal(0, 2.5, n_cohort)     # blood-specific
phenoage = chron + noise + 2.0 * (np.random.normal(0, 1, n_cohort))     # biological age
grimage  = chron + 1.5 * noise + np.random.normal(0, 2.5, n_cohort)     # mortality-tuned (more spread)

# Simulate 10-year mortality: logistic with GrimAge acceleration as driver
grimage_eaa  = grimage - chron
death_prob   = 1 / (1 + np.exp(-(0.05 * chron + 0.15 * grimage_eaa - 3.5)))
died         = np.random.binomial(1, death_prob, n_cohort)

clocks_dict = {'Horvath': horvath, 'Hannum': hannum,
               'PhenoAge': phenoage, 'GrimAge': grimage}
clock_colors = ['#1976D2', '#43A047', '#FB8C00', '#E53935']

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Panel 1: Correlation of each clock with chronological age
mae_vals = {}
r_vals   = {}
for (name, clock), col in zip(clocks_dict.items(), clock_colors):
    r_c, _ = stats.pearsonr(chron, clock)
    mae_c  = np.abs(chron - clock).mean()
    mae_vals[name] = mae_c
    r_vals[name]   = r_c
    axes[0].scatter(chron, clock, alpha=0.4, s=20, color=col, label=f'{name} (r={r_c:.3f})')

axes[0].plot([25, 85], [25, 85], 'k--', lw=1.5)
axes[0].set_xlabel('Chronological Age', fontsize=10)
axes[0].set_ylabel('Clock Age Estimate', fontsize=10)
axes[0].set_title('Clock Accuracy vs Chronological Age', fontsize=10)
axes[0].legend(fontsize=8)

# Panel 2: GrimAge acceleration as predictor of 10-year mortality
grimage_eaa_sorted = np.sort(grimage_eaa)
death_rates = []
bins = np.percentile(grimage_eaa, np.linspace(0, 100, 6))  # quintiles
for lo, hi in zip(bins[:-1], bins[1:]):
    mask = (grimage_eaa >= lo) & (grimage_eaa < hi)
    death_rates.append(died[mask].mean() if mask.sum() > 0 else 0)

bin_centers = [(bins[i] + bins[i+1]) / 2 for i in range(len(bins)-1)]
axes[1].bar(range(5), [d*100 for d in death_rates],
            color=['#BBDEFB','#90CAF9','#64B5F6','#42A5F5','#E53935'],
            edgecolor='white')
axes[1].set_xticks(range(5))
axes[1].set_xticklabels([f'Q{i+1}\n({bin_centers[i]:.1f} yr)' for i in range(5)], fontsize=9)
axes[1].set_ylabel('10-year Mortality Rate (%)', fontsize=10)
axes[1].set_title('GrimAge Acceleration Quintiles\nvs Simulated 10-year Mortality', fontsize=10)

plt.tight_layout()
plt.savefig('clock_comparison_mortality.png', dpi=120, bbox_inches='tight')
plt.show()

print("Clock MAE vs chronological age:")
for name, mae in mae_vals.items():
    print(f"  {name:<12}: {mae:.2f} years")
```

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
