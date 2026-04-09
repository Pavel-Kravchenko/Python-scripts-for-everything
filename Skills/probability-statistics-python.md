---
name: probability-statistics-python
description: Python scipy.stats and statsmodels patterns — distribution fitting, regression diagnostics, bootstrap CI, permutation tests, ANOVA post-hoc, logistic regression, and survival analysis for biological data
---

# Probability & Statistics in Python

## When to Use
- Fitting parametric distributions to gene expression or growth data
- Running regression with full diagnostic output (residuals, influence, VIF)
- Computing confidence intervals via bootstrap when distribution is unknown
- Permutation tests when parametric assumptions are violated
- ANOVA with post-hoc testing (Tukey HSD) in a pure-Python pipeline
- Logistic regression for binary biological outcomes (disease/healthy, response/no-response)
- Survival analysis (Kaplan-Meier, log-rank) without switching to R

See `biostatistics-r` for: test selection decision tree, multiple testing correction, R syntax, NB vs Poisson rationale.

---

## Quick Reference

### scipy.stats Distribution API (all distributions share this interface)
| Method | Purpose | Biological example |
|--------|---------|-------------------|
| `.pdf(x, *params)` | Probability density at x | Likelihood of log2(TPM)=3 under fitted normal |
| `.cdf(x, *params)` | P(X ≤ x) | Fraction of genes below expression threshold |
| `.ppf(q, *params)` | Inverse CDF / quantile | 95th-percentile read depth cutoff |
| `.rvs(size, *params)` | Random samples | Simulate null distribution for permutation |
| `.fit(data)` | MLE parameter estimation | Fit gamma to GC-content across windows |
| `.kstest(data, dist)` | Kolmogorov-Smirnov fit test | Assess whether expression follows lognormal |

Common distributions: `norm`, `lognorm`, `gamma`, `expon`, `nbinom`, `poisson`, `beta`, `t`

### scipy.stats vs statsmodels — When to Use Which
| Task | scipy.stats | statsmodels |
|------|-------------|-------------|
| Single hypothesis test | `ttest_ind`, `mannwhitneyu` | overkill |
| Regression with diagnostics | not available | `OLS`, `GLM` — use this |
| ANOVA + post-hoc | `f_oneway` (stat only) | `anova_lm` + `pairwise_tukeyhsd` |
| Logistic regression | not available | `Logit`, or `GLM(family=Binomial())` |
| Distribution fitting (MLE) | `.fit()` method | not designed for this |
| Power / sample size | limited | `statsmodels.stats.power` |
| Bootstrap CI | manual (use resample) | `bootstrap` in `statsmodels.stats` |
| Survival analysis | not available | use `lifelines` package |

### Regression Diagnostics Checklist
| Assumption | How to Check | Remediation |
|------------|-------------|-------------|
| Linearity | Residual vs fitted plot — no curve | Transform X or Y (log, sqrt) |
| Normality of residuals | Q-Q plot, Shapiro-Wilk on residuals | Large n: CLT covers it; small n: use robust regression |
| Homoscedasticity | Scale-location plot, Breusch-Pagan | WLS or log-transform response |
| Independence | Study design review, DW statistic | Mixed models for repeated measures |
| No multicollinearity | VIF < 5 (warn), < 10 (drop/combine) | PCA on predictors, ridge regression |
| No influential outliers | Cook's D > 4/n flags a point | Investigate then exclude if justified |

### Effect Size Formulas
| Measure | Formula | Use case |
|---------|---------|---------|
| Cohen's d | `(μ₁ − μ₂) / s_pooled` | Two-group continuous (expression levels) |
| Eta-squared η² | `SS_between / SS_total` | ANOVA (proportion of variance explained) |
| Pearson r | `r = t / sqrt(t² + df)` | Convert t-test result to correlation effect |
| Odds ratio | `(a/b) / (c/d)` from 2×2 table | Case-control (disease association) |
| Cramér's V | `sqrt(χ² / (n · (min(r,c)−1)))` | Categorical association (genotype × phenotype) |

Small/medium/large thresholds: d 0.2/0.5/0.8 · η² 0.01/0.06/0.14 · V 0.1/0.3/0.5

### Confidence Interval Methods
| Method | When | scipy/statsmodels call |
|--------|------|----------------------|
| Normal theory (z) | Large n, known σ | `norm.interval(0.95, loc, scale=se)` |
| t-based | Small n, unknown σ | `t.interval(0.95, df, loc, scale=se)` |
| Bootstrap percentile | Any statistic, unknown dist | custom loop — see template below |
| Bootstrap BCa | Skewed statistic | `scipy.stats.bootstrap(..., method='BCa')` |
| Wilson (proportions) | Binomial proportion | `proportion_confint(k, n, method='wilson')` |
| Profile likelihood | GLM parameters | `result.conf_int()` from statsmodels |

---

## Key Patterns

### Distribution Fitting Workflow
1. Fit candidate distribution to data with `.fit()` (returns MLE params)
2. Run KS test — high p-value means no evidence of poor fit (not proof of fit)
3. Visual check with Q-Q plot (systematic deviation = wrong family)
4. Use fitted distribution for downstream inference (quantile thresholds, simulation)

### Assumption Checking Before Parametric Tests
Order matters: check normality first, then variance equality.
- Shapiro-Wilk: reliable for n < 5000; D'Agostino-Pearson for larger n
- Levene's test: preferred over Bartlett (robust to non-normality)
- If normality fails AND n < 30: switch to non-parametric

### Bootstrap Philosophy
Bootstrap resamples the observed data with replacement to approximate the sampling distribution of any statistic — no distributional assumption required. Requires n ≥ 20 for stable estimates; use ≥ 10,000 resamples for CIs.

### Permutation Test Philosophy
Permutation tests the null that group labels are exchangeable. Exact when all permutations are enumerable; Monte Carlo otherwise (shuffle 10,000 times). No distributional assumptions; p-value is exact under the null.

### VIF Calculation for Multicollinearity
VIF_j = 1 / (1 − R²_j) where R²_j is from regressing predictor j on all other predictors. VIF > 10 → serious problem. `variance_inflation_factor` from `statsmodels.stats.outliers_influence`.

---

## Code Templates

### fit_distribution — fit scipy.stats distribution to gene expression
```python
# shared imports for all templates below
import numpy as np, matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm, pandas as pd

def fit_distribution(data: np.ndarray, dist_name: str = "lognorm") -> dict:
    """Fit a scipy.stats distribution to data via MLE. Returns params + KS test."""
    dist = getattr(stats, dist_name)
    params = dist.fit(data)                          # MLE: returns (shape..., loc, scale)
    ks_stat, ks_p = stats.kstest(data, dist_name, args=params)
    return {"dist": dist_name, "params": params, "ks_stat": ks_stat, "ks_p": ks_p}

# Example: TPM values → fit lognormal, get 95th-percentile cutoff
tpm = np.random.lognormal(mean=2.5, sigma=1.2, size=500)
result = fit_distribution(tpm, "lognorm")          # KS p > 0.05 → no evidence against fit
cutoff = stats.lognorm.ppf(0.95, *result["params"])
```

### qq_plot — Q-Q plot for normality assessment
```python
def qq_plot(data: np.ndarray, dist: str = "norm", title: str = "") -> None:
    """Q-Q plot: systematic deviation from the line → wrong distribution family."""
    fig, ax = plt.subplots(figsize=(5, 5))
    (osm, osr), (slope, intercept, r) = stats.probplot(data, dist=dist, plot=None)
    ax.plot(osm, osr, "o", alpha=0.5)
    ax.plot(osm, slope * np.array(osm) + intercept, "r-")
    ax.set(xlabel=f"Theoretical quantiles ({dist})", ylabel="Sample quantiles",
           title=title or f"Q-Q plot (r={r:.3f})")
    plt.tight_layout()

log_expr = np.log2(tpm + 1)
qq_plot(log_expr, title="log2(TPM+1) normality check")
```

### bootstrap_ci — bootstrap confidence interval for any statistic
```python
def bootstrap_ci(data: np.ndarray, stat_fn, n_resamples: int = 10_000,
                 ci: float = 0.95, seed: int = 42) -> tuple[float, float]:
    """Percentile bootstrap CI for any statistic. Use ≥10,000 resamples for stable tails."""
    rng = np.random.default_rng(seed)
    bs = [stat_fn(rng.choice(data, size=len(data), replace=True)) for _ in range(n_resamples)]
    alpha = 1 - ci
    return np.percentile(bs, 100 * alpha / 2), np.percentile(bs, 100 * (1 - alpha / 2))

# Example: 95% CI for median log2(TPM)
lo, hi = bootstrap_ci(log_expr, np.median)
print(f"Median: {np.median(log_expr):.2f}  95% CI [{lo:.2f}, {hi:.2f}]")
```

### regression_diagnostics — OLS with full diagnostics via statsmodels
```python
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.stats.diagnostic import het_breuschpagan

def regression_diagnostics(X, y, feature_names: list[str]):
    """OLS regression with VIF, Breusch-Pagan, and Cook's distance summary."""
    X_const = sm.add_constant(X)           # ALWAYS add; statsmodels does NOT add intercept automatically
    model = sm.OLS(y, X_const).fit()
    print(model.summary())

    # Multicollinearity — skip index 0 (the constant)
    for i, name in enumerate(feature_names, start=1):
        v = variance_inflation_factor(X_const, i)
        flag = "HIGH" if v > 10 else ("warn" if v > 5 else "ok")
        print(f"VIF {name}: {v:.2f} [{flag}]")

    _, bp_p, _, _ = het_breuschpagan(model.resid, X_const)
    print(f"Breusch-Pagan p={bp_p:.4f} ({'heteroscedastic' if bp_p < 0.05 else 'OK'})")

    cooks_d = model.get_influence().cooks_distance[0]
    flagged = np.where(cooks_d > 4 / len(y))[0]
    if len(flagged):
        print(f"Influential points (Cook's D > 4/n): {flagged}")
    return model

# Example: protein abundance ~ mRNA + GC content
rng = np.random.default_rng(0)
mrna, gc = rng.normal(5, 1.5, 120), rng.uniform(0.35, 0.65, 120)
protein = 1.2 * mrna + 0.8 * gc + rng.normal(0, 0.5, 120)
model = regression_diagnostics(np.column_stack([mrna, gc]), protein, ["mRNA", "GC"])
```

### permutation_test — label-shuffling test for two groups
```python
def permutation_test(
    group_a: np.ndarray,
    group_b: np.ndarray,
    stat_fn=None,
    n_perm: int = 10_000,
    seed: int = 42,
) -> tuple[float, float]:
    """Monte Carlo permutation test. Returns observed statistic and two-sided p-value."""
    if stat_fn is None:
        stat_fn = lambda a, b: np.mean(a) - np.mean(b)

    observed = stat_fn(group_a, group_b)
    combined = np.concatenate([group_a, group_b])
    na = len(group_a)
    rng = np.random.default_rng(seed)

    null_dist = []
    for _ in range(n_perm):
        perm = rng.permutation(combined)
        null_dist.append(stat_fn(perm[:na], perm[na:]))
    null_dist = np.array(null_dist)
    p_value = np.mean(np.abs(null_dist) >= np.abs(observed))
    return observed, p_value

# Example: compare mean expression of BRCA1 in tumor vs normal — no normality assumed
tumor_expr  = np.array([8.2, 9.1, 7.8, 10.3, 8.9, 9.5, 7.6])
normal_expr = np.array([5.1, 4.8, 5.5, 6.0, 5.3, 4.9])
obs, p = permutation_test(tumor_expr, normal_expr)
print(f"Observed diff: {obs:.3f}, permutation p={p:.4f}")
```

### anova_with_posthoc — one-way ANOVA with Tukey HSD
```python
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

def anova_with_posthoc(data: pd.DataFrame, value_col: str, group_col: str) -> None:
    """One-way ANOVA via statsmodels formula API, then Tukey HSD for all pairs."""
    formula = f"{value_col} ~ C({group_col})"       # C() marks as categorical
    model = ols(formula, data=data).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)   # Type II SS
    print(anova_table)

    tukey = pairwise_tukeyhsd(data[value_col], data[group_col], alpha=0.05)
    print(tukey.summary())

# Example: log2(expression) of MYC across three tissue types
records = (
    [("liver",  v) for v in [4.1, 3.9, 4.5, 4.2, 4.0]] +
    [("kidney", v) for v in [6.2, 6.8, 6.4, 7.0, 6.5]] +
    [("brain",  v) for v in [5.0, 4.8, 5.2, 5.5, 4.9]]
)
df = pd.DataFrame(records, columns=["tissue", "log2_expr"])
anova_with_posthoc(df, "log2_expr", "tissue")
```

### logistic_regression_bio — binary outcome (disease/healthy)
```python
def logistic_regression_bio(df: pd.DataFrame, outcome: str, predictors: list[str]):
    """Logistic regression via statsmodels. Prints OR with 95% CI per predictor."""
    model = sm.Logit(df[outcome], sm.add_constant(df[predictors])).fit()
    conf = model.conf_int()
    or_table = pd.DataFrame({
        "OR": np.exp(model.params), "CI_low": np.exp(conf[0]),
        "CI_hi": np.exp(conf[1]),   "p": model.pvalues,
    }).drop("const").round(3)
    print(or_table)
    return model

# Example: disease ~ age + SNP dosage + log2(CRP)
rng = np.random.default_rng(1)
n = 200
df = pd.DataFrame({
    "age": rng.normal(55, 10, n), "snp": rng.choice([0,1,2], n, p=[0.5,0.4,0.1]),
    "log_crp": rng.normal(1.5, 0.8, n), "disease": rng.binomial(1, 0.35, n),
})
logistic_regression_bio(df, "disease", ["age", "snp", "log_crp"])
```

### survival_basics — Kaplan-Meier with log-rank test
```python
# pip install lifelines
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

def survival_basics(durations_a, events_a, label_a, durations_b, events_b, label_b):
    """Plot KM curves for two groups and print log-rank p-value."""
    lr = logrank_test(durations_a, durations_b, event_observed_A=events_a, event_observed_B=events_b)
    print(f"Log-rank p = {lr.p_value:.4f}")
    fig, ax = plt.subplots(figsize=(6, 4))
    for dur, ev, lbl in [(durations_a, events_a, label_a), (durations_b, events_b, label_b)]:
        kmf = KaplanMeierFitter()
        kmf.fit(dur, event_observed=ev, label=lbl)
        kmf.plot_survival_function(ax=ax, ci_show=True)
    ax.set(xlabel="Days post-diagnosis", ylabel="Survival probability",
           title=f"KM  log-rank p={lr.p_value:.4f}")
    plt.tight_layout()

# Example: BRCA1-mutant vs wild-type overall survival
rng = np.random.default_rng(7)
survival_basics(
    rng.exponential(800, 60).clip(0, 1825), rng.binomial(1, 0.7, 60), "BRCA1-mut",
    rng.exponential(1400, 80).clip(0, 1825), rng.binomial(1, 0.5, 80), "WT",
)
```

---

## Common Pitfalls

- **scipy one-sided vs two-sided defaults differ by test**: `mannwhitneyu` defaults to `alternative='two-sided'` since scipy 1.7, but always pass `alternative=` explicitly to avoid silent errors across versions.
- **Forgetting `sm.add_constant()`**: statsmodels OLS and Logit do NOT add an intercept by default — omitting it forces regression through the origin, biasing all coefficients.
- **statsmodels formula API vs array API**: `ols("y ~ x", data=df)` handles categorical variables automatically; array API (`OLS(y, X)`) requires manual dummy encoding with `pd.get_dummies`.
- **Bootstrap with too few resamples**: 1,000 resamples gives a rough CI; use ≥ 10,000 for stable 95% CI boundaries, especially in the tails.
- **Conflating statistical significance with effect size**: with n=5,000 RNA-seq samples, a fold change of 1.05× can reach p < 0.001; always report Cohen's d or log2FC alongside the p-value.
- **Parametric tests on raw count data**: Poisson/NB-distributed counts violate normality; use log2(count+1) transformation before t-test or use a GLM with appropriate family (Poisson/NB).
- **KS test as proof of distribution fit**: a high KS p-value means no evidence *against* the fit, not that the distribution is correct — always pair with a Q-Q plot.
- **VIF computed on design matrix including the constant**: `variance_inflation_factor` in statsmodels takes the full design matrix including the intercept column, but VIF is only meaningful for non-constant columns — index from 1, not 0.

---

## Related Skills
- `biostatistics-r` — R syntax, test selection decision tree, multiple testing, NB model for RNA-seq counts
- `numpy-pandas-wrangling` — data preparation before statistical modeling
- `data-visualization-bio` — plotting residuals, volcano plots, survival curves
- `ml-deep-learning-bio` — when to use ML classification instead of logistic regression
