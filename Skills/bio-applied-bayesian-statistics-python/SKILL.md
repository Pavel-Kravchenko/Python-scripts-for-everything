---
name: bio-applied-bayesian-statistics-python
description: "Bayesian statistics with PyMC: prior specification, MCMC sampling, posterior analysis, and hierarchical models for biological data. Use when applying Bayesian inference to experiments."
tool_type: python
primary_tool: NumPy
---

# Bayesian Statistics in Python

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

try:
    import pymc as pm
    import arviz as az
    import bambi as bmb
except ImportError:
    print("Install: pip install pymc arviz bambi")

try:
    from palmerpenguins import load_penguins
    penguins = load_penguins().dropna()
except ImportError:
    rng = np.random.default_rng(42)
    n = 333
    species = rng.choice(["Adelie","Chinstrap","Gentoo"], n)
    bill_length = rng.normal(43, 5, n)
    bill_depth  = rng.normal(17, 2, n)
    flipper     = rng.normal(200, 14, n)
    body_mass   = 100 * flipper + rng.normal(0, 500, n)
    island = rng.choice(["Biscoe","Dream","Torgersen"], n)
    penguins = pd.DataFrame({"species": species, "bill_length_mm": bill_length,
                              "bill_depth_mm": bill_depth, "flipper_length_mm": flipper,
                              "body_mass_g": body_mass, "island": island, "sex": rng.choice(["male","female"],n)})
```

## Frequentist vs Bayesian

**Frequentist:** parameters are fixed unknowns; CIs describe the procedure; p-values test H₀.

**Bayesian:** parameters are uncertain. Posterior ∝ likelihood × prior:
$$p(\theta | y) = \frac{p(y | \theta) \cdot p(\theta)}{p(y)}$$

**Key outputs:** posterior distribution p(θ|y); HDI credible interval (unlike CIs, direct probability statement); posterior predictive for model checking.

```python
df = penguins[["body_mass_g","flipper_length_mm"]].copy()
df = (df - df.mean()) / df.std()
df.columns = ["mass_std","flipper_std"]
X, y = df["flipper_std"].values, df["mass_std"].values

# Frequentist OLS
ols = smf.ols("mass_std ~ flipper_std", data=df).fit()
print(f"β = {ols.params['flipper_std']:.3f}, 95% CI: {ols.conf_int().loc['flipper_std'].values.round(3)}")

# Bayesian PyMC
with pm.Model() as linear_model:
    alpha = pm.Normal("alpha", mu=0, sigma=2)
    beta  = pm.Normal("beta",  mu=0, sigma=1)
    sigma = pm.HalfNormal("sigma", sigma=1)
    mu = alpha + beta * X
    obs = pm.Normal("obs", mu=mu, sigma=sigma, observed=y)
    idata_lin = pm.sample(1000, tune=500, target_accept=0.9, progressbar=False, random_seed=42)

summ = az.summary(idata_lin, var_names=["alpha","beta","sigma"], hdi_prob=0.95)
print(summ[["mean","hdi_2.5%","hdi_97.5%","r_hat"]].round(3))
```

## Prior Specification

| Prior type | When to use | Example |
|---|---|---|
| **Weakly informative** | Little domain knowledge | `Normal(0, 1)` on standardized scale |
| **Informative** | Strong prior knowledge | `Normal(0.5, 0.1)` from previous study |
| **Non-informative (flat)** | Avoid — causes convergence problems | `Uniform(-∞, ∞)` |

**Prior predictive check:** sample from prior and simulate data — verify simulated values are physically plausible before fitting.

```python
with pm.Model() as prior_check_model:
    alpha = pm.Normal("alpha", mu=0, sigma=1)
    beta  = pm.Normal("beta",  mu=0, sigma=0.5)
    sigma = pm.HalfNormal("sigma", sigma=0.5)
    mu = alpha + beta * X
    obs = pm.Normal("obs", mu=mu, sigma=sigma, observed=y)
    prior_pred = pm.sample_prior_predictive(200, random_seed=42)

pp_obs = prior_pred.prior_predictive["obs"].values.reshape(-1, len(y))
```

## Multiple Regression and Collinearity

**VIF > 5–10 indicates collinearity concern:**
$$\text{VIF}_j = \frac{1}{1 - R^2_j}$$

```python
from statsmodels.stats.outliers_influence import variance_inflation_factor

multi_std = (penguins[["body_mass_g","bill_length_mm","bill_depth_mm","flipper_length_mm"]].dropna()
             .pipe(lambda df: (df - df.mean()) / df.std()))

X_mat = sm.add_constant(multi_std[["bill_length_mm","bill_depth_mm","flipper_length_mm"]])
vif = pd.DataFrame({
    "feature": X_mat.columns,
    "VIF": [variance_inflation_factor(X_mat.values, i) for i in range(X_mat.shape[1])]
})
print(vif.to_string(index=False))
```

## Model Comparison: LOO-CV and WAIC

| Method | Notes |
|---|---|
| **LOO-CV** | More robust; prefers fewer parameters |
| **WAIC** | Faster; sensitive to influential observations |

Both computed via Pareto-smoothed importance sampling. `az.compare()` returns table sorted by ELPD — highest (least negative) is preferred.

```python
comp = az.compare({"flipper_only": idata_m1, "flipper+bill": idata_m2}, ic="loo")
print(comp[["elpd_loo","p_loo","d_loo","weight"]].round(2))
az.plot_compare(comp, insample_dev=False)
plt.tight_layout(); plt.show()
```

## Linear Mixed-Effects Models (Bambi)

**Notation:** `y ~ x + (1|group)` — fixed effect + random intercept. `y ~ x + (x|group)` — random slope. Random effects use partial pooling: sparse groups borrow strength from population estimate.

```python
m_me = bmb.Model(
    "body_mass_g ~ flipper_length_mm + (1|species)",
    data=penguins,
    family="gaussian"
)
idata_me = m_me.fit(draws=800, tune=400, target_accept=0.9, progressbar=False, random_seed=42)

summ = az.summary(idata_me, var_names=["flipper_length_mm", "Intercept"], hdi_prob=0.95)
print(summ[["mean","hdi_2.5%","hdi_97.5%","r_hat"]].round(3))

az.plot_forest(idata_me, var_names=["1|species"], combined=True)
plt.title("Random intercepts by species"); plt.tight_layout(); plt.show()
```

## GLM Family Reference

| Family | Link | Use when |
|---|---|---|
| Gaussian | identity | Continuous, symmetric |
| Binomial | logit | Binary outcomes (0/1) |
| Poisson | log | Count data (no overdispersion) |
| Negative Binomial | log | Count data with overdispersion |
| Zero-Inflated Poisson | log | Excess zeros |

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing causes off-by-one errors.
- **Batch effects**: Always check for batch confounding before interpreting biological signal.
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously.
