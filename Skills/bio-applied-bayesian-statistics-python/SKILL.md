---
name: bio-applied-bayesian-statistics-python
description: "Bayesian statistics with PyMC: prior specification, MCMC sampling, posterior analysis, and hierarchical models for biological data. Use when applying Bayesian inference to experiments."
tool_type: python
primary_tool: NumPy
---

# Bayesian Statistics in Python

- Contrast frequentist and Bayesian inference
- Specify priors and run prior predictive checks
- Fit Bayesian linear and generalized linear models with PyMC
- Use Bambi for formula-based mixed-effects models
- Compare models with LOO-CV and WAIC via ArviZ
- Fit Poisson and Negative Binomial models for count data

**Attribution:** *Statistical concepts based on Fränzi Korner-Nievergelt's Applied Statistics course (original R implementation). Python implementation by course authors. Uses public palmerpenguins dataset.*

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
    print(f"pymc {pm.__version__}, arviz {az.__version__}, bambi {bmb.__version__}")
except ImportError:
    print("Install: pip install pymc arviz bambi")

# palmerpenguins dataset (public)
try:
    from palmerpenguins import load_penguins
    penguins = load_penguins().dropna()
except ImportError:
    # fallback: simulate equivalent data
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

print(f"Dataset: {penguins.shape}, species: {penguins['species'].unique()}")
penguins.head()
```python

## Frequentist vs Bayesian Framing

**Frequentist:** Parameters are fixed unknowns. Confidence intervals describe the procedure, not the parameter. p-values test H₀.

**Bayesian:** Parameters are uncertain. We encode prior knowledge as a *prior distribution* p(θ). After observing data y, we update via Bayes' theorem:

$$p(\theta | y) = \frac{p(y | \theta) \cdot p(\theta)}{p(y)}$$

posterior ∝ likelihood × prior

**Key outputs:**
- **Posterior distribution** p(θ|y): full uncertainty over parameter values
- **Credible interval (HDI):** the smallest interval containing X% of the posterior — unlike confidence intervals, we can say "95% probability the parameter is in this range"
- **Posterior predictive:** simulate new data from the model to check fit

```python
# Example: model penguin body mass from flipper length
df = penguins[["body_mass_g","flipper_length_mm"]].copy()
df = (df - df.mean()) / df.std()  # standardize
df.columns = ["mass_std","flipper_std"]

X = df["flipper_std"].values
y = df["mass_std"].values

# Frequentist OLS
ols = smf.ols("mass_std ~ flipper_std", data=df).fit()
print("=== OLS ===")
print(f"β = {ols.params['flipper_std']:.3f}")
print(f"95% CI: {ols.conf_int().loc['flipper_std'].values.round(3)}")
print(f"p-value: {ols.pvalues['flipper_std']:.4f}")

# Bayesian
print("\n=== Bayesian (PyMC) ===")
try:
    with pm.Model() as linear_model:
        alpha = pm.Normal("alpha", mu=0, sigma=2)
        beta  = pm.Normal("beta",  mu=0, sigma=1)
        sigma = pm.HalfNormal("sigma", sigma=1)
        mu = alpha + beta * X
        obs = pm.Normal("obs", mu=mu, sigma=sigma, observed=y)
        idata_lin = pm.sample(1000, tune=500, target_accept=0.9,
                              progressbar=False, random_seed=42)

    summ = az.summary(idata_lin, var_names=["alpha","beta","sigma"], hdi_prob=0.95)
    print(summ[["mean","hdi_2.5%","hdi_97.5%","r_hat"]].round(3))
except Exception as e:
    print(f"PyMC sampling: {e}")
    print("(PyMC requires correct installation; showing structure only)")
```python

## Prior Specification

Choosing priors is the most distinctive Bayesian skill. Priors encode what we know *before* seeing the data.

| Prior type | When to use | Example |
|---|---|---|
| **Weakly informative** | Little domain knowledge; let data dominate | `Normal(0, 1)` on standardized scale |
| **Informative** | Strong prior knowledge (literature, pilot study) | `Normal(0.5, 0.1)` from previous study |
| **Non-informative (flat)** | Avoid — can cause convergence problems | `Uniform(-∞, ∞)` |

**Prior predictive check:** sample from the prior and simulate data. Do the simulated values make sense for your problem? (e.g., negative body masses are impossible)

```python
try:
    with pm.Model() as prior_check_model:
        alpha = pm.Normal("alpha", mu=0, sigma=5)   # wide prior
        beta  = pm.Normal("beta",  mu=0, sigma=3)   # wide prior
        sigma = pm.HalfNormal("sigma", sigma=2)
        mu = alpha + beta * X
        obs = pm.Normal("obs", mu=mu, sigma=sigma, observed=y)
        prior_pred = pm.sample_prior_predictive(200, random_seed=42)

    # Prior predictive samples
    pp_obs = prior_pred.prior_predictive["obs"].values.reshape(-1, len(y))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].hist(pp_obs.flatten(), bins=50, alpha=0.7, color="steelblue")
    axes[0].axvline(y.min(), color="red", lw=2, label=f"Data range")
    axes[0].axvline(y.max(), color="red", lw=2)
    axes[0].set_title("Prior predictive (wide priors)\nbody mass std — does range make sense?")
    axes[0].legend()

    # Now narrow priors
    with pm.Model() as narrow_prior_model:
        alpha = pm.Normal("alpha", mu=0, sigma=1)
        beta  = pm.Normal("beta",  mu=0, sigma=0.5)
        sigma = pm.HalfNormal("sigma", sigma=0.5)
        mu = alpha + beta * X
        obs = pm.Normal("obs", mu=mu, sigma=sigma, observed=y)
        narrow_pred = pm.sample_prior_predictive(200, random_seed=42)

    np_obs = narrow_pred.prior_predictive["obs"].values.reshape(-1, len(y))
    axes[1].hist(np_obs.flatten(), bins=50, alpha=0.7, color="steelblue")
    axes[1].axvline(y.min(), color="red", lw=2); axes[1].axvline(y.max(), color="red", lw=2)
    axes[1].set_title("Prior predictive (weakly informative priors)")
    plt.tight_layout(); plt.show()
    print("Weakly informative priors keep predicted values near the data range.")
except Exception as e:
    print(f"Prior predictive check skipped: {e}")
```python

## Multiple Regression and Collinearity

Adding multiple predictors can improve model fit — but collinearity (predictors correlated with each other) inflates uncertainty and makes interpretation difficult.

**Variance Inflation Factor (VIF):** VIF > 5–10 indicates collinearity concern.
$$\text{VIF}_j = \frac{1}{1 - R^2_j}$$
where R²_j is the R² from regressing predictor j on all other predictors.

```python
from statsmodels.stats.outliers_influence import variance_inflation_factor

multi_df = penguins[["body_mass_g","bill_length_mm","bill_depth_mm","flipper_length_mm"]].dropna()
multi_std = (multi_df - multi_df.mean()) / multi_df.std()

# OLS multiple regression
ols_multi = smf.ols("body_mass_g ~ bill_length_mm + bill_depth_mm + flipper_length_mm",
                     data=penguins).fit()
print(ols_multi.summary().tables[1])

# VIF
X_mat = sm.add_constant(multi_std[["bill_length_mm","bill_depth_mm","flipper_length_mm"]])
vif = pd.DataFrame({
    "feature": X_mat.columns,
    "VIF": [variance_inflation_factor(X_mat.values, i) for i in range(X_mat.shape[1])]
})
print("\nVariance Inflation Factors:")
print(vif.to_string(index=False))
print("\nVIF > 5 suggests collinearity; > 10 is problematic")
```python

## Model Comparison: LOO-CV and WAIC

In Bayesian statistics, we compare models by their **expected log predictive density (ELPD)** — how well each model predicts new data.

| Method | Full name | Notes |
|---|---|---|
| **LOO-CV** | Leave-one-out cross-validation | More robust; prefers fewer parameters |
| **WAIC** | Widely Applicable IC | Faster; sensitive to influential observations |

Both are computed from the posterior using Pareto-smoothed importance sampling (PSIS-LOO).

`az.compare()` returns a table sorted by ELPD — the model with the highest ELPD (least negative) is preferred.

```python
try:
    # Model 1: flipper length only
    with pm.Model() as m1:
        a = pm.Normal("a", 0, 2); b = pm.Normal("b", 0, 1)
        s = pm.HalfNormal("s", 1)
        pm.Normal("y", mu=a + b*X, sigma=s, observed=y)
        idata_m1 = pm.sample(800, tune=400, progressbar=False, random_seed=1)
        pm.sample_posterior_predictive(idata_m1, extend_inferencedata=True)

    # Model 2: flipper + bill length
    X2 = (penguins["bill_length_mm"] - penguins["bill_length_mm"].mean()) / penguins["bill_length_mm"].std()
    with pm.Model() as m2:
        a  = pm.Normal("a", 0, 2)
        b1 = pm.Normal("b1", 0, 1); b2 = pm.Normal("b2", 0, 1)
        s  = pm.HalfNormal("s", 1)
        pm.Normal("y", mu=a + b1*X + b2*X2.values, sigma=s, observed=y)
        idata_m2 = pm.sample(800, tune=400, progressbar=False, random_seed=2)
        pm.sample_posterior_predictive(idata_m2, extend_inferencedata=True)

    comp = az.compare({"flipper_only": idata_m1, "flipper+bill": idata_m2}, ic="loo")
    print("Model comparison (LOO-CV):")
    print(comp[["elpd_loo","p_loo","d_loo","weight"]].round(2))
    az.plot_compare(comp, insample_dev=False)
    plt.title("LOO-CV model comparison"); plt.tight_layout(); plt.show()
except Exception as e:
    print(f"Model comparison skipped: {e}")
    print("Pattern: az.compare({'m1': idata_m1, 'm2': idata_m2}, ic='loo')")
```python

## Linear Mixed-Effects Models

When data has a grouped structure (measurements nested within individuals, sites, or species), **mixed-effects models** partition variance into fixed effects (population-level) and random effects (group-level).

**Notation (Bambi/lme4):**
- `y ~ x` — fixed effect only
- `y ~ x + (1|group)` — fixed effect + random intercept per group
- `y ~ x + (x|group)` — fixed effect + random slope per group

**Partial pooling:** random effects "share information" across groups — groups with few observations borrow strength from the population estimate.

```python
try:
    import bambi as bmb

    # Random intercept per species
    m_me = bmb.Model(
        "body_mass_g ~ flipper_length_mm + (1|species)",
        data=penguins,
        family="gaussian"
    )
    idata_me = m_me.fit(draws=800, tune=400, target_accept=0.9,
                         progressbar=False, random_seed=42)

    print("Mixed-effects model summary:")
    summ = az.summary(idata_me, var_names=["flipper_length_mm", "Intercept"], hdi_prob=0.95)
    print(summ[["mean","hdi_2.5%","hdi_97.5%","r_hat"]].round(3))

    # Plot species random intercepts
    az.plot_forest(idata_me, var_names=["1|species"], combined=True)
    plt.title("Random intercepts by species"); plt.tight_layout(); plt.show()

except Exception as e:
    print(f"Mixed-effects skipped: {e}")
    print("Pattern: bmb.Model('y ~ x + (1|group)', data=df).fit(draws=1000)")
```python

## Generalized Linear Models (GLMs)

Not all outcomes are normally distributed. GLMs extend linear models with:
1. A **link function** connecting the linear predictor to the mean
2. A **family** distribution appropriate for the data type

| Family | Link | Use when | R/Python |
|---|---|---|---|
| Gaussian | identity | Continuous, symmetric | `family="gaussian"` |
| Binomial | logit | Binary outcomes (0/1) | `family="bernoulli"` |
| Poisson | log | Count data (no overdispersion) | `pm.Poisson` |
| Negative Binomial | log | Count data with overdispersion | `pm.NegativeBinomial` |
| Zero-Inflated Poisson | log | Excess zeros | `pm.ZeroInflatedPoisson` |

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
