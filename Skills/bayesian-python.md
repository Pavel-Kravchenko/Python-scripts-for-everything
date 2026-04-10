---
name: bayesian-python
description: Bayesian statistics in Python — linear models with pymc and statsmodels, prior specification, model comparison with WAIC/LOO-CV via arviz, mixed-effects with bambi, GLMs (Poisson, NegBin, Bernoulli)
---

## When to Use

Use this skill when:
- Fitting Bayesian linear or generalized linear models (pymc, bambi)
- Specifying informative or weakly informative priors
- Comparing models with WAIC or LOO-CV
- Fitting linear mixed-effects models (random intercepts/slopes)
- Fitting GLMs: Poisson, Negative Binomial, Bernoulli, Zero-Inflated
- Interpreting credible intervals and posterior predictive checks

## Quick Reference

| Task | Tool | Key Method |
|---|---|---|
| Bayesian linear model | `pymc` | `pm.Model()`, `pm.sample()` |
| Formula-based model | `bambi` | `bmb.Model("y ~ x + (1|group)")` |
| Model comparison | `arviz` | `az.compare({"m1": idata1, "m2": idata2})` |
| LOO-CV | `arviz` | `az.loo(idata)` |
| WAIC | `arviz` | `az.waic(idata)` |
| Posterior summary | `arviz` | `az.summary(idata)` |
| Prior predictive | `pymc` | `pm.sample_prior_predictive()` |
| GLM (Bayesian) | `pymc` | `pm.Poisson`, `pm.NegativeBinomial` |
| Mixed effects | `bambi` | random effects: `(1|group)` syntax |
| VIF (collinearity) | `statsmodels` | `variance_inflation_factor()` |

## Key Patterns

**Pattern 1: Bayesian linear model (pymc)**
```python
import pymc as pm, arviz as az

with pm.Model() as model:
    alpha = pm.Normal("alpha", mu=0, sigma=10)
    beta  = pm.Normal("beta",  mu=0, sigma=1)
    sigma = pm.HalfNormal("sigma", sigma=1)
    mu = alpha + beta * X
    obs = pm.Normal("obs", mu=mu, sigma=sigma, observed=y)
    idata = pm.sample(1000, tune=500, target_accept=0.9, progressbar=False)

az.summary(idata, var_names=["alpha","beta","sigma"])
```

**Pattern 2: Prior predictive check**
```python
with model:
    prior_pred = pm.sample_prior_predictive(100)
# prior_pred.prior_predictive["obs"] shape: (1, 100, n_obs)
az.plot_ppc(prior_pred, group="prior", observed=False)
```

**Pattern 3: Formula model with bambi**
```python
import bambi as bmb

m = bmb.Model("y ~ x1 + x2", data=df, family="gaussian")
idata = m.fit(draws=1000, tune=500, progressbar=False)
az.summary(idata)
```

**Pattern 4: Poisson GLM**
```python
with pm.Model() as poisson_model:
    a = pm.Normal("a", 0, 2)
    b = pm.Normal("b", 0, 1)
    mu = pm.math.exp(a + b * X)
    obs = pm.Poisson("obs", mu=mu, observed=y_count)
    idata = pm.sample(1000, tune=500, progressbar=False)
```

**Pattern 5: Model comparison**
```python
import arviz as az
comp = az.compare({"linear": idata_lin, "quadratic": idata_quad}, ic="loo")
print(comp)  # elpd_loo, p_loo, d_loo, weight columns
```

## Code Templates

**Template 1: Bayesian vs frequentist linear model comparison**
```python
import pymc as pm, arviz as az, statsmodels.formula.api as smf

# Frequentist
ols = smf.ols("y ~ x", data=df).fit()
print(f"OLS: β = {ols.params['x']:.3f}, 95% CI: {ols.conf_int().loc['x'].values}")

# Bayesian
with pm.Model() as bayes_model:
    alpha = pm.Normal("alpha", 0, 10)
    beta  = pm.Normal("beta",  0, 2)
    sigma = pm.HalfNormal("sigma", 2)
    mu = alpha + beta * df["x"].values
    pm.Normal("y_obs", mu=mu, sigma=sigma, observed=df["y"].values)
    idata = pm.sample(1000, tune=500, target_accept=0.9, progressbar=False)

summ = az.summary(idata, var_names=["beta"], hdi_prob=0.95)
print(f"Bayes: β = {summ['mean']['beta']:.3f}, 95% HDI: [{summ['hdi_2.5%']['beta']:.3f}, {summ['hdi_97.5%']['beta']:.3f}]")
```

**Template 2: Linear mixed-effects with bambi**
```python
import bambi as bmb, arviz as az

# Random intercept per group
m = bmb.Model("outcome ~ predictor + (1|group_id)", data=df, family="gaussian")
idata = m.fit(draws=1000, tune=500, target_accept=0.9, progressbar=False)
az.plot_posterior(idata, var_names=["predictor"])
```

## Common Pitfalls

- **Divergences:** increase `target_accept` (0.95) or reparameterize (non-centered); never ignore divergences
- **Slow mixing:** check R-hat > 1.01 and ESS < 100; increase tune steps or simplify model
- **Prior scale matters:** `pm.Normal("b", 0, 100)` is nearly flat — use domain knowledge; `sigma=1` on standardized data is weakly informative
- **LOO vs WAIC:** prefer LOO for models with outliers (more robust); WAIC faster but sensitive to influential points
- **bambi formula syntax:** `(1|group)` = random intercept; `(x|group)` = random slope; multiple random effects allowed
- **Poisson overdispersion:** if variance >> mean, use Negative Binomial instead

## Related Skills

- `probability-statistics-python` — frequentist tests for comparison
- `numerical-methods-bio` — optimization methods underlying MAP estimation
- `gwas-population-genetics` — Bayesian fine-mapping (SuSiE) uses similar posterior framework
