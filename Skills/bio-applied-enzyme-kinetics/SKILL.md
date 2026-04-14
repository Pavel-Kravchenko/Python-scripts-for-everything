---
name: bio-applied-enzyme-kinetics
description: Enzyme kinetics computational patterns — fitting Michaelis-Menten with scipy, bootstrap confidence intervals, inhibition type determination, allosteric cooperativity, and multi-substrate kinetics
tool_type: python
primary_tool: NumPy
---

# Enzyme Kinetics: Computational Patterns

## Fitting Michaelis-Menten Parameters

```python
import numpy as np
from scipy.optimize import curve_fit

def michaelis_menten(S, vmax, km):
    return (vmax * S) / (km + S)

# Substrate concentrations should span 0.1*Km to 10*Km
substrate_uM = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 50.0])

popt, pcov = curve_fit(
    michaelis_menten, substrate_uM, observed_v,
    p0=[max(observed_v), substrate_uM[np.argmin(np.abs(observed_v - max(observed_v)/2))]],
    bounds=([0, 0], [np.inf, np.inf]),
    maxfev=5000
)
vmax_fit, km_fit = popt
perr = np.sqrt(np.diag(pcov))   # standard errors
ci95 = 1.96 * perr               # 95% confidence intervals

kcat = vmax_fit / enzyme_conc_uM                  # turnover number (s⁻¹)
catalytic_efficiency = kcat / km_fit              # M⁻¹s⁻¹ (if units in µM, multiply by 1e6)
```

## Parametric Bootstrap for Confidence Bands

```python
n_boot = 500
residuals = observed_v - michaelis_menten(substrate_uM, *popt)
residual_sd = np.std(residuals)
S_smooth = np.linspace(0, substrate_uM.max() * 1.1, 300)

boot_curves = np.zeros((n_boot, len(S_smooth)))
boot_params = np.zeros((n_boot, 2))

for i in range(n_boot):
    v_boot = michaelis_menten(substrate_uM, *popt) + np.random.normal(0, residual_sd, len(substrate_uM))
    v_boot = np.clip(v_boot, 0, None)
    try:
        p_boot, _ = curve_fit(michaelis_menten, substrate_uM, v_boot,
                               p0=popt, bounds=([0,0],[np.inf,np.inf]), maxfev=2000)
        boot_curves[i] = michaelis_menten(S_smooth, p_boot[0], p_boot[1])
        boot_params[i] = p_boot
    except RuntimeError:
        boot_curves[i] = np.nan

ci_low  = np.nanpercentile(boot_curves, 2.5, axis=0)
ci_high = np.nanpercentile(boot_curves, 97.5, axis=0)
```

## Inhibition Models and Identification

```python
def competitive(S, vmax, km, I, Ki):
    """Km increases; Vmax unchanged. Same lines in Eadie-Hofstee, different x-intercept in LB."""
    return (vmax * S) / (km * (1 + I/Ki) + S)

def uncompetitive(S, vmax, km, I, Ki):
    """Both Vmax and Km decrease by same factor alpha. Parallel lines in Lineweaver-Burk."""
    alpha = 1 + I/Ki
    return (vmax/alpha * S) / (km/alpha + S)

def noncompetitive(S, vmax, km, I, Ki):
    """Vmax decreases; Km unchanged. Lines intersect on x-axis in Lineweaver-Burk."""
    return (vmax/(1 + I/Ki) * S) / (km + S)

def mixed_inhibition(S, vmax, km, I, Ki, alpha):
    """General case: alpha=1 -> noncompetitive; alpha->inf -> competitive."""
    return (vmax * S) / (alpha * km * (1 + I/Ki) + S)
```

### Inhibition Type Identification

| Observation | Inhibition type |
|-------------|----------------|
| Vmax unchanged, Km increases | Competitive |
| Vmax decreases, Km unchanged | Noncompetitive |
| Both Vmax and Km decrease equally | Uncompetitive |
| Both change, different factors | Mixed |
| LB: parallel lines | Uncompetitive |
| LB: lines converge on x-axis | Noncompetitive |
| LB: lines converge left of y-axis | Competitive or mixed |

## Allosteric Cooperativity (Hill Equation)

```python
def hill_equation(S, vmax, K_half, n):
    """n > 1: positive cooperativity; n < 1: negative; n = 1: hyperbolic (MM)."""
    return (vmax * S**n) / (K_half**n + S**n)

popt_hill, _ = curve_fit(hill_equation, substrate_uM, observed_v,
                          p0=[max(observed_v), np.median(substrate_uM), 1.5],
                          bounds=([0, 0, 0.1], [np.inf, np.inf, 10]))
vmax_h, K05, n_hill = popt_hill
# K_half = [S] at half-maximal velocity (equivalent to Km when n=1)
```

## Multi-Substrate Kinetics (Ping-Pong vs Sequential)

```python
# Ping-pong (parallel lines in double reciprocal at varied [B])
def ping_pong(S_A, S_B, vmax, Km_A, Km_B):
    return vmax / (1 + Km_A/S_A + Km_B/S_B)

# Sequential (ternary complex, intersecting lines in double reciprocal)
def sequential(S_A, S_B, vmax, Km_A, Km_B, Ki_A):
    return (vmax * S_A * S_B) / (Ki_A*Km_B + Km_B*S_A + Km_A*S_B + S_A*S_B)
```

## Goodness of Fit

```python
from scipy.stats import chi2

# Residuals
residuals = observed_v - michaelis_menten(substrate_uM, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((observed_v - observed_v.mean())**2)
r_squared = 1 - ss_res / ss_tot

# Chi-squared (if measurement variances are known)
chi2_stat = np.sum((residuals / measurement_sd)**2)
p_val = chi2.sf(chi2_stat, df=len(observed_v) - 2)

print(f"R² = {r_squared:.4f}")
print(f"Vmax = {vmax_fit:.2f} ± {ci95[0]:.2f} µM/s")
print(f"Km   = {km_fit:.2f} ± {ci95[1]:.2f} µM")
```

## Pitfalls

- **Bad initial guesses crash curve_fit**: estimate p0 from the data — Vmax ≈ max(observed_v) × 1.1, Km ≈ [S] where v ≈ Vmax/2; never use default p0=[1,1] for kinetic data
- **Bootstrap assumes residuals are homoscedastic**: enzyme kinetic noise is often proportional to velocity (heteroscedastic); in that case use weighted least squares (`sigma=observed_v`) in curve_fit
- **Hill n is sensitive to data range**: fitting the Hill equation requires data both below and above K_half; sparse data near inflection point gives unstable n estimates
- **Competitive vs mixed inhibition**: mixed inhibition with large alpha looks like competitive — measure at multiple inhibitor concentrations and use Dixon plot (1/v vs [I]) to distinguish
- **kcat requires total enzyme concentration**: convert Vmax to kcat only if [E]_total is accurately known from active-site titration, not total protein concentration (often contains inactive enzyme)
- **pcov returns inf when fit is underdetermined**: too few data points relative to parameters, or parameters are not independently identifiable from the data — add more substrate concentrations
