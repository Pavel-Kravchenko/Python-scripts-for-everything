---
name: bio-applied-enzyme-kinetics
description: "- Describe the enzyme-substrate reaction mechanism and extract initial velocities from progress curves - Derive and fit the Michaelis-Menten equation to obtain Km and Vmax with confidence intervals - "
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/13_Biochemistry_and_Enzyme_Kinetics/02_enzyme_kinetics.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Biochemistry and Enzyme Kinetics

*Source: Course notebook `Tier_3_Applied_Bioinformatics/13_Biochemistry_and_Enzyme_Kinetics/02_enzyme_kinetics.ipynb`*

# Biochemistry and Enzyme Kinetics

## Learning Objectives

By the end of this module you will be able to:

- Describe the enzyme-substrate reaction mechanism and extract initial velocities from progress curves
- Derive and fit the Michaelis-Menten equation to obtain Km and Vmax with confidence intervals
- Apply and critically compare linearization methods: Lineweaver-Burk, Eadie-Hofstee, Hanes-Woolf
- Simulate and distinguish competitive, non-competitive, uncompetitive, and mixed inhibition
- Fit the Hill equation and interpret the cooperativity coefficient for allosteric enzymes
- Build stoichiometric matrices for metabolic pathways and understand flux balance analysis

---
## 1. Introduction to Enzyme Kinetics

Enzymes are biological catalysts that accelerate chemical reactions without being consumed. They lower the **activation energy** of a reaction by binding substrates in an active site and stabilising the transition state.

The overall catalytic cycle:

```
E + S  ⇌  ES  →  E + P
        k₁      k₂
        k₋₁
```

| Symbol | Meaning |
|--------|------------------------------------------|
| E      | Free enzyme |
| S      | Substrate |
| ES     | Enzyme-substrate (Michaelis) complex |
| P      | Product |
| k₁     | Substrate binding rate |
| k₋₁    | Substrate dissociation rate |
| k₂     | Catalytic rate (kcat) |

**Reaction velocity** (v) is the rate of product formation: v = d[P]/dt.

### 1.1 Spectrophotometric Assays

Most enzyme assays follow absorbance over time. At a fixed wavelength, absorbance A is proportional to the concentration of a chromophoric species (substrate or product).

**Beer-Lambert law:**

$$A = \varepsilon \cdot l \cdot c$$

where:
- A = absorbance (dimensionless)
- ε = molar absorption coefficient (M⁻¹ cm⁻¹)
- l = path length (cm), typically 1 cm
- c = concentration (M)

The **initial velocity** (v₀) is measured from the linear portion of the progress curve (before substrate depletion or product inhibition occurs). It is the slope of [P] vs time at t ≈ 0.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.linalg import null_space
import warnings
warnings.filterwarnings('ignore')

rng = np.random.default_rng(seed=42)

# Beer-Lambert: convert absorbance to concentration
epsilon = 6220   # NADH molar absorption coefficient at 340 nm (M^-1 cm^-1)
path_length = 1  # cm

def absorbance_to_concentration(absorbance, eps=epsilon, l=path_length):
    """Convert absorbance to molar concentration via Beer-Lambert."""
    return absorbance / (eps * l)

print(f"A = 0.622 corresponds to [{absorbance_to_concentration(0.622)*1e6:.1f}] µM NADH")
```

### 1.2 Extracting Initial Velocity from a Progress Curve

In practice, the linear region of a progress curve spans the first 5–10% of substrate consumption. We fit a straight line to that region.

```python
# Simulate a product progress curve (molar concentration in µM)
time_s = np.linspace(0, 120, 200)  # seconds

# Simple saturating product accumulation: [P](t) = Pmax * (1 - exp(-t/tau))
P_max_uM = 80.0   # µM total substrate available
tau = 40.0         # time constant (s)
product_uM = P_max_uM * (1 - np.exp(-time_s / tau))
product_uM += rng.normal(0, 0.5, size=len(time_s))  # add measurement noise

# Extract initial velocity: fit linear region (first 10% of Pmax)
linear_mask = product_uM < 0.10 * P_max_uM
if linear_mask.sum() < 3:
    linear_mask = np.arange(len(time_s)) < 15

slope, intercept = np.polyfit(time_s[linear_mask], product_uM[linear_mask], 1)
v0_uM_per_s = slope  # µM/s

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(time_s, product_uM, color='steelblue', lw=1.5, label='Progress curve')
ax.plot(time_s[linear_mask],
        slope * time_s[linear_mask] + intercept,
        color='tomato', lw=2, linestyle='--', label=f'Linear fit (v₀ = {v0_uM_per_s:.2f} µM/s)')
ax.set_xlabel('Time (s)')
ax.set_ylabel('[Product] (µM)')
ax.set_title('Product Progress Curve and Initial Velocity')
ax.legend()
plt.tight_layout()
plt.show()

print(f"Estimated initial velocity v₀ = {v0_uM_per_s:.3f} µM/s")
```

---
## 2. Michaelis-Menten Kinetics

Assuming **quasi-steady-state** for [ES] (d[ES]/dt ≈ 0, Briggs & Haldane 1925), the Michaelis-Menten equation describes reaction velocity as a function of substrate concentration:

$$v = \frac{V_{\max} \cdot [S]}{K_m + [S]}$$

| Parameter | Meaning |
|-----------|------------------------------------------------------------|
| Vmax      | Maximum velocity (all enzyme sites saturated) = kcat·[E]total |
| Km        | Michaelis constant = (k₋₁ + k₂)/k₁; equals [S] at which v = Vmax/2 |
| kcat      | Turnover number = Vmax / [E]total (catalytic rate per enzyme molecule) |

**Km and affinity:** Km = (k₋₁ + k₂)/k₁. Only when k₂ << k₋₁ (rapid equilibrium) does Km ≈ Ks = k₋₁/k₁ (true substrate dissociation constant). In general:
- **Low Km** → half-maximal velocity at low [S] (effective affinity is high)
- **High Km** → requires more substrate to reach half-maximal velocity

**kcat/Km (catalytic efficiency):** The second-order rate constant for substrate capture at low [S]. Approaches the diffusion limit (~10⁸–10⁹ M⁻¹s⁻¹ for evolutionarily optimized enzymes like triosephosphate isomerase).

At low [S]: v ≈ (Vmax/Km)·[S]  — first-order in [S]
At high [S]: v ≈ Vmax           — zero-order; enzyme is saturated

```python
def michaelis_menten(substrate_conc, v_max, k_m):
    """Michaelis-Menten velocity."""
    return (v_max * substrate_conc) / (k_m + substrate_conc)

# True kinetic parameters
true_vmax = 100.0  # µM/s
true_km   = 5.0    # µM

# Substrate concentrations spanning 0.1*Km to 10*Km (good experimental design)
substrate_uM = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 50.0])

# Generate synthetic velocities with 5% Gaussian noise
true_velocities = michaelis_menten(substrate_uM, true_vmax, true_km)
noise_sd = 0.05 * true_velocities
observed_velocities = true_velocities + rng.normal(0, noise_sd)
observed_velocities = np.clip(observed_velocities, 0, None)

print("Substrate [µM] | Observed velocity [µM/s]")
print("-" * 40)
for s, v in zip(substrate_uM, observed_velocities):
    print(f"  {s:6.1f}       |  {v:8.3f}")
```

### 2.1 Fitting Km and Vmax with `scipy.optimize.curve_fit`

`curve_fit` uses nonlinear least squares (Levenberg-Marquardt). It returns the best-fit parameters and the covariance matrix, from which we compute standard errors and 95% confidence intervals.

```python
# Fit Michaelis-Menten model
p0 = [80.0, 3.0]  # initial guesses: [vmax, km]
bounds = ([0, 0], [np.inf, np.inf])

popt, pcov = curve_fit(
    michaelis_menten,
    substrate_uM,
    observed_velocities,
    p0=p0,
    bounds=bounds,
    maxfev=5000
)

fitted_vmax, fitted_km = popt
perr = np.sqrt(np.diag(pcov))  # standard errors

# 95% CI: parameter ± 1.96 * SE
ci_vmax = 1.96 * perr[0]
ci_km   = 1.96 * perr[1]

print(f"Fitted Vmax = {fitted_vmax:.2f} ± {ci_vmax:.2f} µM/s  (true: {true_vmax})")
print(f"Fitted Km   = {fitted_km:.2f}  ± {ci_km:.2f}  µM    (true: {true_km})")
```

```python
# Plot Michaelis-Menten curve with data
s_smooth = np.linspace(0, 55, 300)
v_fit    = michaelis_menten(s_smooth, fitted_vmax, fitted_km)

fig, ax = plt.subplots(figsize=(7, 4))
ax.scatter(substrate_uM, observed_velocities,
           color='steelblue', zorder=5, s=60, label='Observed data')
ax.plot(s_smooth, v_fit, color='tomato', lw=2,
        label=f'MM fit: Vmax={fitted_vmax:.1f}, Km={fitted_km:.2f} µM')
ax.axhline(fitted_vmax, color='gray', linestyle=':', lw=1)
ax.axhline(fitted_vmax / 2, color='gray', linestyle=':', lw=1)
ax.axvline(fitted_km, color='gray', linestyle=':', lw=1)
ax.annotate('Vmax', xy=(55, fitted_vmax), va='bottom', ha='right', color='gray')
ax.annotate('Vmax/2', xy=(55, fitted_vmax / 2), va='bottom', ha='right', color='gray')
ax.annotate('Km', xy=(fitted_km, 2), va='bottom', ha='left', color='gray')
ax.set_xlabel('[S] (µM)')
ax.set_ylabel('v (µM/s)')
ax.set_title('Michaelis-Menten Kinetics')
ax.legend()
plt.tight_layout()
plt.show()
```

### 2.2 Confidence Intervals via Parametric Bootstrap

The covariance matrix from `curve_fit` assumes Gaussian errors and linearity near the optimum. A parametric bootstrap gives a distribution over the fitted curve, making confidence bands visible.

```python
# Parametric bootstrap: resample from fitted model + residual noise
n_bootstrap = 500
residuals = observed_velocities - michaelis_menten(substrate_uM, fitted_vmax, fitted_km)
residual_sd = np.std(residuals)

bootstrap_curves = np.zeros((n_bootstrap, len(s_smooth)))
bootstrap_params = np.zeros((n_bootstrap, 2))

for i in range(n_bootstrap):
    v_boot = michaelis_menten(substrate_uM, fitted_vmax, fitted_km) + rng.normal(0, residual_sd, size=len(substrate_uM))
    v_boot = np.clip(v_boot, 0, None)
    try:
        p_boot, _ = curve_fit(michaelis_menten, substrate_uM, v_boot, p0=popt, bounds=bounds, maxfev=2000)
        bootstrap_curves[i] = michaelis_menten(s_smooth, *p_boot)
        bootstrap_params[i] = p_boot
    except RuntimeError:
        bootstrap_curves[i] = np.nan

ci_low  = np.nanpercentile(bootstrap_curves, 2.5, axis=0)
ci_high = np.nanpercentile(bootstrap_curves, 97.5, axis=0)

fig, ax = plt.subplots(figsize=(7, 4))
ax.fill_between(s_smooth, ci_low, ci_high, alpha=0.25, color='tomato', label='95% CI (bootstrap)')
ax.plot(s_smooth, v_fit, color='tomato', lw=2, label='Best fit')
ax.scatter(substrate_uM, observed_velocities, color='steelblue', zorder=5, s=60, label='Data')
ax.set_xlabel('[S] (µM)')
ax.set_ylabel('v (µM/s)')
ax.set_title('Michaelis-Menten Fit with 95% Confidence Band')
ax.legend()
plt.tight_layout()
plt.show()

print(f"Bootstrap Vmax: {np.nanmean(bootstrap_params[:,0]):.2f} ± {np.nanstd(bootstrap_params[:,0]):.2f} µM/s")
print(f"Bootstrap Km:   {np.nanmean(bootstrap_params[:,1]):.2f} ± {np.nanstd(bootstrap_params[:,1]):.2f} µM")
```

---
## 3. Linearization Methods

Before nonlinear regression was widely available, researchers linearised the Michaelis-Menten equation. Each transformation distorts errors differently.

### 3.1 Lineweaver-Burk (Double Reciprocal)

Take reciprocals of both sides of MM:

$$\frac{1}{v} = \frac{K_m}{V_{\max}} \cdot \frac{1}{[S]} + \frac{1}{V_{\max}}$$

Plot 1/v vs 1/[S]: slope = Km/Vmax, y-intercept = 1/Vmax, x-intercept = -1/Km.

**Problem:** Small velocities (low [S]) are highly amplified in reciprocal space, severely distorting error structure.

### 3.2 Eadie-Hofstee

$$v = V_{\max} - K_m \cdot \frac{v}{[S]}$$

Plot v vs v/[S]: slope = -Km, y-intercept = Vmax.

### 3.3 Hanes-Woolf

$$\frac{[S]}{v} = \frac{K_m}{V_{\max}} + \frac{1}{V_{\max}} \cdot [S]$$

Plot [S]/v vs [S]: slope = 1/Vmax, y-intercept = Km/Vmax. **Most uniform error distribution** of the three.

```python
# Compute all three linearizations
inv_s = 1.0 / substrate_uM
inv_v = 1.0 / observed_velocities
v_over_s = observed_velocities / substrate_uM
s_over_v = substrate_uM / observed_velocities

# Lineweaver-Burk fit
lb_slope, lb_intercept = np.polyfit(inv_s, inv_v, 1)
lb_vmax = 1.0 / lb_intercept
lb_km   = lb_slope * lb_vmax

# Eadie-Hofstee fit
eh_slope, eh_intercept = np.polyfit(v_over_s, observed_velocities, 1)
eh_km   = -eh_slope
eh_vmax = eh_intercept

# Hanes-Woolf fit
hw_slope, hw_intercept = np.polyfit(substrate_uM, s_over_v, 1)
hw_vmax = 1.0 / hw_slope
hw_km   = hw_intercept * hw_vmax

print(f"{'Method':<20} {'Vmax (µM/s)':>12} {'Km (µM)':>10}")
print("-" * 45)
print(f"{'True values':<20} {true_vmax:>12.2f} {true_km:>10.2f}")
print(f"{'Nonlinear (NLS)':<20} {fitted_vmax:>12.2f} {fitted_km:>10.2f}")
print(f"{'Lineweaver-Burk':<20} {lb_vmax:>12.2f} {lb_km:>10.2f}")
print(f"{'Eadie-Hofstee':<20} {eh_vmax:>12.2f} {eh_km:>10.2f}")
print(f"{'Hanes-Woolf':<20} {hw_vmax:>12.2f} {hw_km:>10.2f}")
```

```python
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

# --- Lineweaver-Burk ---
ax = axes[0]
x_lb = np.linspace(min(inv_s) * 0.8, max(inv_s) * 1.1, 100)
ax.scatter(inv_s, inv_v, color='steelblue', zorder=5, s=50)
ax.plot(x_lb, lb_slope * x_lb + lb_intercept, color='tomato', lw=2)
ax.axvline(0, color='k', lw=0.5)
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('1/[S] (µM⁻¹)')
ax.set_ylabel('1/v (s/µM)')
ax.set_title('Lineweaver-Burk')

# --- Eadie-Hofstee ---
ax = axes[1]
x_eh = np.linspace(min(v_over_s) * 0.8, max(v_over_s) * 1.1, 100)
ax.scatter(v_over_s, observed_velocities, color='steelblue', zorder=5, s=50)
ax.plot(x_eh, eh_slope * x_eh + eh_intercept, color='tomato', lw=2)
ax.set_xlabel('v/[S] (s⁻¹)')
ax.set_ylabel('v (µM/s)')
ax.set_title('Eadie-Hofstee')

# --- Hanes-Woolf ---
ax = axes[2]
x_hw = np.linspace(min(substrate_uM) * 0.8, max(substrate_uM) * 1.1, 100)
ax.scatter(substrate_uM, s_over_v, color='steelblue', zorder=5, s=50)
ax.plot(x_hw, hw_slope * x_hw + hw_intercept, color='tomato', lw=2)
ax.set_xlabel('[S] (µM)')
ax.set_ylabel('[S]/v (s)')
ax.set_title('Hanes-Woolf')

plt.suptitle('Linearization Methods Compared', fontsize=13, y=1.02)
plt.tight_layout()
plt.show()
```
