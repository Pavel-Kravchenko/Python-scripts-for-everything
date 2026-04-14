---
name: bio-applied-biochemistry
description: Biochemistry fundamentals for computational biologists — Beer-Lambert law, spectrophotometric assays, Michaelis-Menten kinetics, linearization methods, and inhibition models
tool_type: python
primary_tool: NumPy
---

# Biochemistry: Assays and Enzyme Kinetics

## Beer-Lambert Law

A = ε · l · c

| Symbol | Meaning | Typical value |
|--------|---------|---------------|
| A | Absorbance (dimensionless) | 0.1–1.5 (linear range) |
| ε | Molar absorption coefficient (M⁻¹ cm⁻¹) | NADH at 340 nm: 6,220 |
| l | Path length (cm) | 1 cm (standard cuvette) |
| c | Concentration (M) | — |

```python
# Beer-Lambert: absorbance -> concentration
epsilon_NADH = 6220   # M^-1 cm^-1 at 340 nm
path_length = 1       # cm

def abs_to_conc(A, eps=epsilon_NADH, l=path_length):
    return A / (eps * l)   # returns molar concentration

# A=0.622 -> 100 µM NADH
```

## Extracting Initial Velocity from Progress Curve

Measure the slope in the **linear region** (first ~10% of substrate consumed) before product inhibition or substrate depletion distort the rate.

```python
import numpy as np

# Fit linear region: product_uM < 10% of P_max
linear_mask = product_uM < 0.10 * P_max_uM
slope, intercept = np.polyfit(time_s[linear_mask], product_uM[linear_mask], 1)
v0 = slope   # µM/s — the initial velocity
```

## Michaelis-Menten Equation

v = (Vmax · [S]) / (Km + [S])

| Parameter | Meaning |
|-----------|---------|
| Vmax | Maximum velocity = kcat · [E]_total |
| Km | [S] at half-maximal velocity = (k₋₁ + k₂) / k₁ |
| kcat | Turnover number = Vmax / [E]_total |
| kcat/Km | Catalytic efficiency; diffusion limit ~10⁸–10⁹ M⁻¹s⁻¹ |

Km ≈ Ks (true dissociation constant) only when k₂ << k₋₁ (rapid equilibrium assumption).

```python
from scipy.optimize import curve_fit

def michaelis_menten(S, vmax, km):
    return (vmax * S) / (km + S)

popt, pcov = curve_fit(
    michaelis_menten, substrate_uM, observed_v,
    p0=[80.0, 3.0],
    bounds=([0, 0], [np.inf, np.inf]),
    maxfev=5000
)
vmax_fit, km_fit = popt
se = np.sqrt(np.diag(pcov))    # standard errors
ci_95 = 1.96 * se              # 95% CI
```

## Linearization Methods Comparison

| Method | Plot | Slope | Intercept | Error distortion |
|--------|------|-------|-----------|-----------------|
| Lineweaver-Burk | 1/v vs 1/[S] | Km/Vmax | 1/Vmax | Severe (amplifies low-[S] errors) |
| Eadie-Hofstee | v vs v/[S] | −Km | Vmax | Moderate |
| Hanes-Woolf | [S]/v vs [S] | 1/Vmax | Km/Vmax | Most uniform |

```python
# Lineweaver-Burk
lb_slope, lb_int = np.polyfit(1/substrate_uM, 1/observed_v, 1)
lb_vmax = 1.0 / lb_int
lb_km = lb_slope * lb_vmax

# Hanes-Woolf (best of the three linear methods)
hw_slope, hw_int = np.polyfit(substrate_uM, substrate_uM/observed_v, 1)
hw_vmax = 1.0 / hw_slope
hw_km = hw_int * hw_vmax

# Always prefer nonlinear regression (curve_fit) over these linearizations
```

## Inhibition Models

```python
def competitive_inhibition(S, vmax, km, I, Ki):
    """Competitive: inhibitor raises apparent Km; Vmax unchanged."""
    km_app = km * (1 + I / Ki)
    return (vmax * S) / (km_app + S)

def uncompetitive_inhibition(S, vmax, km, I, Ki):
    """Uncompetitive: both Vmax and Km reduced equally; v/[S] lines parallel in LB."""
    alpha = 1 + I / Ki
    return (vmax * S / alpha) / (km / alpha + S)

def noncompetitive_inhibition(S, vmax, km, I, Ki):
    """Noncompetitive: Vmax reduced; Km unchanged."""
    vmax_app = vmax / (1 + I / Ki)
    return (vmax_app * S) / (km + S)
```

## Experimental Design for Km Estimation

- Use substrate concentrations spanning **0.1·Km to 10·Km** (if Km unknown, pilot with wide range)
- Minimum 8–10 concentrations; include duplicates/triplicates
- Measure initial velocities only (linear phase); avoid >10% substrate depletion
- Include negative control (no enzyme) and blank (no substrate)

## Pitfalls

- **Km ≠ affinity unless rapid equilibrium holds**: Km = (k₋₁ + k₂)/k₁; only when k₂ << k₋₁ does Km ≈ Ks; report Km, not "binding affinity"
- **Lineweaver-Burk distorts errors**: low [S] points (right side of plot) dominate the fit because reciprocal amplifies small-velocity noise — use nonlinear regression as primary method
- **curve_fit requires good initial guesses**: bad p0 causes convergence to wrong local minimum; plot data first, estimate Vmax as observed maximum, Km as [S] at ~v_max/2
- **Absorbance must be in linear range**: Beer-Lambert is only linear for A < 1.5; dilute samples if A > 1.5 (inner filter effect distorts fluorescence assays at A > 0.1)
- **Progress curve non-linearity**: if product inhibition or substrate depletion occurs early, the "linear" region may be only the first few time points — verify with varying enzyme concentration
- **Units**: keep [S] and Km in the same units (µM vs mM mismatch is a common bug); Vmax units must match velocity axis
