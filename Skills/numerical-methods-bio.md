---
name: numerical-methods-bio
description: Interpolation, curve fitting, optimization, and FFT for bioinformatics data analysis
---

# Numerical Methods for Bioinformatics

## When to Use
- Interpolating missing time points in gene expression, protein assay, or pharmacokinetic time series
- Fitting biological models (Michaelis-Menten, Hill, 4PL, two-compartment PK) to experimental measurements and comparing them with AIC
- Computing area under a curve (dose-response AUC, ROC AUC, pharmacokinetic AUC₀₋∞)
- Detecting periodic signals in time series (circadian rhythms, cell cycle oscillations) via FFT power spectrum

## Quick Reference

### Interpolation Methods
| Method | When to use | Stability | scipy / numpy |
|--------|------------|-----------|---------------|
| Lagrange polynomial | Small n, exact nodes, theory | Unstable at high degree on uniform grid | manual (see pattern below) |
| Cubic spline | Smooth biological curves, many nodes | Excellent | `CubicSpline` |
| Linear interp | Quick baseline; no smoothness needed | Good | `interp1d(kind='linear')` |
| Chebyshev nodes | High-degree polynomial, avoid Runge | Excellent | compute with `cos` formula |

Chebyshev node formula for n nodes on [a, b]:
`x_k = (b+a)/2 + (b-a)/2 * cos(pi*(2k-1)/(2n))`, k = 1..n.

### Differentiation Error Orders
| Formula | Error | Use when |
|---------|-------|----------|
| Forward  `(f(x+h) - f(x)) / h`           | O(h)   | One-sided boundary |
| Backward `(f(x) - f(x-h)) / h`           | O(h)   | One-sided boundary |
| Central  `(f(x+h) - f(x-h)) / (2h)`      | O(h²)  | Interior points — always prefer |

Optimal step size with noise epsilon: `h* ~ sqrt(epsilon / M2)`. Never use h smaller than this.

### Integration Rules
| Rule | Error | Notes |
|------|-------|-------|
| Trapezoidal | O(h²) | Works on irregular grids; use for ROC AUC |
| Simpson | O(h⁴) | Requires even number of intervals |
| `scipy.integrate.quad` | adaptive | Use for smooth analytical functions |

### Curve Fitting Functions
| Task | Function | Notes |
|------|----------|-------|
| Linear (overdetermined) | `np.linalg.lstsq` | Returns SVD solution; never form A^T A directly |
| Polynomial | `np.polyfit` | Vandermonde design matrix |
| Nonlinear | `scipy.optimize.curve_fit` | Levenberg-Marquardt; returns popt, pcov |
| MLE / general | `scipy.optimize.minimize` | Use Nelder-Mead (no gradient) or L-BFGS-B |

### Goodness-of-Fit
| Metric | Formula | Interpretation |
|--------|---------|----------------|
| R² | `1 - SS_res / SS_tot` | 1 = perfect; < 0 = worse than mean |
| AIC | `n*ln(SS/n) + 2*p` | Lower is better; penalises parameters |
| chi² / dof | `sum((y - yhat)^2 / sigma^2) / dof` | Should be ~1 for a good fit |

ΔAIC > 2: meaningful difference; ΔAIC > 10: decisive preference.

### FFT Recipes
| Task | Code |
|------|------|
| Forward FFT | `fft_vals = np.fft.rfft(signal)` |
| Frequency axis | `freqs = np.fft.rfftfreq(N, d=dt)` |
| Period axis | `periods = 1.0 / freqs[freqs > 0]` |
| Power spectrum | `power = np.abs(fft_vals)**2` |
| Inverse FFT | `np.fft.irfft(fft_vals, n=N)` |
| 2D FFT | `np.fft.fftshift(np.fft.fft2(image))` |

Nyquist frequency: `f_max = 1 / (2 * dt)`. Any signal above this aliases.

## Key Patterns

### Pattern 1: Cubic spline with periodic boundary (circadian data)
```python
from scipy.interpolate import CubicSpline
import numpy as np

t = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24], dtype=float)
y = np.array([1.0, 3.2, 7.1, 9.4, 6.8, 2.9, 1.1, 0.8, 1.0])
cs = CubicSpline(t, y, bc_type='periodic')
t_fine = np.linspace(0, 24, 500)
peak_t = t_fine[np.argmax(cs(t_fine))]
```

### Pattern 2: Trapezoidal AUC (non-uniform grid)
```python
def trapezoid_auc(x, y):
    return np.sum(0.5*(y[1:] + y[:-1]) * np.diff(x))

# ROC AUC, dose-response AUC, pharmacokinetic AUC all use this
auc = trapezoid_auc(fpr, tpr)
```

### Pattern 3: Michaelis-Menten nonlinear fit with 95% CI
```python
from scipy.optimize import curve_fit

def mm(S, Vmax, Km):
    return Vmax * S / (Km + S)

popt, pcov = curve_fit(mm, S, v, p0=[np.max(v), np.median(S)])
Vmax, Km = popt
Vmax_err, Km_err = np.sqrt(np.diag(pcov)) * 1.96  # 95% CI half-width
```

### Pattern 4: AIC-based model selection
```python
def aic(n_params, y_obs, y_pred):
    n = len(y_obs)
    sse = np.sum((y_obs - y_pred)**2)
    return n * np.log(sse / n) + 2 * n_params

# Compare two models; lower AIC wins
delta_aic = aic(2, v_obs, v_mm) - aic(3, v_obs, v_hill)
```

### Pattern 5: Gradient descent (linear regression)
```python
def gradient_descent(grad_Q, x0, lr=0.01, tol=1e-6, max_iter=5000):
    x = np.asarray(x0, float)
    for _ in range(max_iter):
        g = grad_Q(x)
        x = x - lr * g
        if np.linalg.norm(g) < tol:
            break
    return x
```

### Pattern 6: MLE via Nelder-Mead (negative binomial for RNA-seq)
```python
from scipy.optimize import minimize
from scipy.stats import nbinom

def neg_ll(params):
    r = np.exp(params[0])
    p = 1.0 / (1.0 + np.exp(-params[1]))
    return -np.sum(nbinom.logpmf(counts, r, p))

result = minimize(neg_ll, x0=[np.log(5), 0.0], method='Nelder-Mead')
r_mle = np.exp(result.x[0])
p_mle = 1.0 / (1.0 + np.exp(-result.x[1]))
```

### Pattern 7: FFT periodicity detection
```python
fft_vals = np.fft.rfft(signal)
freqs    = np.fft.rfftfreq(N, d=dt)          # cycles per time unit
power    = np.abs(fft_vals)**2
periods  = np.where(freqs > 0, 1.0/freqs, np.inf)
peak_period = periods[np.argmax(power[1:]) + 1]  # skip DC
```

### Pattern 8: Fourier low-pass filter
```python
fft_sig  = np.fft.rfft(noisy_signal)
freqs    = np.fft.rfftfreq(N, d=dt)
fft_filt = fft_sig.copy()
fft_filt[freqs > cutoff_freq] = 0           # zero out high frequencies
signal_clean = np.fft.irfft(fft_filt, n=N)
```

## Code Templates

### `cubic_spline_interpolate(t, y, t_query, periodic=False)`
```python
from scipy.interpolate import CubicSpline
import numpy as np

def cubic_spline_interpolate(t, y, t_query, periodic=False):
    bc = 'periodic' if periodic else 'not-a-knot'
    cs = CubicSpline(np.asarray(t, float), np.asarray(y, float), bc_type=bc)
    return cs(np.asarray(t_query, float)), cs

# Gene expression time course
t_obs  = [0, 2, 4, 8, 12, 24]
y_obs  = [1.0, 2.3, 4.1, 6.8, 5.2, 3.1]
y_pred, cs = cubic_spline_interpolate(t_obs, y_obs, [6.0, 18.0])
print(y_pred)  # interpolated values at t=6h and t=18h
```

### `curve_fit_with_ci(model, x, y, p0, sigma=None)`
```python
from scipy.optimize import curve_fit
import numpy as np

def curve_fit_with_ci(model, x, y, p0, sigma=None):
    popt, pcov = curve_fit(model, x, y, p0=p0, sigma=sigma,
                           absolute_sigma=(sigma is not None), maxfev=10000)
    perr = np.sqrt(np.diag(pcov))
    return {
        'popt': popt,
        'perr_1sigma': perr,
        'ci_95': perr * 1.96,
        'pcov': pcov
    }

# Michaelis-Menten
S = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
v = np.array([15, 25, 45, 60, 72, 80])
res = curve_fit_with_ci(lambda S, Vm, Km: Vm*S/(Km+S), S, v, p0=[90, 0.5])
print(f"Vmax = {res['popt'][0]:.1f} +/- {res['ci_95'][0]:.1f}")
```

### `fft_power_spectrum(signal, dt)`
```python
import numpy as np

def fft_power_spectrum(signal, dt):
    signal = np.asarray(signal, float)
    N      = len(signal)
    fft_v  = np.fft.rfft(signal)
    freqs  = np.fft.rfftfreq(N, d=dt)
    power  = np.abs(fft_v)**2
    periods = np.where(freqs[1:] > 0, 1.0 / freqs[1:], np.inf)
    peak_idx = np.argmax(power[1:])
    return {
        'freqs': freqs[1:], 'periods': periods,
        'power': power[1:], 'peak_period': periods[peak_idx],
        'peak_power': power[peak_idx+1]
    }

# Circadian clock gene (2h sampling, 72 points = 144h)
dt_h   = 2.0
t_h    = np.arange(72) * dt_h
signal = np.sin(2*np.pi*t_h/24) + 0.3*np.random.randn(72)
result = fft_power_spectrum(signal, dt=dt_h)
print(f"Dominant period: {result['peak_period']:.1f} hours")
```

### `minimize_mle(neg_ll_fn, x0, method='Nelder-Mead')`
```python
from scipy.optimize import minimize
import numpy as np

def minimize_mle(neg_ll_fn, x0, method='Nelder-Mead', tol=1e-8):
    result = minimize(neg_ll_fn, x0=x0, method=method,
                      options={'xatol': tol, 'fatol': tol, 'maxiter': 10000})
    if not result.success:
        raise RuntimeError(f"MLE did not converge: {result.message}")
    return result.x, result.fun

# MLE for normal distribution
data = np.random.normal(5.0, 1.5, 100)
def neg_ll(params):
    mu, log_s = params
    s = np.exp(log_s)
    return 0.5*len(data)*np.log(2*np.pi*s**2) + np.sum((data-mu)**2)/(2*s**2)

params_mle, nll = minimize_mle(neg_ll, x0=[4.0, 0.0])
print(f"MLE: mu={params_mle[0]:.3f}, sigma={np.exp(params_mle[1]):.3f}")
```

## Common Pitfalls

- **Runge phenomenon:** Never use high-degree Lagrange/polynomial interpolation on a uniform grid with more than ~6 points. Use cubic splines or Chebyshev nodes instead.
- **Differentiating noisy data:** h too small amplifies noise. Optimal step is h* ~ sqrt(epsilon/M2). Always use central differences on interior points.
- **Forming A^T A directly:** Condition number squares. Use `np.linalg.lstsq` (SVD-based) rather than solving normal equations with `np.linalg.solve(A.T @ A, A.T @ b)`.
- **Ignoring condition number:** Check `sv[0]/sv[-1]` from `lstsq`. If > 1e6 the fit is numerically unstable; reconsider the design matrix.
- **curve_fit without p0:** Default p0=[1,1,...] can miss the basin of attraction. Use domain knowledge: p0=[max(v), median(S)] for Michaelis-Menten.
- **AIC on different datasets:** AIC comparisons are only valid across models fitted to the same data with the same residual units.
- **DC component in FFT:** Index 0 is always the mean (DC). Skip it when searching for periodic peaks: `np.argmax(power[1:]) + 1`.
- **Aliasing:** Sampling at dt cannot detect periods shorter than 2*dt (Nyquist). If the biological period could be short, reduce sampling interval.
- **sigma parameterisation in MLE:** Always use log(sigma) or log(r) as the free parameter to enforce positivity; avoids optimizer wandering to negative values.
- **Nelder-Mead tolerance:** Default tolerances are loose. Set `xatol=1e-8, fatol=1e-8` for kinetics fitting.

## Related Skills
- `biochemistry-enzymology` — enzyme kinetics curve fitting (Michaelis-Menten, Hill, inhibition models)
- `ml-deep-learning-bio` — gradient descent and optimisation for neural networks
- `probability-statistics-python` — statistical distributions, hypothesis testing, MLE foundations
- `rna-seq-differential-expression` — negative binomial distribution for read count modelling
