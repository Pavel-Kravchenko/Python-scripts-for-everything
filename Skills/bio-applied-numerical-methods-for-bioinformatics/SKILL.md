---
name: bio-applied-numerical-methods-for-bioinformatics
description: - Implement Lagrange and Newton interpolation polynomials and explain the Runge phenomenon - Apply cubic spline interpolation to reconstruct missing time points in biological time series - Compute num
tool_type: python
primary_tool: NumPy
---

# Numerical Methods for Bioinformatics

## Part 1: Interpolation

Interpolation builds a function passing exactly through $(x_1,y_1),\ldots,(x_n,y_n)$. Different from curve fitting — no deviation allowed.

**Weierstrass (1885):** Any continuous function on [a,b] can be approximated to arbitrary precision by a polynomial. But *which* polynomial matters enormously.

### Lagrange Interpolation

$$\pi_n(x) = \sum_i y_i \ell_i(x), \quad \ell_i(x) = \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}$$

Error bound for $|f^{(n)}| \leq M_n$: $|f(x) - \pi_n(x)| = \frac{f^{(n)}(\xi)}{n!} \prod_i (x - x_i)$

```python
def lagrange_interpolate(x_nodes, y_nodes, x_query):
    x_nodes = np.asarray(x_nodes, float); y_nodes = np.asarray(y_nodes, float)
    x_query = np.asarray(x_query, float); n = len(x_nodes)
    result = np.zeros_like(x_query)
    for i in range(n):
        basis = np.ones_like(x_query)
        for j in range(n):
            if j != i: basis *= (x_query - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
        result += y_nodes[i] * basis
    return result

# Mendeleev 1881 solubility data
T_data = np.array([0, 4, 10, 15, 21, 29, 36, 51, 68], dtype=float)
S_data = np.array([66.7, 71.0, 76.3, 80.6, 85.7, 99.4, 99.4, 113.6, 125.1])
T_fine = np.linspace(0, 68, 300)
S_lag = lagrange_interpolate(T_data, S_data, T_fine)
print(f'Solubility at 40C: {lagrange_interpolate(T_data, S_data, [40])[0]:.1f} g')
```

### Runge Phenomenon and Chebyshev Nodes

High-degree Lagrange on a **uniform grid** suffers from oscillations near interval endpoints (Runge 1901: $f(x)=1/(1+25x^2)$, degree 20). Solution: **Chebyshev nodes** minimize $\max_x |\omega_n(x, X)|$:

$$x_k = \frac{b+a}{2} + \frac{b-a}{2} \cos\!\left(\frac{\pi(2k-1)}{2n}\right), \quad k=1,\ldots,n$$

```python
def chebyshev_nodes(n, a=-1.0, b=1.0):
    k = np.arange(1, n+1)
    return 0.5*(b+a) + 0.5*(b-a)*np.cos(np.pi*(2*k-1)/(2*n))

runge_fn = lambda x: 1.0 / (1 + 25*x**2)
x_plot = np.linspace(-1, 1, 400); n_pts = 12
x_unif = np.linspace(-1, 1, n_pts); x_cheb = chebyshev_nodes(n_pts)
y_unif = lagrange_interpolate(x_unif, runge_fn(x_unif), x_plot)
y_cheb = lagrange_interpolate(x_cheb, runge_fn(x_cheb), x_plot)
```

### Newton's Divided Differences

Same polynomial as Lagrange but incremental — adding a new data point requires only one extra term:

$$\pi_n(x) = f[x_1] + f[x_1,x_2](x-x_1) + f[x_1,x_2,x_3](x-x_1)(x-x_2) + \cdots$$

Divided differences defined recursively: $f[x_i,x_{i+1}] = (f(x_{i+1})-f(x_i))/(x_{i+1}-x_i)$

```python
def newton_divided_differences(x, y):
    n = len(x); table = np.zeros((n, n)); table[:, 0] = y
    for j in range(1, n):
        for i in range(n - j):
            table[i, j] = (table[i+1, j-1] - table[i, j-1]) / (x[i+j] - x[i])
    return table[0, :]

def newton_interp(coeffs, x_nodes, x_query):
    x_query = np.asarray(x_query, float); result = np.full_like(x_query, coeffs[-1])
    for i in range(len(coeffs)-2, -1, -1):
        result = result * (x_query - x_nodes[i]) + coeffs[i]
    return result

# Gene expression: interpolate missing time points
t_measured = np.array([0, 2, 4, 8, 12, 24], dtype=float)
expr = np.array([1.0, 2.3, 4.1, 6.8, 5.2, 3.1])
coeffs = newton_divided_differences(t_measured, expr)
missing_t = np.array([6.0, 18.0])
expr_pred = newton_interp(coeffs, t_measured, missing_t)
print(f'Predicted: t=6h: {expr_pred[0]:.2f}, t=18h: {expr_pred[1]:.2f}')
```

### Cubic Spline Interpolation

High-degree polynomials are unstable. **Cubic splines** use piecewise cubics with continuous 1st and 2nd derivatives — solving the variational problem:

$$\min_f \int_a^b [f''(x)]^2 \, dx \quad \text{subject to } f(x_i) = y_i$$

Error bound for $f \in W^4[M_4, I]$: $\|s - f\|_\infty \leq \frac{5}{364} h^4 M_4$ where $h = \max(x_{i+1}-x_i)$.

```python
from scipy.interpolate import CubicSpline, interp1d

# Circadian Per1 expression with periodic boundary
t_sampled = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24], dtype=float)
per1_expr = np.array([1.0, 3.2, 7.1, 9.4, 6.8, 2.9, 1.1, 0.8, 1.0])
cs = CubicSpline(t_sampled, per1_expr, bc_type='periodic')
t_fine = np.linspace(0, 24, 500)
peak_t = t_fine[np.argmax(cs(t_fine))]
print(f'Estimated peak at ZT {peak_t:.1f}h')
```

## Part 2: Numerical Differentiation and Integration

### Finite Difference Formulas

| Formula | Expression | Error |
|---------|-----------|-------|
| Forward | $(f(x+h) - f(x)) / h$ | $O(h)$ |
| Backward | $(f(x) - f(x-h)) / h$ | $O(h)$ |
| Central | $(f(x+h) - f(x-h)) / (2h)$ | $O(h^2)$ |

**Instability with noisy data:** Total error $\approx M_2 h + \varepsilon/h$, minimized at $h^* \sim \sqrt{\varepsilon/M_2}$. Taking $h$ too small amplifies noise — use central differences and moderate $h$.

```python
def forward_diff(f, x, h=1e-5):  return (f(x+h) - f(x)) / h
def central_diff(f, x, h=1e-5):  return (f(x+h) - f(x-h)) / (2*h)

# Drug clearance: dC/dt at interior points via central difference
t_pk = np.array([0, 1, 2, 4, 6, 8, 12, 24], dtype=float)
C_pk = np.array([100, 72, 52, 27, 14, 7.4, 2.0, 0.1])
for i in range(1, len(t_pk)-1):
    dCdt = (C_pk[i+1] - C_pk[i-1]) / (t_pk[i+1] - t_pk[i-1])
    print(f't={t_pk[i]:4.0f}h: dC/dt = {dCdt:6.2f} ng/mL/h')
```

### Numerical Integration

Integration is **stable** — measurement errors average out over the interval.

**Trapezoidal rule** ($f \in W^2$):
$$J_N^T = \frac{b-a}{N}\left[\frac{f(a)+f(b)}{2} + \sum_{k=1}^{N-1} f(x_k)\right], \quad |J - J_N^T| \leq \frac{(b-a)^3}{12N^2} M_2$$

**Simpson's rule** ($f \in W^4$, converges $O(h^4)$ vs trapezoidal $O(h^2)$):
$$J_N^S = \frac{b-a}{6N}\left[f(a) + f(b) + 4\sum_\text{odd} f(x_k) + 2\sum_\text{even} f(x_k)\right], \quad |J - J_N^S| \leq \frac{(b-a)^5}{180N^4} M_4$$

## Pitfalls

- **Coordinate systems**: BED 0-based half-open; VCF/GFF 1-based inclusive — mixing causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) for thousands of features
