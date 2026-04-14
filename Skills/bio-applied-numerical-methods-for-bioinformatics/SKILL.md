---
name: bio-applied-numerical-methods-for-bioinformatics
description: - Implement Lagrange and Newton interpolation polynomials and explain the Runge phenomenon - Apply cubic spline interpolation to reconstruct missing time points in biological time series - Compute num
tool_type: python
primary_tool: NumPy
---

# Numerical Methods for Bioinformatics

## Part 1: Interpolation

Interpolation builds a function that passes exactly through a set of data points (x_1, y_1), ..., (x_n, y_n). This is different from curve fitting, where the function is allowed to deviate from the data to avoid overfitting noise.

**Key theoretical result (Weierstrass, 1885):** Any continuous function on [a, b] can be approximated to arbitrary precision by a polynomial. However, *which* polynomial matters enormously.

### Lagrange Interpolation

Given n distinct nodes x_1, ..., x_n, the **Lagrange interpolation polynomial** of degree n-1 is:

    pi_n(x) = sum_i y_i * l_i(x),   l_i(x) = prod_{j != i} (x - x_j) / (x_i - x_j)

Each basis polynomial l_i equals 1 at x_i and 0 at all other nodes, so pi_n(x_i) = y_i is guaranteed.

**Error bound:** For f in W^n[M_n, I] (i.e., |f^(n)| <= M_n on I = [a,b]):

    |f(x) - pi_n(x)| = f^(n)(xi) / n! * prod_i (x - x_i),   xi in [a, b]

```python
def lagrange_interpolate(x_nodes, y_nodes, x_eval):
    """Evaluate the Lagrange interpolation polynomial at x_eval."""
    x_nodes = np.asarray(x_nodes, float)
    y_nodes = np.asarray(y_nodes, float)
    x_eval  = np.asarray(x_eval, float)
    n = len(x_nodes)
    result = np.zeros_like(x_eval)
    for i in range(n):
        basis = np.ones_like(x_eval)
        for j in range(n):
            if j != i:
                basis *= (x_eval - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
        result += y_nodes[i] * basis
    return result


# Mendeleev 1881 solubility data: NaNO2 dissolved in 100 mL water vs temperature
T_data = np.array([0, 4, 10, 15, 21, 29, 36, 51, 68], dtype=float)
S_data = np.array([66.7, 71.0, 76.3, 80.6, 85.7, 99.4, 99.4, 113.6, 125.1])

T_fine = np.linspace(0, 68, 300)
S_lag  = lagrange_interpolate(T_data, S_data, T_fine)

fig, ax = plt.subplots(figsize=(8, 4))
ax.scatter(T_data, S_data, zorder=5, color='black', label='Mendeleev data')
ax.plot(T_fine, S_lag, label='Lagrange polynomial (degree 8)')
ax.set_xlabel('Temperature (C)')
ax.set_ylabel('NaNO2 per 100 mL water (g)')
ax.set_title('Lagrange interpolation of solubility data')
ax.legend()
plt.tight_layout()
plt.show()
print(f'Interpolated solubility at 40 C: {lagrange_interpolate(T_data, S_data, [40])[0]:.1f} g')
```python

### The Runge Phenomenon and Chebyshev Nodes

High-degree Lagrange interpolation on a **uniform grid** suffers from the **Runge phenomenon**: oscillations near the interval endpoints grow exponentially with the polynomial degree. Runge (1901) demonstrated this on f(x) = 1/(1+25x^2) with a degree-20 polynomial.


    x_k = (b+a)/2 + (b-a)/2 * cos(pi*(2k-1)/(2n)),   k = 1, ..., n

These are zeros of the Chebyshev polynomial T_n(y) = cos(n * arccos(y)), and they minimise max_x |omega_n(x, X)| over all n-point grids.

```python
def chebyshev_nodes(n, a=-1.0, b=1.0):
    """Return n Chebyshev nodes on [a, b]."""
    k = np.arange(1, n+1)
    return 0.5*(b+a) + 0.5*(b-a)*np.cos(np.pi*(2*k-1)/(2*n))


runge_fn = lambda x: 1.0 / (1 + 25*x**2)
x_plot = np.linspace(-1, 1, 400)

n_pts   = 12
x_unif  = np.linspace(-1, 1, n_pts)
x_cheb  = chebyshev_nodes(n_pts)

y_true  = runge_fn(x_plot)
y_unif  = lagrange_interpolate(x_unif, runge_fn(x_unif), x_plot)
y_cheb  = lagrange_interpolate(x_cheb, runge_fn(x_cheb), x_plot)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for ax, y_interp, nodes, title in zip(
    axes,
    [y_unif, y_cheb],
    [x_unif, x_cheb],
    ['Uniform grid (Runge phenomenon)', 'Chebyshev nodes (well-behaved)']
):
    ax.plot(x_plot, y_true, 'k-', lw=2, label='True f(x)')
    ax.plot(x_plot, y_interp, 'b--', label='Interpolant')
    ax.scatter(nodes, runge_fn(nodes), color='red', zorder=5)
    ax.set_ylim(-0.5, 1.5)
    ax.set_title(title)
    ax.legend()
plt.tight_layout()
plt.show()
```python

### Newton's Divided Differences

The same unique interpolation polynomial can be written in **Newton's form**, which makes adding a new data point cheap (just one extra term):

    pi_n(x) = f[x_1] + f[x_1,x_2](x-x_1) + f[x_1,x_2,x_3](x-x_1)(x-x_2) + ...

The **divided differences** are defined recursively:

    f[x_i, x_{i+1}] = (f(x_{i+1}) - f(x_i)) / (x_{i+1} - x_i)
    f[x_i, ..., x_{i+k}] = (f[x_{i+1},...,x_{i+k}] - f[x_i,...,x_{i+k-1}]) / (x_{i+k} - x_i)

On a uniform grid with step delta, f[x_1, x_2] ≈ f'(x_1), showing that Newton's polynomial is a finite-difference Taylor expansion.

```python
def newton_divided_differences(x, y):
    """Build divided difference table. Returns Newton coefficient array."""
    n = len(x)
    table = np.zeros((n, n))
    table[:, 0] = y
    for j in range(1, n):
        for i in range(n - j):
            table[i, j] = (table[i+1, j-1] - table[i, j-1]) / (x[i+j] - x[i])
    return table[0, :]  # coefficients c_0, c_1, ..., c_{n-1}


def newton_eval(coeffs, x_nodes, x_eval):
    """Evaluate Newton interpolation polynomial using Horner-like scheme."""
    x_eval = np.asarray(x_eval, float)
    result  = np.full_like(x_eval, coeffs[-1])
    for i in range(len(coeffs)-2, -1, -1):
        result = result * (x_eval - x_nodes[i]) + coeffs[i]
    return result


# Bio application: gene expression interpolation of missing time points
# Simulated qRT-PCR data: 6 measured time points, 2 to be predicted
t_measured = np.array([0, 2, 4, 8, 12, 24], dtype=float)   # hours
expr       = np.array([1.0, 2.3, 4.1, 6.8, 5.2, 3.1])      # log2 fold change

coeffs  = newton_divided_differences(t_measured, expr)
t_dense = np.linspace(0, 24, 200)
expr_interp = newton_eval(coeffs, t_measured, t_dense)

missing_t  = np.array([6.0, 18.0])
expr_pred  = newton_eval(coeffs, t_measured, missing_t)

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_dense, expr_interp, label='Newton interpolant')
ax.scatter(t_measured, expr, color='black', zorder=5, label='Measured')
ax.scatter(missing_t, expr_pred, color='red', marker='*', s=150, zorder=6,
           label='Predicted missing points')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('log2 fold change')
ax.set_title('Gene expression: Newton interpolation of missing time points')
ax.legend()
plt.tight_layout()
plt.show()
print(f'Predicted at t=6h: {expr_pred[0]:.2f},  t=18h: {expr_pred[1]:.2f}')
```python

### Cubic Spline Interpolation

High-degree polynomial interpolants are unstable. **Splines** avoid this by using low-degree polynomials on each sub-interval, stitched together with continuity conditions.

A **cubic spline** s(x) in S_{3,1}(X, I) is a piecewise cubic with continuous first and second derivatives everywhere. It solves the variational problem:

    min_f  integral_a^b [f''(x)]^2 dx   subject to f(x_i) = y_i


**Error bound:** For f in W^4[M_4, I]:  ||s - f||_inf <= (5/364) * h^4 * M_4,  where h = max(x_{i+1} - x_i).

```python
# Bio application: reconstruct a full circadian time course
# Measured at irregular intervals cubic spline gives a smooth biological curve
t_sampled = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24], dtype=float)  # hours ZT
per1_expr = np.array([1.0, 3.2, 7.1, 9.4, 6.8, 2.9, 1.1, 0.8, 1.0])  # Per1 (AU)

# Periodic boundary condition: same value and derivative at both endpoints
cs = CubicSpline(t_sampled, per1_expr, bc_type='periodic')
t_fine = np.linspace(0, 24, 500)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].scatter(t_sampled, per1_expr, color='black', zorder=5, label='Measured')
axes[0].plot(t_fine, cs(t_fine), label='Cubic spline')
axes[0].set_xlabel('ZT (hours)')
axes[0].set_ylabel('Per1 expression (AU)')
axes[0].set_title('Circadian Per1 expression -- cubic spline')
axes[0].legend()

lin = interp1d(t_sampled, per1_expr, kind='linear')
axes[1].plot(t_fine, cs(t_fine), label='Cubic spline')
axes[1].plot(t_fine, lin(t_fine), '--', label='Linear interpolation')
axes[1].scatter(t_sampled, per1_expr, color='black', zorder=5)
axes[1].set_xlabel('ZT (hours)')
axes[1].set_title('Spline vs linear -- smoothness matters')
axes[1].legend()

plt.tight_layout()
plt.show()

peak_t = t_fine[np.argmax(cs(t_fine))]
print(f'Estimated peak expression at ZT {peak_t:.1f}h')
```python

## Part 2: Numerical Differentiation and Integration

### Finite Difference Formulas

We approximate derivatives by replacing the limiting process with finite differences. Derived from Taylor series:

| Formula | Expression | Error order |
|---------|-----------|-------------|
| Forward difference  | (f(x+h) - f(x)) / h         | O(h)   |
| Backward difference | (f(x) - f(x-h)) / h         | O(h)   |
| Central difference  | (f(x+h) - f(x-h)) / (2h)   | O(h^2) |

**Instability with noisy data:** The total error for a noisy measurement is ~M_2*h + epsilon/h, minimised at h* ~ sqrt(epsilon/M_2). Taking h too small amplifies measurement noise — the classic numerical instability of differentiation.

```python
def forward_diff(f, x, h=1e-5):
    return (f(x + h) - f(x)) / h

def backward_diff(f, x, h=1e-5):
    return (f(x) - f(x - h)) / h

def central_diff(f, x, h=1e-5):
    return (f(x + h) - f(x - h)) / (2*h)


# Error vs step size: central difference is O(h2), much better
f_test       = np.sin
f_prime_true = np.cos(1.0)
h_vals = np.logspace(-8, 0, 200)

err_fwd  = np.abs([forward_diff(f_test, 1.0, h)  - f_prime_true for h in h_vals])
err_cent = np.abs([central_diff(f_test, 1.0, h)  - f_prime_true for h in h_vals])

fig, ax = plt.subplots(figsize=(7, 4))
ax.loglog(h_vals, err_fwd,  label='Forward  O(h)')
ax.loglog(h_vals, err_cent, label='Central  O(h^2)')
ax.loglog(h_vals, h_vals,   'k--', alpha=0.5, label='slope 1')
ax.loglog(h_vals, h_vals**2,'k:',  alpha=0.5, label='slope 2')
ax.set_xlabel('Step size h')
ax.set_ylabel('Absolute error')
ax.set_title('Differentiation error vs step size (sin at x=1)')
ax.legend()
plt.tight_layout()
plt.show()

# Bio application: rate of change of drug plasma concentration
print('\n--- Drug clearance rate dC/dt ---')
t_pk = np.array([0, 1, 2, 4, 6, 8, 12, 24], dtype=float)  # hours
C_pk = np.array([100, 72, 52, 27, 14, 7.4, 2.0, 0.1])      # ng/mL
for i in range(1, len(t_pk)-1):
    dCdt = (C_pk[i+1] - C_pk[i-1]) / (t_pk[i+1] - t_pk[i-1])
    print(f'  t={t_pk[i]:4.0f}h: dC/dt = {dCdt:6.2f} ng/mL/h')
```python

### Numerical Integration: Trapezoidal Rule and Simpson's Rule

Numerical integration is a **stable** operation (unlike differentiation): measurement errors average out over the interval.

**Trapezoidal rule** (piecewise linear, f in W^2):

    J_N^T = [(f(a)+f(b))/2 + sum_{k=1}^{N-1} f(x_k)] * (b-a)/N
    |J - J_N^T| <= (b-a)^3 / (12*N^2) * M_2

where M_2 = max|f''(x)| on [a,b].

**Simpson's rule** (piecewise quadratic, f in W^4):

    J_N^S = (b-a)/(6N) * [f(a) + f(b) + 4*sum_odd f(x_k) + 2*sum_even f(x_k)]
    |J - J_N^S| <= (b-a)^5 / (180*N^4) * M_4

where M_4 = max|f''''(x)| on [a,b].

Simpson's rule converges as O(h^4) vs the trapezoidal O(h^2) — two orders faster.
The factor of 180 (not 12) comes from the derivation: for a single panel of width 2h,

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
