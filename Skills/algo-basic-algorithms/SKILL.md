---
name: algo-basic-algorithms
description: Foundational algorithms — Euclidean GCD, Newton's method for cube root — with complexity analysis and clean implementations.
tool_type: python
primary_tool: Python
---

# Basic Algorithms

## Complexity Table

| Algorithm | Best | Average | Worst | Space |
|-----------|------|---------|-------|-------|
| Euclidean GCD | O(1) | O(log min(a,b)) | O(log min(a,b)) | O(1) |
| Newton's cube root | O(1) | O(log 1/ε) | O(max_iter) | O(1) |

## Euclidean GCD

Key identity: `gcd(a, b) = gcd(b, a % b)`. Base case: `gcd(a, 0) = a`.

```python
def gcd(a: int, b: int) -> int:
    a, b = abs(a), abs(b)
    while b != 0:
        a, b = b, a % b
    return a

def lcm(a: int, b: int) -> int:
    return abs(a * b) // gcd(a, b)
```

Worst case is consecutive Fibonacci numbers (e.g., `gcd(89, 55)`). Each step reduces the problem size by at least half → O(log min(a, b)) steps.

Applications: simplify fractions, compute LCM, RSA key generation, Diophantine equations.

## Newton's Method — Cube Root

Solves `f(x) = x³ - A = 0` via iteration: `x_{n+1} = (2x + A/x²) / 3`.

Convergence is **quadratic** near the root (correct digits roughly double each iteration). Typically 5–15 iterations for ε = 10⁻¹⁰.

```python
def cube_root(a: float, epsilon: float = 1e-10, max_iterations: int = 1000) -> float:
    if a == 0:
        return 0.0
    sign = 1 if a > 0 else -1
    a_abs = abs(a)
    x = a_abs / 2 if a_abs > 1 else a_abs  # initial guess
    for _ in range(max_iterations):
        x_new = (2 * x + a_abs / (x * x)) / 3
        if abs(x_new - x) < epsilon:
            return sign * x_new
        x = x_new
    return sign * x
```

## Pitfalls

- **Euclidean GCD with negatives**: always take `abs()` first; `gcd(-48, 18)` should return 6, not -6.
- **LCM overflow**: in Python integers are arbitrary precision, but in other languages `a * b` overflows before dividing. Compute `(a // gcd(a, b)) * b` instead.
- **Newton's method diverges from bad initial guess**: for cube root, `a/2` works well for `a > 1`; use `a` itself for `0 < a <= 1`.
- **Newton's method doesn't handle `x = 0`**: the update formula divides by `x²`; check `a == 0` as a special case.
- **Using float epsilon for integers**: GCD is exact integer arithmetic — no epsilon needed. Epsilon-based convergence only applies to numerical methods.
