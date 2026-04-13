---
name: algo-basic-algorithms
description: "This notebook covers four fundamental algorithms that form the building blocks of computer science:"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/01_Complexity_Analysis/02_basic_algorithms.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Basic Algorithms: Foundations of Computer Science

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/01_Complexity_Analysis/02_basic_algorithms.ipynb`*

# Basic Algorithms: Foundations of Computer Science

This notebook covers four fundamental algorithms that form the building blocks of computer science:

1. **Euclidean Algorithm** - Finding the Greatest Common Divisor (GCD)
2. **Newton's Method** - Numerical root finding (cube root)
3. **Bubble Sort** - Elementary comparison-based sorting
4. **Caesar Cipher** - Classical substitution cipher

Each algorithm is presented with:
- Theoretical background and historical context
- Step-by-step explanation with visual diagrams
- Complexity analysis
- Clean Python implementation
- Comprehensive examples and edge cases

---

# 1. Euclidean Algorithm (GCD)

## 1.1 Theory

### What Problem Does This Solve?

The **Greatest Common Divisor (GCD)** of two integers is the largest positive integer that divides both numbers without leaving a remainder. The Euclidean algorithm efficiently computes this value.

**Applications:**
- Simplifying fractions (e.g., 8/12 → 2/3)
- Cryptography (RSA algorithm)
- Computing Least Common Multiple: `LCM(a, b) = (a × b) / GCD(a, b)`
- Solving linear Diophantine equations

### Mathematical Background

The algorithm is based on a key mathematical property:

$$\gcd(a, b) = \gcd(b, a \mod b)$$

**Proof:** If $d$ divides both $a$ and $b$, then $d$ also divides $a - kb$ for any integer $k$. 
Since $a \mod b = a - \lfloor a/b \rfloor \cdot b$, any common divisor of $a$ and $b$ is also a common divisor of $b$ and $a \mod b$.

### Historical Context

The Euclidean algorithm appears in **Euclid's Elements** (circa 300 BCE), making it one of the oldest algorithms still in common use today — over 2,300 years old! Euclid described it geometrically: finding the largest square that can tile two rectangles.

## 1.2 How It Works

### Step-by-Step Explanation

1. **Start** with two positive integers `a` and `b`
2. **Divide** `a` by `b` and get the remainder `r`
3. **Replace** `a` with `b`, and `b` with `r`
4. **Repeat** until `b` becomes 0
5. **Return** `a` (which now holds the GCD)

### ASCII Art Diagram

```
gcd(48, 18):
┌─────────────────────────────────────────────────────┐
│  Step 1:  48 = 18 × 2 + 12    →  gcd(18, 12)       │
│           ▲    ▲       ▲                            │
│           a    b    remainder                       │
├─────────────────────────────────────────────────────┤
│  Step 2:  18 = 12 × 1 + 6     →  gcd(12, 6)        │
├─────────────────────────────────────────────────────┤
│  Step 3:  12 = 6  × 2 + 0     →  gcd(6, 0) = 6     │
│                         ▲                           │
│                    remainder = 0, STOP!             │
└─────────────────────────────────────────────────────┘
                         ║
                         ▼
                    Result: 6
```

### Pseudocode

```
FUNCTION gcd(a, b):
    WHILE b ≠ 0:
        temp ← b
        b ← a MOD b
        a ← temp
    RETURN a
```

## 1.3 Complexity Analysis

| Case | Time Complexity | Space Complexity |
|------|-----------------|------------------|
| Best | O(1) | O(1) |
| Average | O(log(min(a,b))) | O(1) |
| Worst | O(log(min(a,b))) | O(1) |

**Notes:**
- The worst case occurs with consecutive Fibonacci numbers (e.g., gcd(F_n, F_{n-1}))
- Each step reduces the larger number by at least half
- The number of steps is at most 5 times the number of digits in the smaller number (Lamé's theorem)

## 1.4 Implementation

```python
def gcd(a: int, b: int) -> int:
    """
    Compute the Greatest Common Divisor using the Euclidean algorithm.
    
    The algorithm repeatedly replaces the larger number with the remainder
    of dividing the larger by the smaller until one number becomes zero.
    
    Args:
        a: First positive integer
        b: Second positive integer
    
    Returns:
        The greatest common divisor of a and b
    
    Complexity:
        Time: O(log(min(a, b)))
        Space: O(1)
    
    Examples:
        >>> gcd(48, 18)
        6
        >>> gcd(17, 13)
        1
    """
    # Ensure we work with absolute values
    a, b = abs(a), abs(b)
    
    # The core Euclidean algorithm
    while b != 0:
        a, b = b, a % b  # Simultaneous assignment: a becomes b, b becomes remainder
    
    return a
```

```python
def gcd_verbose(a: int, b: int) -> int:
    """
    Verbose version of GCD that shows each step of the computation.
    
    Args:
        a: First positive integer
        b: Second positive integer
    
    Returns:
        The greatest common divisor of a and b
    """
    a, b = abs(a), abs(b)
    step = 1
    
    print(f"Computing gcd({a}, {b})")
    print("=" * 50)
    
    while b != 0:
        quotient = a // b
        remainder = a % b
        print(f"Step {step}: {a} = {b} × {quotient} + {remainder}")
        print(f"         → gcd({b}, {remainder})")
        a, b = b, remainder
        step += 1
    
    print("=" * 50)
    print(f"Result: GCD = {a}")
    return a
```

```python
# Basic example with step-by-step walkthrough
print("Example 1: gcd(48, 18)")
print()
result = gcd_verbose(48, 18)
```

```python
# Example with coprime numbers (GCD = 1)
print("Example 2: Coprime numbers gcd(17, 13)")
print()
result = gcd_verbose(17, 13)
```

```python
# Worst case: Fibonacci numbers
print("Example 3: Worst case - Fibonacci numbers gcd(89, 55)")
print()
result = gcd_verbose(89, 55)
```

```python
# Edge cases
print("Edge Cases:")
print("=" * 50)

# One number is zero
print(f"gcd(42, 0) = {gcd(42, 0)}  # GCD with 0 is the other number")

# Same numbers
print(f"gcd(15, 15) = {gcd(15, 15)}  # GCD of same number is itself")

# One divides the other
print(f"gcd(100, 25) = {gcd(100, 25)}  # 25 divides 100 evenly")

# Negative numbers
print(f"gcd(-48, 18) = {gcd(-48, 18)}  # Works with negatives")

# Large numbers
print(f"gcd(123456789, 987654321) = {gcd(123456789, 987654321)}")
```

---

# 2. Newton's Method (Cube Root)

## 2.1 Theory

### What Problem Does This Solve?

**Newton's Method** (also called Newton-Raphson method) is an iterative numerical technique for finding successively better approximations to the roots of a function. Here we apply it to find cube roots.

Given a number $A$, we want to find $x$ such that:
$$x^3 = A \quad \Rightarrow \quad x = \sqrt[3]{A}$$

**Applications:**
- Scientific computing (roots, optimization)
- Computer graphics (fast inverse square root)
- Machine learning (optimization algorithms)
- Financial modeling

### Mathematical Background

To find $\sqrt[3]{A}$, we solve $f(x) = x^3 - A = 0$.

Newton's iteration formula:
$$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$$

For $f(x) = x^3 - A$, we have $f'(x) = 3x^2$, giving us:

$$x_{n+1} = x_n - \frac{x_n^3 - A}{3x_n^2} = \frac{2x_n + A/x_n^2}{3}$$

### Historical Context

The method is named after **Isaac Newton** (1643-1727), though similar techniques were known earlier. Newton described it in *De analysi* (1669). Joseph Raphson published a simplified version in 1690. The method showcases the power of calculus in numerical computation.

## 2.2 How It Works

### Step-by-Step Explanation

1. **Choose** an initial guess $x_0$ (often $A/2$ or $A$ itself)
2. **Iterate** using the formula: $x_{n+1} = \frac{2x_n + A/x_n^2}{3}$
3. **Check** if $|x_{n+1} - x_n| < \epsilon$ (tolerance)
4. **Return** $x_{n+1}$ when converged

### ASCII Art Diagram

```
Finding ∛27 (should converge to 3.0):
┌──────────────────────────────────────────────────────────────┐
│  Initial guess: x₀ = 27/2 = 13.5                            │
├──────────────────────────────────────────────────────────────┤
│  Iteration 1:                                                │
│    x₁ = (2×13.5 + 27/13.5²) / 3                             │
│       = (27 + 0.148) / 3 = 9.049...                         │
├──────────────────────────────────────────────────────────────┤
│  Iteration 2:                                                │
│    x₂ = (2×9.049 + 27/9.049²) / 3 = 6.361...                │
├──────────────────────────────────────────────────────────────┤
│  Iteration 3: x₃ = 4.573...                                  │
│  Iteration 4: x₄ = 3.480...                                  │
│  Iteration 5: x₅ = 3.061...                                  │
│  Iteration 6: x₆ = 3.001...                                  │
│  Iteration 7: x₇ = 3.000...                                  │
├──────────────────────────────────────────────────────────────┤
│  Converged! |x₇ - x₆| < ε                                    │
└──────────────────────────────────────────────────────────────┘
                              ║
                              ▼
                         Result: 3.0

Geometric interpretation:
                     y
                     │     ╭─── f(x) = x³ - 27
                     │    ╱
                     │   ╱
                     │  ╱
            ─────────┼─●────────── x
                     │  3 (root)
                     │
```

### Pseudocode

```
FUNCTION cube_root(A, epsilon, max_iterations):
    x ← A / 2                          # Initial guess
    FOR i FROM 1 TO max_iterations:
        x_new ← (2 * x + A / x²) / 3   # Newton's formula
        IF |x_new - x| < epsilon:      # Convergence check
            RETURN x_new
        x ← x_new
    RETURN x                           # Best approximation
```

## 2.3 Complexity Analysis

| Case | Time Complexity | Space Complexity |
|------|-----------------|------------------|
| Best | O(1) | O(1) |
| Average | O(log(1/ε)) | O(1) |
| Worst | O(max_iterations) | O(1) |

**Notes:**
- Newton's method has **quadratic convergence** near the root: the number of correct digits roughly doubles each iteration
- For cube root, typically converges in 5-15 iterations for ε = 10⁻¹⁰
- Convergence depends on the quality of the initial guess

## 2.4 Implementation

```python
def cube_root(a: float, epsilon: float = 1e-10, max_iterations: int = 1000) -> float:
    """
    Compute the cube root of a number using Newton's method.
    
    The method iteratively improves an estimate using the formula:
    x_new = (2 * x + a / x^2) / 3
    
    This is derived from Newton's iteration for f(x) = x³ - a.
    
    Args:
        a: The number to find the cube root of (can be negative)
        epsilon: Convergence tolerance (default: 1e-10)
        max_iterations: Maximum number of iterations (default: 1000)
    
    Returns:
        The cube root of a
    
    Complexity:
        Time: O(log(1/epsilon)) due to quadratic convergence
        Space: O(1)
    
    Examples:
        >>> cube_root(27)
        3.0
        >>> cube_root(-8)
        -2.0
    """
    # Handle special cases
    if a == 0:
        return 0.0
    
    # Handle negative numbers: ∛(-a) = -∛a
    sign = 1 if a > 0 else -1
    a = abs(a)
    
    # Initial guess: a/2 works well for most cases
    x = a / 2 if a > 1 else a
    
    for _ in range(max_iterations):
        # Newton's iteration formula for cube root
        x_new = (2 * x + a / (x * x)) / 3
        
        # Check for convergence
        if abs(x_new - x) < epsilon:
            return sign * x_new
        
        x = x_new
    
    return sign * x
```

```python
def cube_root_verbose(a: float, epsilon: float = 1e-10, max_iterations: int = 20) -> float:
    """
    Verbose version of cube_root that shows convergence progress.
    
    Args:
        a: The number to find the cube root of
        epsilon: Convergence tolerance
        max_iterations: Maximum number of iterations to display
    
    Returns:
        The cube root of a
    """
    if a == 0:
        print("∛0 = 0")
        return 0.0
    
    sign = 1 if a > 0 else -1
    a_abs = abs(a)
    
    x = a_abs / 2 if a_abs > 1 else a_abs
    
    print(f"Computing ∛{a}")
    print("=" * 60)
    print(f"Initial guess: x₀ = {x:.10f}")
    print()
    
    for i in range(max_iterations):
        x_new = (2 * x + a_abs / (x * x)) / 3
        error = abs(x_new - x)
        
        print(f"Iteration {i+1:2d}: x = {sign * x_new:15.10f}  |  error = {error:.2e}")
        
        if error < epsilon:
            print()
            print("=" * 60)
            print(f"Converged after {i+1} iterations!")
            print(f"Result: ∛{a} ≈ {sign * x_new:.10f}")
            print(f"Verification: ({sign * x_new:.10f})³ = {(sign * x_new)**3:.10f}")
            return sign * x_new
        
        x = x_new
    
    print(f"\nDid not converge within {max_iterations} iterations")
    return sign * x
```
