---
name: foundations-r-hypothesis-testing-and-nonparametrics
description: "**Tier 0 -- Computational Foundations | Module 8, Notebook 1**"
tool_type: python
source_notebook: "Tier_0_Computational_Foundations/08_Advanced_R_Statistics/01_r_hypothesis_testing_and_nonparametrics.ipynb"
---

# R Hypothesis Testing and Nonparametric Methods

*Source: Course notebook `Tier_0_Computational_Foundations/08_Advanced_R_Statistics/01_r_hypothesis_testing_and_nonparametrics.ipynb`*

# R Hypothesis Testing and Nonparametric Methods

**Tier 0 -- Computational Foundations | Module 8, Notebook 1**

---

## Overview

This notebook covers the core nonparametric and exact hypothesis tests used in biological and medical research. These methods make minimal distributional assumptions and are therefore essential when sample sizes are small, distributions are skewed, or normality cannot be verified.

### Learning Objectives

By the end of this notebook you will be able to:
- Apply R's `d/p/q/r` distribution function convention to any distribution
- Run exact and asymptotic binomial tests and interpret p-values
- Perform the sign test, Wilcoxon signed-rank test, and Wilcoxon rank-sum / Mann-Whitney U test in R
- Compute Hodges-Lehmann estimates and robust confidence intervals
- Apply the Kruskal-Wallis test for multi-group comparisons and Dunn post-hoc tests
- Construct power functions and determine required sample sizes

## How to use this notebook
1. All cells use the R kernel (IRkernel). Ensure R is installed with base stats package (included by default).
2. Run each section sequentially — later sections reference datasets created earlier.
3. For each test, note: (a) the null hypothesis being tested, (b) the assumptions, and (c) what the output numbers mean.
4. When comparing parametric and non-parametric tests on the same data, think about which is more appropriate for that data shape.

## Common stumbling points

- **`wilcox.test` vs. Mann-Whitney**: In R, `wilcox.test` with two unpaired samples performs the Mann-Whitney U test. With `paired=TRUE`, it performs the Wilcoxon signed-rank test. The function name is confusing — know which one you need.
- **`p.adjust` method names**: Use `"BH"` for Benjamini-Hochberg FDR, `"bonferroni"` for Bonferroni. Passing `"fdr"` will not work — R uses `"BH"` as the method name.
- **`t.test` default is two-sided**: Use `alternative="greater"` or `alternative="less"` for one-sided tests. In bioinformatics, you should almost always use two-sided unless you have a pre-specified directional hypothesis.
- **ANOVA tells you groups differ, not which ones**: A significant ANOVA F-test only tells you that *at least one* group is different. Use post-hoc tests (TukeyHSD, pairwise.t.test with correction) to find which pairs differ.

---

## 1. R's Distribution Functions

R was built for statistical computing. One of its most practical features is a **consistent naming convention** for distribution functions: every distribution is accessible through four functions that share a common name prefix. Once you learn this convention for one distribution, you know it for all of them.

### 1.1 The d/p/q/r Distribution Function Convention

Every probability distribution in R is accessed through four functions sharing a common prefix:

| Prefix | Returns | Example |
|--------|---------|--------|
| `d` | Probability density (or mass) at *x* | `dnorm(x, mean, sd)` |
| `p` | Cumulative probability P(X ≤ x) | `pnorm(q, mean, sd)` |
| `q` | Quantile — inverse CDF | `qnorm(p, mean, sd)` |
| `r` | Random sample | `rnorm(n, mean, sd)` |

The suffix names the distribution: `norm`, `binom`, `t`, `chisq`, `f`, `exp`, `pois`, `unif`, `nbinom`.

The `lower.tail` argument controls which tail is used:
- `lower.tail = TRUE` (default): P(X ≤ x)
- `lower.tail = FALSE`: P(X > x)

```python
# P(X > 54) for Binomial(100, 0.25) — strictly greater than 54
pbinom(54, size = 100, prob = 0.25, lower.tail = FALSE)

# P(X > 69) for Binomial(100, 0.85)
pbinom(69, size = 100, prob = 0.85, lower.tail = FALSE)
```

```python
# Exponential distribution: P(X <= 12) where mean = 80 min
# rate = lambda = 1 / mean
pexp(12, rate = 1/80)

# Normal distribution: P(X > 252) where mean = 262.5, sd = 12
sigma <- 12
mu    <- 262.5
pnorm(252, mu, sigma, lower.tail = FALSE)
```

```python
# Sample size estimation with normal quantile
# How many samples needed so that the 99.5% CI for the difference of two means
# has half-width <= 0.5, given sigma = 12?
sigma <- 12
n_2 <- (sigma / 0.5 * qnorm(0.995, mean = 0, sd = 1))^2
ceiling(n_2)   # ceiling() rounds up to the nearest integer
```

---

## 2. Exact Binomial Test

The exact binomial test evaluates whether the probability of "success" in a sequence of Bernoulli trials equals a specified value $p_0$.

**Hypotheses:**
- $H_0$: $p = p_0$
- $H_1$: $p > p_0$ (right-sided), $p < p_0$ (left-sided), or $p \neq p_0$ (two-sided)

**Test statistic:** $B = \sum_{i=1}^n \mathbb{1}[X_i = \text{success}]$ — total number of successes.

Under $H_0$, $B \sim \text{Binomial}(n, p_0)$.  
The p-value is computed exactly from the binomial CDF.

```python
# Example: In a survey of 1002 patients, 701 preferred treatment A over B.
# Test H0: p = 0.5 against H1: p > 0.5 (more patients prefer A).
n   <- 1002
b   <- 701
p_0 <- 0.5

my_test <- binom.test(b, n, p_0, alternative = "greater")
my_test

# Access just the p-value:
my_test$p.value
```

```python
# Manual p-value using pbinom (equivalent to binom.test for right-sided alternative)
# P(B >= b) = P(B > b-1) with lower.tail = FALSE
pbinom(b - 1, n, p_0, lower.tail = FALSE)
```

```python
# Asymptotic (normal approximation) binomial test
# Standardized statistic B* ~ N(0,1) under H0
B_star <- (b - n * p_0) / sqrt(n * p_0 * (1 - p_0))

# Critical value for right-sided test at alpha = 0.05
qnorm(0.95, mean = 0, sd = 1)

# p-value from the asymptotic test
pnorm(B_star, mean = 0, sd = 1, lower.tail = FALSE)
```

### 2.1 Finding the Critical Region

For the exact test, the critical region is the set of values of $B$ that lead to rejection. With a **left-sided** alternative at $\alpha = 0.05$:

```python
# Left-sided exact binomial critical region
# H0: p = 0.3, H1: p < 0.3, n = 50
n   <- 50
p_0 <- 0.3

# qbinom with left-sided alternative returns the value ONE ABOVE the critical boundary
qb <- qbinom(0.05, n, p_0, lower.tail = TRUE)
qb

# Verify: P(B <= qb-1) < 0.05 and P(B <= qb) > 0.05
pbinom(qb - 1, n, p_0, lower.tail = TRUE)
pbinom(qb,     n, p_0, lower.tail = TRUE)
```

```python
# Asymptotic critical region using normal quantile
C_krit <- qnorm(0.05) * sqrt(n * p_0 * (1 - p_0)) + n * p_0
floor(C_krit)   # floor() rounds down
```

---

## 3. Power Analysis for the Binomial Test

The **power function** $\beta(p)$ gives the probability of rejecting $H_0$ when the true probability of success is $p$. It is defined using the critical region found above.

```python
# Power function for left-sided binomial test (H0: p=0.3, critical region: B <= 9)
# n = 50, critical boundary at qb - 1 = 9
n   <- 50
p_0 <- 0.3
qb  <- qbinom(0.05, n, p_0, lower.tail = TRUE)

power <- function(p) { pbinom(qb - 1, n, p, lower.tail = TRUE) }

# Plot power as a function of p
plot(power, 0, 0.5,
     xlab = "True probability p",
     ylab = "Power  beta(p)",
     main = "Power Function — Exact Binomial Test")

# Power at a specific alternative (H1: p = 0.25)
power(0.25)

# Type II error probability (beta): probability of failing to reject H0 when H1 is true
cat("Type II error at p=0.25:", 1 - power(0.25), "\n")
```

### 3.1 Sample Size Determination

How many trials are needed so that the test has:
- Type I error probability $\leq 5\%$
- Type II error probability $\leq 2.5\%$ when $H_1: p = p_1$

```python
# Required sample size as a function of the alternative p1
# H0: p = 0.3, left-sided alternative, alpha = 0.05, beta = 0.025
p_0 <- 0.3

m <- function(p) {
  ceiling(
    ((qnorm(0.975) * sqrt(p * (1 - p)) - sqrt(p_0 * (1 - p_0)) * qnorm(0.05)) /
     (p_0 - p))^2
  )
}

# Required n when the true p is 0.25
cat("Required n for p1=0.25:", m(0.25), "\n")

# Required n when the true p is 0.20
cat("Required n for p1=0.20:", m(0.20), "\n")
```

---

## 4. Sign Test

The **sign test** is the simplest nonparametric test for paired data. It tests whether the median difference between paired observations is zero.

**Idea:** Let $D_i = Y_i - X_i$ (after minus before). Under $H_0$, each $D_i$ is equally likely to be positive or negative, so $B = \sum \mathbb{1}[D_i > 0] \sim \text{Binomial}(n, 0.5)$.

The sign test is **identical** to a binomial test with $p_0 = 0.5$.

```python
# Load weight data: patients measured before and after a diet intervention
# Columns: weight_before, weight_after
data <- read.table("../Assets/data/r_statistics/Weight.txt",
                   header = TRUE, sep = " ", dec = ".")

weight_before <- data[[1]]
weight_after  <- data[[2]]

head(data)
```

```python
# H0: median difference = 0, H1: median difference < 0 (weight decreases after diet)
b <- sum(weight_after > weight_before)   # count of positive differences
n <- length(weight_before)

cat("Positive differences (weight increased):", b, "of", n, "\n")

# Exact sign test = exact binomial test with p0 = 0.5
binom.test(b, n, p = 0.5, alternative = "less")
```

```python
# Equivalent manual calculation using pbinom
pbinom(b, n, 0.5, lower.tail = TRUE)

# Asymptotic sign test: standardized B*
B_star <- (b - n * 0.5) / sqrt(n * 0.5 * 0.5)
cat("Asymptotic statistic B*:", B_star, "\n")

# Normal critical value for left-sided test at alpha = 0.05
qnorm(0.05, mean = 0, sd = 1, lower.tail = TRUE)

# Asymptotic p-value
pnorm(B_star, mean = 0, sd = 1, lower.tail = TRUE)
```

---

## 5. Wilcoxon Signed-Rank Test

The **Wilcoxon signed-rank test** is more powerful than the sign test because it uses both the sign *and* the magnitude of each difference $D_i = Y_i - X_i$.

**Procedure:**
1. Compute $D_i = Y_i - X_i$ for each pair; discard ties ($D_i = 0$).
2. Rank the absolute differences $|D_i|$ from 1 to *n*.
3. Compute $W^+ = $ sum of ranks where $D_i > 0$, $W^- = $ sum of ranks where $D_i < 0$.
4. The test statistic is $V = \min(W^+, W^-)$.

**Assumption:** The differences $D_i$ are symmetrically distributed around a common median $\theta$.

**Null hypothesis:** $\theta = 0$.

```python
# Install exactRankTests for exact p-values with ties or small n
# install.packages("exactRankTests")
library(exactRankTests)

# Wilcoxon signed-rank test on weight data (paired = TRUE)
wilcox.test(weight_after, weight_before,
            paired = TRUE, alternative = "less", exact = TRUE)
```

```python
# wilcox.exact from exactRankTests — preferred for exact p-values with potential ties
wilcox.exact(weight_after, weight_before,
             paired = TRUE, alternative = "less", exact = TRUE)
```

### 5.1 Manual Step-by-Step Wilcoxon Sign-Rank Calculation

The following recreates the Wilcoxon signed-rank test by hand using the journey time example from Practicum 2: minutes of commute before and after a new road was opened.

```python
# Journey times (minutes) before and after road opening
x <- c(125.3, 101.0, 117.2, 133.7,  96.4, 124.5, 118.7, 106.2, 116.3, 120.2, 125.0, 128.8)
y <- c(127.3, 120.2, 126.2, 125.4, 115.1, 118.5, 135.5, 118.2, 122.9, 120.1, 120.8, 130.7)

# Built-in test for comparison
wilcox.test(x, y, paired = TRUE)
```

```python
# Manual calculation
df <- data.frame(x = x, y = y)
df$diff <- with(df, y - x)

# Assign ranks to |diff| (sorted by absolute value)
df[order(abs(df$diff)), 'rk'] <- 1:nrow(df)

# Signed ranks
df$signrk <- df$rk * sign(df$diff)
df

# W+ = sum of positive signed ranks, W- = |sum of negative signed ranks|
r <- df$signrk
W_plus  <- sum(r[r > 0])
W_minus <- sum(-r[r < 0])
cat("W+ =", W_plus, " W- =", W_minus, "\n")

# Asymptotic approximation for V = min(W+, W-)
n_eff <- nrow(df)
V <- min(W_plus, W_minus)
z <- (V - 0) / sqrt(n_eff * (n_eff + 1) * (2 * n_eff + 1) / 6)
cat("Two-sided asymptotic p-value:", 2 * pnorm(-abs(z)), "\n")
```

---

## 6. Wilcoxon Rank-Sum Test (Mann-Whitney U)

The **Wilcoxon rank-sum test** (equivalent to Mann-Whitney U) compares the location of two **independent** samples without assuming normality.

**Null hypothesis:** The two populations have the same distribution (or equivalently, $P(X > Y) = 0.5$).

Use `paired = FALSE` for independent samples.

```python
# Load blood pressure data for female and male patients
female_data <- read.table("../Assets/data/r_statistics/female.csv",
                          header = TRUE, sep = ";", dec = ",")
male_data   <- read.table("../Assets/data/r_statistics/male.csv",
                          header = TRUE, sep = ";", dec = ",")

female <- female_data[[1]]
male   <- male_data[[1]]

cat("Female n =", length(female), "  Male n =", length(male), "\n")
cat("Female median:", median(female), "  Male median:", median(male), "\n")
```
