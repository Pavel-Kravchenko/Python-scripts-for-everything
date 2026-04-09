---
name: advanced-r-statistics
description: Nonparametric tests, parametric inference, regression, correlation, and normality diagnostics in R for biological and medical data
---

# Advanced R Statistical Computing

## When to Use
- Choosing between parametric and nonparametric tests when normality is uncertain
- Comparing paired biological measurements (before/after, treated/control) with small samples
- Building linear regression models and checking whether residual assumptions hold
- Testing association between categorical variables in contingency tables (chi-squared, prop.test)

---

## Quick Reference

### Test Selection: Parametric vs Nonparametric
| Situation | Parametric (assumes normality) | Nonparametric alternative |
|-----------|-------------------------------|--------------------------|
| One sample vs constant | `t.test(x, mu=mu0)` | `binom.test` (sign test) |
| Paired samples | `t.test(paired=TRUE)` | `wilcox.test(paired=TRUE)` |
| Two independent samples | `t.test(var.equal=?)` | `wilcox.test(paired=FALSE)` |
| 3+ groups | `lm` + `anova` | `kruskal.test` |
| Post-hoc after 3+ groups | `TukeyHSD` | `posthoc.kruskal.dunn.test` (PMCMR) |

### Confidence Interval Types
| Parameter | Quantile | Degrees of freedom |
|-----------|----------|--------------------|
| Mean (σ unknown) | `qt(1-α/2, df=n-1)` | n-1 |
| Variance σ² | `qchisq(1-α/2, df=n-1)` | n-1 |
| Proportion (large n) | `qnorm(1-α/2)` | — |

### Correlation Method Selection
| Condition | Method | Function |
|-----------|--------|----------|
| Both variables normal, linear relationship | Pearson | `cor.test(..., method="pearson")` |
| Non-normal or monotone (not linear) | Spearman | `cor.test(..., method="spearman")` |
| Small n, ordinal data | Kendall | `cor.test(..., method="kendall")` |

---

## Key Patterns

### 1. Exact Binomial Test
```r
# H0: p = 0.5, H1: p > 0.5; b successes in n trials
binom.test(b, n, p = 0.5, alternative = "greater")

# Manual p-value (right-sided): P(B >= b)
pbinom(b - 1, n, 0.5, lower.tail = FALSE)

# Asymptotic statistic B*
B_star <- (b - n * 0.5) / sqrt(n * 0.5 * 0.5)
pnorm(B_star, lower.tail = FALSE)
```

### 2. Power Function for Binomial Test
```r
# Find critical region boundary (left-sided, alpha=0.05)
qb <- qbinom(0.05, n, p_0, lower.tail = TRUE)
# qbinom returns one above boundary for left-sided; use qb-1

power <- function(p) { pbinom(qb - 1, n, p, lower.tail = TRUE) }
plot(power, 0, p_0)
1 - power(p_alternative)  # Type II error
```

### 3. Wilcoxon Tests (both variants)
```r
library(exactRankTests)  # exact p-values with ties or small n

# Paired (signed-rank) — before/after comparisons
wilcox.exact(after, before, paired = TRUE,  alternative = "less", conf.int = TRUE)

# Independent (rank-sum / Mann-Whitney U)
wilcox.exact(group1, group2, paired = FALSE, alternative = "two.sided", conf.int = TRUE)
```

### 4. Hodges-Lehmann Robust Estimate
```r
# conf.int = TRUE adds the H-L estimate and CI to the output
result <- wilcox.test(x, y, paired = FALSE, conf.int = TRUE)
result$estimate   # H-L estimate of location shift
result$conf.int   # confidence interval
```

### 5. Kruskal-Wallis + Dunn Post-Hoc
```r
kruskal.test(list(group1, group2, group3))

library(PMCMR)
posthoc.kruskal.dunn.test(data_frame, p.adjust.method = "bonferroni")
```

### 6. t-Test Workflow (Independent Samples)
```r
# Step 1: test variance equality
var.test(x, y, ratio = 1, alternative = "two.sided")

# Step 2: use result to set var.equal
t.test(x, y, alternative = "greater", mu = 0,
       paired = FALSE, var.equal = FALSE)  # Welch if variances unequal
```

### 7. Normality Assessment
```r
shapiro.test(x)          # p < 0.05 → evidence against normality

library(nortest)
pearson.test(x, adjust = FALSE)  # conservative (n-1 df)
pearson.test(x, adjust = TRUE)   # corrected (n-3 df); truth is between two

qqnorm(x); qqline(x, col = "red")
```

### 8. Pearson Correlation + Fisher z CI
```r
cor.test(x, y, method = "pearson", conf.level = 0.95)

# Manual Fisher z CI
r     <- cor(x, y)
z1    <- 0.5 * log((1 + r) / (1 - r))
zq    <- qnorm(0.975)
z_lo  <- z1 - zq / sqrt(n - 3)
z_hi  <- z1 + zq / sqrt(n - 3)
r_lo  <- tanh(z_lo)
r_hi  <- tanh(z_hi)
```

### 9. Linear Regression Workflow
```r
fit <- lm(y ~ x1 + x2, data = df)
summary(fit)          # coefficients, R², F-test
confint(fit)          # CI for each coefficient
shapiro.test(fit$residuals)   # check normality assumption

# Polynomial: degree 2 and degree 7
lm(y ~ x + I(x^2))
lm(y ~ poly(x, 7))   # orthogonal polynomial basis
```

### 10. ANOVA vs Kruskal-Wallis Comparison
```r
# Parametric ANOVA
anova(lm(response ~ group))

# Nonparametric alternative (no normality required)
kruskal.test(response, group)

# If ANOVA residuals are non-normal (shapiro p < 0.05), prefer Kruskal-Wallis
```

---

## Code Templates

### Confidence Interval for Mean (Student-t)
```r
n        <- length(x)
xbar     <- mean(x)
s        <- sd(x)
alpha    <- 0.05
tq       <- qt(1 - alpha/2, df = n - 1)
ci_left  <- xbar - tq * s / sqrt(n)
ci_right <- xbar + tq * s / sqrt(n)
cat(sprintf("95%% CI: [%.3f, %.3f]\n", ci_left, ci_right))
```

### Confidence Interval for Variance (Chi-squared)
```r
# s2 = sample variance; n = sample size
ci_lo <- (n - 1) * s2 / qchisq(1 - alpha/2, df = n - 1)
ci_hi <- (n - 1) * s2 / qchisq(    alpha/2, df = n - 1)
```

### CLT / LLN Simulation
```r
N <- 10000
sizes <- c(8, 16, 30, 50, 100)
for (n in sizes) {
    M   <- matrix(runif(n * N), nrow = N, ncol = n)
    x   <- apply(M, 1, mean)
    cat(n, "-> sd(xbar) =", sd(x), " theory:", 1/(sqrt(12)*sqrt(n)), "\n")
}
```

### Contingency Table from Data
```r
# From two factor columns
tbl <- xtabs(Freq ~ VarA + VarB, data = df)
ftable(tbl)
chisq.test(tbl)

# Manually specified
tbl2 <- matrix(c(60, 440, 42, 558), nrow = 2, byrow = TRUE)
rownames(tbl2) <- c("GroupA", "GroupB")
colnames(tbl2) <- c("Success", "Failure")
prop.test(c(60, 42), c(500, 600))   # equivalent two-proportion test
```

---

## Common Pitfalls

- **Left-sided `qbinom` is off-by-one**: `qbinom(0.05, n, p)` returns the value *above* the critical boundary for a left-sided test. Subtract 1: use `qb - 1` as the boundary.
- **`lower.tail = FALSE` is strict**: `pbinom(b, n, p, lower.tail=FALSE)` computes P(X > b), not P(X >= b). Use `pbinom(b-1, ...)` for ≥.
- **Paired vs independent**: `paired=TRUE` requires vectors of the same length in matching order. Mispairing inflates variance and loses power.
- **`wilcox.test` vs `wilcox.exact`**: `wilcox.test` with `exact=TRUE` can be slow or incorrect with ties. Use `exactRankTests::wilcox.exact` for small n or tied data.
- **`sd()` in R is the corrected (n-1) version**: matches the `iq_s` convention used in seminars. The formula gives $s$, not $\sigma$.
- **ANOVA F-test requires normality of residuals**: always run `shapiro.test(fit$residuals)` after `lm`. If violated, switch to `kruskal.test`.
- **`adjust=FALSE` vs `adjust=TRUE` in `pearson.test`**: the true p-value lies *between* the two results. Report both bounds.
- **Pearson r formula in R**: `cor(x, y)` uses n-1. The manual formula from seminars uses `(n-1)/n` factor to compensate for `sd()`.

---

## Related Skills
- `biostatistics-r` — R basics, distribution function convention, test selection, multiple testing correction
- `probability-statistics-python` — Python equivalents using scipy.stats and statsmodels
