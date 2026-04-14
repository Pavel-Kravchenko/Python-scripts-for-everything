---
name: foundations-r-regression-correlation-and-diagnostics
description: "R regression, correlation, and diagnostics: confidence intervals, t-tests, normality testing, Pearson/Spearman/Kendall correlation, Fisher z-transform CI."
tool_type: python
primary_tool: Python
---

# R Regression, Correlation, and Diagnostics

## Pitfalls

- **R-squared alone is misleading**: High R² does not mean the model is correct or causal. Always check residual diagnostics.
- **Residuals vs. Fitted**: The single most important diagnostic. A pattern (curve, funnel) = model misspecification or violated assumptions. Random scatter around zero = good.
- **Pearson vs. Spearman**: Pearson measures linear association (requires approximately normal data). Spearman measures monotonic association via ranks — use for skewed data (counts, survival times) or when linearity is uncertain.
- **Confounders in multiple regression**: Adding a confounder changes all other coefficients. Each coefficient means "effect of this variable holding all others constant."
- **var.test before t.test**: Always test variance equality before choosing `var.equal=TRUE/FALSE` in a two-sample t-test.
- **Bonferroni is conservative**: For many correlated tests (genomic), prefer Benjamini-Hochberg FDR.

## CI Selection Decision Table

| Situation | Function | Distribution |
|---|---|---|
| Mean, σ unknown (usual case) | `qt(1-α/2, df=n-1)` | Student-t |
| Variance σ² | `qchisq(...)` | Chi-squared |
| Count / Poisson mean, large n | `qnorm(1-α/2)` | Normal asymptotic |
| Proportion, large n | `qnorm(1-α/2)` | Normal asymptotic |

## Key Patterns

### Student-t CI for the mean
```r
n <- length(x); xbar <- mean(x); s <- sd(x); alpha <- 0.05
q <- qt(1 - alpha/2, df = n - 1)
ci <- c(xbar - q * s / sqrt(n), xbar + q * s / sqrt(n))
t.test(x, conf.level = 1 - alpha)   # also gives CI
```

### Chi-squared CI for variance
```r
# H0: sigma^2 = sigma0^2
T_stat <- (n - 1) * var(x) / sigma0^2
crit   <- qchisq(1 - alpha, df = n - 1)
ci_var <- c((n-1)*var(x)/qchisq(1-alpha/2, df=n-1),
            (n-1)*var(x)/qchisq(alpha/2,   df=n-1))
```

### Two-sample t-test workflow
```r
var.test(x, y, ratio = 1, alternative = "two.sided")  # F-test first
t.test(x, y, alternative = "greater", paired = FALSE, var.equal = FALSE)  # Welch

# Nonparametric alternative
library(exactRankTests)
wilcox.exact(x, y, paired = FALSE, alternative = "greater", exact = TRUE)
```

### Normality testing
```r
shapiro.test(x)                           # most powerful for n <= 5000
library(nortest)
pearson.test(x, adjust = FALSE)           # chi-sq GoF; true p-value between adjust=F and T
qqnorm(x); qqline(x, col = "red")        # visual check
```

### Correlation
```r
# Choose method based on normality
shapiro.test(x); shapiro.test(y)
cor.test(x, y, method = "pearson",  conf.level = 0.95)   # Fisher z CI included
cor.test(x, y, method = "spearman")
cor.test(x, y, method = "kendall")
```

### Fisher z-transform CI for Pearson r (manual)
```r
r      <- cor(x, y)
z1     <- 0.5 * log((1 + r) / (1 - r))
norm_q <- qnorm(1 - alpha/2)
z_lo   <- z1 - norm_q / sqrt(n - 3)
z_hi   <- z1 + norm_q / sqrt(n - 3)
r_ci   <- (exp(2 * c(z_lo, z_hi)) - 1) / (exp(2 * c(z_lo, z_hi)) + 1)
```

### When Pearson fails — scenario guide
```r
# Quadratic relationship: Pearson near 0, Spearman detects monotone dependence
y3 <- x^2 + rnorm(n, 0, sd)
cor.test(x, y3, method = "pearson")$estimate   # low
cor.test(x, y3, method = "spearman")$estimate  # higher
```
