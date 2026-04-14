---
name: foundations-r-hypothesis-testing-and-nonparametrics
description: "R hypothesis testing and nonparametric methods: binomial test, sign test, Wilcoxon, Kruskal-Wallis, and power analysis. Use when choosing and applying statistical tests in R."
tool_type: python
primary_tool: Python
---

# R Hypothesis Testing and Nonparametric Methods

## Pitfalls

- **`wilcox.test` naming:** With two unpaired samples it performs Mann-Whitney U. With `paired=TRUE` it performs Wilcoxon signed-rank. The name is ambiguous — always specify `paired=`.
- **`p.adjust` method names:** Use `"BH"` for Benjamini-Hochberg FDR. Passing `"fdr"` fails silently (not a valid method name).
- **`t.test` default is two-sided:** Use `alternative="greater"` or `alternative="less"` only when you have a pre-specified directional hypothesis.
- **ANOVA ≠ pairwise differences:** Significant F-test only means at least one group differs. Use `TukeyHSD` or `pairwise.t.test` with correction to find which pairs.
- **Assumption violations:** Check normality and homoscedasticity before parametric tests; use nonparametric alternatives when assumptions fail.
- **Bootstrap ≠ probability:** A bootstrap value of 95 means 95% of replicates recovered that result — not a 95% probability of being correct.

## R Distribution Function Convention

| Prefix | Returns | Example |
|---|---|---|
| `d` | Density/mass at x | `dnorm(x, mean, sd)` |
| `p` | CDF: P(X ≤ x) | `pnorm(q, mean, sd)` |
| `q` | Quantile (inverse CDF) | `qnorm(p, mean, sd)` |
| `r` | Random sample | `rnorm(n, mean, sd)` |

Suffixes: `norm`, `binom`, `t`, `chisq`, `f`, `pois`, `exp`, `unif`, `nbinom`.
`lower.tail=FALSE` gives P(X > x).

```r
# P(X > 54) for Binomial(100, 0.25)
pbinom(54, size=100, prob=0.25, lower.tail=FALSE)

# Normal: P(X > 252) where mean=262.5, sd=12
pnorm(252, mean=262.5, sd=12, lower.tail=FALSE)

# Sample size for 99.5% CI half-width = 0.5, sigma = 12
ceiling((12 / 0.5 * qnorm(0.995))^2)
```

## Exact Binomial Test

H0: p = p0. Test statistic B ~ Binomial(n, p0). Use for: proportion of successes.

```r
# 701 of 1002 patients prefer treatment A — test H0: p=0.5 (right-sided)
my_test <- binom.test(701, 1002, 0.5, alternative="greater")
my_test$p.value

# Equivalent manual calculation
pbinom(700, 1002, 0.5, lower.tail=FALSE)

# Asymptotic (normal approximation)
B_star <- (701 - 1002 * 0.5) / sqrt(1002 * 0.5 * 0.5)
pnorm(B_star, lower.tail=FALSE)
```

### Power analysis
```r
n   <- 50; p_0 <- 0.3
qb  <- qbinom(0.05, n, p_0, lower.tail=TRUE)
power <- function(p) pbinom(qb - 1, n, p, lower.tail=TRUE)

# Required sample size for alpha=0.05, beta=0.025, p1=0.25
m <- function(p) ceiling(
  ((qnorm(0.975)*sqrt(p*(1-p)) - sqrt(p_0*(1-p_0))*qnorm(0.05)) / (p_0 - p))^2
)
m(0.25)
```

## Sign Test

Paired data, tests H0: median difference = 0. Identical to binomial test with p0 = 0.5.

```r
b <- sum(weight_after > weight_before)
n <- length(weight_before)
binom.test(b, n, p=0.5, alternative="less")
```

## Wilcoxon Signed-Rank Test (paired)

More powerful than sign test — uses sign and magnitude of differences. Assumption: differences are symmetrically distributed.

```r
library(exactRankTests)
wilcox.exact(weight_after, weight_before, paired=TRUE, alternative="less")

# Manual: compute differences, rank by |diff|, sum positive ranks
df$diff   <- df$y - df$x
df$rk     <- rank(abs(df$diff))
df$signrk <- df$rk * sign(df$diff)
W_plus  <- sum(df$signrk[df$signrk > 0])
W_minus <- sum(-df$signrk[df$signrk < 0])
```

## Wilcoxon Rank-Sum / Mann-Whitney U (independent)

H0: two populations have the same distribution.

```r
wilcox.test(female, male, paired=FALSE, alternative="two.sided")
```

## Test Selection Guide

| Data | Groups | Parametric | Nonparametric |
|---|---|---|---|
| Proportions | 1 | `binom.test` | — |
| Paired continuous | 2 | `t.test(paired=T)` | `wilcox.test(paired=T)` |
| Independent continuous | 2 | `t.test` | `wilcox.test` |
| Independent continuous | ≥3 | `aov` + `TukeyHSD` | `kruskal.test` |
