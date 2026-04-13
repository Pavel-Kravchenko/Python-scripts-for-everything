---
name: foundations-r-regression-correlation-and-diagnostics
description: "R regression, correlation, and diagnostics: Pearson/Spearman correlation, linear regression, residual analysis, and model diagnostics. Use when modeling relationships in biological data."
tool_type: python
source_notebook: "Tier_0_Computational_Foundations/08_Advanced_R_Statistics/02_r_regression_correlation_and_diagnostics.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# R Regression, Correlation, and Diagnostics

*Source: Course notebook `Tier_0_Computational_Foundations/08_Advanced_R_Statistics/02_r_regression_correlation_and_diagnostics.ipynb`*


**Tier 0 -- Computational Foundations | Module 8, Notebook 2**

---

## Overview

This notebook covers parametric inference, normality testing, correlation analysis, contingency table methods, linear regression, ANOVA, and simulation-based demonstrations of the Law of Large Numbers and Central Limit Theorem — all implemented in R.

### Learning Objectives

By the end of this notebook you will be able to:
- Construct Student-t, chi-squared, and normal asymptotic confidence intervals
- Run paired and unpaired t-tests and test variance equality with `var.test`
- Assess normality using Shapiro-Wilk, Pearson chi-squared, and Q-Q plots
- Compute Pearson, Spearman, and Kendall correlation with Fisher z-transform CIs
- Build contingency tables and apply chi-squared independence and homogeneity tests
- Fit simple and multiple linear regression models, evaluate residual diagnostics, and select polynomial degree
- Compare ANOVA and Kruskal-Wallis on the same dataset
- Demonstrate the Law of Large Numbers and Central Limit Theorem by simulation

## How to use this notebook
1. Run cells top-to-bottom — datasets are simulated at the start of each section and reused in subsequent cells.
2. When you run a `lm()` call, always run the diagnostic plots (`plot(model)`) immediately after. A model without diagnostics is incomplete.
3. For each regression output, locate the coefficient estimates, standard errors, and p-values — these are the quantities you will report in a paper.
4. The multiple regression section introduces confounders: a concept critical to valid bioinformatics analysis.

## Common stumbling points

- **R-squared is not the whole story**: A high R² does not mean the model is correct or that the relationship is causal. It just means the model explains variance. Check diagnostics and think about biology.
- **Residuals vs. Fitted plot**: The single most important diagnostic. A clear pattern in residuals (curve, funnel shape) means your model is misspecified or assumptions are violated. Random scatter around zero is what you want.
- **Pearson vs. Spearman correlation**: Pearson measures linear association and requires roughly normal data. Spearman measures monotonic association using ranks and is appropriate for skewed data (read counts, survival times) or when you just need to know "do they go up and down together."
- **Confounders in multiple regression**: Adding a confounder variable to a model changes the coefficients of all other variables. This is not a bug — it is the model correctly partitioning variance. A coefficient's meaning is "the effect of this variable, *holding all other variables constant*."

---

## 1. Confidence Intervals

A **confidence interval** (CI) gives a range of plausible values for a population parameter based on sample data. A 95% CI means: if we repeated the experiment many times and computed a CI each time, 95% of those intervals would contain the true parameter.

In R, the correct CI depends on what you know:
- **Student-t CI**: sigma (population SD) is unknown — use the t-distribution with n-1 degrees of freedom. This is the most common case.
- **Chi-squared CI**: for a variance (sigma²) when data is normal.
- **Normal asymptotic CI**: for counts (Poisson mean) or proportions (binomial parameter) with large n.

```python
# Load IQ data
iq_data  <- read.table("../Assets/data/r_statistics/IQ.txt",
                       header = TRUE, sep = " ", dec = ".")
iq_data  <- iq_data[[1]]

n        <- length(iq_data)
iq_mean  <- mean(iq_data)           # sample mean
iq_s     <- sd(iq_data)             # corrected sample standard deviation
alpha    <- 0.10

cat("n =", n, "  mean =", round(iq_mean, 2), "  s =", round(iq_s, 2), "\n")
```python

```python
# Student-t confidence interval for the mean (sigma unknown)
student_q <- qt(1 - alpha/2, df = n - 1)

ci_left  <- iq_mean - student_q * iq_s / sqrt(n)
ci_right <- iq_mean + student_q * iq_s / sqrt(n)
cat(sprintf("90%% two-sided CI for mean: [%.2f, %.2f]\n", ci_left, ci_right))

# One-sided (right) CI: mean > lower bound
ci_right_one <- iq_mean + qt(1 - alpha, df = n - 1) * iq_s / sqrt(n)
cat(sprintf("90%% right-sided CI:  mean > %.2f\n",
            iq_mean - qt(1 - alpha, df = n - 1) * iq_s / sqrt(n)))
```python

```python
# t-test also returns a CI:
t.test(iq_data, alternative = "greater", mu = 110, conf.level = 1 - alpha)

# Manual t-statistic for H0: mu = 110
t_stat   <- (iq_mean - 110) * sqrt(n) / iq_s
iq_crit  <- qt(1 - alpha, df = n - 1)
cat(sprintf("t = %.3f  critical value = %.3f  reject H0: %s\n",
            t_stat, iq_crit, t_stat > iq_crit))
```python

```python
# Chi-squared confidence interval for variance
# Example: n=30 product quality measurements, sample variance = 0.3
n       <- 30
var_obs <- 0.3    # this is s^2, not sd
alpha   <- 0.05

# Test H0: sigma^2 = 0.2 against H1: sigma^2 > 0.2
T_stat  <- (n - 1) * var_obs / 0.2
crit    <- qchisq(1 - alpha, df = n - 1)
cat(sprintf("chi^2 statistic = %.2f  critical value = %.2f  reject H0: %s\n",
            T_stat, crit, T_stat > crit))

# Two-sided CI for sigma^2
ci_var_left  <- (n - 1) * var_obs / qchisq(1 - alpha/2, df = n - 1)
ci_var_right <- (n - 1) * var_obs / qchisq(alpha/2,     df = n - 1)
cat(sprintf("95%% CI for sigma^2: [%.4f, %.4f]\n", ci_var_left, ci_var_right))
```python

```python
# Normal asymptotic CI for a count parameter (Poisson mean)
# Example: store receives on average 256 customers; n=100 days observed
n         <- 100
mean_shop <- 256
alpha     <- 0.01

norm_q    <- qnorm(1 - alpha/2)
ci_left   <- mean_shop - norm_q * sqrt(mean_shop) / sqrt(n)
ci_right  <- mean_shop + norm_q * sqrt(mean_shop) / sqrt(n)
cat(sprintf("99%% CI for Poisson mean: [%.2f, %.2f]\n", ci_left, ci_right))
```python

```python
# Normal asymptotic CI for proportion (binomial parameter)
# Example: 180 of 300 patients improved
n      <- 300
scs    <- 180
alpha  <- 0.05
p_hat  <- scs / n

norm_q <- qnorm(1 - alpha/2)
ci_left  <- p_hat - norm_q * sqrt(p_hat * (1 - p_hat) / n)
ci_right <- p_hat + norm_q * sqrt(p_hat * (1 - p_hat) / n)
cat(sprintf("95%% CI for proportion: [%.4f, %.4f]\n", ci_left, ci_right))
```python

---

## 2. Parametric t-Tests

When the data can be assumed normally distributed, the **t-test** is more powerful than nonparametric alternatives. Use `var.test` to check variance equality before choosing between equal- and unequal-variance t-tests.

```python
# Paired t-test: weight before and after diet intervention
# (same data as Notebook 1 — comparing parametric and nonparametric results)
weight_data <- read.table("../Assets/data/r_statistics/Weight.txt",
                          header = TRUE, sep = " ", dec = ".")
before <- weight_data[[1]]
after  <- weight_data[[2]]

# H1: weight decreases after diet (after < before)
t.test(after, before, alternative = "less", mu = 0, paired = TRUE)
```python

```python
# Independent two-sample t-test: compare two towns
towns_data <- read.table("../Assets/data/r_statistics/town.txt",
                         header = TRUE, sep = " ", dec = ".")
skov    <- towns_data[[1]]
kastrul <- towns_data[[2]]

# Step 1: Test equality of variances (F-test)
var.test(skov, kastrul, ratio = 1, alternative = "two.sided")
```python

```python
# Step 2: If variances are unequal, use var.equal = FALSE (Welch's t-test)
t.test(skov, kastrul,
       alternative = "greater", mu = 0,
       paired = FALSE, var.equal = FALSE)
```python

```python
# Comparison: nonparametric Wilcoxon rank-sum on the same data
library(exactRankTests)
wilcox.exact(skov, kastrul, paired = FALSE, alternative = "greater", exact = TRUE)
```python

---

## 3. Normality Testing

Before applying parametric tests, check whether the data are consistent with a normal distribution.

| Method | Function | Notes |
|--------|----------|-------|
| Histogram | `hist` | Visual only |
| Q-Q plot | `qqnorm` + `qqline` | Points should fall on the line |
| Shapiro-Wilk | `shapiro.test` | Most powerful for $n \leq 5000$ |
| Pearson chi-squared GoF | `pearson.test` (nortest) | Requires choosing bin count |

A **low p-value** from `shapiro.test` is evidence *against* normality.

```python
# Load weather data: temperature measurements
weather_data <- read.table("../Assets/data/r_statistics/weather.txt",
                           header = TRUE, sep = "\t", dec = ".")

# Extract summer temperatures (months 6-8)
temp_summer <- as.numeric(weather_data$Temp[
    weather_data$Month >= 6 & weather_data$Month <= 8
])

cat("Summer temperature observations: n =", length(temp_summer), "\n")
cat("Mean:", round(mean(temp_summer), 1), " SD:", round(sd(temp_summer), 1), "\n")
```python

```python
# Visual normality checks
par(mfrow = c(1, 2))

hist(temp_summer, breaks = 15,
     main = "Summer Temperatures",
     xlab = "Temperature (°C)", col = "lightblue")

qqnorm(temp_summer, main = "Q-Q Plot: Summer Temperatures")
qqline(temp_summer, col = "red", lwd = 2)

par(mfrow = c(1, 1))
```python

```python
# Shapiro-Wilk test
shapiro.test(temp_summer)
```python

```python
# Pearson chi-squared goodness-of-fit test for normality
# install.packages("nortest")
library(nortest)

# adjust=FALSE: n-1 degrees of freedom (conservative)
# adjust=TRUE:  n-3 degrees of freedom (corrects for estimated mean and sd)
# True p-value lies between the two
pearson.test(temp_summer, adjust = FALSE)
pearson.test(temp_summer, adjust = TRUE)
```python

```python
# Second dataset: lake depth measurements
lake_data   <- read.table("../Assets/data/r_statistics/Lake.txt",
                          header = TRUE, sep = " ", dec = ".")
lake_column <- lake_data[[1]]

par(mfrow = c(1, 2))
hist(lake_column, breaks = 12, main = "Lake Depths",
     xlab = "Depth (m)", col = "lightgreen")
qqnorm(lake_column, main = "Q-Q Plot: Lake Depths")
qqline(lake_column, col = "red", lwd = 2)
par(mfrow = c(1, 1))

shapiro.test(lake_column)
pearson.test(lake_column, adjust = TRUE)
```python

---

## 4. Correlation Analysis

### Pearson Correlation

The Pearson coefficient $r$ measures linear association. It assumes both variables are normally distributed (for inference). The `cor.test` function tests $H_0: \rho = 0$ and returns a confidence interval via Fisher's z-transform.

### Spearman and Kendall Rank Correlations

Rank-based correlations detect **monotone** (not just linear) relationships and do not require normality.

```python
# Load country data: two continuous variables per country
corr_data <- read.table("../Assets/data/r_statistics/country.txt",
                        header = TRUE, sep = ";", dec = ",")

# Note: column 1 contains row numbers — skip it
x1 <- corr_data[[2]]
y1 <- corr_data[[3]]

plot(x1, y1,
     xlab = colnames(corr_data)[2],
     ylab = colnames(corr_data)[3],
     main = "Scatter Plot", pch = 20)
```python

```python
# Check normality before choosing correlation method
shapiro.test(x1)
shapiro.test(y1)
```python

```python
# Pearson correlation with 95% CI
cor.test(x1, y1, alternative = "two.sided", method = "pearson", conf.level = 0.95)
```python

```python
# Spearman rank correlation
cor.test(x1, y1, alternative = "two.sided", method = "spearman")

# Kendall tau
cor.test(x1, y1, alternative = "two.sided", method = "kendall")
```python

### 4.1 Fisher z-Transform CI for Pearson Correlation

The Fisher z-transform $z = \frac{1}{2}\ln\frac{1+r}{1-r}$ is approximately $N(\rho_z, 1/(n-3))$ for large $n$. This stabilizes the variance and allows direct CI construction.

```python
n  <- length(x1)

# Manual Pearson r
r  <- (mean(x1 * y1) - mean(x1) * mean(y1)) / (sd(x1) * sd(y1) * (n - 1) / n)
cat("Pearson r (manual):", round(r, 4), "\n")

# Fisher z-transform
z1       <- 0.5 * log((1 + r) / (1 - r))
alpha    <- 0.05
norm_q   <- qnorm(1 - alpha/2)

# CI bounds on the z scale
z1_left  <- z1 - norm_q / sqrt(n - 3)
z1_right <- z1 + norm_q / sqrt(n - 3)

# Back-transform to r scale
r_left   <- (exp(2 * z1_left)  - 1) / (exp(2 * z1_left)  + 1)
r_right  <- (exp(2 * z1_right) - 1) / (exp(2 * z1_right) + 1)

cat(sprintf("95%% CI for rho: [%.4f, %.4f]\n", r_left, r_right))
```python

```python
# Demonstration: different relationship shapes affect Pearson vs Spearman

# Scenario A: Pearson picks up linear relationship well
x2 <- abs(x1) + rnorm(length(x1), 0, 0.03)
y2 <- abs(y1) + rnorm(length(y1), 0, 0.02)
cat("Pearson r (abs transformation):",
    round(cor.test(x2, y2, method="pearson")$estimate, 3), "\n")
cat("Spearman rho (abs transformation):",
    round(cor.test(x2, y2, method="spearman")$estimate, 3), "\n\n")

# Scenario B: Quadratic relationship — Spearman detects monotone dependence
x3 <- x1 + rnorm(length(x1), 0, 0.03)
y3 <- x1^2 + rnorm(length(x1), 0, 0.02)
cat("Pearson r (quadratic):",
    round(cor.test(x3, y3, method="pearson")$estimate, 3), "\n")
cat("Spearman rho (quadratic):",
    round(cor.test(x3, y3, method="spearman")$estimate, 3), "\n")
```python

## Common Pitfalls

- **Assumption violations**: Check normality and homoscedasticity before parametric tests; use nonparametric alternatives when assumptions fail
- **Multiple comparisons**: Bonferroni is conservative; use Benjamini-Hochberg FDR for large-scale testing
- **Correlation ≠ causation**: Statistical association does not imply biological mechanism
