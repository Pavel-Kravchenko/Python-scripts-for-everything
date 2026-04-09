# 0.08 Advanced R Statistics

**Tier 0: Computational Foundations**

Two hands-on notebooks covering the full range of classical hypothesis testing and regression methods implemented in R — from exact nonparametric tests and robust estimation to linear regression, ANOVA, and simulation-based demonstrations of fundamental statistical theorems. All examples use biological and medical datasets.

## Topics Covered

- R distribution function convention: `d/p/q/r` prefixes for `norm`, `binom`, `t`, `chisq`, `exp`, and more
- Exact binomial test, asymptotic approximation, power function construction, sample size determination
- Sign test, Wilcoxon signed-rank test, and manual rank-statistic calculation for paired data
- Wilcoxon rank-sum / Mann-Whitney U test for independent samples; manual U-statistic derivation
- Hodges-Lehmann robust location estimation with confidence intervals via `wilcox.test(conf.int=TRUE)`
- Kruskal-Wallis multi-group nonparametric test and Dunn post-hoc pairwise comparisons
- Student-t, chi-squared, and normal asymptotic confidence intervals for mean, variance, and proportion
- Paired and unpaired t-tests; F-test for variance equality before choosing test variant
- Normality assessment: Shapiro-Wilk, Pearson chi-squared goodness-of-fit, Q-Q plots
- Pearson, Spearman, and Kendall correlation with Fisher z-transform confidence intervals
- Chi-squared tests of independence and homogeneity; `xtabs`, `ftable`, `prop.test`
- Simple and multiple linear regression with `lm`; polynomial degree selection; residual diagnostics
- ANOVA via `lm` + `anova` compared directly to Kruskal-Wallis on the same dataset
- Law of Large Numbers and Central Limit Theorem demonstrated by simulation; CI coverage experiment

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_r_hypothesis_testing_and_nonparametrics.ipynb](01_r_hypothesis_testing_and_nonparametrics.ipynb) | Exact binomial, sign test, Wilcoxon signed-rank and rank-sum, Hodges-Lehmann estimation, Kruskal-Wallis with Dunn post-hoc, and power analysis |
| [02_r_regression_correlation_and_diagnostics.ipynb](02_r_regression_correlation_and_diagnostics.ipynb) | Confidence intervals, parametric t-tests, normality testing, correlation analysis, chi-squared contingency table methods, linear and polynomial regression, ANOVA, and LLN/CLT simulation |

## Prerequisites

- [0.06 Biostatistics Fundamentals](../06_Biostatistics/)
- [0.07 Probability and Statistics Python](../07_Probability_and_Statistics_Python/)

---

[← Previous: 0.07 Probability and Statistics Python](../07_Probability_and_Statistics_Python/) | [Course Overview](../../README.md)
