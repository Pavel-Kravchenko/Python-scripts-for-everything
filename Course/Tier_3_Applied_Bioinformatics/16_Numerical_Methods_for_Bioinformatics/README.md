# 3.16 Numerical Methods for Bioinformatics
**Tier 3: Applied Bioinformatics**

Numerical methods are the computational backbone of quantitative biology. This module translates the classical mathematical apparatus — polynomial interpolation, finite differences, quadrature, least squares, iterative optimisation, and the Fourier transform — into practical Python tools for biological data analysis. The material is grounded in the lecture course taught at Moscow State University (Faculty of Bioengineering and Bioinformatics), where these techniques were motivated by Mendeleev's 1881 solubility experiments and extended to modern problems such as circadian rhythm detection, enzyme kinetics fitting, RNA-seq count modelling, and X-ray crystallography.

## Topics Covered

- **Interpolation:** Lagrange polynomial, Newton divided differences, Runge phenomenon, Chebyshev nodes, cubic splines — reconstructing missing time points in gene expression data
- **Numerical differentiation:** forward, backward, and central finite differences; error analysis; instability with noisy measurements; drug clearance rate calculation
- **Numerical integration:** trapezoidal rule, Simpson's rule, error bounds; area under dose-response curves, ROC-AUC, pharmacokinetic AUC
- **Linear least squares:** design matrix, normal equations, SVD via `numpy.linalg.lstsq`, condition number; polynomial regression on biological data
- **Nonlinear curve fitting:** `scipy.optimize.curve_fit` (Levenberg-Marquardt); Michaelis-Menten, Hill equation, four-parameter logistic (4PL) models
- **Goodness-of-fit:** R², chi-squared test, AIC/AICc; model selection between kinetic models
- **Optimisation:** gradient descent from scratch (with convergence analysis); coordinate descent; Nelder-Mead; `scipy.optimize.minimize`; maximum likelihood estimation for normal and negative binomial distributions
- **Fourier transform:** continuous FT, Fourier series, convolution theorem, uncertainty principle; DFT definition; Nyquist criterion and aliasing; FFT via `numpy.fft`; power spectrum analysis; low-pass Fourier filtering; 2D FFT for crystallography

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_numerical_methods_for_bioinformatics.ipynb](01_numerical_methods_for_bioinformatics.ipynb) | Complete coverage of interpolation, differentiation, integration, least squares, nonlinear fitting, optimisation, and FFT — all with biological examples |

## Prerequisites

- Tier 1: Python fundamentals, NumPy arrays, Matplotlib
- Tier 1: Pandas data wrangling
- Basic calculus concepts (derivatives, integrals, Taylor series)
- Familiarity with enzyme kinetics notation is helpful but not required (see module 3.13)

---
[← Previous: 3.15 Population Genetics](../15_Population_Genetics_and_Molecular_Evolution/) | [Course Overview](../../README.md) | [Next: 3.17 Genome Assembly →](../17_Genome_Assembly_and_Advanced_NGS/)
