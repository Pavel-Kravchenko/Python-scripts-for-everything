---
name: bio-applied-population-genetics
description: Population genetics — Hardy-Weinberg equilibrium, Wright-Fisher drift simulation, selection models, molecular clock, dN/dS, Tajima's D, Fst, and linkage disequilibrium
tool_type: python
primary_tool: Matplotlib
---

# Population Genetics and Molecular Evolution

## Core Evolutionary Forces

| Force | Effect on allele freq | Effect on diversity |
|-------|-----------------------|--------------------|
| Genetic drift | Random walk (stronger in small N) | Reduces diversity |
| Natural selection | Directional shift | Reduces (purifying/directional) or maintains (balancing) |
| Mutation | Slow increase of new alleles | Increases diversity |
| Gene flow (migration) | Homogenizes populations | Increases local diversity |
| Non-random mating | No direct effect on p/q | Excess homozygotes (inbreeding) |
| Population stratification | No effect within strata | Wahlund effect: excess homozygotes in pooled sample |

## Hardy-Weinberg Equilibrium

HWE predicts genotype frequencies from allele frequencies in an ideal population:
p² (AA) + 2pq (Aa) + q² (aa) = 1

```python
from scipy import stats
import numpy as np

# Allele frequencies from genotype counts
n = obs_MM + obs_MN + obs_NN
p = (2 * obs_MM + obs_MN) / (2 * n)
q = 1 - p

# Expected counts
exp_MM, exp_MN, exp_NN = p**2 * n, 2*p*q * n, q**2 * n

# Chi-squared test (df=1 for biallelic locus)
chi2 = sum((o - e)**2 / e for o, e in
           zip([obs_MM, obs_MN, obs_NN], [exp_MM, exp_MN, exp_NN]))
p_val = stats.chi2.sf(chi2, df=1)

# Inbreeding coefficient
F = 1 - (obs_MN / n) / (2 * p * q)
# F > 0: excess homozygotes; F < 0: excess heterozygotes
```

## Wright-Fisher Drift Simulation

```python
def wright_fisher_trajectory(N, p0, n_generations, rng):
    """Simulate one WF trajectory. Returns allele frequency array."""
    freq = np.empty(n_generations + 1)
    freq[0] = p0
    p = p0
    for g in range(n_generations):
        k = rng.binomial(2 * N, p)
        p = k / (2 * N)
        freq[g + 1] = p
        if p == 0.0 or p == 1.0:
            freq[g + 2:] = p   # absorbing state
            break
    return freq

rng = np.random.default_rng(42)
# P(fixation) of neutral allele = p0 (Kimura 1962)
# Mean time to fixation/loss ~ -4N[p*ln(p) + q*ln(q)] generations
```

### Population Size Effects

| N | Drift strength | Time to fixation |
|---|---------------|-----------------|
| Small (10–50) | Strong | Few generations |
| Large (1000+) | Weak | ~4N generations |
| Bottleneck | Extreme during event | Permanent diversity loss |

## Selection Models

```python
# Directional selection (additive): one allele increases fitness
def selection_trajectory(N, p0, s, h, n_gen, rng):
    """s = selection coeff; h = dominance (0.5 = additive, 1 = dominant)."""
    p = p0
    traj = [p]
    for _ in range(n_gen):
        # Mean fitness
        w_bar = p**2*(1+s) + 2*p*(1-p)*(1+h*s) + (1-p)**2
        # Allele freq after selection
        p_sel = (p**2*(1+s) + p*(1-p)*(1+h*s)) / w_bar
        # Drift
        p = rng.binomial(2*N, p_sel) / (2*N)
        traj.append(p)
        if p in (0.0, 1.0): break
    return traj
```

## Molecular Clock and dN/dS

```python
# Jukes-Cantor correction for multiple hits
def jukes_cantor(p_diff):
    """Correct observed proportion of differences for multiple hits."""
    if p_diff >= 0.75: return float('inf')
    return -0.75 * np.log(1 - 4 * p_diff / 3)

# Divergence time estimate
# t = d / (2 * rate)  where d = JC-corrected distance, rate = subs/site/year

# dN/dS (Nei-Gojobori method)
# dN/dS < 1 = purifying selection (most genes)
# dN/dS ≈ 1 = neutral evolution
# dN/dS > 1 = positive/adaptive selection (rare)
```

## Summary Statistics

```python
# Tajima's D (neutrality test)
# D > 0: excess intermediate-frequency variants (balancing selection / population contraction)
# D < 0: excess rare variants (selective sweep / population expansion)
# |D| > 2 is typically significant

# Fst (population differentiation)
# Fst = (Ht - Hs) / Ht
# Fst = 0: no differentiation; Fst = 1: complete differentiation
def fst(p1, p2, n1, n2):
    p_bar = (n1*p1 + n2*p2) / (n1 + n2)
    ht = 2 * p_bar * (1 - p_bar)
    hs = (2*p1*(1-p1)*n1 + 2*p2*(1-p2)*n2) / (n1 + n2)
    return (ht - hs) / ht if ht > 0 else 0.0

# Linkage disequilibrium (D')
# r² = D² / (p_A * q_A * p_B * q_B)
# r² > 0.8 = strong LD (blocks)
```

## Pitfalls

- **HWE departure does not identify the cause**: chi-squared only flags departure; genotyping errors, selection, inbreeding, and stratification all produce the same test result — need additional tests to distinguish
- **Wahlund effect**: pooling samples from differentiated subpopulations creates apparent HWE departure (excess homozygotes) without any selection or inbreeding — always stratify by population before HWE testing
- **Wright-Fisher assumes non-overlapping generations**: does not model age structure; use coalescent models or continuous-time approximations for organisms with overlapping generations
- **Jukes-Cantor correction breaks down at p > 0.5**: multiple substitutions per site saturate; use more sophisticated models (HKY, GTR) for highly diverged sequences
- **dN/dS > 1 requires averaging over many sites**: a single gene rarely achieves dN/dS > 1 genome-wide; use codon-level tests (PAML codeml, HyPhy) to detect episodic selection at specific sites
- **Tajima's D is sensitive to demographic history**: population expansion gives D < 0 (mimics sweeps); always interpret with a demographic null model
