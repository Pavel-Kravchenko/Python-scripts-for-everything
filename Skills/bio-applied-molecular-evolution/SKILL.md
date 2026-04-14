---
name: bio-applied-molecular-evolution
description: Population genetics and molecular evolution — Hardy-Weinberg, Wright-Fisher drift, selection models, dN/dS, Tajima's D, Fst, and the neutral theory
tool_type: python
primary_tool: Matplotlib
---

# Population Genetics and Molecular Evolution

## Core Forces Summary

| Force | Effect on allele freq | Effect on diversity |
|-------|-----------------------|---------------------|
| Genetic drift | Random walk → fixation or loss | Reduces |
| Natural selection (purifying) | Removes deleterious alleles | Reduces at site |
| Natural selection (positive) | Drives beneficial allele to fixation | Selective sweep reduces nearby diversity |
| Balancing selection | Maintains multiple alleles | Increases |
| Mutation | Creates new alleles | Increases |
| Gene flow / migration | Homogenizes populations | Reduces Fst |

## Hardy-Weinberg Equilibrium

$$p^2 + 2pq + q^2 = 1 \quad \text{(biallelic locus)}$$

HWE is a null model. Departures signal: selection, inbreeding, drift, admixture, or genotyping error.

```python
import numpy as np
from scipy import stats

def hwe_test(obs_AA, obs_Aa, obs_aa):
    n = obs_AA + obs_Aa + obs_aa
    p = (2*obs_AA + obs_Aa) / (2*n)
    q = 1 - p
    exp = [p**2*n, 2*p*q*n, q**2*n]
    obs = [obs_AA, obs_Aa, obs_aa]
    chi2 = sum((o-e)**2/e for o, e in zip(obs, exp))
    p_val = stats.chi2.sf(chi2, df=1)   # 1 df for biallelic locus
    F = 1 - (obs_Aa/n) / (2*p*q)       # inbreeding coefficient
    return chi2, p_val, F

# F ≈ 0: no inbreeding; F > 0: excess homozygotes; F < 0: excess heterozygotes
```

## Wright-Fisher Drift

$$k_{t+1} \sim \text{Binomial}(2N,\; p_t), \quad \text{Var}(p_{t+1}) = \frac{p_t(1-p_t)}{2N}$$

- Absorbing states: p=0 (loss) and p=1 (fixation)
- P(fixation of neutral allele) = p₀
- Expected time to fixation ≈ −4N[p₀ ln p₀ + (1−p₀) ln(1−p₀)] / p₀

```python
def wright_fisher_trajectory(N, p0, n_gen, rng):
    freq = np.empty(n_gen + 1)
    freq[0] = p = p0
    for g in range(n_gen):
        p = rng.binomial(2*N, p) / (2*N)
        freq[g+1] = p
        if p in (0.0, 1.0):
            freq[g+2:] = p
            break
    return freq

def wf_with_bottleneck(N_normal, N_bottle, bottle_start, bottle_dur, p0, n_gen, rng):
    freq = np.empty(n_gen + 1); freq[0] = p = p0
    for g in range(n_gen):
        N = N_bottle if bottle_start <= g < bottle_start + bottle_dur else N_normal
        p = rng.binomial(2*N, p) / (2*N); freq[g+1] = p
        if p in (0.0, 1.0): freq[g+2:] = p; break
    return freq
```

## Effective Population Size (Ne)

Ne << census N due to:
- Unequal sex ratio: `Ne = 4*Nm*Nf / (Nm + Nf)`
- Variance in reproductive success
- Bottlenecks: Ne ≈ harmonic mean of per-generation sizes
- Human Ne ≈ 10,000 despite 8 billion census size

## Molecular Clock and Jukes-Cantor

```python
def jukes_cantor_distance(p_observed):
    """p: proportion of observed differences. Returns JC-corrected distance."""
    if p_observed >= 0.75:
        return float('inf')
    return -0.75 * np.log(1 - 4*p_observed/3)

# Divergence time: t = d / (2 * mu)  where mu = substitution rate per site per generation
```

## dN/dS (Nei-Gojobori)

- dN/dS < 1: purifying selection (most protein-coding genes)
- dN/dS ≈ 1: neutral evolution
- dN/dS > 1: positive (adaptive) selection

## Summary Statistics

| Statistic | Formula (simplified) | Interpretation |
|-----------|---------------------|----------------|
| Tajima's D | (π − θ_W) / std | D<0: excess rare variants (selection/expansion); D>0: balancing selection/bottleneck |
| McDonald-Kreitman NI | (Dn/Ds) / (Pn/Ps) | NI<1: positive selection; NI>1: slightly deleterious |
| Fst | (Ht − Hs) / Ht | 0=no differentiation; 1=complete differentiation |
| LD (r²) | corr(hap_A, hap_B)² | Decays with recombination and time |

```python
def tajimas_d(sequences):
    """sequences: list of same-length strings."""
    n = len(sequences)
    L = len(sequences[0])
    S = sum(1 for i in range(L) if len(set(s[i] for s in sequences)) > 1)
    pi = sum(sum(a!=b for a,b in zip(sequences[i],sequences[j]))
             for i in range(n) for j in range(i+1,n)) / (L * n*(n-1)/2)
    a1 = sum(1/i for i in range(1,n))
    theta_w = S / a1 / L
    # simplified D (full formula needs variance term)
    return (pi - theta_w)

def fst(allele_freq_pop1, allele_freq_pop2):
    """Two-population Fst for a biallelic locus."""
    ht = 2 * ((allele_freq_pop1+allele_freq_pop2)/2) * (1-(allele_freq_pop1+allele_freq_pop2)/2)
    hs = (2*allele_freq_pop1*(1-allele_freq_pop1) + 2*allele_freq_pop2*(1-allele_freq_pop2)) / 2
    return (ht - hs) / ht if ht > 0 else 0.0
```

## Pitfalls

- **HWE departures from genotyping error**: batch-specific or low-quality genotypes produce systematic HWE failures — always check per-batch before biological interpretation
- **Ne ≠ census N**: using census size in drift calculations severely underestimates drift; human Ne ≈ 10,000
- **Tajima's D is confounded**: demographic history (expansion, bottleneck) produces the same signature as selection — always interpret alongside a demographic null model
- **JC correction breaks down**: the Jukes-Cantor formula assumes equal base frequencies and single substitutions; fails for p > 0.6 (returns inf) and for highly diverged sequences
- **dN/dS at single sites is noisy**: use branch-site or codon models (PAML, HyPhy) rather than site-by-site ratios
- **Multiple testing**: apply Benjamini-Hochberg when scanning many loci for selection signals
