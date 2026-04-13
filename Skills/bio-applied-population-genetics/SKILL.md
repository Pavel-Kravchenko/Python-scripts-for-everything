---
name: bio-applied-population-genetics
description: "Split from `01_population_genetics_and_molecular_evolution.ipynb` to keep this topic self-contained."
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/15_Population_Genetics_and_Molecular_Evolution/01_population_genetics.ipynb"
---

# Population Genetics

*Source: Course notebook `Tier_3_Applied_Bioinformatics/15_Population_Genetics_and_Molecular_Evolution/01_population_genetics.ipynb`*

# Population Genetics

Split from `01_population_genetics_and_molecular_evolution.ipynb` to keep this topic self-contained.

**Navigation:** [Topic overview](./01_population_genetics_and_molecular_evolution.ipynb) · [Next: Molecular Evolution](./02_molecular_evolution.ipynb)

# Population Genetics and Molecular Evolution

**Tier 3 -- Applied Bioinformatics**

Population genetics explains how allele frequencies shift across generations under the competing forces of mutation, genetic drift, natural selection, and gene flow. Molecular evolution extends these ideas to DNA sequences, asking which substitutions are neutral, purifying, or adaptive.

**Learning objectives:**
- Compute allele frequencies and test Hardy-Weinberg equilibrium with chi-squared
- Simulate Wright-Fisher genetic drift and understand the effect of population size
- Model directional, balancing, purifying, and frequency-dependent selection
- Apply the molecular clock and Jukes-Cantor correction to estimate divergence times
- Calculate dN/dS ratios with the Nei-Gojobori method and interpret selection pressure
- Compute Tajima's D, McDonald-Kreitman neutrality index, Fst, and linkage disequilibrium

**Prerequisites:** Tier 2 Python, NumPy, basic probability  
**Libraries:** `numpy`, `pandas`, `matplotlib`, `scipy`

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from itertools import combinations

%matplotlib inline
plt.rcParams['figure.figsize'] = (12, 5)
plt.rcParams['font.size'] = 12
np.random.seed(42)
```

---

## 1. Allele Frequencies and Hardy-Weinberg Equilibrium

The **Hardy-Weinberg principle** (1908) states that, in an infinitely large, randomly mating population with no mutation, selection, or migration, allele frequencies remain constant across generations and genotype frequencies satisfy:

$$p^2 + 2pq + q^2 = 1$$

where *p* = frequency of allele A and *q* = 1 − *p* = frequency of allele a.

HWE is a null model: departures signal interesting biology.

```python
# Observed genotype counts from a sample of 500 individuals
# Locus: MN blood group (codominant: MM, MN, NN)
observed_MM = 233
observed_MN = 200
observed_NN = 67
n_individuals = observed_MM + observed_MN + observed_NN

# Allele frequencies
n_M_alleles = 2 * observed_MM + observed_MN
n_N_alleles = 2 * observed_NN + observed_MN
total_alleles = 2 * n_individuals

p_M = n_M_alleles / total_alleles   # frequency of M
q_N = n_N_alleles / total_alleles   # frequency of N

print(f"Sample size (n): {n_individuals} individuals")
print(f"Allele counts  -- M: {n_M_alleles}, N: {n_N_alleles}")
print(f"Allele frequencies -- p(M) = {p_M:.4f}, q(N) = {q_N:.4f}")
print(f"p + q = {p_M + q_N:.4f}  (must equal 1)")
```

```python
# HWE expected genotype frequencies and bar chart comparison
expected_MM = p_M**2 * n_individuals
expected_MN = 2 * p_M * q_N * n_individuals
expected_NN = q_N**2 * n_individuals

genotypes  = ['MM', 'MN', 'NN']
obs_counts = [observed_MM, observed_MN, observed_NN]
exp_counts = [expected_MM, expected_MN, expected_NN]

x = np.arange(len(genotypes))
width = 0.35
fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(x - width/2, obs_counts, width, label='Observed', color='steelblue')
ax.bar(x + width/2, exp_counts, width, label='HWE Expected', color='coral', alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(genotypes)
ax.set_ylabel('Count')
ax.set_title('MN Blood Group: Observed vs HWE Expected')
ax.legend()
plt.tight_layout()
plt.show()

print(f"  MM: observed={observed_MM}, expected={expected_MM:.1f}")
print(f"  MN: observed={observed_MN}, expected={expected_MN:.1f}")
print(f"  NN: observed={observed_NN}, expected={expected_NN:.1f}")
```

```python
# Chi-squared test for HWE departure
# df = (number of genotype classes - 1) - (independent allele freqs estimated)
# For biallelic locus: df = 3 - 1 - 1 = 1
observed_arr = np.array([observed_MM, observed_MN, observed_NN], dtype=float)
expected_arr = np.array([expected_MM, expected_MN, expected_NN], dtype=float)

chi2_stat = np.sum((observed_arr - expected_arr)**2 / expected_arr)
df = 1
p_value = stats.chi2.sf(chi2_stat, df)   # survival function = 1 - CDF

print(f"Chi-squared statistic: {chi2_stat:.4f}")
print(f"Degrees of freedom:    {df}")
print(f"p-value:               {p_value:.4f}")
print()
if p_value > 0.05:
    print("Conclusion: No significant departure from HWE (p > 0.05).")
else:
    print("Conclusion: Significant departure from HWE (p <= 0.05).")
```

### Causes of HWE Departure

| Force | Effect on allele freq | Effect on genotype freq |
|---|---|---|
| Natural selection | Changes p and q | Excess/deficit of certain genotypes |
| Genetic drift | Changes p and q (random) | Excess homozygotes (inbreeding-like) |
| Gene flow (migration) | Introduces new alleles | Shifts all genotype freqs |
| Mutation | Slow, negligible per generation | Tiny shift toward new alleles |
| Non-random mating | No direct effect on p, q | Excess homozygotes (inbreeding) or excess heterozygotes |
| Population stratification | No effect within strata | Wahlund effect: excess homozygotes in pooled sample |

The **inbreeding coefficient F** measures the fractional reduction in heterozygosity:
$$F = 1 - \frac{H_{\text{obs}}}{2pq}$$

```python
# Inbreeding coefficient F and multi-allelic HWE extension
H_obs = observed_MN / n_individuals
H_exp = 2 * p_M * q_N
F_inbreeding = 1 - H_obs / H_exp

print(f"Inbreeding coefficient F = {F_inbreeding:.4f}")
print(f"  H_obs={H_obs:.4f}, H_exp={H_exp:.4f}")
print("  F ≈ 0: no evidence of inbreeding or outbreeding.")
print()

# Multi-allelic HWE for a 3-allele ABO-like locus
p_A, p_B, p_O = 0.28, 0.06, 0.66
allele_freqs = {'A': p_A, 'B': p_B, 'O': p_O}
alleles = list(allele_freqs.keys())

print("Expected genotype frequencies under HWE (3-allele ABO-like locus):")
total = 0.0
for i, a1 in enumerate(alleles):
    for a2 in alleles[i:]:
        freq = allele_freqs[a1]**2 if a1 == a2 else 2*allele_freqs[a1]*allele_freqs[a2]
        total += freq
        print(f"  {a1}{a2}: {freq:.4f}")
print(f"  Sum: {total:.4f}")
```

---

## 2. Genetic Drift

**Genetic drift** is the random change in allele frequency caused by sampling a finite number of gametes each generation. It is the dominant evolutionary force in small populations.

### Wright-Fisher Model

In a diploid population of size *N*, the next generation's count of allele A is drawn from a binomial distribution:

$$k_{t+1} \sim \text{Binomial}(2N,\; p_t)$$

Variance per generation: $\text{Var}(p_{t+1}) = p_t(1-p_t)/(2N)$

The allele frequency performs a random walk bounded by absorbing states at 0 (loss) and 1 (fixation). Smaller N means larger variance and faster drift.

### Absorbing States: Fixation and Loss

Once an allele reaches frequency 0 or 1, it stays there forever (in the absence of mutation).
These are called **absorbing states** of the Wright-Fisher Markov chain.

- **Fixation** (p = 1): the allele has replaced all others — it is now the only allele at this locus
- **Loss** (p = 0): the allele has been eliminated from the population

For a neutral allele starting at frequency p₀ in a population of size N:
- P(fixation) = p₀
- P(loss) = 1 − p₀
- Expected time to fixation (given fixation) ≈ −4N × [p₀ ln(p₀) + (1-p₀) ln(1-p₀)] / p₀

```python
def wright_fisher_trajectory(N, p0, n_generations, rng):
    """Simulate one Wright-Fisher trajectory.
    Returns array of allele frequencies, length n_generations+1."""
    freq = np.empty(n_generations + 1)
    freq[0] = p0
    p = p0
    for g in range(n_generations):
        k = rng.binomial(2 * N, p)
        p = k / (2 * N)
        freq[g + 1] = p
        if p == 0.0 or p == 1.0:
            freq[g + 2:] = p   # absorbing state: stay at 0 or 1
            break
    return freq

rng = np.random.default_rng(42)

# Single trajectory for three different population sizes
n_gen = 200
p_initial = 0.5
pop_sizes = [20, 100, 1000]
colors = ['firebrick', 'steelblue', 'seagreen']

fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
for ax, N, color in zip(axes, pop_sizes, colors):
    traj = wright_fisher_trajectory(N, p_initial, n_gen, rng)
    ax.plot(traj, color=color, lw=1.5)
    ax.axhline(0.0, ls='--', color='black', lw=0.8)
    ax.axhline(1.0, ls='--', color='black', lw=0.8)
    ax.set_title(f'N = {N}')
    ax.set_xlabel('Generation')
    ax.set_ylim(-0.05, 1.05)
axes[0].set_ylabel('Allele frequency (p)')
fig.suptitle('Wright-Fisher Drift: Single Trajectory per Population Size')
plt.tight_layout()
plt.show()
```

```python
# Many replicates: visualise spread, fixation, and loss
# Under neutral drift, P(fixation) = p0 (Kimura 1962)
n_replicates = 50
n_gen = 150
N_small = 50
p0 = 0.3

trajectories = np.array([
    wright_fisher_trajectory(N_small, p0, n_gen, rng)
    for _ in range(n_replicates)
])

n_fixed = np.sum(trajectories[:, -1] == 1.0)
n_lost  = np.sum(trajectories[:, -1] == 0.0)
n_poly  = n_replicates - n_fixed - n_lost

fig, ax = plt.subplots(figsize=(12, 5))
for traj in trajectories:
    final = traj[-1]
    col = 'steelblue' if final == 1.0 else ('tomato' if final == 0.0 else 'lightgray')
    ax.plot(traj, color=col, alpha=0.6, lw=0.9)

ax.axhline(p0, ls=':', color='black', lw=1.2, label=f'Initial p = {p0}')
ax.set_xlabel('Generation')
ax.set_ylabel('Allele frequency')
ax.set_title(f'N={N_small}, p0={p0}: {n_replicates} replicates | '
             f'Fixed={n_fixed} (blue), Lost={n_lost} (red), Polymorphic={n_poly} (gray)')
ax.legend()
plt.tight_layout()
plt.show()
print(f"Theoretical P(fixation) = p0 = {p0:.2f}")
print(f"Observed fixation fraction   = {n_fixed/n_replicates:.2f}")
```

```python
# Time to fixation or loss as a function of N
# Theory: E[T_fix or loss] ≈ -4N[p*ln(p) + (1-p)*ln(1-p)]

def mean_time_to_fixation_loss(N, p0, n_replicates, rng, max_gen=10_000):
    times = []
    for _ in range(n_replicates):
        p = p0
        for g in range(1, max_gen + 1):
            p = rng.binomial(2 * N, p) / (2 * N)
            if p == 0.0 or p == 1.0:
                times.append(g)
                break
    return np.mean(times) if times else np.nan

pop_sizes_test = [10, 25, 50, 100, 200]
mean_times = [mean_time_to_fixation_loss(N, 0.5, 200, rng) for N in pop_sizes_test]

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(pop_sizes_test, mean_times, 'o-', color='steelblue', lw=2, label='Simulated')
theory_x = np.linspace(10, 200, 100)
theory_y = -4 * theory_x * (0.5 * np.log(0.5) + 0.5 * np.log(0.5))
ax.plot(theory_x, theory_y, '--', color='coral', label='Theory: −4N[p ln p + q ln q]')
ax.set_xlabel('Population size (N)')
ax.set_ylabel('Mean generations to fixation/loss')
ax.set_title('Time to Fixation or Loss vs Population Size (p₀ = 0.5)')
ax.legend()
plt.tight_layout()
plt.show()
```

```python
# Bottleneck: simulating a transient reduction in Ne

def wf_with_bottleneck(N_normal, N_bottleneck, bottleneck_start, bottleneck_duration,
                        p0, n_gen, rng):
    """Wright-Fisher drift with a population size bottleneck."""
    freq = np.empty(n_gen + 1)
    freq[0] = p0
    p = p0
    for g in range(n_gen):
        in_bottleneck = bottleneck_start <= g < bottleneck_start + bottleneck_duration
        N = N_bottleneck if in_bottleneck else N_normal
        p = rng.binomial(2 * N, p) / (2 * N)
        freq[g + 1] = p
        if p == 0.0 or p == 1.0:
            freq[g + 2:] = p
            break
    return freq

n_gen = 200
n_reps = 30
bottle_trajs = [wf_with_bottleneck(500, 5, 50, 10, 0.5, n_gen, rng) for _ in range(n_reps)]
normal_trajs = [wright_fisher_trajectory(500, 0.5, n_gen, rng) for _ in range(n_reps)]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
for traj in normal_trajs:
    ax1.plot(traj, color='steelblue', alpha=0.4, lw=0.8)
ax1.set_title('No bottleneck (N=500 throughout)')

for traj in bottle_trajs:
    ax2.plot(traj, color='tomato', alpha=0.4, lw=0.8)
ax2.axvspan(50, 60, color='orange', alpha=0.2, label='Bottleneck (N=5, 10 gen)')
ax2.set_title('Bottleneck at generation 50-60 (N drops to 5)')
ax2.legend()

for ax in (ax1, ax2):
    ax.set_xlabel('Generation')
    ax.set_ylabel('Allele frequency')
plt.suptitle('Bottleneck Effects on Drift', y=1.01)
plt.tight_layout()
plt.show()
```
