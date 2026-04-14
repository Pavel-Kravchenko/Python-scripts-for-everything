---
name: bio-applied-phylodynamics
description: Viral phylodynamics — molecular clocks, root-to-tip regression, time-scaled phylogenies, Bayesian skyline plots, and phylogeography tool selection
tool_type: python
primary_tool: Matplotlib
---

# Viral Phylodynamics and Molecular Epidemiology

## Substitution Rates (Common Viruses)

| Virus | Rate (subs/site/year) | Notes |
|-------|-----------------------|-------|
| SARS-CoV-2 | ~1 × 10⁻³ | ~1–2 SNPs per 2 weeks |
| Influenza A (HA) | 3–5 × 10⁻³ | High immune selection |
| HIV-1 | 2–4 × 10⁻³ | Highest among RNA viruses |
| Hepatitis C | 1–2 × 10⁻³ | Chronic, slow |
| Measles | ~5 × 10⁻⁴ | More conserved |

## Tool Selection

| Tool | Approach | Speed | Use Case |
|------|----------|-------|----------|
| **TreeTime** | ML + least squares | Fast | Real-time surveillance (Nextstrain) |
| **BEAST2** | Bayesian MCMC | Slow | Rigorous inference, posterior distributions |
| **IQ-TREE + LSD2** | Least-squares dating | Fast | Large datasets |
| **Augur** (Nextstrain) | TreeTime wrapper | Fast | Integrated pipeline |

## Clock Models

- **Strict clock**: single global rate — too restrictive for long time scales or cross-species comparisons
- **Relaxed clock (UCLN)**: each branch draws rate from lognormal distribution; standard in BEAST2; captures rate heterogeneity

## Root-to-Tip Regression

Validates clock signal before full Bayesian analysis. R² > 0.70 is generally sufficient.

```python
from scipy import stats
import numpy as np

# root_to_tip = genetic distance from root; sampling_dates = decimal years
slope, intercept, r, p, se = stats.linregress(sampling_dates, root_to_tip)

# Estimated clock rate
print(f"Rate: {slope:.2e} subs/site/year")
print(f"R² = {r**2:.3f}  (>0.70 = good clock signal)")
print(f"Expected SNPs/2 weeks: {slope * (14/365) * genome_len:.1f}")
```

## Nextstrain Augur Workflow

```bash
# Filter and subsample
augur filter --sequences seqs.fasta --metadata metadata.tsv \
    --min-date 2020-01-01 --subsample-max-sequences 300 \
    --output filtered.fasta

# Align to reference
augur align --sequences filtered.fasta \
    --reference-sequence reference.gb --output aligned.fasta

# Build ML tree
augur tree --alignment aligned.fasta --output tree_raw.nwk

# Time-scale
augur refine --tree tree_raw.nwk --alignment aligned.fasta \
    --metadata metadata.tsv --timetree --coalescent opt \
    --output-tree timetree.nwk --output-node-data branch_lengths.json

# Export for Auspice
augur export v2 --tree timetree.nwk --output auspice.json
```

## TMRCA and Effective Population Size

```python
# Simulate time-scaled phylogeny (TreeTime-style)
mrca_date = 2019.8
clock_rate = 1e-3
genome_len = 29903

n_snps_from_root = np.random.poisson(
    clock_rate * (sampling_dates - mrca_date) * genome_len
)
```

## Bayesian Skyline Plot Interpretation

The BSP estimates viral effective population size (Ne) over time from coalescent events in a Bayesian MCMC framework (BEAST2).

| Ne Trend | Epidemiological Meaning |
|----------|------------------------|
| Rising Ne | Exponential epidemic growth |
| Plateau | Endemic equilibrium |
| Sharp decline | Intervention, seasonal end, or population immunity |

N_e is a **genetic diversity metric** (coalescent rate), not actual viral particle count.

## Phylogeography Approaches

| Method | Tool | Use Case |
|--------|------|----------|
| Discrete (DTA) | BEAST2 | Location states per node, directional migration rates |
| Continuous | BEAST2 | Lat/lon Brownian motion, diffusion coefficients |
| Parsimony/ML | Nextstrain | Fast, visual — animated Auspice maps |

Key metrics: migration rate matrix (transitions/year), Bayes factor (evidence for routes), geographic diffusion coefficient (km²/year).

## Pitfalls

- **R² > 0.70 is necessary but not sufficient**: a high R² from root-to-tip regression can result from non-clock-like structure (e.g., recombination, strong selection) — inspect the residuals
- **Strict clock for long time scales**: rates vary across lineages; use relaxed clock (UCLN) in BEAST2 when analyzing datasets spanning years or different host species
- **Ne ≠ actual viral population size**: BSP estimates the *effective* population size from coalescent timing; large Ne in highly sampled outbreaks reflects sampling density, not census size
- **Sampling bias distorts phylogeography**: oversampling one location creates spurious migration routes toward that location — always subsample with `augur filter`
- **BEAST2 MCMC convergence**: ESS > 200 for all parameters is required; check in Tracer before interpreting posteriors; insufficient chain length is the most common mistake
- **Recombinant viruses violate clock assumptions**: HIV, coronaviruses have recombination — screen with RDP4 or ClonalFrameML before applying molecular clock models
