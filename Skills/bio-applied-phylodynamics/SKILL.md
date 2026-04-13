---
name: bio-applied-phylodynamics
description: "**Tier 3 — Applied Bioinformatics | Module 37 · Notebook 2**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/37_Virology_Bioinformatics/02_phylodynamics.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Viral Phylodynamics and Molecular Epidemiology

*Source: Course notebook `Tier_3_Applied_Bioinformatics/37_Virology_Bioinformatics/02_phylodynamics.ipynb`*

# Viral Phylodynamics and Molecular Epidemiology

**Tier 3 — Applied Bioinformatics | Module 37 · Notebook 2**

*Prerequisites: Notebook 1 (Viral Genome Assembly), Module 06 (Phylogenetics)*

---

**By the end of this notebook you will be able to:**
1. Build timed phylogenetic trees using Nextstrain / Augur pipeline
2. Estimate effective population size through time with coalescent models (BEAST2)
3. Identify transmission clusters and superspreader events
4. Perform ancestral state reconstruction for geographic spread
5. Interpret R_e (effective reproduction number) from genomic data



**Key resources:**
- [Nextstrain documentation](https://docs.nextstrain.org/)
- [BEAST2 tutorials](https://www.beast2.org/tutorials/)
- [Treetime documentation](https://treetime.readthedocs.io/)
- [Pango lineage nomenclature](https://cov-lineages.org/)

## 1. Molecular Clocks and Evolutionary Rates

A **molecular clock** assumes that genetic changes accumulate at a roughly constant rate over time. This allows us to date evolutionary events from sequence data alone.

### Clock rate for common viruses

| Virus | Rate (subs/site/year) | Notes |
|---|---|---|
| SARS-CoV-2 | ~1 x 10^-3 | ~1-2 SNPs per 2 weeks |
| Influenza A (HA) | 3-5 x 10^-3 | High immune selection |
| HIV-1 | 2-4 x 10^-3 | Highest among RNA viruses |
| Hepatitis C | 1-2 x 10^-3 | Chronic, slow evolution |
| Measles | ~5 x 10^-4 | More conserved |

### Clock models

**Strict clock**: single global substitution rate for all lineages — too restrictive for long time scales.

**Relaxed clock (uncorrelated lognormal, UCLN)**:
- Each branch draws its rate from a lognormal distribution
- Captures rate heterogeneity across lineages
- Standard in BEAST2, supported in TreeTime

### Root-to-tip regression
Before fitting a full Bayesian model, a simple linear regression of genetic distance (root-to-tip) vs. sampling date validates the molecular clock signal. R² > 0.70 is generally considered sufficient.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

np.random.seed(42)

# ----- Simulate: Root-to-tip regression (clock signal) -----
n_sequences = 80
sampling_dates = np.sort(np.random.uniform(2020.0, 2022.0, n_sequences))

# True clock rate: 1e-3 subs/site/year (SARS-CoV-2 like)
clock_rate = 1e-3
genome_len = 29903
root_to_tip = clock_rate * (sampling_dates - 2020.0) + np.random.normal(0, 5e-5, n_sequences)
root_to_tip = np.maximum(root_to_tip, 0)

slope, intercept, r, p, se = stats.linregress(sampling_dates, root_to_tip)
x_fit = np.linspace(2020, 2022, 100)
y_fit = slope * x_fit + intercept

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Molecular Clock Analysis', fontsize=13, fontweight='bold')

axes[0].scatter(sampling_dates, root_to_tip, alpha=0.6, s=30, c='steelblue', label='Sequences')
axes[0].plot(x_fit, y_fit, 'r-', linewidth=2,
             label=f'Rate = {slope*1e3:.3f}e-3 subs/site/yr\nR2 = {r**2:.3f}')
axes[0].set_xlabel('Sampling Date (decimal year)')
axes[0].set_ylabel('Root-to-tip Distance')
axes[0].set_title('Root-to-Tip Regression (Clock Signal)')
axes[0].legend(fontsize=9)

np.random.seed(1)
strict_rates = np.random.normal(slope, slope*0.05, 1000)
relaxed_rates = np.random.lognormal(np.log(slope), 0.3, 1000)

axes[1].hist(strict_rates * 1e3, bins=40, alpha=0.6, color='steelblue', label='Strict clock', density=True)
axes[1].hist(relaxed_rates * 1e3, bins=40, alpha=0.6, color='orange', label='Relaxed clock (UCLN)', density=True)
axes[1].axvline(slope*1e3, color='red', linestyle='--', label=f'Estimated rate: {slope*1e3:.3f}')
axes[1].set_xlabel('Substitution Rate (x1e-3 subs/site/year)')
axes[1].set_ylabel('Density')
axes[1].set_title('Strict vs Relaxed Clock Rate Distributions')
axes[1].legend(fontsize=9)

plt.tight_layout()
plt.show()

print(f"Clock rate estimate: {slope:.2e} subs/site/year")
print(f"R2 = {r**2:.3f} (>0.70 indicates good clock signal)")
print(f"Expected SNPs per 2 weeks: {slope * (14/365) * genome_len:.1f}")
```

## 2. Time-Scaled Phylogenies

A **time-scaled phylogeny** (timetree) places sequences in calendar time, enabling:
- Dating of common ancestors (TMRCA: time to most recent common ancestor)
- Estimating epidemic origins
- Identifying transmission chains

### Key tools

| Tool | Approach | Speed | Use case |
|---|---|---|---|
| **TreeTime** | ML + least squares | Fast | Real-time surveillance (Nextstrain) |
| **BEAST2** | Bayesian MCMC | Slow | Rigorous inference, posterior distributions |
| **IQ-TREE + LSD2** | Least squares dating | Fast | Large datasets |
| **Augur** (Nextstrain) | TreeTime wrapper | Fast | Integrated pipeline |

### Nextstrain Augur workflow
```bash
# 1. Filter and subsample sequences
augur filter --sequences seqs.fasta --metadata metadata.tsv \
    --min-date 2020-01-01 --subsample-max-sequences 300 --output filtered.fasta

# 2. Align to reference
augur align --sequences filtered.fasta --reference-sequence reference.gb --output aligned.fasta

# 3. Build ML tree and time-scale
augur tree --alignment aligned.fasta --output tree_raw.nwk
augur refine --tree tree_raw.nwk --alignment aligned.fasta \
    --metadata metadata.tsv --timetree --coalescent opt \
    --output-tree timetree.nwk --output-node-data branch_lengths.json

# 4. Export for Auspice visualization
augur export v2 --tree timetree.nwk --output auspice.json
```

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

np.random.seed(7)
n = 25
sampling_dates = np.sort(np.random.uniform(2020.0, 2021.5, n))
mrca_date = 2019.8
clock_rate = 1e-3
genome_len = 29903
n_snps_from_root = np.random.poisson(clock_rate * (sampling_dates - mrca_date) * genome_len)

countries = np.random.choice(['England', 'USA', 'Germany', 'India', 'South_Africa'],
                               n, p=[0.3, 0.3, 0.15, 0.15, 0.1])
clades = np.where(sampling_dates < 2020.6, 'Clade_A',
         np.where(sampling_dates < 2021.0, 'Clade_B', 'Clade_C'))

fig, axes = plt.subplots(1, 2, figsize=(14, 8))

country_colors = {'England': '#E74C3C', 'USA': '#3498DB', 'Germany': '#2ECC71',
                  'India': '#F39C12', 'South_Africa': '#9B59B6'}

for i, (date, snps, country) in enumerate(zip(sampling_dates, n_snps_from_root, countries)):
    divergence = snps / genome_len
    axes[0].plot([0, divergence], [i, i], color='gray', alpha=0.3, linewidth=0.8)
    axes[0].scatter(divergence, i, c=country_colors[country], s=40, zorder=5)

axes[0].set_xlabel('Divergence from root (substitutions/site)')
axes[0].set_ylabel('Sequences')
axes[0].set_title('Divergence Tree')
axes[0].set_yticks([])
handles = [mpatches.Patch(color=v, label=k) for k, v in country_colors.items()]
axes[0].legend(handles=handles, fontsize=8, title='Country')

for i, (date, country, clade) in enumerate(zip(sampling_dates, countries, clades)):
    internal_date = date - np.random.uniform(0.05, 0.4)
    axes[1].plot([internal_date, date], [i, i], color=country_colors[country], linewidth=2)
    axes[1].scatter(date, i, c=country_colors[country], s=40, zorder=5)

axes[1].axvline(mrca_date, color='black', linestyle='--', linewidth=1.5, label=f'MRCA: {mrca_date}')
axes[1].set_xlabel('Calendar Date (year)')
axes[1].set_ylabel('Sequences')
axes[1].set_title('Time-Scaled Phylogeny (TreeTime-style)')
axes[1].set_yticks([])
axes[1].legend(fontsize=9)

plt.tight_layout()
plt.show()
print("TMRCA estimate:", mrca_date)
print("Clade distribution:")
print(pd.Series(clades).value_counts().to_string())
```

## 3. Bayesian Skyline Plot — Effective Population Size

The **Bayesian Skyline Plot (BSP)** estimates how the effective viral population size (N_e) changed through time.

### Interpretation
- **Rising N_e**: exponential epidemic growth
- **Plateau**: endemic equilibrium
- **Sharp decline**: intervention effect, seasonal end, or population immunity
- **N_e != actual viral particles** — it is a genetic diversity metric reflecting the coalescent rate

### BEAST2 Bayesian skyline
In BEAST2, the skyline model uses the timing of coalescent events (when lineages merge in the tree) to estimate N_e(t). More coalescent events per unit time → smaller N_e → faster epidemic growth.

The output is a posterior distribution of N_e at each time point, summarized as median + 95% HPD interval.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

np.random.seed(5)

days = np.linspace(0, 730, 500)
dates = 2020.0 + days / 365

def ne_trajectory(t):
    wave1 = 500 * np.exp(0.015 * t) * np.exp(-0.008 * np.maximum(t - 90, 0))
    wave2 = 800 * np.exp(-0.5 * ((t - 280) / 60)**2)
    wave3 = 2000 * np.exp(-0.5 * ((t - 500) / 80)**2)
    return wave1 + wave2 + wave3 + 50

ne_mean = gaussian_filter1d(ne_trajectory(days), sigma=10)
ne_lower = gaussian_filter1d(ne_mean * np.exp(np.random.normal(-0.3, 0.05, len(days))), sigma=10)
ne_upper = gaussian_filter1d(ne_mean * np.exp(np.random.normal(0.3, 0.05, len(days))), sigma=10)
reported_cases = gaussian_filter1d(ne_mean * 50 * np.exp(np.random.normal(0, 0.2, len(days))), sigma=8)

fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
fig.suptitle('Bayesian Skyline Plot -- Viral Effective Population Size', fontsize=13, fontweight='bold')

axes[0].fill_between(dates, ne_lower, ne_upper, alpha=0.3, color='steelblue', label='95% HPD')
axes[0].plot(dates, ne_mean, 'b-', linewidth=2, label='Median N_e')
axes[0].set_ylabel('Effective Population Size (N_e)')
axes[0].set_yscale('log')
axes[0].set_title('Estimated Viral Effective Population Size Over Time')
axes[0].legend()

for t, label, color in [(90/365+2020, 'Wave 1 Peak', 'red'),
                         (280/365+2020, 'Wave 2 Peak', 'orange'),
                         (500/365+2020, 'Wave 3 Variant', 'purple')]:
    axes[0].axvline(t, color=color, linestyle='--', alpha=0.6, label=label)

axes[0].legend(fontsize=8)

axes[1].fill_between(dates, reported_cases, alpha=0.5, color='salmon')
axes[1].plot(dates, reported_cases, 'r-', linewidth=1)
axes[1].set_xlabel('Date (year)')
axes[1].set_ylabel('Reported Cases (scaled)')
axes[1].set_title('Epidemiological Case Data (for comparison)')

plt.tight_layout()
plt.show()
print("Rising Ne = epidemic growth; declining Ne = transmission reduction")
```

## 4. Phylogeography — Geographic Spread

**Phylogeography** integrates phylogenetics with geographic information to reconstruct how a virus spread across regions.

### Approaches

**Discrete phylogeography (BEAST2 DTA)**
- Each node gets a discrete location state
- Infers most likely ancestral location via symmetric or asymmetric substitution model
- Counts directional transitions between locations

**Continuous phylogeography**
- Lat/lon coordinates as continuous traits (Brownian motion)
- Kernel density estimation of spread velocity

**Nextstrain geographic reconstruction**
- Parsimony or ML location assignment to internal nodes
- Visualized as animated map in Auspice

### Key metrics
- **Migration rate matrix**: transitions per year between locations
- **Bayes factor**: evidence for significant migration routes
- **Geographic diffusion coefficient**: km²/year (continuous model)

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.random.seed(9)

cities = {
    'Wuhan': (30.6, 114.3), 'Singapore': (1.3, 103.8),
    'Milan': (45.5, 9.2), 'New_York': (40.7, -74.0),
    'London': (51.5, -0.1), 'Mumbai': (19.1, 72.9),
    'Nairobi': (-1.3, 36.8), 'Sydney': (-33.9, 151.2),
    'Toronto': (43.7, -79.4), 'Sao_Paulo': (-23.5, -46.6),
}

city_names = list(cities.keys())
coords = np.array(list(cities.values()))

emergence_week = {
    'Wuhan': 0, 'Singapore': 2, 'Milan': 5, 'New_York': 6,
    'London': 6, 'Mumbai': 8, 'Nairobi': 9, 'Sydney': 4,
    'Toronto': 7, 'Sao_Paulo': 9
}

migrations = [
    ('Wuhan', 'Singapore', 0.8), ('Wuhan', 'Milan', 0.6),
    ('Wuhan', 'Sydney', 0.5), ('Milan', 'New_York', 0.4),
    ('Milan', 'London', 0.7), ('London', 'Toronto', 0.5),
    ('Mumbai', 'Nairobi', 0.3), ('New_York', 'Sao_Paulo', 0.4),
]

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

emergence_vals = np.array([emergence_week[c] for c in city_names])
sc = axes[0].scatter(coords[:, 1], coords[:, 0], c=emergence_vals, s=120,
                cmap='viridis_r', zorder=5, edgecolors='black', linewidths=0.5)
plt.colorbar(sc, ax=axes[0], label='Emergence week')
for name, (lat, lon) in cities.items():
    axes[0].text(lon + 2, lat + 2, name.replace('_', ' '), fontsize=6.5)

for src, dst, intensity in migrations:
    sc2 = cities[src]; dc = cities[dst]
    axes[0].annotate('', xy=(dc[1], dc[0]), xytext=(sc2[1], sc2[0]),
                arrowprops=dict(arrowstyle='->', color='red', alpha=intensity, lw=1.5*intensity))

axes[0].set_xlim(-150, 170); axes[0].set_ylim(-60, 75)
axes[0].set_facecolor('#d4e6f1')
axes[0].set_title('Phylogeographic Spread', fontsize=11, fontweight='bold')
axes[0].set_xlabel('Longitude'); axes[0].set_ylabel('Latitude')

for i, city in enumerate(city_names):
    ew = emergence_week[city]
    t_local = np.arange(0, 10)
    cases = np.exp(0.5 * t_local) * 100 * np.random.uniform(0.5, 2.0)
    axes[1].plot(t_local + ew, cases / cases.max(), alpha=0.6, linewidth=1.5,
                 label=city.replace('_', ' '))

axes[1].set_xlabel('Epidemic Week (global)')
axes[1].set_ylabel('Normalized Local Incidence')
axes[1].set_title('Epidemic Curves by Region')
axes[1].legend(fontsize=7, ncol=2)
axes[1].set_xlim(0, 20)

plt.tight_layout()
plt.show()
```
