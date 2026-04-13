---
name: bio-applied-variant-surveillance
description: "**Tier 3 — Applied Bioinformatics | Module 37 · Notebook 3**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/37_Virology_Bioinformatics/03_variant_surveillance.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Variant Surveillance and Wastewater Epidemiology

*Source: Course notebook `Tier_3_Applied_Bioinformatics/37_Virology_Bioinformatics/03_variant_surveillance.ipynb`*

# Variant Surveillance and Wastewater Epidemiology

**Tier 3 — Applied Bioinformatics | Module 37 · Notebook 3**

*Prerequisites: Notebook 2 (Phylodynamics)*

---

**By the end of this notebook you will be able to:**
1. Build automated surveillance pipelines using Nextclade and Freyja
2. Deconvolve variant proportions from wastewater sequencing with Freyja
3. Track variant emergence and displacement dynamics over time
4. Predict functional effects of variant mutations (spike protein, protease)
5. Submit sequences to GISAID / GenBank with automated metadata



**Key resources:**
- [Freyja documentation](https://github.com/andersen-lab/Freyja)
- [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/)
- [GISAID EpiCoV](https://www.gisaid.org/)
- [outbreak.info](https://outbreak.info/)

## 1. Lineage Classification — Pango Nomenclature

The **Pango** (Phylogenetic Assignment of Named Global Outbreak) system provides a hierarchical naming scheme for SARS-CoV-2 lineages.

### How Pango works
1. **Pangolin tool** aligns a consensus sequence to reference + known lineage sequences
2. Assigns to the most probable Pango lineage using probabilistic model or scorpio (constellation-based)
3. Outputs lineage + confidence score + WHO designation if applicable

### Lineage naming rules
- Root: A, B (early Wuhan sequences)
- Each sub-lineage adds a letter or number (B.1 → B.1.1 → B.1.1.7 = Alpha)
- After 3 sub-levels: recombinants get X prefix (XBB, XBF)
- After 3+sub-levels: new alias issued (BA, BQ, XBB)

### WHO variant categories
| Category | Criteria | Examples |
|---|---|---|
| **VOC** (Concern) | Increased transmissibility + immune escape + severity | Alpha, Delta, Omicron |
| **VOI** (Interest) | Some evidence of phenotypic change | Lambda, Mu |
| **VUM** (Monitoring) | Under investigation | BA.2.86, JN.1 |

```bash
# Assign lineages with Pangolin
pangolin sequences.fasta --outfile lineages.csv --threads 4

# Assign clades with Nextclade
nextclade run --input-fasta sequences.fasta \
    --input-dataset sars-cov-2 --output-tsv nextclade.tsv
```

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

np.random.seed(42)

# Simulate lineage frequency tracking over time (SARS-CoV-2-like)
weeks = np.arange(0, 52)  # one year of weekly data
dates_str = [f"2022-W{w:02d}" for w in range(1, 53)]

# Lineage names and their emergence/dominance profiles
lineages = {
    'BA.1 (Omicron)':   (0, 15, 3),    # (peak_week, peak_freq, spread)
    'BA.2':              (8, 35, 4),
    'BA.2.12.1':         (18, 20, 3),
    'BA.4':              (20, 15, 4),
    'BA.5':              (24, 40, 5),
    'BQ.1':              (35, 25, 4),
    'XBB.1.5':           (40, 30, 5),
    'Other':             None
}

def lineage_freq(t, peak, max_freq, spread):
    return max_freq * np.exp(-0.5 * ((t - peak) / spread)**2)

freqs = {}
for lin, params in lineages.items():
    if params is not None:
        freqs[lin] = lineage_freq(weeks, *params)
    
# Normalize to sum to 100% per week
freq_matrix = np.array(list(freqs.values()))
row_sums = freq_matrix.sum(axis=0)
row_sums = np.maximum(row_sums, 0.1)
freq_matrix_norm = freq_matrix / row_sums * 100
other_freq = np.maximum(100 - freq_matrix_norm.sum(axis=0), 0)
freq_matrix_norm = np.vstack([freq_matrix_norm, other_freq])
lin_labels = list(freqs.keys()) + ['Other']

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('SARS-CoV-2 Lineage Surveillance (simulated)', fontsize=13, fontweight='bold')

# Stacked area plot
colors = cm.tab10(np.linspace(0, 1, len(lin_labels)))
axes[0].stackplot(weeks, freq_matrix_norm, labels=lin_labels, colors=colors, alpha=0.85)
axes[0].set_xlabel('Epidemic week')
axes[0].set_ylabel('Frequency (%)')
axes[0].set_title('Variant Frequency Over Time')
axes[0].legend(loc='upper left', fontsize=7, ncol=1, bbox_to_anchor=(1.01, 1))
axes[0].set_ylim(0, 100)

# Dominant lineage per week
dominant = np.argmax(freq_matrix_norm, axis=0)
dom_colors = [colors[d] for d in dominant]
axes[1].bar(weeks, [100]*len(weeks), color=dom_colors, alpha=0.8)
axes[1].set_xlabel('Epidemic week')
axes[1].set_ylabel('Dominant lineage')
axes[1].set_title('Weekly Dominant Lineage')
axes[1].set_yticks([])
# Add lineage legend
from matplotlib.patches import Patch
handles = [Patch(color=colors[i], label=lin_labels[i]) for i in range(len(lin_labels))]
axes[1].legend(handles=handles, fontsize=7, ncol=1, bbox_to_anchor=(1.01, 1), loc='upper left')

plt.tight_layout()
plt.show()

# Summary table
peak_weeks = {lin_labels[i]: weeks[np.argmax(freq_matrix_norm[i])] for i in range(len(lin_labels)-1)}
print("Estimated peak weeks:")
for lin, pw in peak_weeks.items():
    print(f"  {lin}: week {pw}")
```

## 2. Mutation Tracking and Spike Protein Analysis

Tracking mutations at key functional sites enables:
- Early detection of immune escape
- Vaccine effectiveness monitoring
- Treatment resistance surveillance

### Key functional regions of Spike protein

| Region | Positions | Function |
|---|---|---|
| N-terminal domain (NTD) | 13-305 | Antibody binding site |
| Receptor-binding domain (RBD) | 319-541 | ACE2 contact; most mutations here |
| Furin cleavage site | 681-685 | PRRAR motif; P681H/R increases fitness |
| Fusion peptide | 788-806 | Membrane fusion mechanism |

### Mutation effect classification
- **Immune escape**: reduces neutralization by antibodies (E484K, K417N, F486V)
- **ACE2 affinity**: increases binding to host receptor (N501Y, Q493R)
- **Transmissibility**: increases viral fitness without known mechanism (D614G, P681H)

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.random.seed(7)

# Spike protein mutations of concern in SARS-CoV-2
spike_mutations = {
    'D614G':     {'position': 614, 'effect': 'transmission', 'freq_increase': 0.95},
    'N501Y':     {'position': 501, 'effect': 'ACE2_affinity', 'freq_increase': 0.85},
    'E484K':     {'position': 484, 'effect': 'immune_escape', 'freq_increase': 0.60},
    'K417N':     {'position': 417, 'effect': 'immune_escape', 'freq_increase': 0.70},
    'L452R':     {'position': 452, 'effect': 'both', 'freq_increase': 0.75},
    'P681H':     {'position': 681, 'effect': 'furin_cleavage', 'freq_increase': 0.88},
    'F486V':     {'position': 486, 'effect': 'immune_escape', 'freq_increase': 0.45},
    'R346T':     {'position': 346, 'effect': 'immune_escape', 'freq_increase': 0.40},
    'Q493R':     {'position': 493, 'effect': 'ACE2_affinity', 'freq_increase': 0.55},
    'G339D':     {'position': 339, 'effect': 'structural', 'freq_increase': 0.30},
}

df_mut = pd.DataFrame(spike_mutations).T.reset_index()
df_mut.columns = ['Mutation', 'Position', 'Effect', 'Freq']
df_mut['Position'] = df_mut['Position'].astype(int)
df_mut['Freq'] = df_mut['Freq'].astype(float)

# RBD domain: positions 319-541
rbd_mask = (df_mut['Position'] >= 319) & (df_mut['Position'] <= 541)

# Simulate time series for each mutation
weeks = np.arange(0, 52)
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle('Spike Protein Mutation Tracking', fontsize=13, fontweight='bold')

# 1. Mutation positions on spike (linear map)
effect_colors = {'transmission': 'steelblue', 'ACE2_affinity': 'green',
                 'immune_escape': 'red', 'furin_cleavage': 'orange',
                 'both': 'purple', 'structural': 'gray'}
spike_len = 1273  # AA
axes[0,0].barh([0]*spike_len, [1]*spike_len, left=range(spike_len), height=0.3,
               color='lightgray', alpha=0.5)
# RBD domain
axes[0,0].barh(0, 541-319, left=319, height=0.3, color='lightyellow',
               alpha=0.8, label='RBD (319-541)')
for _, row in df_mut.iterrows():
    axes[0,0].scatter(row['Position'], 0, s=150, c=effect_colors[row['Effect']],
                      zorder=5, edgecolors='black', linewidths=0.5)
    axes[0,0].text(row['Position'], 0.25, row['Mutation'], ha='center',
                   va='bottom', fontsize=6.5, rotation=45)
axes[0,0].set_xlim(0, 1273)
axes[0,0].set_yticks([])
axes[0,0].set_xlabel('Spike protein position (aa)')
axes[0,0].set_title('Mutations on Spike Protein')
from matplotlib.patches import Patch
handles = [Patch(color=v, label=k.replace('_',' ')) for k,v in effect_colors.items()]
axes[0,0].legend(handles=handles, fontsize=7, loc='upper right')

# 2. Mutation frequency barplot
bar_colors = [effect_colors[e] for e in df_mut['Effect']]
axes[0,1].barh(df_mut['Mutation'], df_mut['Freq']*100, color=bar_colors, alpha=0.8)
axes[0,1].axvline(50, color='black', linestyle='--', label='50% prevalence')
axes[0,1].set_xlabel('Current Prevalence (%)')
axes[0,1].set_title('Mutation Prevalence in Circulating Sequences')
axes[0,1].legend()

# 3. Temporal trajectory for top mutations
for _, row in df_mut.iterrows():
    emergence_week = np.random.randint(0, 20)
    growth_rate = np.random.uniform(0.1, 0.4)
    freq_traj = row['Freq'] / (1 + np.exp(-growth_rate * (weeks - emergence_week - 15)))
    freq_traj += np.random.normal(0, 0.02, len(weeks))
    freq_traj = np.clip(freq_traj, 0, 1)
    axes[1,0].plot(weeks, freq_traj * 100, label=row['Mutation'],
                   color=effect_colors[row['Effect']], alpha=0.7)
axes[1,0].axhline(50, color='black', linestyle='--', alpha=0.5)
axes[1,0].set_xlabel('Epidemic week')
axes[1,0].set_ylabel('Frequency (%)')
axes[1,0].set_title('Mutation Frequency Trajectories')
axes[1,0].legend(fontsize=7, ncol=2)

# 4. Co-occurrence matrix
n = len(df_mut)
co_matrix = np.random.uniform(0.3, 0.9, (n, n))
np.fill_diagonal(co_matrix, 1.0)
co_matrix = (co_matrix + co_matrix.T) / 2
im = axes[1,1].imshow(co_matrix, cmap='RdYlGn', vmin=0, vmax=1)
axes[1,1].set_xticks(range(n))
axes[1,1].set_yticks(range(n))
axes[1,1].set_xticklabels(df_mut['Mutation'], rotation=45, ha='right', fontsize=7)
axes[1,1].set_yticklabels(df_mut['Mutation'], fontsize=7)
axes[1,1].set_title('Mutation Co-occurrence Matrix')
plt.colorbar(im, ax=axes[1,1], label='Co-occurrence freq')

plt.tight_layout()
plt.show()
print(f"RBD mutations: {rbd_mask.sum()}/{len(df_mut)} mutations in receptor binding domain")
```

## 3. Wastewater Genomic Surveillance

**Wastewater-based epidemiology (WBE)** detects SARS-CoV-2 (and other pathogens) in sewage:
- Captures infections regardless of clinical testing rates
- Provides 3-7 days **lead time** before clinical case surge
- Population-level signal: ~1 infected person detectable per 100,000 catchment

### Freyja lineage deconvolution
Freyja uses a **mixture model** to estimate variant proportions from metagenomics:
1. Call variants from wastewater BAM with `freyja variants`
2. Deconvolve using curated **UShER barcodes** (lineage-defining mutations)
3. Output: fraction of each lineage in the sample

```bash
# Step 1: Variant calling from wastewater BAM
freyja variants wastewater.bam --variants variants.tsv --depths depths.tsv

# Step 2: Lineage deconvolution
freyja demix variants.tsv depths.tsv --output lineages.csv

# Step 3: Aggregate across time
freyja aggregate --inputdir ./samples/ --output aggregated.tsv
freyja plot aggregated.tsv --output lineage_plot.pdf
```

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

np.random.seed(13)

# Simulate wastewater surveillance data
n_weeks = 60
weeks = np.arange(n_weeks)
dates = pd.date_range('2022-01-01', periods=n_weeks, freq='W')

# True variant frequencies in community (from clinical)
ba2_true = np.clip(0.9 * np.exp(-0.5 * ((weeks - 10) / 8)**2), 0, 1)
ba5_true = np.clip(0.8 * np.exp(-0.5 * ((weeks - 30) / 10)**2), 0, 1)
xbb_true = np.clip(0.6 * np.exp(-0.5 * ((weeks - 50) / 8)**2), 0, 1)
other_true = np.maximum(0, 1 - ba2_true - ba5_true - xbb_true)

# Wastewater estimates (Freyja deconvolution with uncertainty)
noise = 0.05
ba2_ww = np.clip(ba2_true + np.random.normal(0, noise, n_weeks), 0, 1)
ba5_ww = np.clip(ba5_true + np.random.normal(0, noise, n_weeks), 0, 1)
xbb_ww = np.clip(xbb_true + np.random.normal(0, noise, n_weeks), 0, 1)

# Viral load signal (proxy for incidence)
viral_load = (gaussian_filter1d(ba2_true + ba5_true + xbb_true, sigma=2) *
              1e6 * np.exp(np.random.normal(0, 0.3, n_weeks)))
viral_load_smooth = gaussian_filter1d(viral_load, sigma=3)

# Clinical cases (lag 5 days behind wastewater)
cases = np.roll(viral_load / 1000, 5) * np.exp(np.random.normal(0, 0.2, n_weeks))
cases[cases < 0] = 0

fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
fig.suptitle('Wastewater-Based Epidemiology (Freyja-style)', fontsize=13, fontweight='bold')

# Panel 1: Wastewater viral load
axes[0].fill_between(range(n_weeks), viral_load, alpha=0.3, color='steelblue')
axes[0].plot(range(n_weeks), viral_load_smooth, 'b-', linewidth=2, label='SARS-CoV-2 RNA (WW)')
axes[0].set_ylabel('Viral load (copies/L)')
axes[0].set_yscale('log')
axes[0].set_title('Wastewater SARS-CoV-2 Viral Load')
axes[0].legend()

# Panel 2: Variant deconvolution from wastewater (Freyja)
axes[1].stackplot(range(n_weeks),
                  [ba2_ww * 100, ba5_ww * 100, xbb_ww * 100,
                   np.maximum(0, 100 - ba2_ww*100 - ba5_ww*100 - xbb_ww*100)],
                  labels=['BA.2', 'BA.5', 'XBB.1.5', 'Other'],
                  colors=['steelblue', 'orange', 'green', 'lightgray'], alpha=0.85)
axes[1].set_ylabel('Lineage frequency (%)')
axes[1].set_title('Lineage Deconvolution (Freyja) from Wastewater Samples')
axes[1].legend(loc='upper left', fontsize=9)
axes[1].set_ylim(0, 100)

# Panel 3: Wastewater vs clinical cases (lead indicator)
ax2 = axes[2].twinx()
axes[2].bar(range(n_weeks), cases, color='salmon', alpha=0.6, label='Clinical cases (reported)')
ax2.plot(range(n_weeks), viral_load_smooth, 'b-', linewidth=2, label='Wastewater signal (lead)')
axes[2].set_xlabel('Epidemic week')
axes[2].set_ylabel('Clinical cases')
ax2.set_ylabel('WW viral load (copies/L)')
ax2.set_yscale('log')
axes[2].set_title('Wastewater as Early Warning Signal (leads cases by ~5 days)')
lines1, labels1 = axes[2].get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
axes[2].legend(lines1 + lines2, labels1 + labels2, fontsize=9)

plt.tight_layout()
plt.show()

# Correlation analysis
from scipy.stats import pearsonr
# Shift clinical cases back 5 weeks to align with WW
ww_signal = viral_load_smooth[:-5]
cli_shifted = cases[5:]
r, p = pearsonr(ww_signal, cli_shifted)
print(f"Pearson r (WW vs clinical, 5-week lead): {r:.3f} (p={p:.2e})")
print("Freyja command: freyja variants <bam> --variants <vcf> --depths <depths>")
print("                freyja demix <variants> <depths> --output <lineages.csv>")
```
