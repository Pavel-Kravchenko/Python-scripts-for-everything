---
name: python-bio-data-visualization
description: "Data visualization for bioinformatics: matplotlib, seaborn, volcano plots, heatmaps, and genome browser tracks. Use when creating publication-quality scientific figures."
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/18_Data_Visualization/01_data_visualization.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Data Visualization for Bioinformatics

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/18_Data_Visualization/01_data_visualization.ipynb`*


**A comprehensive guide to creating publication-quality biological figures**

---

## Learning Objectives

By the end of this notebook, you will be able to:

1. **Matplotlib fundamentals** -- figures, axes, and all major plot types
2. **Customization** -- colors, labels, legends, annotations, styles
3. **Seaborn** -- statistical plots, heatmaps, pair plots, violin/box plots
4. **Plotly** -- interactive plots for exploratory analysis
5. **Bioinformatics visualizations** -- volcano plots, GC landscapes, expression heatmaps, PCA, and more
6. **Publication-quality figures** -- fonts, DPI, colorblind-friendly palettes, multi-panel layouts

---

| Section | Topic |
|---------|-------|
| 1 | Matplotlib Fundamentals |
| 2 | Plot Customization |
| 3 | Seaborn Statistical Plots |
| 4 | Bio: GC Content Sliding Window |
| 5 | Bio: Volcano Plot |
| 6 | Bio: Gene Expression Heatmap |
| 7 | Bio: Sequence Length Distribution |
| 8 | Bio: Quality Score Distribution |
| 9 | Bio: Genome Coverage Plot |
| 10 | Bio: PCA of Samples |
| 11 | Bio: Phylogenetic Dendrogram |
| 12 | Plotly Interactive Plots |
| 13 | Publication-Quality Figures |
| 14 | Exercises |

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. Matplotlib's two interfaces**
`plt.plot(...)` (pyplot interface) acts on the current active axes. `fig, ax = plt.subplots(); ax.plot(...)` (OO interface) is explicit and required for multi-panel figures. Use the OO interface for anything beyond a single plot.

**2. `savefig` must come before `plt.show()`**
`plt.show()` clears the figure. If you call it first, `savefig` writes a blank file. Always save before showing.

**3. Seaborn expects long-form data**
Most seaborn functions expect one observation per row. If your data is wide (one column per sample), melt it with `pd.melt()` first.

**4. Color choices**
Default matplotlib colors are not colorblind-friendly. Use `seaborn`'s `colorblind` palette for categorical data. For expression heatmaps, use `RdBu_r` or `viridis` (perceptually uniform), not `jet`.

**5. Log axes and zeros**
`plt.yscale('log')` fails on data containing zeros. For RNA-seq data, log-transform values before plotting (`log2(count + 1)`) rather than using a log-scale axis on raw counts.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

np.random.seed(42)
plt.rcParams['figure.dpi'] = 100
```python

### 1.1 The Figure and Axes Model

```python
Figure (the canvas)
 +-- Axes (a single plot area)
      +-- Title, xlabel, ylabel
      +-- x-axis, y-axis (ticks, limits)
      +-- plotted data (lines, bars, points)
```python

```python
# Create a figure with one axes
fig, ax = plt.subplots(figsize=(8, 4))

# Simulated growth curve of bacterial culture (OD600 over time)
time_hours = np.linspace(0, 24, 200)
od600 = 0.05 * np.exp(0.3 * time_hours) / (1 + 0.05 * (np.exp(0.3 * time_hours) - 1) / 2.0)

ax.plot(time_hours, od600, color='darkgreen', linewidth=2)
ax.set_xlabel('Time (hours)', fontsize=12)
ax.set_ylabel('OD$_{600}$', fontsize=12)
ax.set_title('Bacterial Growth Curve (Logistic Model)', fontsize=14)
ax.set_xlim(0, 24)
ax.set_ylim(0, 2.2)

# Annotate growth phases
ax.axvspan(0, 3, alpha=0.15, color='gray', label='Lag')
ax.axvspan(3, 14, alpha=0.15, color='green', label='Exponential')
ax.axvspan(14, 24, alpha=0.15, color='orange', label='Stationary')
ax.legend(loc='upper left', fontsize=10)

plt.tight_layout()
plt.show()
```python

### 1.2 Core Plot Types

```python
fig, axes = plt.subplots(2, 3, figsize=(15, 8))

# --- Line plot: enzyme kinetics (Michaelis-Menten) ---
substrate = np.linspace(0, 100, 200)
vmax, km = 120, 15
velocity = vmax * substrate / (km + substrate)
axes[0, 0].plot(substrate, velocity, 'b-', linewidth=2)
axes[0, 0].axhline(vmax, color='red', linestyle='--', alpha=0.6, label=f'V_max = {vmax}')
axes[0, 0].axvline(km, color='gray', linestyle=':', alpha=0.6, label=f'K_m = {km}')
axes[0, 0].set_xlabel('[S] (mM)')
axes[0, 0].set_ylabel('Velocity (nmol/min)')
axes[0, 0].set_title('Line: Michaelis-Menten')
axes[0, 0].legend(fontsize=8)

# --- Scatter plot: gene expression correlation ---
gene_a = np.random.lognormal(3, 1.5, 200)
gene_b = gene_a * 0.7 + np.random.lognormal(2, 1, 200)
axes[0, 1].scatter(np.log2(gene_a + 1), np.log2(gene_b + 1),
                    alpha=0.5, s=15, c='steelblue', edgecolors='none')
axes[0, 1].set_xlabel('Gene A (log2 FPKM)')
axes[0, 1].set_ylabel('Gene B (log2 FPKM)')
axes[0, 1].set_title('Scatter: Expression Correlation')

# --- Bar plot: nucleotide composition ---
nucleotides = ['A', 'T', 'G', 'C']
human_freq = [29.3, 29.3, 20.7, 20.7]
ecoli_freq = [24.7, 23.6, 25.7, 26.0]
x_pos = np.arange(len(nucleotides))
width = 0.35
axes[0, 2].bar(x_pos - width/2, human_freq, width, label='Human', color='#2196F3')
axes[0, 2].bar(x_pos + width/2, ecoli_freq, width, label='E. coli', color='#FF9800')
axes[0, 2].set_xticks(x_pos)
axes[0, 2].set_xticklabels(nucleotides)
axes[0, 2].set_ylabel('Frequency (%)')
axes[0, 2].set_title('Bar: Nucleotide Composition')
axes[0, 2].legend(fontsize=8)

# --- Histogram: protein lengths ---
protein_lengths = np.concatenate([
    np.random.lognormal(5.5, 0.8, 800),
    np.random.lognormal(7, 0.5, 200)
])
axes[1, 0].hist(protein_lengths, bins=60, color='mediumpurple',
                edgecolor='white', linewidth=0.5, range=(0, 3000))
axes[1, 0].set_xlabel('Protein Length (aa)')
axes[1, 0].set_ylabel('Count')
axes[1, 0].set_title('Histogram: Protein Lengths')

# --- Pie chart: genome composition ---
labels = ['Exons', 'Introns', 'Intergenic', 'Repeats']
sizes = [1.5, 25, 28.5, 45]
colors = ['#4CAF50', '#2196F3', '#FFC107', '#F44336']
axes[1, 1].pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
               startangle=90, textprops={'fontsize': 9})
axes[1, 1].set_title('Pie: Human Genome Composition')

# --- Stem plot: SNP positions ---
snp_positions = np.sort(np.random.choice(range(1, 1001), 20, replace=False))
snp_effects = np.random.choice([-1, 0, 1], 20, p=[0.2, 0.6, 0.2])
colors_snp = ['red' if e == -1 else ('green' if e == 1 else 'gray') for e in snp_effects]
markerline, stemlines, baseline = axes[1, 2].stem(snp_positions, snp_effects)
plt.setp(stemlines, linewidth=1)
plt.setp(markerline, markersize=4)
axes[1, 2].set_xlabel('Genomic Position')
axes[1, 2].set_ylabel('Effect')
axes[1, 2].set_title('Stem: SNP Effects')
axes[1, 2].set_yticks([-1, 0, 1])
axes[1, 2].set_yticklabels(['Deleterious', 'Neutral', 'Beneficial'])

plt.tight_layout()
plt.show()
```python

### 1.3 Subplots and Layout

Multiple approaches for arranging subplots:
- `plt.subplots(nrows, ncols)` -- regular grid
- `fig.add_gridspec()` -- flexible grid with spanning

```python
# GridSpec for flexible layout: one wide plot on top, two below
fig = plt.figure(figsize=(12, 7))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.35, wspace=0.3)

ax_top = fig.add_subplot(gs[0, :])
ax_bl = fig.add_subplot(gs[1, 0])
ax_br = fig.add_subplot(gs[1, 1])

# Top: simulated expression timecourse for 3 genes
t = np.linspace(0, 48, 100)
for i, (name, phase) in enumerate([('p53', 0), ('MDM2', 4), ('CDKN1A', 8)]):
    signal = 2 * np.sin(2 * np.pi * t / 24 - phase) * np.exp(-t / 80) + 5
    ax_top.plot(t, signal, linewidth=2, label=name)
ax_top.set_xlabel('Time (hours)')
ax_top.set_ylabel('Expression (log2)')
ax_top.set_title('Gene Expression Timecourse After DNA Damage')
ax_top.legend()

# Bottom-left: bar chart of pathway enrichment
pathways = ['Apoptosis', 'Cell Cycle', 'DNA Repair', 'p53 Signaling', 'Autophagy']
neg_log_p = [8.5, 7.2, 6.1, 10.3, 3.2]
ax_bl.barh(pathways, neg_log_p, color=plt.cm.viridis(np.linspace(0.2, 0.8, 5)))
ax_bl.set_xlabel('-log$_{10}$(p-value)')
ax_bl.set_title('Pathway Enrichment')
ax_bl.axvline(x=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
ax_bl.legend(fontsize=8)

# Bottom-right: scatter of fold-change vs base mean
base_mean = np.random.lognormal(5, 2, 500)
log2fc = np.random.normal(0, 0.5, 500)
log2fc[:30] = np.random.normal(2, 0.3, 30)
ax_br.scatter(np.log10(base_mean), log2fc, s=10, alpha=0.5, c='gray')
ax_br.scatter(np.log10(base_mean[:30]), log2fc[:30], s=15, alpha=0.8, c='red', label='Significant')
ax_br.set_xlabel('log$_{10}$(Base Mean)')
ax_br.set_ylabel('log$_2$(Fold Change)')
ax_br.set_title('MA Plot')
ax_br.axhline(0, color='black', linewidth=0.8)
ax_br.legend(fontsize=8)

plt.show()
```python

---
## 2. Plot Customization

Making plots informative and visually appealing.

```python
# Demonstrate customization on a single figure
fig, ax = plt.subplots(figsize=(9, 5))

# Simulated RT-qPCR data: expression of 5 genes in control vs treatment
genes = ['GAPDH', 'TP53', 'BRCA1', 'MYC', 'BCL2']
control = [1.0, 1.0, 1.0, 1.0, 1.0]
treatment = [1.05, 3.2, 0.4, 5.1, 0.25]
treatment_err = [0.1, 0.5, 0.08, 0.8, 0.05]

x = np.arange(len(genes))
width = 0.35

bars1 = ax.bar(x - width/2, control, width, label='Control',
               color='#90CAF9', edgecolor='#1565C0', linewidth=1.2)
bars2 = ax.bar(x + width/2, treatment, width, yerr=treatment_err,
               label='Treatment (6h cisplatin)',
               color='#EF9A9A', edgecolor='#C62828', linewidth=1.2,
               capsize=4)

# Customization
ax.set_xticks(x)
ax.set_xticklabels(genes, fontsize=12, fontweight='bold')
ax.set_ylabel('Relative Expression\n(fold change)', fontsize=12)
ax.set_title('RT-qPCR: Gene Expression After Cisplatin Treatment', fontsize=14, pad=15)
ax.legend(fontsize=11, frameon=True, fancybox=True, shadow=True)
ax.axhline(1, color='gray', linestyle='--', linewidth=0.8, zorder=0)

# Add significance stars
for i, (c, t) in enumerate(zip(control, treatment)):
    if abs(t - c) > 0.5:
        y_pos = max(c, t + treatment_err[i]) + 0.3
        ax.text(i + width/2, y_pos, '*', ha='center', fontsize=16, color='red')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(0, 7)

plt.tight_layout()
plt.show()
```python

### 2.1 Colors, Colormaps, and Colorblind-Friendly Palettes

```python
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# 1. Named color palettes for categorical data
palettes = {
    'colorblind': sns.color_palette('colorblind', 6),
    'Set2': sns.color_palette('Set2', 6),
    'deep': sns.color_palette('deep', 6),
}
for row, (name, pal) in enumerate(palettes.items()):
    for i, c in enumerate(pal):
        axes[0].add_patch(plt.Rectangle((i, row), 0.9, 0.9, color=c))
axes[0].set_xlim(-0.2, 6.2)
axes[0].set_ylim(-0.2, 3.2)
axes[0].set_yticks([0.45, 1.45, 2.45])
axes[0].set_yticklabels(list(palettes.keys()))
axes[0].set_title('Categorical Palettes')
axes[0].set_xticks([])

# 2. Sequential colormap (great for heatmaps)
gradient = np.linspace(0, 1, 256).reshape(1, -1)
for row, cmap_name in enumerate(['viridis', 'plasma', 'cividis']):
    axes[1].imshow(gradient, aspect='auto', cmap=cmap_name,
                   extent=[0, 10, row, row+0.8])
axes[1].set_yticks([0.4, 1.4, 2.4])
axes[1].set_yticklabels(['viridis', 'plasma', 'cividis'])
axes[1].set_title('Sequential Colormaps (perceptually uniform)')
axes[1].set_xticks([])
axes[1].set_ylim(-0.2, 3)

# 3. Diverging colormap (centered data like fold-change)
diverging = np.linspace(-1, 1, 256).reshape(1, -1)
for row, cmap_name in enumerate(['RdBu_r', 'coolwarm', 'PiYG']):
    axes[2].imshow(diverging, aspect='auto', cmap=cmap_name,
                   extent=[0, 10, row, row+0.8])
axes[2].set_yticks([0.4, 1.4, 2.4])
axes[2].set_yticklabels(['RdBu_r', 'coolwarm', 'PiYG'])
axes[2].set_title('Diverging Colormaps (centered at 0)')
axes[2].set_xticks([])
axes[2].set_ylim(-0.2, 3)

plt.tight_layout()
plt.show()

print('Tip: viridis, cividis, and the colorblind palette are safe for color vision deficiency.')
```python

---
## 3. Seaborn Statistical Plots

Seaborn is built on Matplotlib and provides a higher-level interface for statistical graphics. It integrates closely with Pandas DataFrames.

## Common Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
