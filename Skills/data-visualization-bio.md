---
name: data-visualization-bio
description: Publication-quality biological plots — volcano plots, MA plots, heatmaps, genome tracks, PCA, quality scores, survival curves with matplotlib and seaborn
---

# Data Visualization for Bioinformatics

## When to Use
- Differential expression results → volcano plot or MA plot
- Multi-sample transcriptomics → clustered heatmap or PCA
- Genomic region → coverage/track plot or GC sliding window
- Sequencing QC → per-base quality (FastQC-style)
- Phylogenetics → dendrogram
- Clinical/time-to-event data → Kaplan-Meier survival curve
- Statistical distributions → seaborn box/violin/KDE plots

## Quick Reference

| Goal | Function |
|------|----------|
| Multiple panels | `fig, axes = plt.subplots(nrows, ncols)` |
| Flexible layout | `gs = fig.add_gridspec(r, c, height_ratios=[...])` |
| Clustered heatmap | `sns.clustermap(df, cmap='RdBu_r', center=0)` |
| Pair plot | `sns.pairplot(df, hue='condition')` |
| Save publication | `fig.savefig('fig.svg')` or `fig.savefig('fig.png', dpi=300)` |
| Temporary settings | `with plt.rc_context(pub_params): ...` |

## Key Patterns

### Imports
```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA

plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette('colorblind')
```

### Color Palettes
- **Categorical (colorblind-safe):** `sns.color_palette('colorblind', n)` — always prefer over default
- **Sequential (heatmaps):** `viridis`, `plasma`, `cividis` — perceptually uniform, print safe
- **Diverging (fold-change, z-score):** `RdBu_r`, `coolwarm`, `PiYG` — center at 0
- **Manual colorblind-safe:** `#2196F3` (blue), `#E53935` (red), `#4CAF50` (green), `#BDBDBD` (gray)

### Axes Cleanup
```python
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
```

### Significance Annotation
```python
for i, (ctrl, treat, err) in enumerate(zip(control, treatment, errs)):
    if abs(treat - ctrl) > threshold:
        ax.text(i + offset, max(ctrl, treat + err) + pad, '*',
                ha='center', fontsize=16, color='red')
```

## Code Templates

### Publication RC Parameters
```python
pub_params = {
    'figure.figsize': (7, 5),        # Nature single column ~89mm
    'figure.dpi': 150,
    'font.size': 10,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'axes.linewidth': 1.0,
    'lines.linewidth': 1.5,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
}

with plt.rc_context(pub_params):
    fig, ax = plt.subplots()
    # ... plot code ...
    fig.savefig('Figure1.svg')   # vector, editable in Illustrator
    fig.savefig('Figure1.pdf')   # vector, for LaTeX
    fig.savefig('Figure1.png', dpi=300)
```

### Multi-Panel Layout with Panel Labels
```python
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)
ax_a = fig.add_subplot(gs[0, 0])
ax_top = fig.add_subplot(gs[0, :])   # span full row

# Add panel letter (A, B, C...)
ax_a.text(-0.15, 1.05, 'A', transform=ax_a.transAxes,
          fontsize=16, fontweight='bold', va='bottom')
```

### Volcano Plot
```python
fc_thresh, p_thresh = 1.0, 0.05
df['neg_log10_p'] = -np.log10(df['pvalue'])
df['category'] = 'Not Significant'
df.loc[(df['log2FC'] > fc_thresh) & (df['pvalue'] < p_thresh), 'category'] = 'Up'
df.loc[(df['log2FC'] < -fc_thresh) & (df['pvalue'] < p_thresh), 'category'] = 'Down'

color_map = {'Not Significant': '#BDBDBD', 'Up': '#E53935', 'Down': '#1E88E5'}
fig, ax = plt.subplots(figsize=(10, 8))
for cat in ['Not Significant', 'Up', 'Down']:
    sub = df[df['category'] == cat]
    ax.scatter(sub['log2FC'], sub['neg_log10_p'],
               c=color_map[cat], s=12, alpha=0.6, edgecolors='none', label=f'{cat} ({len(sub)})')

ax.axhline(-np.log10(p_thresh), color='black', linestyle='--', linewidth=0.8)
ax.axvline(fc_thresh, color='black', linestyle='--', linewidth=0.8)
ax.axvline(-fc_thresh, color='black', linestyle='--', linewidth=0.8)

# Label top hits
for _, row in df.nlargest(8, 'neg_log10_p').iterrows():
    ax.annotate(row['gene'], xy=(row['log2FC'], row['neg_log10_p']),
                xytext=(5, 5), textcoords='offset points', fontsize=8, fontweight='bold',
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

ax.set_xlabel('log$_2$(Fold Change)', fontsize=13)
ax.set_ylabel('-log$_{10}$(p-value)', fontsize=13)
ax.legend(title='Category', fontsize=10)
```

### MA Plot
```python
A = np.log2(base_mean + 1)   # average expression (log2 scale, pseudocount for zeros)
M = log2fc                # fold change

fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(A[~sig], M[~sig], s=8, alpha=0.3, c='gray', label='Not significant')
ax.scatter(A[sig & (M > 0)], M[sig & (M > 0)], s=20, alpha=0.7, c='red', label='Up')
ax.scatter(A[sig & (M < 0)], M[sig & (M < 0)], s=20, alpha=0.7, c='blue', label='Down')

# Optional lowess-style trend line
z = np.polyfit(A, M, 3)
x_smooth = np.linspace(A.min(), A.max(), 200)
ax.plot(x_smooth, np.poly1d(z)(x_smooth), 'k-', linewidth=1.5, alpha=0.6)

ax.axhline(0, color='black', linestyle='--', linewidth=0.8)
ax.set_xlabel('A: log$_2$(Base Mean Expression)', fontsize=12)
ax.set_ylabel('M: log$_2$(Fold Change)', fontsize=12)
```

### Clustered Heatmap
```python
# Z-score normalize rows before clustering
expr_z = (expr - expr.mean(axis=1, keepdims=True)) / expr.std(axis=1, keepdims=True)
df_z = pd.DataFrame(expr_z, index=gene_labels, columns=samples)

# Condition color bar
condition_colors = ['#1E88E5'] * n_ctrl + ['#E53935'] * n_treat
col_colors = pd.Series(condition_colors, index=samples, name='Condition')

g = sns.clustermap(df_z,
                   cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
                   figsize=(10, 12),
                   col_colors=col_colors,
                   dendrogram_ratio=(0.12, 0.08),
                   cbar_pos=(0.02, 0.82, 0.03, 0.12),
                   linewidths=0.3, linecolor='white')

# Add legend for condition colors
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#1E88E5', label='Control'),
                   Patch(facecolor='#E53935', label='Drug')]
g.ax_heatmap.legend(handles=legend_elements, loc='lower right',
                    bbox_to_anchor=(1.3, 1.02), fontsize=9)
```

### PCA Plot with Confidence Ellipses
```python
from sklearn.decomposition import PCA
from matplotlib.patches import Ellipse

pca = PCA(n_components=2)
pcs = pca.fit_transform(log_expr.T)   # samples as rows

palette = {'WT': '#2196F3', 'KO': '#E53935', 'Drug': '#4CAF50'}
markers = {'WT': 'o', 'KO': 's', 'Drug': '^'}

fig, ax = plt.subplots(figsize=(9, 7))
for cond, color in palette.items():
    sub = pca_df[pca_df['Condition'] == cond]
    ax.scatter(sub['PC1'], sub['PC2'], c=color, marker=markers[cond],
               s=120, label=cond, edgecolors='black', linewidth=0.8)
    ellipse = Ellipse((sub['PC1'].mean(), sub['PC2'].mean()),
                      width=sub['PC1'].std()*3, height=sub['PC2'].std()*3,
                      fill=True, alpha=0.1, color=color)
    ax.add_patch(ellipse)

ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
```

### GC Content Sliding Window
```python
def sliding_gc(sequence, window_size=100, step=10):
    positions, gc_values = [], []
    for i in range(0, len(sequence) - window_size + 1, step):
        w = sequence[i:i + window_size]
        gc = (w.count('G') + w.count('C')) / window_size * 100
        positions.append(i + window_size // 2)
        gc_values.append(gc)
    return np.array(positions), np.array(gc_values)

fig, ax = plt.subplots(figsize=(14, 4))
ax.fill_between(positions, gc_values, 50,
                where=(gc_values >= 50), color='#E53935', alpha=0.4, label='GC-rich (>50%)')
ax.fill_between(positions, gc_values, 50,
                where=(gc_values < 50), color='#1E88E5', alpha=0.4, label='AT-rich (<50%)')
ax.plot(positions, gc_values, 'k-', linewidth=0.8)
ax.axhline(50, color='gray', linestyle='--', linewidth=1)
```

### Genome Coverage + Gene Track
```python
fig, (ax_cov, ax_gene) = plt.subplots(2, 1, figsize=(14, 6),
                                       height_ratios=[4, 1],
                                       sharex=True, gridspec_kw={'hspace': 0.05})
# Smooth coverage
smoothed = np.convolve(base_coverage, np.ones(20)/20, mode='same')
ax_cov.fill_between(positions, smoothed, alpha=0.3, color='steelblue')
ax_cov.plot(positions, smoothed, color='steelblue', linewidth=0.8)
ax_cov.axhline(expected_cov, color='black', linestyle='--', linewidth=0.8, alpha=0.5)

# Highlight CNV regions
ax_cov.axvspan(del_start, del_end, alpha=0.15, color='red')
ax_cov.axvspan(dup_start, dup_end, alpha=0.15, color='blue')

# Gene arrows
for gene in genes:
    ax_gene.annotate('', xy=(gene['end'], 0.5), xytext=(gene['start'], 0.5),
                     arrowprops=dict(arrowstyle='->', color='steelblue', lw=8))
    ax_gene.text((gene['start']+gene['end'])/2, 0.1, gene['name'],
                 ha='center', fontsize=9, fontweight='bold')
ax_gene.set_yticks([])
```

### Per-Base Sequencing Quality (FastQC-style)
```python
fig, ax = plt.subplots(figsize=(14, 5))
bp = ax.boxplot([quality_array[i, :] for i in range(0, read_length, step)],
                positions=base_positions[::step], widths=2,
                patch_artist=True, showfliers=False,
                medianprops=dict(color='red', linewidth=1.5))

for i, box in enumerate(bp['boxes']):
    median = np.median(box_data[i])
    box.set_facecolor('#A5D6A7' if median >= 30 else ('#FFF176' if median >= 20 else '#EF9A9A'))
    box.set_alpha(0.8)

ax.axhline(30, color='green', linestyle='--', linewidth=1, label='Q30')
ax.axhline(20, color='orange', linestyle='--', linewidth=1, label='Q20')
ax.axhspan(30, 42, alpha=0.05, color='green')
ax.axhspan(20, 30, alpha=0.05, color='yellow')
ax.axhspan(0, 20, alpha=0.05, color='red')
ax.set_ylim(0, 42)
```

### Kaplan-Meier Survival Curve
```python
def kaplan_meier(times, events):
    """times and events are parallel lists; events: 1=death, 0=censored."""
    pairs = sorted(zip(times, events))
    km_times, km_survival = [0], [1.0]
    censored_times, censored_surv = [], []
    survival = 1.0
    n = len(pairs)
    i = 0
    while i < n:
        t, e = pairs[i]
        at_risk = n - i
        if e == 1:
            # count all deaths at this time point
            deaths = sum(1 for tt, ee in pairs[i:] if tt == t and ee == 1)
            survival *= (at_risk - deaths) / at_risk
            km_times.append(t); km_survival.append(survival)
            i += sum(1 for tt, _ in pairs[i:] if tt == t)
        else:
            censored_times.append(t); censored_surv.append(survival)
            i += 1
    return km_times, km_survival, censored_times, censored_surv

fig, ax = plt.subplots(figsize=(10, 6))
ax.step(km_t_a, km_s_a, where='post', color='#2196F3', linewidth=2, label='Group A')
ax.step(km_t_b, km_s_b, where='post', color='#E53935', linewidth=2, label='Group B')
ax.plot(cens_t_a, cens_s_a, '|', color='#2196F3', markersize=10, markeredgewidth=2)
ax.plot(cens_t_b, cens_s_b, '|', color='#E53935', markersize=10, markeredgewidth=2)
ax.set_ylim(0, 1.05); ax.set_ylabel('Survival Probability')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
```

### Seaborn Statistical Plots
```python
# Box + violin + strip + KDE (pick one per use case)
sns.boxplot(data=df, x='Cell_Type', y='expression', hue='Condition', palette='Set2', ax=ax)
sns.violinplot(data=df, x='Cell_Type', y='expression', hue='Condition',
               split=True, inner='quartile', palette='muted', ax=ax)
sns.stripplot(data=df, x='Cell_Type', y='expression', hue='Condition',
              dodge=True, alpha=0.6, size=4, palette='deep', ax=ax)
sns.kdeplot(data=subset, x='expression', fill=True, alpha=0.3, label='group', ax=ax)

# Pair plot
g = sns.pairplot(df[['Condition', 'gene1', 'gene2', 'gene3']], hue='Condition',
                 palette='colorblind', plot_kws={'alpha': 0.5, 's': 20})
```

## Common Pitfalls
- **Row normalization before heatmap:** always z-score rows (`axis=1`); raw counts produce unreadable heatmaps
- **log p-values of 0:** guard with `np.clip(pvalues, 1e-300, 1)` before `-log10()`
- **PCA direction:** pass samples as rows to `fit_transform` (`expr.T` if genes×samples)
- **GridSpec shared axis:** use `sharex=True` + `gridspec_kw={'hspace': 0.05}` for genome tracks
- **clustermap layout:** `cbar_pos`, `dendrogram_ratio` must be tuned per figure size; use `fig.suptitle(..., y=1.01)` to avoid overlap
- **Font substitution:** specify `['Arial', 'Helvetica', 'DejaVu Sans']` in order; journals require Arial/Helvetica
- **Saving figures:** call `fig.savefig()` inside `plt.rc_context(pub_params)` block; use `.svg` for editable vector output

## Related Skills
- `rnaseq` — DESeq2/edgeR results feeding into volcano/MA plots
- `numpy-pandas-wrangling` — data wrangling to prepare DataFrames for seaborn
- `ml-deep-learning-bio` — dimensionality reduction and unsupervised analysis
