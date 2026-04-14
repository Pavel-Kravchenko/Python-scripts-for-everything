---
name: python-bio-data-visualization
description: "Matplotlib and seaborn for bioinformatics figures: volcano plots, heatmaps, MA plots, genome tracks. Critical gotchas for publication-quality figures."
tool_type: python
primary_tool: NumPy
---

# Data Visualization for Bioinformatics

## Critical Gotchas

- **`savefig` before `plt.show()`**: `show()` clears the figure; if you call it first, `savefig` writes a blank file.
- **Log axes and zeros**: `plt.yscale('log')` fails on zeros. Log-transform before plotting: `np.log2(counts + 1)`, not log-scale axis on raw counts.
- **Seaborn expects long-form data**: if data is wide (one column per sample), melt first: `pd.melt(df, id_vars=['gene'], value_vars=samples)`.
- **OO interface for multi-panel**: use `fig, ax = plt.subplots()` + `ax.plot()` for anything with multiple panels. `plt.plot()` acts on the current active axes and breaks in loops.
- **Colormaps**: use `RdBu_r` or `viridis` for expression heatmaps, never `jet`. Use seaborn `colorblind` palette for categorical data.

## Plot Type Reference

| Plot | Use case | Key params |
|------|----------|-----------|
| `ax.scatter` | expression correlation, UMAP | `alpha`, `s`, `c` |
| `ax.bar` / `barh` | enrichment, composition | `yerr` for error bars |
| `ax.hist` | quality score distribution | `bins`, `range` |
| `sns.heatmap` | expression matrix | `cmap='RdBu_r'`, `center=0` |
| Volcano | DE results | log2FC vs −log10(padj) |
| MA plot | DE mean vs fold-change | log10(baseMean) vs log2FC |

## Common Figure Patterns

```python
import matplotlib.pyplot as plt
import numpy as np

# Always use OO interface
fig, ax = plt.subplots(figsize=(8, 5))

# Multi-panel with GridSpec
fig = plt.figure(figsize=(12, 7))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)
ax_top = fig.add_subplot(gs[0, :])   # spans both columns
ax_bl  = fig.add_subplot(gs[1, 0])
ax_br  = fig.add_subplot(gs[1, 1])
```

## Volcano Plot

```python
def plot_volcano(log2fc, neg_log10_padj, ax=None, fc_thresh=1.0, padj_thresh=0.05):
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))
    sig = (np.abs(log2fc) > fc_thresh) & (neg_log10_padj > -np.log10(padj_thresh))
    ax.scatter(log2fc[~sig], neg_log10_padj[~sig], s=5, alpha=0.4, c='gray')
    ax.scatter(log2fc[sig],  neg_log10_padj[sig],  s=8, alpha=0.7, c='red')
    ax.axvline(-fc_thresh, color='gray', linestyle='--', linewidth=0.8)
    ax.axvline( fc_thresh, color='gray', linestyle='--', linewidth=0.8)
    ax.axhline(-np.log10(padj_thresh), color='gray', linestyle='--', linewidth=0.8)
    ax.set_xlabel('log$_2$ Fold Change')
    ax.set_ylabel('-log$_{10}$(padj)')
    return ax
```

## Heatmap Pattern

```python
import seaborn as sns

# Z-score rows before plotting expression heatmaps
from scipy.stats import zscore
mat_z = zscore(mat, axis=1)  # per-gene z-score across samples

g = sns.clustermap(
    mat_z,
    cmap='RdBu_r', center=0, vmin=-3, vmax=3,
    row_cluster=True, col_cluster=True,
    figsize=(10, 8), yticklabels=False
)
```

## Colormap Quick Reference

| Data type | Recommended colormap |
|-----------|---------------------|
| Expression fold-change (diverging) | `RdBu_r`, `coolwarm` |
| Expression level (sequential) | `viridis`, `YlOrRd` |
| Categorical (cell types) | `sns.color_palette('colorblind')` |
| P-value / significance | `plasma` or custom threshold-based |

## Saving Figures

```python
fig.savefig('figure.pdf', dpi=300, bbox_inches='tight')   # vector for publication
fig.savefig('figure.png', dpi=150, bbox_inches='tight')   # raster for display
plt.close(fig)  # free memory in loops
```

## Pitfalls

- **Legend outside plot area gets clipped on save**: use `bbox_inches='tight'` in `savefig`
- **Tick labels overlap**: rotate with `ax.set_xticklabels(labels, rotation=45, ha='right')`
- **seaborn changes global rcParams**: import order matters; set `plt.rcParams` after seaborn imports
- **`clustermap` returns a `ClusterGrid`, not a `Figure`**: access figure via `g.fig`, axes via `g.ax_heatmap`
