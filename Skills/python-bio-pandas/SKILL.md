---
name: python-bio-pandas
description: NumPy arrays and Pandas DataFrames for bioinformatics — vectorized ops, broadcasting, CPM/RPKM normalization, PWM scoring, loc/iloc, groupby, and annotation merges.
tool_type: python
primary_tool: Pandas
---

## Complicated Moments

**`loc` vs `iloc` vs `[]`**: `df['col']` selects a column. `df.loc[row_label, col_label]` selects by label. `df.iloc[row_int, col_int]` selects by integer position. After filtering, integer positions no longer match labels — always know which accessor you need.

**Chain indexing creates copies unpredictably**: `df[df['gc'] > 0.5]['length'] = 100` may silently fail. Always use `df.loc[mask, 'length'] = 100`.

**`groupby` + `transform` vs `agg`**: `agg` reduces each group to one row. `transform` returns a same-shape array with each row filled with its group statistic — ideal for group-wise normalization.

**Left join for annotation merges**: inner join silently drops genes missing from the annotation table. Use left join and inspect the resulting NaNs.

## NumPy Core Patterns

```python
import numpy as np

# Gene lengths (vectorized, no loop)
lengths = stops - starts + 1

# Log-transform counts
log_counts = np.log2(counts + 1)   # +1 avoids log(0)

# RPKM normalization
def compute_rpkm(counts, lengths):
    rpm  = counts * (1_000_000 / np.sum(counts))
    return rpm * (1_000 / lengths)

# CPM via broadcasting: (4 genes × 3 samples) / (1 × 3) row vector
sample_totals = expr.sum(axis=0, keepdims=True)   # shape (1, 3)
cpm = expr / sample_totals * 1_000_000

# Z-score each gene (row) across samples
gene_means = expr.mean(axis=1, keepdims=True)     # shape (4, 1)
gene_stds  = expr.std(axis=1, keepdims=True)
z_scores   = (expr - gene_means) / gene_stds
```

### Indexing
```python
arr[0]          # first element
arr[-1]         # last element
arr[2:5]        # slice
arr[::2]        # every other

mat[0]          # entire row 0
mat[:, 1]       # entire column 1
mat[2, 2]       # single element
mat[:2, :2]     # submatrix

# Boolean indexing
reads[reads > 50]
reads[reads != 0]
```

### Sliding-Window GC Content
```python
def sliding_gc(sequence, window=50):
    gc_binary = np.array([1 if n in 'GC' else 0 for n in sequence])
    cumsum = np.insert(np.cumsum(gc_binary), 0, 0)
    return (cumsum[window:] - cumsum[:-window]) / window * 100
```

### PWM Scoring
```python
BASES = list('ATGC')

def build_pwm(sequences):
    counts = np.zeros((4, len(sequences[0])))
    for seq in sequences:
        for i, nuc in enumerate(seq):
            counts[BASES.index(nuc), i] += 1
    ppm = counts / len(sequences)
    # avoid log(0): add pseudocount before calling build_pwm, or:
    freq = np.where(ppm > 0, ppm, 1e-9)
    return np.log2(freq / 0.25)   # log-odds PWM

def score_sequence(pwm, sequence):
    return sum(pwm[BASES.index(nuc), i] for i, nuc in enumerate(sequence))
```

## Pandas Core Patterns

```python
import pandas as pd

# Quick overview
df.shape, df.dtypes, df.describe(), df.info()

# loc vs iloc
df.loc[0, 'gene']      # by label
df.iloc[0, 0]          # by integer position
df.loc[mask, 'col']    # boolean mask selection (safe for assignment)
```

### Annotation Merge
```python
# Always use left join — inner join silently drops unmatched genes
merged = expr_df.merge(annotation_df, on='gene_id', how='left')
print(merged['gene_name'].isna().sum(), 'genes missing annotation')
```

### GroupBy Patterns
```python
# Aggregate: one row per group
summary = df.groupby('condition', as_index=False).agg(
    mean_expr=('expression', 'mean'),
    n=('expression', 'count'),
)

# Transform: same shape as input, filled with group statistic
df['expr_zscore'] = df.groupby('condition')['expression'].transform(
    lambda x: (x - x.mean()) / x.std()
)
```

### Reading Expression Data
```python
# Count matrix: genes × samples
counts = pd.read_csv('counts.csv', index_col=0)    # genes as index

# Long format from melted matrix
long = counts.reset_index().melt(
    id_vars='gene_id', var_name='sample', value_name='count'
)
```

## Pitfalls

- **Chain indexing**: `df[mask]['col'] = val` sets a copy silently. Use `df.loc[mask, 'col'] = val`.
- **`groupby` index**: by default the grouping key becomes the index; use `as_index=False` to keep it as a column.
- **Integer index after filtering**: after `df = df[df['qc_pass']]`, `df.iloc[0]` is the first remaining row but `df.loc[0]` still refers to original label 0 (may raise KeyError). Reset index with `df.reset_index(drop=True)` if you need position-based access.
- **`np.log2` on zero counts**: always add a pseudocount (`+ 1` or `+ 0.5`) before log-transforming read counts.
- **Broadcasting axis confusion**: `expr.mean(axis=0)` averages across genes (per sample); `axis=1` averages across samples (per gene). Wrong axis is a silent error.
- **`np.corrcoef` returns a matrix**: `corrcoef(a, b)[0, 1]` is the Pearson r, not `corrcoef(a, b)`.
