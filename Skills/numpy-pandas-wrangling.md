---
name: numpy-pandas-wrangling
description: NumPy arrays and Pandas DataFrames for expression matrices, annotation merges, and tidy biological data
---

# NumPy, Pandas & Data Wrangling for Bioinformatics

## When to Use
- Building or normalizing expression matrices (CPM, RPKM, Z-score)
- Merging expression data with gene annotations or sample metadata
- Reshaping between wide (matrix) and long (tidy) format for plotting/stats
- Cleaning messy annotation files (gene names, chromosomes, missing values)
- Scoring sequences against a PWM or computing sliding-window statistics

---

## Quick Reference

### NumPy array creation
```python
expr = np.array([[120, 135], [45, 50]], dtype=float)  # expression matrix
np.zeros((n_genes, n_samples))
np.full((3, 3), np.nan)
np.arange(0, 1000, 100)        # genomic positions
np.linspace(0, 1, 5)
```

### Key NumPy axes
- `axis=0` — operates down rows (per-column result, e.g., per-sample total)
- `axis=1` — operates across columns (per-row result, e.g., per-gene mean)
- `keepdims=True` — preserves shape for broadcasting

### Pandas I/O
```python
pd.read_csv('expr.csv', index_col='gene_id')
pd.read_csv('data.tsv', sep='\t')
pd.read_csv('genes.gtf', sep='\t', comment='#', header=None)
df.to_csv('out.csv', index=False)
```

---

## Key Patterns

### CPM normalization (broadcasting)
```python
sample_totals = expr.sum(axis=0, keepdims=True)   # shape (1, n_samples)
cpm = expr / sample_totals * 1_000_000
```

### RPKM normalization
```python
def compute_rpkm(counts, lengths):
    rpm = counts * (1_000_000 / np.sum(counts))
    return rpm * (1_000 / lengths)
```

### Z-score per gene across samples
```python
gene_means = expr.mean(axis=1, keepdims=True)
gene_stds  = expr.std(axis=1, keepdims=True)
z_scores   = (expr - gene_means) / gene_stds
```

### Log2 fold change
```python
ctrl_mean  = expr_df[ctrl_cols].mean(axis=1)
treat_mean = expr_df[treat_cols].mean(axis=1)
log2fc     = np.log2(treat_mean / ctrl_mean)
```

### Boolean masking
```python
reads[reads > 50]                    # 1D filter
expression[0]                        # row 0 (Gene A)
expression[:, 1]                     # column 1 (Sample 2)
df[(df['chr'] == '17') & (df['gc'] > 0.45)]   # multi-condition
df.query('chr == "17" and gc > 0.45')          # equivalent, more readable
```

### GroupBy aggregations
```python
df.groupby('biotype')['expression'].mean()
df.groupby('chromosome')['expression'].agg(['count', 'mean', 'median', 'std'])
# Average technical replicates keyed by gene_id
df.groupby('gene_id', as_index=False).agg({'sample_1': 'mean', 'gene_type': 'first'})
```

### Merge patterns
```python
# Inner: genes present in both tables
pd.merge(expr, annot, on='gene_id', how='inner')
# Left: keep all expression rows, NaN for unmatched annotations
pd.merge(expr, annot, on='gene_id', how='left')
# Two-step: expression -> id_map -> pathways
expr.merge(id_map, on='ensembl_id').merge(pathways, on='gene_name', how='left')
```

---

## Code Templates

### Position Weight Matrix (PWM) from aligned sequences
```python
def build_pwm(sequences):
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    counts = np.zeros((4, len(sequences[0])))
    for seq in sequences:
        for i, nuc in enumerate(seq):
            counts[mapping[nuc], i] += 1
    return counts / len(sequences)   # shape (4, L), rows = A/T/G/C

def score_sequence(pwm, sequence):
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    score = sum(np.log2(pwm[mapping[n], i] / 0.25) if pwm[mapping[n], i] > 0 else -10
                for i, n in enumerate(sequence))
    return score
```

### Sliding-window GC content (O(n) via cumsum)
```python
def sliding_gc(sequence, window=50):
    gc = np.array([1 if n in 'GC' else 0 for n in sequence])
    cumsum = np.insert(np.cumsum(gc), 0, 0)
    return (cumsum[window:] - cumsum[:-window]) / window * 100
```

### Wide -> Long -> annotated pipeline
```python
# Melt expression matrix
long = wide_expr.melt(
    id_vars=['gene', 'gene_type'],
    value_vars=['ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'],
    var_name='sample',
    value_name='expression',
)
long['condition'] = long['sample'].str.split('_').str[0]

# Merge sample metadata
long = long.merge(sample_metadata, on='sample')

# Compute per-gene per-condition mean, then fold change
summary = long.groupby(['gene', 'condition'])['expression'].mean().unstack()
summary['log2FC'] = np.log2(summary['treatment'] / summary['control'])
```

### Long -> Wide (pivot back)
```python
wide = long.pivot_table(index=['gene', 'gene_type'], columns='sample', values='expression')
wide.columns.name = None
wide = wide.reset_index()
```

### pipe()-based processing pipeline
```python
def remove_low_expression(df, threshold=3.0):
    cols = [c for c in df.columns if c.startswith(('ctrl', 'treat'))]
    return df[df[cols].max(axis=1) >= threshold]

def add_fold_change(df):
    df = df.copy()
    df['log2FC'] = np.log2(df[treat_cols].mean(axis=1) / df[ctrl_cols].mean(axis=1))
    return df

result = df.pipe(remove_low_expression, threshold=2.5).pipe(add_fold_change)
```

### group-wise Z-score normalization (transform)
```python
long['zscore'] = long.groupby('gene_type')['expression'].transform(
    lambda x: (x - x.mean()) / x.std()
)
```

### Cleaning annotation tables
```python
df['gene_name'] = df['gene_name'].str.strip().str.upper()
df['chromosome'] = df['chromosome'].str.lower()
df = df.drop_duplicates(subset='gene_name', keep='first')
# Parse GTF attribute string
df['ensembl_id'] = df['raw_id'].str.extract(r'gene_id "(ENSG\d+)"')
df['gene_name']  = df['raw_id'].str.extract(r'gene_name "(\w+)"')
df['biotype']    = df['raw_id'].str.extract(r'gene_biotype "([^"]+)"')
# Cytogenetic band parsing
chromosomes = bands.str.extract(r'(\d+)[pq]')[0]
arms        = bands.str.extract(r'\d+([pq])')[0]
```

### Missing value imputation strategies
```python
# Drop (small random missingness)
df.dropna()

# Fill with column mean (sample-wise imputation)
df[sample_cols] = df[sample_cols].fillna(df[sample_cols].mean())

# Fill with row median (gene-wise imputation — preferred for expression)
row_medians = df[sample_cols].median(axis=1)
for col in sample_cols:
    df[col] = df[col].fillna(row_medians)

# Fill count data where NaN means zero reads
df[sample_cols] = df[sample_cols].fillna(0)

# Clinical categoricals
df['tumor_stage'] = df['tumor_stage'].fillna('Unknown')
df['age'] = df['age'].fillna(df['age'].median())
```

### Type conversion for clinical tables
```python
df['age']            = pd.to_numeric(df['age'], errors='coerce')   # 'unknown' -> NaN
df['mutation_count'] = pd.to_numeric(df['mutation_count'], errors='coerce')
df['alive']          = df['alive'].map({'True': True, 'False': False})
stage_order = ['I', 'II', 'III', 'IV']
df['tumor_stage'] = pd.Categorical(df['tumor_stage'], categories=stage_order, ordered=True)
# Now supports: df[df['tumor_stage'] >= 'III']
```

---

## Common Pitfalls

- **Log of zero**: always use `np.log2(counts + 1)` (pseudocount) with raw counts.
- **Fold change direction**: `log2(treat / ctrl)`, not `ctrl / treat`.
- **Boolean combine**: use `&`/`|`/`~` with parentheses, not `and`/`or`/`not`.
- **Chained assignment**: use `df = df.copy()` before mutating inside a function; avoid `df[col][mask] = value`.
- **`loc` vs `iloc`**: `loc` is label-based (inclusive on both ends for slices); `iloc` is integer-position (exclusive end).
- **`apply` on columns vs rows**: `axis=0` = column-wise (default), `axis=1` = row-wise.
- **Missing annotation keys after merge**: use `how='left'` to keep expression rows that lack annotation; inspect with `merged[merged['gene_name'].isna()]`.
- **Duplicate gene symbols after merge**: check with `df.duplicated(subset='gene_id')` before merging; average with `groupby(...).mean()`.
- **`pivot_table` MultiIndex columns**: flatten with `df.columns.name = None` after reset.
- **Wide-format for NumPy, long-format for seaborn/statsmodels**: convert before plotting.

---

## Related Skills
- `rnaseq` — downstream from the expression matrix built here
- `data-visualization-bio` — seaborn/matplotlib require long-format DataFrames
- `biopython-databases` — BioPython complements NumPy for sequence-level work
