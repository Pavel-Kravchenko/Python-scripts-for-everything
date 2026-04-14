---
name: python-bio-data-wrangling
description: Pandas patterns for bio data — missing values, duplicates, type conversion, wide/long reshaping, string parsing, apply/transform/pipe.
tool_type: python
primary_tool: NumPy
---

## Complicated Moments

**`fillna` is not in-place by default**: `df['col'].fillna(0)` returns a new Series. Use `df['col'] = df['col'].fillna(0)` or `df.fillna({'col': 0}, inplace=True)`.

**`melt` vs `stack`**: `melt` converts specific value columns to rows (wide → long), keeping id columns. `stack` pivots the innermost column level into the row index. For expression data reshaping, `melt` is usually what you want.

**`groupby` and the index**: `groupby('gene_id').agg(...)` moves `gene_id` into the index by default. Use `as_index=False` to keep it as a regular column.

**String ops on NaN columns**: `df['col'].str.upper()` returns `NaN` for missing values silently. Always check for unexpected NaNs after parsing annotation columns.

## Missing Values

```python
# Strategy selection for expression data:
# Drop:           missingness is small and random (failed library preps)
# Fill with mean: need complete matrix for PCA/clustering, low missingness
# Fill with 0:    count data where NaN = "no reads detected"

sample_cols = ['s1', 's2', 's3']

# Column-mean imputation (sample-wise)
df[sample_cols] = df[sample_cols].fillna(df[sample_cols].mean())

# Row-mean imputation (gene-wise)
row_means = df[sample_cols].mean(axis=1)
for col in sample_cols:
    df[col] = df[col].fillna(row_means)
```

## Duplicates

```python
# Keep first occurrence
deduped = df.drop_duplicates(subset='gene_id', keep='first')

# Average technical replicates (most correct for expression data)
averaged = df.groupby('gene_id', as_index=False).agg({
    's1': 'mean', 's2': 'mean', 's3': 'mean',
    'gene_type': 'first',   # keep first annotation
})
```

## Type Conversion

```python
# Clinical/annotation tables often load everything as str
clinical['age'] = pd.to_numeric(clinical['age'], errors='coerce')      # 'unknown' → NaN
clinical['mut_count'] = pd.to_numeric(clinical['mut_count'], errors='coerce')
clinical['alive'] = clinical['alive'].map({'True': True, 'False': False})

# Ordered categorical — enables comparison operators
stage_order = ['I', 'II', 'III', 'IV']
clinical['tumor_stage'] = pd.Categorical(
    clinical['tumor_stage'], categories=stage_order, ordered=True
)
clinical[clinical['tumor_stage'] >= 'III']   # works with ordered categorical
```

## Wide ↔ Long Reshaping

```
Wide (one row per gene):  gene | s1 | s2 | s3
Long (one row per obs):   gene | sample | expression
```

Wide → long for seaborn/ggplot/statistical models. Wide for NumPy matrix ops and heatmaps.

```python
# Wide → Long
long = df.melt(
    id_vars=['gene', 'gene_type'],
    value_vars=['ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'],
    var_name='sample',
    value_name='expression',
)
long['condition'] = long['sample'].str.split('_').str[0]   # 'ctrl' / 'treat'

# Long → Wide
wide = long.pivot_table(
    index=['gene', 'gene_type'],
    columns='sample',
    values='expression',
).reset_index()
wide.columns.name = None   # flatten MultiIndex column names

# stack/unstack for MultiIndex DataFrames
stacked = df.set_index(['gene', 'gene_type'])[sample_cols].stack()
stacked.name = 'expression'
unstacked = stacked.unstack()
```

## String Parsing (GTF-style attributes)

```python
# Parse GTF attribute column: gene_id "ENSG..."; gene_name "TP53"; ...
df['ensembl_id'] = df['raw_id'].str.extract(r'gene_id "(ENSG\d+)"')
df['gene_name']  = df['raw_id'].str.extract(r'gene_name "(\w+)"')
df['biotype']    = df['raw_id'].str.extract(r'gene_biotype "([^"]+)"')

# Cytogenetic band parsing: '17p13.1' → chrom='17', arm='p'
chromosomes = bands.str.extract(r'(\d+)[pq]')[0]
arms        = bands.str.extract(r'\d+([pq])')[0]
```

## Apply / Transform / Pipe

| Method | Returns | Use case |
|--------|---------|----------|
| `apply(fn, axis=1)` | Anything | Row-wise custom logic |
| `apply(fn, axis=0)` | Anything | Column-wise (e.g. z-score) |
| `transform(fn)` | Same shape as input | Group-wise normalization |
| `pipe(fn)` | DataFrame | Chaining processing steps |

```python
# Row-wise classification
def classify_expr(row):
    m = row[sample_cols].mean()
    return 'high' if m > 6 else 'medium' if m > 4 else 'low'

df['category'] = df.apply(classify_expr, axis=1)

# Column-wise z-score normalization
def zscore(col):
    return (col - col.mean()) / col.std()

normalized = df[sample_cols].apply(zscore)

# Group-wise normalization with transform (keeps original shape)
df['expr_zscore'] = df.groupby('condition')['expression'].transform(zscore)
```

## Pitfalls

- **Chain indexing**: `df[mask]['col'] = value` silently sets a copy. Always use `df.loc[mask, 'col'] = value`.
- **`groupby` moves key into index**: use `as_index=False` to keep it as a column, or call `.reset_index()` after.
- **`pd.to_numeric(errors='coerce')`**: silently converts unparseable strings to NaN; always inspect `.isna().sum()` afterward.
- **`melt` column order**: `value_vars` defaults to all non-id columns if omitted; be explicit to avoid including unexpected columns.
- **Off-by-one in coordinates**: BED is 0-based half-open; VCF/GFF are 1-based. Converting between them without adjusting causes systematic errors.
- **String ops on mixed-case gene IDs**: normalize with `.str.upper()` before any merge or filter to avoid missed matches.
