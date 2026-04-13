---
name: python-bio-data-wrangling
description: "- Handle missing values and duplicates in biological datasets - Convert data types and clean messy annotation tables - Reshape data between wide and long formats (melt, pivot, stack/unstack) - Use str"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/17_Data_Wrangling/01_data_wrangling.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 17: Data Wrangling with Pandas

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/17_Data_Wrangling/01_data_wrangling.ipynb`*

# Module 17: Data Wrangling with Pandas

**Estimated time: 90-120 minutes**

---

## Learning Objectives

- Handle missing values and duplicates in biological datasets
- Convert data types and clean messy annotation tables
- Reshape data between wide and long formats (melt, pivot, stack/unstack)
- Use string operations to extract information from annotation columns
- Apply custom transformations with `apply`, `transform`, and `pipe`
- Merge gene annotations with expression data

### Why Data Wrangling?

Real biological data is messy. Gene names have inconsistent capitalization, clinical tables have missing values, expression matrices arrive in the wrong shape, and annotation files mix numeric IDs with free text. Before any analysis, you must clean and reshape the data. This module teaches the Pandas tools for that.

```
Raw Data        Clean          Reshape         Analyze
 Missing   -->  Fill/Drop  -->  Melt/Pivot -->  Ready for
 Dupes          Fix types       Merge           statistics
 Messy          Strings         Stack
```

---

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. `fillna` is not in-place by default**
`df['col'].fillna(0)` returns a new Series — it does not modify `df`. Use `df['col'] = df['col'].fillna(0)` or `df.fillna({'col': 0}, inplace=True)`.

**2. `melt` vs `stack`**
`melt` converts specific value columns to rows (wide → long), keeping id columns unchanged. `stack` pivots the innermost column level into the row index. For typical expression data reshaping, `melt` is what you want.

**3. `groupby` and the index**
By default, `groupby('gene_id').agg(...)` moves `gene_id` into the index. Use `as_index=False` to keep it as a regular column, avoiding surprises when merging downstream.

**4. String operations on columns with NaN**
`df['col'].str.upper()` silently returns `NaN` for missing values instead of raising an error. Always check for unexpected NaNs after string operations on annotation columns.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(42)
```

```python
# A messy gene expression dataset with missing values and duplicates
messy_expr = pd.DataFrame({
    'gene_id':  ['TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS', 'TP53', 'PIK3CA'],
    'sample_1': [5.2,    np.nan,   7.1,    2.4,   4.5,    5.3,    6.3],
    'sample_2': [5.4,    3.8,      7.3,    np.nan, 4.2,   5.1,    6.1],
    'sample_3': [5.0,    3.5,      6.9,    2.6,   4.8,    5.5,    np.nan],
    'gene_type': ['TSG', 'TSG', 'Oncogene', 'Oncogene', 'Oncogene', 'TSG', 'Oncogene'],
})

print("Messy expression data:")
print(messy_expr)
print(f"\nMissing values per column:\n{messy_expr.isna().sum()}")
print(f"\nTotal missing values: {messy_expr.isna().sum().sum()}")
```

```python
# Strategy 1: Drop rows with any missing value
dropped = messy_expr.dropna()
print(f"After dropna(): {len(dropped)} rows (from {len(messy_expr)})")
print(dropped)
```

```python
# Strategy 2: Fill with column mean (common for expression data)
sample_cols = ['sample_1', 'sample_2', 'sample_3']
filled = messy_expr.copy()
filled[sample_cols] = filled[sample_cols].fillna(filled[sample_cols].mean())

print("After filling NaN with column means:")
print(filled)
```

```python
# Strategy 3: Fill with row mean (gene-wise imputation)
filled_row = messy_expr.copy()
row_means = filled_row[sample_cols].mean(axis=1)
for col in sample_cols:
    filled_row[col] = filled_row[col].fillna(row_means)

print("After filling NaN with row means:")
print(filled_row)
```

### Biological context: when to drop vs. fill

- **Drop** when the fraction of missing data is small and the missingness is random (e.g., a few failed library preps). Dropping preserves unbiased statistics.
- **Fill with mean/median** when you need a complete matrix for downstream analysis (e.g., clustering, PCA) and missingness is low.
- **Fill with 0** for count data when NaN genuinely means "no reads detected."
- Never silently fill without documenting which strategy you used.

---

## Part 2: Handling Duplicates

```python
# Detect duplicates
print("Duplicate gene_ids:")
print(filled[filled.duplicated(subset='gene_id', keep=False)])
```

```python
# Strategy 1: Keep the first occurrence
deduped = filled.drop_duplicates(subset='gene_id', keep='first')
print(f"After drop_duplicates (keep='first'): {len(deduped)} rows")
print(deduped)
```

```python
# Strategy 2: Average the duplicates (common for technical replicates)
averaged = filled.groupby('gene_id', as_index=False).agg({
    'sample_1': 'mean',
    'sample_2': 'mean',
    'sample_3': 'mean',
    'gene_type': 'first',  # keep the first annotation
})
print(f"After averaging duplicates: {len(averaged)} rows")
print(averaged)
```

---

## Part 3: Type Conversion

Data loaded from CSV files often has columns stored as the wrong type (e.g., chromosome numbers as strings, p-values as text). Use `astype()` and `pd.to_numeric()` to fix this.

```python
# Messy clinical metadata
clinical = pd.DataFrame({
    'patient_id': ['P001', 'P002', 'P003', 'P004', 'P005'],
    'age': ['45', '52', '38', 'unknown', '61'],           # should be numeric
    'tumor_stage': ['II', 'III', 'I', 'II', 'IV'],
    'mutation_count': ['12', '8', '25', '3', 'N/A'],       # should be numeric
    'alive': ['True', 'True', 'False', 'True', 'False'],   # should be boolean
})

print("Original dtypes:")
print(clinical.dtypes)
print()
print(clinical)
```

```python
# Convert with error handling
clinical['age'] = pd.to_numeric(clinical['age'], errors='coerce')  # 'unknown' -> NaN
clinical['mutation_count'] = pd.to_numeric(clinical['mutation_count'], errors='coerce')
clinical['alive'] = clinical['alive'].map({'True': True, 'False': False})

print("Fixed dtypes:")
print(clinical.dtypes)
print()
print(clinical)
```

```python
# Categorical dtype -- saves memory and enables ordering
stage_order = ['I', 'II', 'III', 'IV']
clinical['tumor_stage'] = pd.Categorical(clinical['tumor_stage'],
                                          categories=stage_order,
                                          ordered=True)

print("Tumor stage is now categorical:")
print(clinical['tumor_stage'])
print(f"\nPatients with stage >= III:")
print(clinical[clinical['tumor_stage'] >= 'III'][['patient_id', 'tumor_stage']])
```

---

## Part 4: Reshaping Data -- Melt and Pivot

Biological data often needs to be converted between **wide** and **long** format.

```
WIDE (one row per gene):
  gene  | sample_1 | sample_2 | sample_3
  TP53  |   5.2    |   5.4    |   5.0

          melt()  <-->  pivot_table()

LONG (one row per measurement):
  gene  | sample   | expression
  TP53  | sample_1 |   5.2
  TP53  | sample_2 |   5.4
  TP53  | sample_3 |   5.0
```

- **Wide format** is natural for matrix computations (NumPy, heatmaps).
- **Long format** is required by most plotting libraries (seaborn, ggplot) and statistical models.

```python
# Start with a clean wide expression table
wide_expr = pd.DataFrame({
    'gene':     ['TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS'],
    'gene_type': ['TSG', 'TSG', 'Oncogene', 'Oncogene', 'Oncogene'],
    'ctrl_1':   [5.2, 3.8, 7.1, 2.4, 4.5],
    'ctrl_2':   [5.4, 3.6, 7.3, 2.6, 4.2],
    'treat_1':  [3.1, 5.9, 7.0, 4.8, 2.1],
    'treat_2':  [2.9, 6.2, 6.8, 5.1, 1.9],
})

print("Wide format:")
print(wide_expr)
```

```python
# Melt: wide -> long
long_expr = wide_expr.melt(
    id_vars=['gene', 'gene_type'],
    value_vars=['ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'],
    var_name='sample',
    value_name='expression',
)

print("Long format (melted):")
print(long_expr)
```

```python
# Add a condition column derived from the sample name
long_expr['condition'] = long_expr['sample'].str.split('_').str[0]

print("With condition column:")
print(long_expr.head(8))
```

```python
# Pivot: long -> wide
wide_again = long_expr.pivot_table(
    index=['gene', 'gene_type'],
    columns='sample',
    values='expression',
).reset_index()

# Flatten the MultiIndex column names
wide_again.columns.name = None

print("Pivoted back to wide:")
print(wide_again)
```

### Stack and Unstack

`stack()` and `unstack()` work with hierarchical (MultiIndex) indices. They are useful when your data already has a MultiIndex.

```python
# Create a MultiIndex DataFrame
multi_df = wide_expr.set_index(['gene', 'gene_type'])[['ctrl_1', 'ctrl_2', 'treat_1', 'treat_2']]
print("MultiIndex DataFrame:")
print(multi_df)

# Stack: columns -> rows (wide -> long)
stacked = multi_df.stack()
stacked.name = 'expression'
print("\nStacked (long):")
print(stacked.head(8))

# Unstack: rows -> columns (long -> wide)
unstacked = stacked.unstack()
print("\nUnstacked (wide again):")
print(unstacked)
```

---

## Part 5: String Operations in Pandas

Gene annotation files often contain free-text columns that need parsing. Pandas provides vectorized string methods through the `.str` accessor.

```python
# Messy gene annotation table
annotations = pd.DataFrame({
    'raw_id': [
        'gene_id "ENSG00000141510"; gene_name "TP53"; gene_biotype "protein_coding";',
        'gene_id "ENSG00000012048"; gene_name "BRCA1"; gene_biotype "protein_coding";',
        'gene_id "ENSG00000146648"; gene_name "EGFR"; gene_biotype "protein_coding";',
        'gene_id "ENSG00000136997"; gene_name "MYC"; gene_biotype "protein_coding";',
        'gene_id "ENSG00000133703"; gene_name "KRAS"; gene_biotype "protein_coding";',
        'gene_id "ENSG00000228630"; gene_name "HOTAIR"; gene_biotype "lncRNA";',
    ]
})

print("Raw GTF-like attribute column:")
print(annotations)
```

```python
# Extract Ensembl ID using str.extract with regex
annotations['ensembl_id'] = annotations['raw_id'].str.extract(r'gene_id "(ENSG\d+)"')

# Extract gene name
annotations['gene_name'] = annotations['raw_id'].str.extract(r'gene_name "(\w+)"')

# Extract biotype
annotations['biotype'] = annotations['raw_id'].str.extract(r'gene_biotype "([^"]+)"')

print("Parsed annotations:")
print(annotations[['ensembl_id', 'gene_name', 'biotype']])
```

```python
# Common string operations
genes = pd.Series(['brca1', 'TP53', 'egfr', 'Myc', 'KRAS'])

print("Original:   ", genes.tolist())
print("Upper:      ", genes.str.upper().tolist())
print("Lower:      ", genes.str.lower().tolist())
print("Starts 'T': ", genes.str.upper().str.startswith('T').tolist())
print("Contains 'A':", genes.str.upper().str.contains('A').tolist())
print("Length:     ", genes.str.len().tolist())
```

```python
# Splitting strings: extract chromosome arm from cytogenetic band
bands = pd.Series(['17p13.1', '17q21.31', '7p11.2', '8q24.21', '12p12.1'])

# Split on the arm letter (p or q)
chromosomes = bands.str.extract(r'(\d+)[pq]')[0]
arms = bands.str.extract(r'\d+([pq])')[0]

band_df = pd.DataFrame({
    'band': bands,
    'chromosome': chromosomes,
    'arm': arms,
})
print(band_df)
```

---

## Part 6: Apply, Transform, and Pipe

These three methods let you apply custom functions to DataFrames.

| Method | Input | Returns | Use case |
|--------|-------|---------|----------|
| `apply()` | Row or column | Anything (scalar, Series) | General-purpose transformation |
| `transform()` | Row or column | Same-shape output | Group-wise normalization |
| `pipe()` | Entire DataFrame | DataFrame | Chaining processing steps |

```python
# Apply a function to each row
expr_data = pd.DataFrame({
    'gene':     ['TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS'],
    'ctrl_1':   [5.2, 3.8, 7.1, 2.4, 4.5],
    'ctrl_2':   [5.4, 3.6, 7.3, 2.6, 4.2],
    'treat_1':  [3.1, 5.9, 7.0, 4.8, 2.1],
    'treat_2':  [2.9, 6.2, 6.8, 5.1, 1.9],
})

sample_cols = ['ctrl_1', 'ctrl_2', 'treat_1', 'treat_2']


def classify_expression(row):
    """Classify a gene based on its mean expression."""
    mean_val = row[sample_cols].mean()
    if mean_val > 6:
        return 'high'
    elif mean_val > 4:
        return 'medium'
    else:
        return 'low'


expr_data['category'] = expr_data.apply(classify_expression, axis=1)
print(expr_data[['gene', 'category']])
```

```python
# Apply a function column-wise: z-score normalize each sample
def zscore(col):
    return (col - col.mean()) / col.std()

normalized = expr_data[sample_cols].apply(zscore)
normalized.insert(0, 'gene', expr_data['gene'])

print("Z-score normalized:")
print(normalized.round(2))
```
