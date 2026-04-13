---
name: python-bio-numpy
description: "NumPy for bioinformatics: arrays, vectorized operations, broadcasting, linear algebra, and statistical functions for biological data. Use when performing numerical computation on biological datasets."
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/16_NumPy_and_Pandas/01_numpy.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Numpy

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/16_NumPy_and_Pandas/01_numpy.ipynb`*


Split from `01_numpy_and_pandas.ipynb` to keep this topic self-contained.

**Navigation:** [Topic overview](./01_numpy_and_pandas.ipynb) · [Next: Pandas](./02_pandas.ipynb)

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. Vectorized operations are element-wise**
`arr * 2` doubles every element; `arr1 + arr2` adds matching elements. No loop needed. This is much faster than `[x * 2 for x in arr]` because the loop runs in compiled C code.

**2. Broadcasting aligns from the right**
When array shapes don't match, NumPy aligns from the trailing dimension and broadcasts size-1 dimensions outward. A `(4, 3)` minus a `(3,)` works (broadcasts across rows). A `(4, 3)` minus a `(4,)` fails — reshape to `(4, 1)` first.

**3. Views vs copies**
`arr[2:5]` is a *view* — modifying it modifies `arr`. Fancy indexing `arr[[0, 2, 4]]` and boolean indexing `arr[arr > 0]` both return *copies*. Use `.copy()` to be explicit.

**4. `axis` in aggregations**
`arr.mean(axis=0)` collapses rows → per-column mean. `arr.mean(axis=1)` collapses columns → per-row mean. In a genes × samples matrix: `axis=0` gives per-sample statistics; `axis=1` gives per-gene statistics.

```python
import numpy as np
import pandas as pd

# Reproducibility
np.random.seed(42)
```python

---

## Part 1: NumPy Fundamentals

NumPy's core object is the **ndarray** -- an n-dimensional array of elements that all share the same type. Unlike Python lists, arrays support vectorized operations: you write the operation once and it applies to every element without an explicit loop.

### 1.1 Creating Arrays

```python
# From a Python list
gc_values = np.array([0.42, 0.51, 0.48, 0.55, 0.44])
print("GC content per gene:", gc_values)
print("Type:", type(gc_values))
print("Dtype:", gc_values.dtype)
```python

```python
# Convenience constructors
zeros = np.zeros(10)                  # 10 zeros
ones = np.ones((3, 4))               # 3x4 matrix of ones
positions = np.arange(0, 1000, 100)  # like range(), but returns array
fractions = np.linspace(0, 1, 5)     # 5 evenly spaced values in [0, 1]
empty_matrix = np.full((3, 3), np.nan)  # 3x3 filled with NaN

print("zeros:", zeros)
print("ones shape:", ones.shape)
print("positions:", positions)
print("fractions:", fractions)
print("empty_matrix:\n", empty_matrix)
```python

```python
# 2D array -- a gene expression matrix (rows=genes, columns=samples)
expression = np.array([
    [120, 135, 128],   # Gene A across 3 samples
    [ 45,  50,  48],   # Gene B
    [300, 280, 310],   # Gene C
])

print("Shape:", expression.shape)    # (rows, cols)
print("Dimensions:", expression.ndim)
print("Total elements:", expression.size)
print("Matrix:\n", expression)
```python

### 1.2 Indexing and Slicing

NumPy arrays use comma-separated indices for multiple dimensions: `arr[row, col]`. Slices, integer arrays, and boolean arrays all work as indices.

```python
# 1D indexing
reads = np.array([63, 0, 17, 250, 89, 42, 130])

print("First element:", reads[0])
print("Last element:", reads[-1])
print("Slice [2:5]:", reads[2:5])
print("Every other:", reads[::2])
```python

```python
# 2D indexing on the expression matrix
print("Row 0 (Gene A):", expression[0])        # entire row
print("Column 1 (Sample 2):", expression[:, 1])  # entire column
print("Gene C, Sample 3:", expression[2, 2])    # single element
print("Genes A-B, Samples 1-2:\n", expression[:2, :2])  # submatrix
```python

```python
# Boolean indexing -- filter without loops
print("Reads > 50:", reads[reads > 50])
print("Non-zero reads:", reads[reads != 0])

# Boolean mask
high_expr_mask = reads > 100
print("Mask:", high_expr_mask)
print("High-count reads:", reads[high_expr_mask])
```python

### 1.3 Vectorized Operations

Every arithmetic operation on NumPy arrays is **vectorized** -- it applies element-by-element with no Python loop. This is both faster and more readable.

```python
# Gene coordinates: start and stop positions
starts = np.array([11869, 14404, 17369, 29554, 30366, 34554, 52473])
stops  = np.array([14409, 29570, 17436, 31109, 30503, 36081, 53312])

# Calculate gene lengths -- vectorized, no loop needed
lengths = stops - starts + 1
print("Gene lengths:", lengths)

# Log-transform expression counts
counts = np.array([63, 119, 17, 250, 89, 42, 130], dtype=float)
log_counts = np.log2(counts + 1)  # +1 to avoid log(0)
print("Log2 counts:", np.round(log_counts, 2))
```python

```python
# RPKM normalization (Reads Per Kilobase per Million mapped reads)
def compute_rpkm(counts, lengths):
    """Normalize raw read counts to RPKM."""
    rpm = counts * (1_000_000 / np.sum(counts))   # per million
    rpkm = rpm * (1_000 / lengths)                 # per kilobase
    return rpkm

rpkm = compute_rpkm(counts, lengths)
print("RPKM values:", np.round(rpkm, 2))
```python

### 1.4 Broadcasting

When arrays have different but compatible shapes, NumPy **broadcasts** the smaller array to match. This lets you write concise code for operations across rows or columns.

```python
# Expression matrix: 4 genes x 3 samples
expr = np.array([
    [120, 135, 128],
    [ 45,  50,  48],
    [300, 280, 310],
    [ 10,  12,   9],
], dtype=float)

# Normalize each sample (column) to its total -- broadcasting divides
# a (4,3) matrix by a (1,3) row vector
sample_totals = expr.sum(axis=0, keepdims=True)  # shape (1, 3)
cpm = expr / sample_totals * 1_000_000

print("Sample totals:", sample_totals.flatten())
print("CPM (counts per million):\n", np.round(cpm, 1))
```python

```python
# Z-score each gene (row) across samples
gene_means = expr.mean(axis=1, keepdims=True)  # shape (4, 1)
gene_stds  = expr.std(axis=1, keepdims=True)

z_scores = (expr - gene_means) / gene_stds
print("Z-scores:\n", np.round(z_scores, 2))
```python

### 1.5 Statistical Functions

```python
# Simulated gene expression for 100 genes
np.random.seed(42)
data = np.random.lognormal(mean=5, sigma=1, size=100)

print(f"Mean:   {np.mean(data):.2f}")
print(f"Median: {np.median(data):.2f}")
print(f"Std:    {np.std(data):.2f}")
print(f"Min:    {np.min(data):.2f}")
print(f"Max:    {np.max(data):.2f}")
print(f"25th percentile: {np.percentile(data, 25):.2f}")
print(f"75th percentile: {np.percentile(data, 75):.2f}")
```python

```python
# axis parameter: compute along rows or columns
print("Per-gene mean (across samples):", np.round(expr.mean(axis=1), 1))
print("Per-sample mean (across genes):", np.round(expr.mean(axis=0), 1))
print("Grand mean:", np.round(expr.mean(), 1))
```python

### 1.6 Reshape and Flatten

```python
# A 2x6 matrix reshaped into 3x4, 6x2, or flattened
matrix = np.arange(12).reshape(2, 6)
print("Original (2x6):\n", matrix)
print("Reshaped (3x4):\n", matrix.reshape(3, 4))
print("Flattened:", matrix.flatten())
```python

### 1.7 Linear Algebra Basics

```python
# Dot product: correlation-like score between two gene profiles
gene_a = np.array([1.2, 3.4, 2.1, 5.0])
gene_b = np.array([1.1, 3.5, 2.0, 4.8])

dot = np.dot(gene_a, gene_b)
print(f"Dot product: {dot:.2f}")

# Correlation coefficient
corr = np.corrcoef(gene_a, gene_b)
print(f"Pearson r: {corr[0, 1]:.4f}")
```python

```python
# Matrix multiplication: transform a count matrix
# rows = genes, cols = samples
count_matrix = np.array([[10, 20], [30, 40], [50, 60]])

# A 2x2 normalization/rotation matrix
transform = np.array([[0.5, 0.5], [0.5, -0.5]])

result = count_matrix @ transform  # or np.matmul(count_matrix, transform)
print("Original (3x2):\n", count_matrix)
print("After transform (3x2):\n", result)
```python

### 1.8 Biological Application: Position Weight Matrix

A PWM represents the frequency of each nucleotide at each position in a set of aligned binding-site sequences. It is stored as a 4 x L NumPy array (rows = A, T, G, C; columns = positions).

```python
def build_pwm(sequences):
    """Build a Position Weight Matrix from aligned sequences."""
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    seq_len = len(sequences[0])
    counts = np.zeros((4, seq_len))
    for seq in sequences:
        for i, nuc in enumerate(seq):
            counts[mapping[nuc], i] += 1
    return counts / len(sequences)


def score_sequence(pwm, sequence):
    """Score a sequence against a PWM (sum of log-odds)."""
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    score = 0.0
    for i, nuc in enumerate(sequence):
        freq = pwm[mapping[nuc], i]
        score += np.log2(freq / 0.25) if freq > 0 else -10
    return score


# TATA-box motif examples
tata_seqs = [
    "TATAAAG",
    "TATAAAT",
    "TATAAAA",
    "TATAAAT",
    "TATAAAG",
    "TATAAAC",
]

pwm = build_pwm(tata_seqs)
print("PWM (rows: A, T, G, C):")
for i, nuc in enumerate("ATGC"):
    print(f"  {nuc}: {np.round(pwm[i], 2)}")

# Score a candidate sequence
print(f"\nScore 'TATAAAA': {score_sequence(pwm, 'TATAAAA'):.2f}")
print(f"Score 'GCGCGCG': {score_sequence(pwm, 'GCGCGCG'):.2f}")
```python

### 1.9 Biological Application: Sliding Window GC Content

```python
def sliding_gc(sequence, window=50):
    """Calculate GC content in a sliding window using NumPy cumsum."""
    gc_binary = np.array([1 if n in 'GC' else 0 for n in sequence])
    cumsum = np.insert(np.cumsum(gc_binary), 0, 0)
    window_sums = cumsum[window:] - cumsum[:-window]
    return window_sums / window * 100


# A synthetic sequence with variable GC content
np.random.seed(7)
low_gc = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 100, p=[0.35, 0.35, 0.15, 0.15]))
high_gc = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 100, p=[0.15, 0.15, 0.35, 0.35]))
seq = low_gc + high_gc

gc = sliding_gc(seq, window=20)
print(f"Sequence length: {len(seq)}")
print(f"GC values computed: {len(gc)}")
print(f"Mean GC in first half:  {gc[:90].mean():.1f}%")
print(f"Mean GC in second half: {gc[90:].mean():.1f}%")
```python

---

## Part 2: Pandas Fundamentals

Pandas extends NumPy with **labeled axes** (row index and column names) and rich I/O. Its two main structures are:

| Structure | Description |
|-----------|-------------|
| **Series** | A 1D labeled array (like a column in a spreadsheet) |
| **DataFrame** | A 2D labeled table (rows and columns) |

```python
# From a list (default integer index)
gc_series = pd.Series([0.42, 0.51, 0.48, 0.55, 0.44],
                      name='gc_content')
print(gc_series)
print("\nType:", type(gc_series))
```python

```python
# From a dictionary -- keys become the index
complement = pd.Series({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'})
print(complement)
print("\nComplement of G:", complement['G'])
print("Values:", complement.values)
```python

```python
# Vectorized operations on Series
gene_lengths = pd.Series({'BRCA1': 7088, 'TP53': 2512, 'EGFR': 5616, 'MYC': 2357})
print("Lengths in kilobases:")
print(gene_lengths / 1000)
```python

### 2.2 Creating DataFrames

```python
# From a dictionary -- each key becomes a column
genes_df = pd.DataFrame({
    'gene':       ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
    'chromosome':  ['17',    '17',   '7',    '8',   '12'],
    'length_bp':  [7088,    2512,   5616,   2357,  5764],
    'gc_content': [0.423,   0.512,  0.487,  0.551, 0.448],
    'biotype':    ['protein_coding'] * 5,
})

print(genes_df)
print(f"\nShape: {genes_df.shape}")
print(f"Columns: {genes_df.columns.tolist()}")
```python

```python
# Quick overview methods
print("--- dtypes ---")
print(genes_df.dtypes)
print("\n--- describe ---")
print(genes_df.describe())
print("\n--- info ---")
genes_df.info()
```python

### 2.3 Indexing with loc and iloc

| Accessor | Lookup by | Example |
|----------|-----------|--------|
| `.loc[]` | **label** (index/column name) | `df.loc[0, 'gene']` |
| `.iloc[]` | **integer position** | `df.iloc[0, 0]` |

```python
# Column access
print("Single column (Series):")
print(genes_df['gene'])

print("\nMultiple columns (DataFrame):")
print(genes_df[['gene', 'gc_content']])
```python

```python
# loc -- label-based
print("Row 0:")
print(genes_df.loc[0])

print(f"\ngenes_df.loc[0, 'gene'] = {genes_df.loc[0, 'gene']}")

print("\nRows 1-3, columns gene and length_bp:")
print(genes_df.loc[1:3, ['gene', 'length_bp']])
```python

```python
# iloc -- integer position
print("First 2 rows, first 3 columns:")
print(genes_df.iloc[:2, :3])

print(f"\nLast row, last column: {genes_df.iloc[-1, -1]}")
```python

### 2.4 Filtering Rows

## Common Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
