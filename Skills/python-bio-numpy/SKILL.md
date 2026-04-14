---
name: python-bio-numpy
description: NumPy for bioinformatics — vectorized operations, broadcasting for expression matrices, PWM construction, sliding-window GC with cumsum, and normalization recipes
tool_type: python
primary_tool: NumPy
---

# NumPy for Bioinformatics

## Key Concepts

- **Vectorized operations** run in compiled C — `arr * 2` is much faster than `[x * 2 for x in arr]`.
- **Broadcasting aligns from the right.** `(4,3) - (3,)` works (broadcasts across rows). `(4,3) - (4,)` fails — reshape to `(4,1)` first.
- **Views vs copies:** `arr[2:5]` is a view (modifies original). `arr[[0,2,4]]` and `arr[arr > 0]` return copies.
- **`axis` in aggregations:** `axis=0` collapses rows (per-column stats); `axis=1` collapses columns (per-row stats). In genes x samples: `axis=0` = per-sample, `axis=1` = per-gene.

## Bio Recipes

### RPKM Normalization

```python
def compute_rpkm(counts, lengths):
    rpm = counts * (1_000_000 / np.sum(counts))
    return rpm * (1_000 / lengths)
```

### CPM via Broadcasting

```python
expr = np.array([[120, 135, 128], [45, 50, 48], [300, 280, 310]], dtype=float)
sample_totals = expr.sum(axis=0, keepdims=True)   # shape (1, 3)
cpm = expr / sample_totals * 1_000_000
```

### Z-score Per Gene (row) Across Samples

```python
gene_means = expr.mean(axis=1, keepdims=True)
gene_stds  = expr.std(axis=1, keepdims=True)
z_scores = (expr - gene_means) / gene_stds
```

### Position Weight Matrix (PWM)

```python
def build_pwm(sequences):
    """Build PWM from aligned sequences. Returns (4, seq_len) frequency matrix."""
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    seq_len = len(sequences[0])
    counts = np.zeros((4, seq_len))
    for seq in sequences:
        for i, nuc in enumerate(seq):
            counts[mapping[nuc], i] += 1
    return counts / len(sequences)

def score_sequence(pwm, sequence):
    """Score sequence against PWM (sum of log-odds vs uniform background)."""
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    return sum(np.log2(pwm[mapping[nuc], i] / 0.25) if pwm[mapping[nuc], i] > 0 else -10
               for i, nuc in enumerate(sequence))
```

### Sliding Window GC with Cumulative Sum (O(n))

```python
def sliding_gc(sequence, window=50):
    """O(n) sliding window GC content using cumsum trick."""
    gc_binary = np.array([1 if n in 'GC' else 0 for n in sequence])
    cumsum = np.insert(np.cumsum(gc_binary), 0, 0)
    window_sums = cumsum[window:] - cumsum[:-window]
    return window_sums / window * 100
```

### Log-transform Expression Counts

```python
log_counts = np.log2(counts + 1)  # +1 pseudocount to avoid log(0)
```

### Correlation Between Gene Profiles

```python
corr = np.corrcoef(gene_a, gene_b)[0, 1]  # Pearson r
```
