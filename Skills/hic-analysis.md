---
name: hic-analysis
description: Hi-C data analysis — loading cooler files, contact decay, A/B compartments, TAD boundaries, pileup plots using cooler and cooltools
---

## When to Use

Use this skill when:
- Analyzing 3D genome organization from Hi-C or Micro-C experiments
- Working with contact matrices in `.cool` or `.mcool` format
- Computing contact decay (expected) curves to normalize distance effects
- Detecting A/B compartments via eigenvector decomposition
- Identifying TAD boundaries using insulation score
- Creating pileup/aggregate plots around genomic features (CTCF, loop anchors)

## Quick Reference

| Tool / Method | Purpose | Notes |
|---|---|---|
| `cooler.Cooler(path)` | Open a .cool file | Use `::resolutionN` suffix for .mcool |
| `clr.info` | Inspect resolution, genome, shape | dict with metadata |
| `clr.chromnames` | List chromosome names | Check `chr1` vs `1` format |
| `clr.bins()` | Bin table (chrom, start, end, weight) | `weight` = balancing factor |
| `clr.matrix(balance=True).fetch(region)` | Fetch balanced submatrix | Returns numpy ndarray |
| `cooltools.expected_cis()` | Contact decay by genomic distance | Returns DataFrame |
| `cooltools.eigs_cis()` | Eigenvector decomposition (A/B) | E1 sign needs GC-content flip |
| `cooltools.saddle()` | Saddle plot (compartment strength) | Needs expected + eigenvectors |
| `cooltools.insulation()` | Insulation score + TAD boundaries | window parameter in bp |
| `cooltools.pileup()` | Aggregate/snipping around features | Needs snipper + features BED |
| Common resolutions | 1 kb, 5 kb, 10 kb, 25 kb, 50 kb, 100 kb | TADs: 25–50 kb; compartments: 100 kb |

## Key Patterns

**Pattern 1: Load a cooler file and inspect it**
```python
import cooler

clr = cooler.Cooler("path/to/file.cool")
print(clr.info)           # resolution, genome assembly, shape
print(clr.chromnames)     # ['chr1', 'chr2', ...] or ['1', '2', ...]
print(clr.bins().head())  # chrom, start, end, weight columns
```

**Pattern 2: Fetch a contact matrix for a region**
```python
region = "chr1:0-5000000"
matrix = clr.matrix(balance=True).fetch(region)
import numpy as np
matrix_log = np.log1p(matrix)  # log-transform for visualization
```

**Pattern 3: Compute contact decay (expected) curve**
```python
import cooltools

view_df = cooler.util.make_chromarms(clr.chromsizes, mid_point_flag=False)
view_df.columns = ["chrom", "start", "end", "name"]
expected = cooltools.expected_cis(clr, view_df=view_df)
```

**Pattern 4: A/B compartment eigenvectors**
```python
eigvals, eigvecs = cooltools.eigs_cis(clr, view_df=view_df, n_eigs=3)
# eigvecs has columns E1, E2, E3 per genomic bin
# Flip E1 sign so positive = A compartment (GC-rich)
```

**Pattern 5: Insulation score and TAD boundaries**
```python
insulation = cooltools.insulation(clr, [200_000, 400_000], view_df=view_df)
# boundary_strength_200000 column indicates TAD boundaries
boundaries = insulation[insulation["is_boundary_200000"]]
```

## Code Templates

**Template 1: Contact matrix heatmap**
```python
import cooler, numpy as np, matplotlib.pyplot as plt

clr = cooler.Cooler("sample.cool")
region = "chr1:0-10000000"
mat = clr.matrix(balance=True).fetch(region)
mat_log = np.log1p(np.nan_to_num(mat))

fig, ax = plt.subplots(figsize=(8, 8))
ax.matshow(mat_log, cmap="YlOrRd", origin="upper")
ax.set_title(f"Hi-C contact matrix — {region}")
plt.tight_layout()
plt.show()
```

**Template 2: Contact decay curve**
```python
import cooltools, matplotlib.pyplot as plt

expected = cooltools.expected_cis(clr, view_df=view_df)
fig, ax = plt.subplots()
ax.loglog(expected["dist"] * clr.binsize,
          expected["balanced.avg"],
          color="steelblue")
ax.set_xlabel("Genomic distance (bp)")
ax.set_ylabel("Average contact frequency")
ax.set_title("Contact decay (P(s) curve)")
plt.tight_layout()
plt.show()
```

**Template 3: Saddle plot**
```python
import cooltools

eigvals, eigvecs = cooltools.eigs_cis(clr, view_df=view_df, n_eigs=3)
q_lo, q_hi = 0.025, 0.975
saddle = cooltools.saddle(
    clr,
    expected,
    eigvecs[["chrom", "start", "end", "E1"]],
    view_df=view_df,
    n_bins=50,
    qrange=(q_lo, q_hi)
)
```

## Common Pitfalls

- **`chr` prefix mismatch:** cooler stores chroms exactly as in the assembly; always check `clr.chromnames` before fetching — `"chr1:0-1000000"` vs `"1:0-1000000"` both work but must match
- **NaN after balancing:** bins with weight=NaN are masked; `np.nan_to_num(matrix, nan=0)` before visualization
- **E1 sign ambiguity:** eigenvector sign is arbitrary; flip so active compartment (high GC) = positive E1
- **Resolution vs. window:** insulation window must be larger than the bin size; `window ≥ 5 × binsize` is a good rule
- **mcool files:** use `cooler.Cooler("file.mcool::resolutions/10000")` to open a specific resolution
- **Memory:** fetching large chromosomes at high resolution into dense arrays can exceed RAM — use `balance=False` or restrict to a region

## Related Skills

- `ngs-variant-calling` — upstream FASTQ → BAM pipeline that produces Hi-C alignments
- `motif-discovery` — used to define genomic features for pileup analysis (CTCF motif sites)
- `structural-bioinformatics` — complementary 3D genome perspective at the protein-DNA level
