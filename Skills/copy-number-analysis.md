---
name: copy-number-analysis
description: DNA copy number analysis — read depth normalization, log2 ratio computation, CBS segmentation, CN state calling, genome-wide visualization, gene annotation with pyranges, tumor purity concepts
---

## When to Use

Use this skill when:
- Computing log2 copy number ratios from read depth arrays
- Segmenting a genome into regions of uniform copy number (CBS algorithm)
- Calling amplifications, deletions, and LOH events
- Annotating CNV segments with overlapping genes (pyranges)
- Working with TCGA CNV segment files or GISTIC2 outputs
- Accounting for tumor purity in somatic copy number analysis

## Quick Reference

| Concept | Formula | Notes |
|---|---|---|
| Log2 ratio | `log2(observed_depth / expected_depth)` | expected = genome-wide median |
| CN=0 (homo del) | log2 ≈ -∞ (observed ≈ 0) | add pseudocount before log |
| CN=1 (loss) | log2 ≈ -0.58 (= log2(0.75) if 50% purity) | purity-dependent |
| CN=2 (diploid) | log2 ≈ 0 | baseline |
| CN=3 (gain) | log2 ≈ +0.58 | |
| CN=4 (amp) | log2 ≈ +1.0 | |
| CN=6 (high amp) | log2 ≈ +1.58 | |
| Purity correction | observed_lr = purity × lr_tumor + (1-purity) × 0 | dilution by normal cells |
| GC normalization | rolling median normalization within windows | correct systematic depth bias |
| GISTIC2 threshold | amplification: log2 ≥ 0.1; deep deletion: log2 ≤ -1.2 | TCGA convention |

| Tool | Purpose |
|---|---|
| `scipy.ndimage.median_filter` | Rolling median normalization |
| `scipy.signal.find_peaks` | Detect breakpoints in log2 ratio |
| Custom CBS function | Binary segmentation on log2 array |
| `pyranges` | Gene-level annotation of segments |
| TCGA GDC / cBioPortal | Somatic CNV segment files |

## Key Patterns

**Pattern 1: Normalize read depth (GC-correction)**
```python
import numpy as np
from scipy.ndimage import median_filter

def normalize_depth(depth, window=50):
    """Rolling median normalization — corrects GC bias."""
    smoothed = median_filter(depth.astype(float), size=window, mode="reflect")
    # Normalize by genome-wide median
    return depth / np.median(depth)

log2_ratio = np.log2(normalize_depth(depth) + 1e-4)  # pseudocount for CN=0
```

**Pattern 2: Simulate read depth with CNV events**
```python
def simulate_read_depth(n_bins=500, base_depth=50, seed=42):
    rng = np.random.default_rng(seed)
    depth = rng.poisson(base_depth, n_bins).astype(float)
    return depth

chr1 = simulate_read_depth(500)
chr1[100:150] *= 3   # amplification (CN=6 in diploid context)
chr1[300:380] *= 0.5 # heterozygous deletion (CN=1)
```

**Pattern 3: CBS segmentation (Circular Binary Segmentation)**
```python
def cbs_segment(log2_ratios, min_segment=10, alpha=0.01):
    """Simplified CBS: recursively split at maximum-deviation breakpoint."""
    segments = []
    def _split(lo, hi):
        if hi - lo < min_segment * 2:
            segments.append((lo, hi, log2_ratios[lo:hi].mean()))
            return
        seg = log2_ratios[lo:hi]
        n = len(seg)
        # t-statistic at each potential split point
        cumsum = np.cumsum(seg)
        i_range = np.arange(min_segment, n - min_segment)
        t_stats = [(cumsum[i]/i - (cumsum[-1]-cumsum[i])/(n-i))**2 *
                   i*(n-i)/n for i in i_range]
        best = np.argmax(t_stats) + min_segment
        from scipy.stats import ttest_ind
        _, pval = ttest_ind(seg[:best], seg[best:], equal_var=False)
        if pval < alpha:
            _split(lo, lo + best)
            _split(lo + best, hi)
        else:
            segments.append((lo, hi, seg.mean()))
    _split(0, len(log2_ratios))
    return sorted(segments)
```

**Pattern 4: CN state calling from log2 ratio**
```python
def call_cn_state(mean_log2, ploidy=2, thresholds=None):
    """Map log2 ratio to integer copy number state."""
    if thresholds is None:
        # (log2_upper_bound, CN_state)
        thresholds = [(-1.5, 0), (-0.4, 1), (0.25, 2), (0.75, 3), (1.3, 4)]
    estimated_cn = ploidy * (2 ** mean_log2)
    for upper, cn in thresholds:
        if mean_log2 < upper:
            return cn
    return max(4, round(estimated_cn))

# Annotate segments
results = []
for s, e, mean_lr in segments:
    cn = call_cn_state(mean_lr)
    label = {0: "homozygous_del", 1: "loss", 2: "neutral",
             3: "gain", 4: "amplification"}.get(cn, f"CN{cn}")
    results.append({"start": s, "end": e, "log2": mean_lr, "CN": cn, "type": label})
```

**Pattern 5: Tumor purity correction**
```python
def purity_corrected_log2(observed_log2, purity):
    """
    Corrects observed bulk log2 ratio for tumor purity.
    observed_lr = purity * tumor_lr + (1 - purity) * 0 (diploid normal)
    """
    if purity <= 0 or purity > 1:
        raise ValueError("purity must be in (0, 1]")
    return observed_log2 / purity

# Example: true amplification (log2=1.58) at 60% purity
# observed = 0.6 * 1.58 + 0.4 * 0 = 0.95
observed = 0.95
true_lr = purity_corrected_log2(observed, purity=0.6)
print(f"True log2 ratio: {true_lr:.2f}")  # ~1.58
```

## Code Templates

**Template 1: Genome-wide CNV profile**
```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Build per-chromosome segments and plot genome-wide
chromosomes = {"chr1": chr1_log2, "chr2": chr2_log2, "chr3": chr3_log2}
fig, ax = plt.subplots(figsize=(14, 4))
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
offset = 0; xticks = []; xticklabels = []

for i, (name, lr) in enumerate(chromosomes.items()):
    segs = cbs_segment(lr)
    x = np.arange(len(lr)) + offset
    ax.scatter(x, lr, c=colors[i % 3], s=2, alpha=0.4, label=name)
    for s, e, m in segs:
        ax.hlines(m, s + offset, e + offset, colors="black", lw=2)
    xticks.append(offset + len(lr) / 2)
    xticklabels.append(name)
    offset += len(lr) + 20

ax.axhline(0, color="gray", lw=0.8, ls="--")
ax.axhline(np.log2(3/2), color="orange", lw=0.5, ls=":")   # gain threshold
ax.axhline(np.log2(1/2), color="red", lw=0.5, ls=":")      # loss threshold
ax.set_xticks(xticks); ax.set_xticklabels(xticklabels)
ax.set_xlabel("Genomic position (bins)"); ax.set_ylabel("Log₂ ratio")
ax.set_title("Genome-wide copy number profile")
plt.tight_layout()
```

**Template 2: Gene-level annotation with pyranges**
```python
import pyranges as pr
import pandas as pd

# Gene annotation BED (from Ensembl/UCSC)
genes = pr.read_bed("hg38_genes.bed")  # chr, start, end, name

# CNV segments as pyranges
seg_df = pd.DataFrame(segs, columns=["Chromosome", "Start", "End", "log2"])
seg_df["Chromosome"] = "chr1"  # set per chromosome
cnv = pr.PyRanges(seg_df)

# Find genes overlapping each segment
overlaps = cnv.join(genes)
print(overlaps.df[["Chromosome", "Start", "End", "log2", "Name_b"]])
```

**Template 3: Load TCGA CNV segment file**
```python
import pandas as pd

# TCGA segment file format: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean
seg = pd.read_csv("TCGA_BRCA_segments.txt", sep="\t")
seg.columns = ["sample", "chrom", "start", "end", "num_probes", "log2"]

# Filter for high-confidence segments
seg = seg[seg["num_probes"] >= 10]

# Call type
def classify_segment(log2):
    if log2 < -1.2: return "deep_deletion"
    if log2 < -0.4: return "shallow_deletion"
    if log2 > 1.0:  return "high_amplification"
    if log2 > 0.1:  return "amplification"
    return "neutral"

seg["type"] = seg["log2"].apply(classify_segment)
print(seg["type"].value_counts())
```

## Common Pitfalls

- **Log2 of zero depth:** always add pseudocount (`+ 1e-4` or `+ 0.01`) before `log2`; CN=0 gives -inf otherwise
- **GC bias:** rolling-median normalization corrects local depth biases; skipping it creates wave artifacts across chromosomes
- **Purity dilution:** a 50% purity tumor with CN=4 shows log2 ≈ 0.5 in bulk — always interpret log2 relative to expected purity
- **Min segment length:** CBS `min_segment` must be ≥ 5–10 bins; too short produces noisy micro-segments
- **GISTIC vs raw thresholds:** TCGA GISTIC uses log2 ≥ 0.1 (any gain) vs. log2 ≥ 1.0 (high amplification); know which convention you're applying
- **Bin size matters:** 50 bp bins → too noisy for segmentation; 1–10 kb bins give stable log2 ratios for most applications

## Related Skills

- `ngs-variant-calling` — upstream BAM file generation; read depth computed from aligned reads
- `gwas-population-genetics` — germline CNVs use similar calling but with expected diploid baseline
- `rnaseq-metagenomics` — copy number affects RNA expression (gene dosage)
- `tf-footprinting-atac` — ATAC-seq workflow shares read depth QC patterns
