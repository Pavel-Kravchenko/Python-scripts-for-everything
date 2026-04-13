---
name: bio-applied-copy-number-analysis
description: "**By the end you will be able to:** - Understand copy number variation (CNV) concepts - Compute read depth across genomic windows - Normalize read depth and detect segments - Implement Circular Binary"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/21_Copy_Number_Analysis/21_copy_number_analysis.ipynb"
primary_tool: Pandas
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 21: DNA Copy Number Analysis

*Source: Course notebook `Tier_3_Applied_Bioinformatics/21_Copy_Number_Analysis/21_copy_number_analysis.ipynb`*

**Tier 3 — Applied Bioinformatics | Module 21**
Prerequisites: Module 01 (NGS Fundamentals), Module 02 (Variant Calling and SNP Analysis).

---

**By the end you will be able to:**
- Understand copy number variation (CNV) concepts
- Compute read depth across genomic windows
- Normalize read depth and detect segments
- Implement Circular Binary Segmentation (CBS) concept
- Annotate genes in CNV segments
- Visualize genome-wide copy number profiles

**Attribution:** *Patterns inspired by NGSchool 2023 DNA copy number practical. Uses simulated read depth data.*

## Background: Copy Number Variation

**Copy number variation (CNV)** refers to regions of the genome where the number of copies differs from the diploid expectation of 2.

| Type | Copy number | Log2 ratio | Effect |
|---|---|---|---|
| Homozygous deletion | 0 | −∞ (→ −∞ after clipping) | Gene loss; common in tumor suppressors |
| Heterozygous deletion | 1 | log2(1/2) = −1.0 | Loss of one allele (LOH if other allele mutated) |
| Normal diploid | 2 | log2(2/2) = 0.0 | Baseline |
| Single-copy gain | 3 | log2(3/2) ≈ +0.585 | Increased gene dosage |
| Amplification (high) | 4+ | log2(4/2) = +1.0 etc. | Common in oncogenes (ERBB2, MYCN) |

**Log2 ratio:**
The fundamental measurement is the log2 ratio of observed depth to expected depth:

$$\text{log2 ratio} = \log_2\left(\frac{\text{depth}_{\text{tumor}}}{\text{depth}_{\text{normal}}}\right) = \log_2\left(\frac{\text{CN}}{2}\right)$$

**Detection from sequencing:**
1. Compute read depth (coverage) in fixed-width genomic windows (typically 1 kb – 10 kb)
2. Normalize for GC bias (GC content correlates with sequencing efficiency), mappability (low-complexity regions get fewer reads), and global depth scaling
3. Segment the normalized log2-ratio profile into regions of uniform copy number using CBS
4. Call integer copy number states from segment means

**Circular Binary Segmentation (CBS; Olshen et al. 2004):** recursively splits chromosomes at the point that maximizes the absolute difference in mean log2 ratio between the two resulting segments. Splitting continues until no significant breakpoints remain. The significance threshold uses a permutation test.

**Tumour purity:** In cancer samples, tumour DNA is mixed with normal stromal cells. At 50% purity, the effective copy number is a mixture of tumour CN and normal CN=2: effective_CN = 0.5×1 + 0.5×2 = 1.5. The log2 ratio is log2(1.5/2) = log2(0.75) ≈ −0.415, not −1.0. This attenuated signal is why purity must be accounted for when calling CN states. Tools like CNVkit, Sequenza, and PURPLE jointly estimate purity and ploidy from the data.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import find_peaks
import warnings
warnings.filterwarnings("ignore")
plt.rcParams.update({"figure.dpi":120,"axes.spines.top":False,"axes.spines.right":False})
print("Libraries loaded.")
```python

```python
rng = np.random.default_rng(42)

# Simulate 3 chromosomes with known CNV events
def simulate_read_depth(n_bins=500, base_depth=50, seed=42):
    rng = np.random.default_rng(seed)
    depth = rng.poisson(base_depth, n_bins).astype(float)
    return depth

# Chromosome 1: two amplifications and one deletion
chr1 = simulate_read_depth(500)
chr1[100:150] *= 3   # amplification (CN=6)
chr1[300:380] *= 0.5 # heterozygous deletion (CN=1)

# Chromosome 2: large deletion
chr2 = simulate_read_depth(400)
chr2[200:350] = rng.poisson(25, 150)  # hemizygous deletion

# Chromosome 3: normal
chr3 = simulate_read_depth(300)

print(f"chr1: {len(chr1)} bins, chr2: {len(chr2)} bins, chr3: {len(chr3)} bins")
print(f"Chr1 mean depth (amplified region 100-150): {chr1[100:150].mean():.1f}")
print(f"Chr1 mean depth (normal region 0-100): {chr1[:100].mean():.1f}")
```python

```python
def normalize_depth(depth, window=50):
    """
    Normalize read depth by dividing each bin by the local rolling median.
    The rolling median captures regional GC-content bias and mappability effects.
    Dividing by the local (not global) median corrects for systematic regional biases
    while preserving CNV signal. A final global median rescaling centers the log2 ratio at 0.
    """
    from scipy.ndimage import median_filter
    local_median = median_filter(depth.astype(float), size=window, mode="reflect")
    # Normalize by local median to correct regional GC/mappability bias
    norm = depth / np.where(local_median > 0.1, local_median, 0.1)
    # Rescale so that diploid regions have mean ~1.0
    norm = norm / np.median(norm)
    return norm

chr1_norm = normalize_depth(chr1)
chr2_norm = normalize_depth(chr2)
chr3_norm = normalize_depth(chr3)

# Log2 ratio (log2(observed/expected))
# Expected = 1.0 (diploid)
chr1_log2 = np.log2(np.clip(chr1_norm, 0.01, None))
chr2_log2 = np.log2(np.clip(chr2_norm, 0.01, None))
chr3_log2 = np.log2(np.clip(chr3_norm, 0.01, None))

fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=False)
for i, (name, lr) in enumerate(zip(["chr1","chr2","chr3"],
                                    [chr1_log2, chr2_log2, chr3_log2])):
    axes[i].plot(lr, lw=0.5, alpha=0.6, color="gray")
    axes[i].axhline(0, color="black", lw=0.8, ls="--", label="CN=2 (diploid)")
    axes[i].axhline(np.log2(3/2), color="red", lw=0.5, ls=":", alpha=0.5, label="CN=3")
    axes[i].axhline(np.log2(1/2), color="blue", lw=0.5, ls=":", alpha=0.5, label="CN=1")
    axes[i].set_ylabel("log2 ratio"); axes[i].set_title(name)
    if i == 0: axes[i].legend(frameon=False, ncol=3, fontsize=8)
plt.tight_layout(); plt.show()
```python

```python
def cbs_segment(log2_ratios, min_segment=10, alpha=0.01):
    """
    Simplified CBS: split array at point of maximum absolute mean difference.
    Returns list of (start, end, mean_log2) tuples.
    """
    n = len(log2_ratios)
    segments = []

    def _split(lo, hi):
        if hi - lo < min_segment * 2:
            segments.append((lo, hi, log2_ratios[lo:hi].mean()))
            return
        # Find split point that maximizes |mean_left - mean_right|
        best_delta, best_k = 0, None
        for k in range(lo + min_segment, hi - min_segment):
            left_mean  = log2_ratios[lo:k].mean()
            right_mean = log2_ratios[k:hi].mean()
            delta = abs(left_mean - right_mean)
            if delta > best_delta:
                best_delta, best_k = delta, k
        if best_k is None or best_delta < 0.3:
            segments.append((lo, hi, log2_ratios[lo:hi].mean()))
            return
        _split(lo, best_k)
        _split(best_k, hi)

    _split(0, n)
    return sorted(segments)

segs = cbs_segment(chr1_log2)
print(f"Segments detected in chr1: {len(segs)}")
for s, e, m in segs:
    cn_est = 2 * 2**m
    print(f"  bins {s:3d}–{e:3d}: log2={m:+.3f}, est. CN={cn_est:.1f}")
```python

```python
# Combine chromosomes for genome-wide view
all_lr = np.concatenate([chr1_log2, chr2_log2, chr3_log2])
offsets = [0, len(chr1_log2), len(chr1_log2)+len(chr2_log2)]
chr_names = ["chr1","chr2","chr3"]

fig, ax = plt.subplots(figsize=(14, 4))
colors = ["#1f77b4","#ff7f0e","#2ca02c"]
xticks, xticklabels = [], []
for i, (name, offset, lr) in enumerate(zip(chr_names, offsets,
                                            [chr1_log2, chr2_log2, chr3_log2])):
    x = np.arange(len(lr)) + offset
    ax.scatter(x, lr, s=1, alpha=0.4, color=colors[i])
    # Segmentation
    segs = cbs_segment(lr)
    for s, e, m in segs:
        ax.hlines(m, s+offset, e+offset, colors="red", lw=2, alpha=0.9)
    xticks.append(offset + len(lr)//2)
    xticklabels.append(name)

ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.7)
ax.axhline(np.log2(3/2), color="orange", lw=0.5, ls=":", alpha=0.6, label="Gain (CN=3)")
ax.axhline(np.log2(1/2), color="dodgerblue", lw=0.5, ls=":", alpha=0.6, label="Loss (CN=1)")
ax.set_xticks(xticks); ax.set_xticklabels(xticklabels)
ax.set_ylabel("log2(depth/expected)")
ax.set_title("Genome-wide copy number profile (red = CBS segments)")
ax.legend(frameon=False)
plt.tight_layout(); plt.show()
```python

```python
def call_cn_state(mean_log2, ploidy=2, thresholds=None):
    """Map log2 ratio to integer copy number state using log2 cutoff thresholds.

    Parameters
    ----------
    mean_log2 : float
        Mean log2 ratio for the segment (log2(observed/expected)).
    ploidy : int
        Expected ploidy of the sample (2 = diploid).
    thresholds : dict, optional
        Dict mapping integer CN state to the *lower* log2 boundary for that state.
        E.g. {0: -inf, 1: -1.5, 2: -0.4, 3: 0.25, 4: 0.75, 5: 1.3}.
        If None, uses defaults appropriate for a diploid sample.

    Returns
    -------
    int : called copy number state
    """
    if thresholds is None:
        # Lower log2 boundary for each CN state (diploid, ploidy=2)
        # CN=0: log2(0/2) → -inf; CN=1: log2(1/2) = -1.0; CN=2: log2(2/2) = 0.0
        # CN=3: log2(3/2) ≈ +0.585; CN=4: log2(4/2) = +1.0
        # Boundaries are midpoints between adjacent theoretical values:
        # (-inf, -1.5) → CN=0; (-1.5, -0.4) → CN=1; (-0.4, +0.35) → CN=2
        # (+0.35, +0.80) → CN=3; (+0.80, +1.25) → CN=4; (+1.25, +inf) → CN≥5
        thresholds = {
            0: float('-inf'),  # lower boundary for CN=0
            1: -1.5,           # lower boundary for CN=1
            2: -0.4,           # lower boundary for CN=2 (includes diploid baseline)
            3:  0.35,          # lower boundary for CN=3
            4:  0.80,          # lower boundary for CN=4
            5:  1.25,          # lower boundary for CN≥5
        }

    # Walk from highest CN state down; return first state whose lower boundary ≤ mean_log2
    for cn in sorted(thresholds.keys(), reverse=True):
        if mean_log2 >= thresholds[cn]:
            return cn
    return 0  # below the lowest threshold → homozygous deletion


segs_chr1 = cbs_segment(chr1_log2)
results = []
for s, e, m in segs_chr1:
    cn = call_cn_state(m)
    # Also compute the naive rounded estimate for comparison
    cn_naive = max(0, round(2 * (2 ** m)))
    results.append({
        "chrom": "chr1", "start_bin": s, "end_bin": e,
        "n_bins": e - s, "log2_ratio": round(m, 3), "CN_call": cn, "CN_naive": cn_naive
    })

df = pd.DataFrame(results)
print("Chr1 copy number segments:")
print(df.to_string(index=False))
print()
print("Note: CN_call uses log2 thresholds (robust to noise); CN_naive rounds 2*2^log2.")
print("Both agree for clean segments; thresholds are more robust at segment boundaries.")
print()
print("Expected log2 values by copy number (diploid, ploidy=2):")
for cn in range(6):
    import math
    if cn == 0:
        print(f"  CN={cn}: log2(0/2) = -inf  (homozygous deletion)")
    else:
        print(f"  CN={cn}: log2({cn}/2) = {math.log2(cn/2):+.3f}")
```python

## Gene-level Annotation

After calling segments, the next step is to determine which genes overlap each CNV. In practice:

```python
# Conceptual pattern (requires pybedtools or pyranges)
import pyranges as pr

# Gene annotation BED (from Ensembl/UCSC)
genes = pr.read_bed("hg38_genes.bed")

# CNV segments
cnv_df = pd.DataFrame({"Chromosome": ["chr1","chr1"], "Start": [2500000, 7500000],
                        "End": [3750000, 9500000], "CN": [1, 4]})
cnvs = pr.from_dict(cnv_df)

# Overlap
overlap = cnvs.join(genes)
print(overlap.df[["CN","Name"]].head())
```python

**Downstream analysis:**
- Amplified oncogenes: look up in COSMIC Cancer Gene Census
- Deleted tumor suppressors: TP53, BRCA1/2, RB1 are common targets
- Focal vs. arm-level events: focal (<3 Mb) are more often driver events

## Copy Number in Cancer

Copy number data from TCGA is publicly available for 33 cancer types. The TCGA CNV pipeline uses GISTIC2 to identify recurrently amplified/deleted regions.

**Access:** `https://gdc.cancer.gov/` (portal) or `TCGAbiolinks` (R) or `xenahubs.net` (Python).

```python
# Conceptual: load TCGA BRCA copy number (segment files)
import pandas as pd
# Download: https://portal.gdc.cancer.gov/exploration?filters=...
brca_cnv = pd.read_csv("TCGA-BRCA.cnv.tsv", sep="\t")
# Columns: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean
print(brca_cnv.head())
```python

Copy number analysis extends the variant calling concepts in Module 02. For somatic SV detection, see Module 17 (Genome Assembly and Advanced NGS).

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
