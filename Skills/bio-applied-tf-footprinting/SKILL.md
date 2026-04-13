---
name: bio-applied-tf-footprinting
description: "*Prerequisites: Modules 01–22. Module 01 (NGS), Module 15 (Motif Discovery) strongly recommended.*"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/23_TF_Footprinting/23_tf_footprinting.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 23: TF Footprinting & Chromatin Accessibility

*Source: Course notebook `Tier_3_Applied_Bioinformatics/23_TF_Footprinting/23_tf_footprinting.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 23**

*Prerequisites: Modules 01–22. Module 01 (NGS), Module 15 (Motif Discovery) strongly recommended.*

---

Transcription factor footprinting identifies where specific TFs are bound in the genome by analyzing Tn5 insertion patterns from ATAC-seq data. This module covers the complete analysis from fragment size QC to footprint scoring and accumulation plots.

**By the end of this module you will be able to:**
- Understand the ATAC-seq Tn5 insertion bias concept
- Compute expected vs. observed insertion profiles around motif sites
- Calculate footprint scores
- Perform genomic interval operations with pybedtools
- Build accumulation (aggregate) plots around genomic features

**Attribution:** *Code patterns inspired by TotipotencyLab chromatin analysis workflows, Max Planck Institute of Biochemistry. Uses simulated ATAC-seq data and public ENCODE ChIP-seq patterns.*

## Background: ATAC-seq and TF Footprinting

**ATAC-seq** (Assay for Transposase-Accessible Chromatin with sequencing) uses Tn5 transposase to cut and tag open chromatin. Nucleosome-free regions (NFRs) are cut frequently; nucleosome-occupied DNA is protected.

**Fragment size distribution:**
- < 150 bp: nucleosome-free region (NFR) — most informative for TF binding
- ~200 bp: mononucleosome
- ~400 bp: dinucleosome
- etc.

**Footprinting concept:** When a TF is bound, it physically blocks Tn5 from cutting at the exact binding site. This creates a "footprint" — a local depletion of cuts at the TF motif center, flanked by elevated cuts at the accessible chromatin on either side.

| Feature | TF-bound sites | Unbound accessible sites |
|---|---|---|
| Cuts at motif center | Depleted (footprint) | Elevated |
| Cuts flanking motif | Elevated | Uniform |
| Footprint score | High | Low |

**Analysis workflow:**
1. Map ATAC-seq reads → BAM
2. Extract Tn5 insertion positions (5' end of each read, shifted +4/-5 bp for Tn5 offset)
3. Define TF motif sites (from PWM scan or ChIP-seq peaks)
4. Compute insertion profiles around each motif site
5. Average across all sites → aggregate profile
6. Compute footprint score = flanking signal / central signal

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({
    "figure.dpi": 120,
    "axes.spines.top": False,
    "axes.spines.right": False,
})
print("Libraries loaded.")

try:
    import pybedtools
    print(f"pybedtools available: {pybedtools.__version__}")
    PYBEDTOOLS_AVAILABLE = True
except ImportError:
    print("pybedtools not installed — interval operations shown as patterns")
    PYBEDTOOLS_AVAILABLE = False
```python

## 1. Fragment Size Distribution

The fragment size distribution is a key ATAC-seq quality metric. A high-quality ATAC-seq library shows a characteristic nucleosomal ladder pattern.

We simulate fragment sizes here; in practice, these are computed from a BAM file using `pysam`.

```python
rng = np.random.default_rng(42)

def simulate_atac_fragments(n=50_000):
    """Simulate ATAC-seq fragment sizes with nucleosomal ladder."""
    # NFR: ~100–150 bp (40% of fragments)
    nfr = rng.normal(120, 25, int(n * 0.40))
    nfr = nfr[(nfr > 50) & (nfr < 200)]

    # Mono-nucleosome: ~180–250 bp (35%)
    mono = rng.normal(200, 20, int(n * 0.35))
    mono = mono[(mono > 150) & (mono < 300)]

    # Di-nucleosome: ~350–450 bp (15%)
    di = rng.normal(400, 30, int(n * 0.15))
    di = di[(di > 300) & (di < 500)]

    # Tri+: ~550+ bp (10%)
    tri = rng.normal(600, 50, int(n * 0.10))
    tri = tri[(tri > 480) & (tri < 800)]

    return np.concatenate([nfr, mono, di, tri]).astype(int)

fragment_sizes = simulate_atac_fragments()
print(f"Total fragments: {len(fragment_sizes):,}")
print(f"NFR (<150 bp): {(fragment_sizes < 150).mean():.1%}")
print(f"Mono (150–300 bp): {((fragment_sizes >= 150) & (fragment_sizes < 300)).mean():.1%}")

fig, ax = plt.subplots(figsize=(10, 4))
ax.hist(fragment_sizes, bins=range(50, 801, 10), color="steelblue", alpha=0.8, density=True)
for pos, label in [(150, "NFR"), (200, "Mono"), (400, "Di"), (600, "Tri+")]:
    ax.axvline(pos, color="red", lw=1, ls="--", alpha=0.7)
ax.set_xlabel("Fragment size (bp)")
ax.set_ylabel("Density")
ax.set_title("ATAC-seq fragment size distribution — nucleosomal ladder")
ax.text(105, ax.get_ylim()[1]*0.9, "NFR", ha="center", fontsize=9)
ax.text(200, ax.get_ylim()[1]*0.85, "Mono", ha="center", fontsize=9)
ax.text(400, ax.get_ylim()[1]*0.7, "Di", ha="center", fontsize=9)
plt.tight_layout()
plt.show()
```python

## 2. Tn5 Insertion Profiles Around Motif Sites

The footprinting analysis centers on the *insertion profile*: how many Tn5 cuts occur at each position relative to a TF motif center.

**Tn5 offset correction:** The Tn5 enzyme cuts as a dimer; the actual insertion is offset from the read start by +4 bp (forward strand) and -5 bp (reverse strand). This correction is essential for footprinting.

We simulate insertion profiles for:
1. **Bound sites** — TF is present → footprint (depletion at center)
2. **Unbound accessible sites** — no TF → uniform or slightly elevated cuts

```python
WINDOW = 200   # ± 200 bp around motif center
positions = np.arange(-WINDOW, WINDOW + 1)

def simulate_insertion_profile(n_sites=500, footprint_depth=0.4, noise_level=0.15):
    """
    Simulate Tn5 insertion profile around TF motif sites.
    footprint_depth: how much the center is depleted (0 = no footprint, 1 = complete depletion)
    """
    rng = np.random.default_rng(42)
    profile = np.ones(len(positions)) * 50  # baseline

    # Flanking peaks at ±8 bp (nucleosome positioning)
    for flank in [-8, 8]:
        idx = WINDOW + flank
        profile[idx-3:idx+3] += 30

    # Central depletion (footprint)
    center_depletion = profile[WINDOW-5:WINDOW+5].mean() * footprint_depth
    profile[WINDOW-5:WINDOW+5] -= center_depletion

    # Add Poisson noise
    profile = rng.poisson(np.maximum(profile, 1)).astype(float)

    # Smooth
    return gaussian_filter1d(profile, sigma=2)

# Bound sites (strong footprint)
profile_bound = simulate_insertion_profile(footprint_depth=0.55)
# Unbound sites (no footprint)
profile_unbound = simulate_insertion_profile(footprint_depth=0.0, noise_level=0.2)

fig, axes = plt.subplots(1, 2, figsize=(14, 4))
for ax, profile, title, color in [
    (axes[0], profile_bound,   "TF-bound sites (footprint visible)", "steelblue"),
    (axes[1], profile_unbound, "Unbound accessible sites",           "salmon"),
]:
    ax.plot(positions, profile, color=color, lw=2)
    ax.axvline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.axvspan(-8, 8, alpha=0.1, color="red", label="Motif window")
    ax.set_xlabel("Position relative to motif center (bp)")
    ax.set_ylabel("Tn5 insertions (smoothed)")
    ax.set_title(title)
    ax.legend(frameon=False)

plt.tight_layout()
plt.show()
```python

## 3. Footprint Score

The footprint score quantifies how pronounced the depletion is at the motif center relative to its flanks:

$$\text{Footprint Score} = \frac{\text{mean(flanking signal)}}{\text{mean(central signal) + ε}}$$

Higher score = stronger footprint = more likely TF binding.

A score > 1 means the center is depleted relative to flanks. A well-footprinted TF typically has scores of 1.5–3.

```python
def footprint_score(profile, center_window=8, flank_window=(20, 60)):
    """
    Compute footprint score.
    center_window: bp around motif center considered as footprint
    flank_window: (inner, outer) bp range for flanking signal
    """
    mid = len(profile) // 2
    central = profile[mid - center_window : mid + center_window].mean()
    left_flank  = profile[mid - flank_window[1] : mid - flank_window[0]].mean()
    right_flank = profile[mid + flank_window[0] : mid + flank_window[1]].mean()
    flanking = (left_flank + right_flank) / 2
    return flanking / (central + 1e-6)

score_bound   = footprint_score(profile_bound)
score_unbound = footprint_score(profile_unbound)

print(f"Footprint score — bound sites:   {score_bound:.3f}")
print(f"Footprint score — unbound sites: {score_unbound:.3f}")
print(f"Ratio: {score_bound / score_unbound:.2f}x")

# Score distribution across many simulated site types
rng = np.random.default_rng(0)
scores_bound_dist   = [footprint_score(simulate_insertion_profile(footprint_depth=rng.uniform(0.4,0.7))) for _ in range(100)]
scores_unbound_dist = [footprint_score(simulate_insertion_profile(footprint_depth=rng.uniform(0.0,0.1))) for _ in range(100)]

fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(scores_bound_dist,   bins=20, alpha=0.7, color="steelblue", label="Bound sites")
ax.hist(scores_unbound_dist, bins=20, alpha=0.7, color="salmon",    label="Unbound sites")
ax.set_xlabel("Footprint score"); ax.set_ylabel("Count")
ax.set_title("Footprint score distributions")
ax.legend(frameon=False)
plt.tight_layout()
plt.show()

stat, p = mannwhitneyu(scores_bound_dist, scores_unbound_dist, alternative="greater")
print(f"\nMann-Whitney U test: p = {p:.2e}  (bound > unbound)")
```python

## 4. Genomic Interval Operations with pybedtools

Footprinting workflows require intersecting multiple genomic coordinate files:
- ATAC-seq peaks (BED format)
- TF motif sites (from PWM scan or ChIP-seq)
- Gene annotations (GTF/BED)

`pybedtools` wraps BEDTools and provides a Python interface for these operations.

```python
# pybedtools patterns (shown with synthetic data)
print("=== pybedtools Genomic Interval Patterns ===\n")

if PYBEDTOOLS_AVAILABLE:
    import pybedtools

    # Create synthetic ATAC peaks and motif sites
    atac_peaks_str = """chr1\t1000\t2000\tpeak1\t500\t.
chr1\t5000\t6000\tpeak2\t300\t.
chr1\t10000\t11000\tpeak3\t200\t.
chr2\t2000\t3000\tpeak4\t400\t."""

    motif_sites_str = """chr1\t1400\t1420\tCTCF_1\t800\t+
chr1\t1600\t1620\tCTCF_2\t700\t-
chr1\t5200\t5220\tCTCF_3\t600\t+
chr2\t2500\t2520\tCTCF_4\t900\t+
chr3\t1000\t1020\tCTCF_5\t400\t+"""  # chr3 has no peaks

    peaks = pybedtools.BedTool(atac_peaks_str, from_string=True)
    motifs = pybedtools.BedTool(motif_sites_str, from_string=True)

    # Intersection: motifs overlapping ATAC peaks
    motifs_in_peaks = motifs.intersect(peaks, u=True)
    print("Motif sites overlapping ATAC peaks:")
    for interval in motifs_in_peaks:
        print(f"  {interval.chrom}:{interval.start}-{interval.end} ({interval.name})")

    # Sites NOT in peaks (negative control)
    motifs_not_in_peaks = motifs.intersect(peaks, v=True)
    print(f"\nMotif sites NOT in peaks (background): {len(motifs_not_in_peaks)}")

    # Extend sites to ±200 bp window for footprinting
    motifs_extended = motifs_in_peaks.slop(b=200, g={"chr1": 250_000_000, "chr2": 243_000_000})
    print(f"\nExtended to ±200 bp window: {len(motifs_extended)} sites")

else:
    print("""pybedtools patterns (install with: pip install pybedtools)

# Load BED files
peaks  = pybedtools.BedTool("atac_peaks.bed")
motifs = pybedtools.BedTool("ctcf_motifs.bed")

# Intersection: motifs in open chromatin
motifs_in_peaks = motifs.intersect(peaks, u=True)

# Extend to footprinting window
motifs_extended = motifs_in_peaks.slop(b=200, genome="hg38")

# Subtract blacklist regions
blacklist = pybedtools.BedTool("hg38_blacklist.bed")
clean = motifs_extended.subtract(blacklist)

# Closest gene TSS
tss = pybedtools.BedTool("hg38_tss.bed")
nearest = motifs_in_peaks.closest(tss, d=True)
""")
```python

## 5. Accumulation Plots (Signal Around Genomic Features)

An **accumulation plot** (also called a meta-profile or aggregate plot) shows the average signal across many genomic loci. This reveals systematic patterns that are too weak to see at individual sites.

**Applications:**
- Average ChIP-seq signal around TSS regions
- ATAC-seq signal around TF binding sites
- H3K27ac signal around enhancers

**Workflow:**
1. Define anchor points (e.g. motif centers)
2. Extract signal in a window around each anchor
3. Average signal across all anchors
4. Plot mean ± SE

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
