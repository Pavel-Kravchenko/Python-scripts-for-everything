---
name: bio-applied-tf-footprinting
description: TF footprinting from ATAC-seq — Tn5 insertion profiles, footprint score calculation, pybedtools interval operations, accumulation plots.
tool_type: python
primary_tool: NumPy
---

# TF Footprinting & Chromatin Accessibility

## Background

**ATAC-seq**: Tn5 transposase cuts and tags open chromatin. Nucleosome-free regions (NFRs) are cut frequently; nucleosome-occupied DNA is protected.

**Footprinting**: A bound TF blocks Tn5 at the binding site → local depletion of cuts at the motif center flanked by elevated cuts on either side.

| Feature | TF-bound sites | Unbound accessible sites |
|---|---|---|
| Cuts at motif center | Depleted (footprint) | Elevated |
| Cuts flanking motif | Elevated | Uniform |
| Footprint score | High (>1.5) | ~1.0 |

## Pitfalls

- **Tn5 offset correction is mandatory**: The Tn5 dimer inserts with a +4 bp (forward strand) / -5 bp (reverse strand) offset from the read 5' end. Skip this and footprint centers will be shifted.
- **Fragment size selection**: Use only NFR fragments (<150 bp) for footprinting. Mono/dinucleosomal fragments dilute the footprint signal.
- **Low-occupancy TFs don't footprint**: Only highly occupied (>50% of sites bound) TFs produce detectable footprints. CTCF, cohesin, and pioneer factors work well; transient binders (MYC) rarely do.
- **Sequence bias**: Tn5 has a strong insertion sequence preference (~10 bp periodic). Always compare to expected (background) insertion profile, not just observed.
- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors.
- **Batch effects**: Always check for batch confounding before interpreting biological signal.

## Fragment Size Distribution (QC)

Expected ladder for a high-quality ATAC-seq library:

| Peak | Size | Interpretation |
|---|---|---|
| NFR | <150 bp | Nucleosome-free — most informative for TF binding |
| Mono | ~200 bp | Mononucleosome |
| Di | ~400 bp | Dinucleosome |
| Tri+ | ~600 bp | Trinucleosome and higher |

NFR fraction >40% of total fragments = good library quality.

## Analysis Workflow

```
1. Map ATAC-seq reads → BAM
2. Extract Tn5 insertion positions (5' end, +4/-5 bp offset)
3. Define TF motif sites (PWM scan or ChIP-seq peaks)
4. Compute insertion profiles in ±200 bp window around each site
5. Average across all sites → aggregate profile
6. Compute footprint score = flanking signal / central signal
```

## Footprint Score

$$\text{Footprint Score} = \frac{\text{mean(flanking signal)}}{\text{mean(central signal) + ε}}$$

Score >1 = center depleted relative to flanks. Well-footprinted TFs: 1.5–3.0.

```python
from scipy.ndimage import gaussian_filter1d

def footprint_score(profile, center_window=8, flank_window=(20, 60)):
    """
    center_window: bp around motif center counted as footprint
    flank_window: (inner, outer) bp range for flanking signal
    """
    mid = len(profile) // 2
    central     = profile[mid - center_window : mid + center_window].mean()
    left_flank  = profile[mid - flank_window[1] : mid - flank_window[0]].mean()
    right_flank = profile[mid + flank_window[0] : mid + flank_window[1]].mean()
    flanking = (left_flank + right_flank) / 2
    return flanking / (central + 1e-6)

# Smooth profile before scoring
profile_smooth = gaussian_filter1d(raw_profile, sigma=2)
score = footprint_score(profile_smooth)
```

## pybedtools Patterns

```python
import pybedtools

peaks  = pybedtools.BedTool("atac_peaks.bed")
motifs = pybedtools.BedTool("ctcf_motifs.bed")

# Motifs overlapping ATAC peaks (accessible sites)
motifs_in_peaks = motifs.intersect(peaks, u=True)

# Negative control: sites NOT in peaks
background = motifs.intersect(peaks, v=True)

# Extend to ±200 bp window for footprinting
motifs_extended = motifs_in_peaks.slop(b=200, genome="hg38")

# Remove ENCODE blacklist regions
blacklist = pybedtools.BedTool("hg38_blacklist.bed")
clean = motifs_extended.subtract(blacklist)

# Closest gene TSS
tss = pybedtools.BedTool("hg38_tss.bed")
nearest = motifs_in_peaks.closest(tss, d=True)
```

Install: `pip install pybedtools` (requires BEDTools in PATH)

## Accumulation Plot (Meta-Profile)

```python
import numpy as np

def accumulation_plot(bam, sites, window=200):
    """
    Compute average Tn5 insertion profile across all sites.
    bam: pysam.AlignmentFile
    sites: list of (chrom, center) tuples
    Returns: array of shape (2*window+1,)
    """
    import pysam
    profile = np.zeros(2 * window + 1)
    for chrom, center in sites:
        for read in bam.fetch(chrom, max(0, center - window), center + window):
            if read.is_unmapped: continue
            # Tn5 offset correction
            if not read.is_reverse:
                pos = read.reference_start + 4
            else:
                pos = read.reference_end - 5 - 1
            offset = pos - center
            if -window <= offset <= window:
                profile[offset + window] += 1
    return profile / max(len(sites), 1)
```
