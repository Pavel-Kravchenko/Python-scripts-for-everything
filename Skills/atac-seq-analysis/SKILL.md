---
name: atac-seq-analysis
description: ATAC-seq quality control and accessibility analysis.
tool_type: python
primary_tool: NumPy
---

## When to Use

Use this atomic skill for focused work on **atac-seq-analysis** without bundling unrelated topics.

## Quick Reference

This skill was split from `tf-footprinting-atac.md` to keep topics independent and self-contained.

## Core Patterns

Use the parent material below as the source reference, then keep implementations specific to this topic.

## Source Reference (from merged skill)

name: tf-footprinting-atac
description: ATAC-seq analysis and TF footprinting — fragment size QC, Tn5 insertion profiles, footprint score, pybedtools interval operations, accumulation/meta-profile plots, NFR fraction

## When to Use

Use this skill when:
- Performing ATAC-seq quality control (fragment size distribution, NFR fraction)
- Computing Tn5 insertion profiles around transcription factor motif sites
- Calculating footprint scores to detect TF occupancy from chromatin accessibility
- Using pybedtools to intersect ATAC peaks with motif sites or gene annotations
- Building accumulation (meta-profile) plots averaged across many genomic loci
- Integrating ATAC-seq with motif discovery (Module 15 → Module 23 workflow)

## Quick Reference

| Metric | Threshold | Notes |
|---|---|---|
| NFR fraction (< 150 bp) | ≥ 40% | Key ATAC-seq quality indicator |
| Mono-nucleosome fraction | ~35% (180–250 bp) | characteristic ladder band |
| Di-nucleosome fraction | ~15% (350–450 bp) | |
| Tn5 forward offset | +4 bp | shift forward read starts +4 |
| Tn5 reverse offset | -5 bp | shift reverse read starts -5 |
| Footprint score > 1.5 | Strong footprint | flanking/central insertion ratio |
| Center window for score | ± 8 bp from motif center | |
| Flanking window for score | 20–60 bp from motif center | |
| FRiP (Fraction of Reads in Peaks) | ≥ 0.20 | enrichment quality check |

| Tool | Purpose |
|---|---|
| `pysam` | Read BAM file, compute Tn5 insertions per strand |
| `pybedtools` | Intersect BED files (peaks ∩ motifs, slop, subtract) |
| `scipy.ndimage.gaussian_filter1d` | Smooth insertion profiles |
| `scipy.stats.mannwhitneyu` | Test footprint depth difference across TF families |
| `numpy` | Accumulation / meta-profile averaging |

## Key Patterns

**Pattern 1: Fragment size distribution and NFR fraction**
```python
import numpy as np
import matplotlib.pyplot as plt

def classify_fragments(sizes):
    """Classify ATAC fragments by nucleosomal content."""
    nfr   = sizes[(sizes > 50) & (sizes < 150)]
    mono  = sizes[(sizes >= 150) & (sizes < 300)]
    di    = sizes[(sizes >= 300) & (sizes < 500)]
    tri   = sizes[sizes >= 500]
    total = len(sizes)
    return {
        "NFR":  len(nfr)/total,
        "Mono": len(mono)/total,
        "Di":   len(di)/total,
        "Tri+": len(tri)/total,
    }

fracs = classify_fragments(fragment_sizes)
print(f"NFR fraction: {fracs['NFR']:.1%}")  # target ≥ 40%

# Plot fragment size histogram
plt.hist(fragment_sizes, bins=200, range=(0, 800), density=True, color="steelblue")
plt.xlabel("Fragment size (bp)"); plt.ylabel("Density")
plt.title("ATAC-seq fragment size distribution (nucleosomal ladder)")
```

**Pattern 2: Tn5 offset correction**
```python
def tn5_corrected_insertions(bam_file, chrom, start, end):
    """
    Compute per-base Tn5 insertion counts with offset correction.
    Forward reads: insert at read_start + 4
    Reverse reads: insert at read_end - 5 (= read_start + len - 5)
    Requires pysam.
    """
    import pysam
    insertions = np.zeros(end - start)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary:
                continue
            if not read.is_reverse:
                pos = read.reference_start + 4  # forward +4
            else:
                pos = read.reference_end - 5    # reverse -5
            if start <= pos < end:
                insertions[pos - start] += 1
    return insertions
```

**Pattern 3: Simulate Tn5 insertion profile (without BAM)**
```python
WINDOW = 200   # ± 200 bp around motif center

def simulate_insertion_profile(n_sites=500, footprint_depth=0.4, noise=0.15):
    """Simulate average Tn5 insertion profile with footprint depletion."""
    rng = np.random.default_rng(42)
    positions = np.arange(-WINDOW, WINDOW + 1)
    profile = np.ones(len(positions)) * 50

    # Flanking nucleosome positioning peaks at ±8 bp
    for flank in [-8, 8]:
        idx = WINDOW + flank
        profile[idx-5:idx+5] += 40

    # Central depletion (footprint)
    profile[WINDOW-8:WINDOW+8] *= (1 - footprint_depth)

    # Add Poisson noise over n_sites
    noisy = rng.poisson(profile, (n_sites, len(positions)))
    return noisy.mean(axis=0), noisy.std(axis=0) / np.sqrt(n_sites)
```

**Pattern 4: Footprint score**
```python
def footprint_score(profile, center_window=8, flank_window=(20, 60)):
    """
    Footprint score = flanking_mean / (center_mean + ε)
    Higher score = stronger TF occupancy signal.
    """
    mid = len(profile) // 2
    central = profile[mid - center_window : mid + center_window].mean()
    left  = profile[mid - flank_window[1] : mid - flank_window[0]].mean()
    right = profile[mid + flank_window[0] : mid + flank_window[1]].mean()
    flanking = (left + right) / 2
    return flanking / (central + 1e-9)
```

**Pattern 5: pybedtools interval operations**
```python
import pybedtools

peaks  = pybedtools.BedTool("atac_peaks.bed")
motifs = pybedtools.BedTool("ctcf_motifs.bed")

# Motif sites overlapping ATAC peaks
motifs_in_peaks = motifs.intersect(peaks, u=True)

# Extend motif sites by 200 bp for footprint window
motif_windows = motifs_in_peaks.slop(b=200, genome="hg38")

# Remove blacklisted regions
blacklist = pybedtools.BedTool("hg38_blacklist.bed")
clean = motif_windows.subtract(blacklist)

# Save as FASTA for motif scanning
# clean.sequence(fihg38.fa, fowindows.fa)
print(f"Motifs in peaks: {len(motifs_in_peaks)}")
```

## Code Templates

**Template 1: Accumulation (meta-profile) plot**
```python
import numpy as np, matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def accumulation_plot(profiles, positions, title="Meta-profile", sigma=5):
    """
    profiles: (n_sites, window_size) array of per-site insertion counts
    positions: 1D array of genomic offsets (e.g. -200 to +200)
    """
    mean_profile = profiles.mean(axis=0)
    sem_profile  = profiles.std(axis=0) / np.sqrt(len(profiles))
    smoothed     = gaussian_filter1d(mean_profile, sigma=sigma)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(positions, smoothed, color="steelblue", lw=2)
    ax.fill_between(positions,
                    smoothed - sem_profile,
                    smoothed + sem_profile,
                    alpha=0.3, color="steelblue")
    ax.axvline(0, color="red", lw=1, ls="--", label="Motif center")
    ax.set_xlabel("Position relative to motif (bp)")
    ax.set_ylabel("Mean Tn5 insertions")
    ax.set_title(title); ax.legend(frameon=False)
    plt.tight_layout()
    return fig, ax

fig, ax = accumulation_plot(all_profiles, np.arange(-WINDOW, WINDOW+1))
```

**Template 2: Full footprinting QC pipeline**
```python
import numpy as np
from scipy.ndimage import gaussian_filter1d

def footprinting_pipeline(bam_file, peaks_bed, motifs_bed, window=200):
    """End-to-end footprinting: load peaks → filter motifs → extract profiles → score."""
    import pybedtools

    # Step 1: motifs in peaks
    peaks  = pybedtools.BedTool(peaks_bed)
    motifs = pybedtools.BedTool(motifs_bed).intersect(peaks, u=True)

    # Step 2: extract insertion profiles per motif site
    profiles = []
    for site in motifs:
        chrom, center = site.chrom, (site.start + site.end) // 2
        profile = tn5_corrected_insertions(
            bam_file, chrom, center - window, center + window)
        profiles.append(profile)

    if not profiles:
        return None, None

    profiles = np.array(profiles)

    # Step 3: compute aggregate profile and footprint score
    mean_p = gaussian_filter1d(profiles.mean(axis=0), sigma=3)
    score  = footprint_score(mean_p)
    return mean_p, score

# mean_profile, fp_score  footprinting_pipeline(sample.bam, peaks.bed, ctcf.bed)
```

## Pitfalls

- **Missing Tn5 offset:** uncorrected reads shift the apparent footprint ±4–5 bp; always apply +4/-5 before aggregating
- **Not filtering NFR vs. nucleosomal reads:** footprinting uses only NFR (< 150 bp) fragments — nucleosomal reads dilute the signal
- **Genome blacklist regions:** always subtract ENCODE blacklist before analysis; these regions have artifactually high depth
- **Motif strand:** CTCF and other factors have asymmetric binding; flip reverse-strand sites before averaging profiles
- **Footprint score vs. occupancy:** high score only means accessibility asymmetry; validate with ChIP-seq if possible
- **pybedtools genome arg:** `slop` requires genome file or dict `{'chr1': (0, 248956422), ...}`; use `pybedtools.genome_registry.hg38` for convenience

## Related Skills

- `motif-discovery` — PWM scanning to generate motif sites BED file for footprinting
- `hic-analysis` — pileup/aggregate plots share the meta-profile concept
- `ngs-variant-calling` — upstream BAM alignment from ATAC-seq reads
- `copy-number-analysis` — ATAC-seq depth normalization uses similar rolling-median approach


## Related Skills

- `atac-seq-analysis` (this file)
- `tf-footprinting-atac` (legacy merged skill)
