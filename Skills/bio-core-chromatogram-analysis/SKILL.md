---
name: bio-core-chromatogram-analysis
description: "- Explain how Sanger (dideoxy chain termination) sequencing works - Read `.ab1` chromatogram files with BioPython - Plot and interpret trace data with matplotlib - Work with Phred quality scores and b"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/09_Chromatogram_Analysis/01_chromatogram_analysis.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Chromatogram Analysis: Sanger Sequencing in Practice

*Source: Course notebook `Tier_2_Core_Bioinformatics/09_Chromatogram_Analysis/01_chromatogram_analysis.ipynb`*


---

## Learning Objectives

By the end of this notebook you will be able to:

- Explain how Sanger (dideoxy chain termination) sequencing works
- Read `.ab1` chromatogram files with BioPython
- Plot and interpret trace data with matplotlib
- Work with Phred quality scores and base calling
- Identify mixed bases, heterozygous positions, and low-quality regions
- Understand assembly of multiple Sanger reads
- Know when to use Sanger sequencing vs NGS in a modern lab

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **Phred scores are log-scale**: Q20 = 1 error per 100 bases (99% accuracy); Q30 = 1 per 1,000 (99.9%); Q40 = 1 per 10,000 (99.99%). Most analyses require Q20 as a minimum threshold, and Q30 for high-confidence variant calling. The difference between Q20 and Q30 is enormous in practice: at Q20 you expect 7 errors per 700 bp read; at Q30 only 0.7.
- **Sanger read quality profile**: The first ~20-30 bases have poor quality due to primer signal and dye blobs. Quality is highest from ~bases 30-700. After ~800 bp, signal decay causes quality to drop sharply. Always trim both ends before downstream analysis.
- **Double peaks indicate heterozygosity (or contamination)**: Two overlapping peaks of comparable height at a position indicate that two different bases are present. This could be a genuine heterozygous SNP, a mixed PCR product, two overlapping clones, or a frameshift that creates two reading frames. Check the surrounding context to distinguish.
- **ab1 vs. abi vs. scf**: The `.ab1` format (Applied Biosystems) stores the raw fluorescence traces, base calls, quality scores, and instrument metadata. `.abi` is an alternative extension for the same format. `.scf` (Standard Chromatogram Format) is an older format from MRC-LMB. BioPython reads all three via `SeqIO.read(file, 'abi')`.
- **Why 4 fluorescent dyes?**: Each dideoxynucleotide (ddA, ddC, ddG, ddT) is labeled with a different fluorescent dye (e.g., ddA=green, ddC=blue, ddG=black, ddT=red in ABI nomenclature). Capillary electrophoresis separates fragments by size, and a laser excites the dye as it passes — this is how the instrument determines the sequence.

## Environment check (run this first)

```python
# Environment check
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt

print("Imports ready.")

# Phred quality scale reference
print("\nPhred quality scale:")
print(f"{'Q score':<10} {'Error prob':<15} {'Accuracy':<12} {'Threshold'}")
print("-" * 55)
for q in [10, 20, 30, 40]:
    prob = 10 ** (-q/10)
    acc = (1 - prob) * 100
    threshold = "Minimum" if q == 20 else ("Standard" if q == 30 else ("High" if q == 40 else "Poor"))
    print(f"Q{q:<9} {prob:<15.4f} {acc:<12.3f}% {threshold}")

print("\nNote: Real .ab1 files are needed for Sections 2-3.")
print("Sections 3-10 use synthetic chromatogram data generated in-notebook.")
```python

```python
def read_ab1(filepath):
    """
    Read an .ab1 chromatogram file and print a summary.
    Returns the SeqRecord.
    """
    record = SeqIO.read(filepath, "abi")

    print(f"File:            {filepath}")
    print(f"Sample ID:       {record.id}")
    print(f"Sequence length: {len(record.seq)} bp")
    print(f"First 60 bases:  {record.seq[:60]}")

    # Quality scores
    quals = record.letter_annotations["phred_quality"]
    print(f"Mean Phred:      {np.mean(quals):.1f}")
    print(f"Median Phred:    {np.median(quals):.1f}")

    return record

# ----- Uncomment when you have a real file: -----
# record = read_ab1(DATA_DIR + "sample_forward.ab1")

# For this tutorial we will also build synthetic examples
# so the notebook runs without data files.
```python

### 2.2 What is Inside the Record?

The `SeqRecord` returned by BioPython stores traces in `record.annotations["abif_raw"]`.
Key tags:

| Tag | Content |
|-----|---------|
| `DATA9`-`DATA12` | Processed trace channels (order depends on dye set) |
| `DATA1`-`DATA4` | Raw trace channels |
| `PLOC1` / `PLOC2` | Peak locations (scan numbers) |
| `PBAS1` / `PBAS2` | Called bases (raw / edited) |
| `PCON1` / `PCON2` | Quality values |
| `FWO_1` | Base order for the four channels (e.g. `b'GATC'`) |

```python
def inspect_ab1_channels(record):
    """Show metadata stored in an .ab1 record."""
    raw = record.annotations.get("abif_raw", {})

    # Base order tells us which channel corresponds to which nucleotide
    base_order = raw.get("FWO_1", b"GATC").decode()
    print(f"Channel base order: {base_order}")

    # Processed trace data
    channels = {}
    for i, base in enumerate(base_order):
        tag = f"DATA{9 + i}"  # DATA9, DATA10, DATA11, DATA12
        trace = raw.get(tag, ())
        channels[base] = trace
        print(f"  Channel {base} ({tag}): {len(trace)} data points")

    # Peak locations
    peak_locs = raw.get("PLOC1", ())
    print(f"  Peak locations: {len(peak_locs)} called bases")

    return channels, peak_locs

# Uncomment with real data:
# channels, peaks = inspect_ab1_channels(record)
```python

---

## 3. Chromatogram Visualization

Visualizing the four fluorescence traces is the most informative way to assess sequencing
quality. Good chromatograms have:
- Well-separated, evenly-spaced peaks
- Low baseline noise
- Consistent peak heights

We will first build a synthetic chromatogram for demonstration, then show how to
plot a real `.ab1` file.

```python
def make_synthetic_trace(sequence, peak_spacing=15, peak_width=4.0, noise_level=0.05):
    """
    Generate synthetic chromatogram traces for demonstration.

    Returns:
        traces: dict with keys A, T, G, C -> numpy arrays
        peak_positions: list of int scan positions for each base
    """
    n_bases = len(sequence)
    n_points = (n_bases + 2) * peak_spacing
    x = np.arange(n_points)

    traces = {b: np.zeros(n_points) for b in "ATGC"}
    peak_positions = []

    rng = np.random.default_rng(42)

    for i, base in enumerate(sequence.upper()):
        center = (i + 1) * peak_spacing
        peak_positions.append(center)
        # Gaussian peak for the called base
        amplitude = 0.8 + 0.2 * rng.random()
        if base in traces:
            traces[base] += amplitude * np.exp(-0.5 * ((x - center) / peak_width) ** 2)

    # Add background noise
    for b in traces:
        traces[b] += noise_level * rng.random(n_points)
        traces[b] = np.clip(traces[b], 0, None)

    return traces, peak_positions


demo_seq = "ATGCTAGCGATCGATCGATCGATCGATCGA"
traces, peaks = make_synthetic_trace(demo_seq)

print(f"Sequence:    {demo_seq}")
print(f"Trace length: {len(traces['A'])} scan points")
print(f"Bases called: {len(peaks)}")
```python

```python
# Standard colors used for chromatograms
TRACE_COLORS = {"A": "green", "C": "blue", "G": "black", "T": "red"}


def plot_chromatogram(traces, peak_positions, sequence, start=0, end=None, title="Chromatogram"):
    """
    Plot a four-color chromatogram trace.

    Parameters:
        traces: dict of base -> numpy array of fluorescence values
        peak_positions: list of scan positions per called base
        sequence: string of called bases
        start, end: base index range to display (default: full sequence)
    """
    if end is None:
        end = len(sequence)

    # Determine scan range
    scan_start = max(0, peak_positions[start] - 10)
    scan_end = min(len(traces["A"]), peak_positions[min(end, len(peak_positions)) - 1] + 10)

    fig, ax = plt.subplots(figsize=(14, 3.5))

    x = np.arange(scan_start, scan_end)
    for base, color in TRACE_COLORS.items():
        ax.plot(x, traces[base][scan_start:scan_end], color=color, linewidth=0.8, label=base)

    # Label called bases at peak positions
    for i in range(start, min(end, len(peak_positions))):
        pos = peak_positions[i]
        base = sequence[i]
        color = TRACE_COLORS.get(base, "gray")
        ax.text(pos, -0.07 * ax.get_ylim()[1] if ax.get_ylim()[1] else -0.05,
                base, ha="center", va="top", fontsize=7, fontweight="bold", color=color)

    ax.set_xlabel("Scan position")
    ax.set_ylabel("Fluorescence intensity")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8, ncol=4)
    ax.set_xlim(scan_start, scan_end)
    plt.tight_layout()
    plt.show()


plot_chromatogram(traces, peaks, demo_seq, start=0, end=30, title="Synthetic Chromatogram")
```python

```python
def plot_ab1_chromatogram(record, start=0, end=None, title=None):
    """
    Plot chromatogram directly from a BioPython .ab1 SeqRecord.
    """
    raw = record.annotations["abif_raw"]
    base_order = raw.get("FWO_1", b"GATC").decode()

    # Extract processed traces
    channels = {}
    for i, base in enumerate(base_order):
        channels[base] = np.array(raw[f"DATA{9 + i}"])

    peak_locs = list(raw["PLOC1"])
    sequence = str(record.seq)

    if title is None:
        title = f"Chromatogram: {record.id}"

    plot_chromatogram(channels, peak_locs, sequence, start=start,
                      end=end if end else len(sequence), title=title)

# Uncomment with real data:
# plot_ab1_chromatogram(record, start=50, end=120)
```python

---

## 4. Phred Quality Scores and Base Calling

### 4.1 Phred Scores

Phred quality scores quantify the probability that a base call is **wrong**:

$$Q = -10 \times \log_{10}(P_{\text{error}})$$

| Phred Score | Error Probability | Accuracy |
|:-----------:|:-----------------:|:--------:|
| 10 | 1 in 10 | 90% |
| 20 | 1 in 100 | 99% |
| 30 | 1 in 1,000 | 99.9% |
| 40 | 1 in 10,000 | 99.99% |
| 50 | 1 in 100,000 | 99.999% |

**Practical guidelines:**
- Q < 20: unreliable -- consider trimming or re-sequencing
- Q 20-30: acceptable for most applications
- Q >= 30: high confidence

### 4.2 How Base Calling Works

Base-calling software (e.g. KB Basecaller, `phred`) performs:
1. **Baseline correction** -- subtract background fluorescence
2. **Color deconvolution** -- separate overlapping dye spectra ("spectral calibration")
3. **Peak detection** -- locate local maxima
4. **Spacing model** -- predict expected peak positions to handle compressions/expansions
5. **Quality assignment** -- evaluate peak shape, resolution, signal-to-noise to assign Phred scores

```python
def make_synthetic_qualities(n_bases, good_start=20, good_end=None, rng_seed=42):
    """
    Simulate realistic Phred quality profile for a Sanger read.
    Qualities ramp up at the start, plateau in the middle, and decay at the end.
    """
    if good_end is None:
        good_end = n_bases - 40
    rng = np.random.default_rng(rng_seed)

    quals = np.zeros(n_bases)
    for i in range(n_bases):
        if i < good_start:
            # Ramp-up region (primer / dye blobs)
            base_q = 8 + (30 - 8) * (i / good_start)
        elif i < good_end:
            # High-quality plateau
            base_q = 38 + 7 * rng.random()
        else:
            # Decay region
            frac = (i - good_end) / (n_bases - good_end)
            base_q = 38 * (1 - frac) + 5 * frac
        quals[i] = max(0, base_q + rng.normal(0, 3))

    return np.round(quals).astype(int)


# Simulate a 750-base read
synth_quals = make_synthetic_qualities(750)
print(f"Simulated {len(synth_quals)} quality scores")
print(f"Mean Q: {np.mean(synth_quals):.1f}  |  Median Q: {np.median(synth_quals):.0f}")
print(f"Min Q: {np.min(synth_quals)}  |  Max Q: {np.max(synth_quals)}")
```python

```python
def plot_quality_scores(qualities, title="Phred Quality Scores"):
    """
    Plot quality scores along the read with color-coded bars.
    Green >= Q30, orange Q20-29, red < Q20.
    """
    positions = np.arange(1, len(qualities) + 1)
    colors = ["red" if q < 20 else "orange" if q < 30 else "green" for q in qualities]

    fig, ax = plt.subplots(figsize=(13, 3.5))
    ax.bar(positions, qualities, color=colors, width=1.0, edgecolor="none")
    ax.axhline(y=20, color="orange", linestyle="--", alpha=0.7, label="Q20 threshold")
    ax.axhline(y=30, color="green", linestyle="--", alpha=0.7, label="Q30 threshold")

    ax.set_xlabel("Base position")
    ax.set_ylabel("Phred quality")
    ax.set_title(title)
    ax.set_xlim(0, len(qualities) + 1)
    ax.set_ylim(0, max(qualities) + 5)
    ax.legend(loc="lower right")
    plt.tight_layout()
    plt.show()

    # Summary
    low = sum(1 for q in qualities if q < 20)
    mid = sum(1 for q in qualities if 20 <= q < 30)
    high = sum(1 for q in qualities if q >= 30)
    total = len(qualities)
    print(f"Quality distribution:")
    print(f"  Low  (Q < 20):  {low:4d} bases ({100 * low / total:.1f}%)")
    print(f"  Med  (Q 20-29): {mid:4d} bases ({100 * mid / total:.1f}%)")
    print(f"  High (Q >= 30): {high:4d} bases ({100 * high / total:.1f}%)")


plot_quality_scores(synth_quals, "Synthetic Read -- Quality Profile")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
