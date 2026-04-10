# Course Extension Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the bioinformatics course with 10 new notebooks, 8 skill files, a new Tier 5, and updated READMEs — executed by 7 parallel Sonnet agents.

**Architecture:** Seven independent parallel tracks (Tasks 1–6 run in parallel; Task 7 runs after all content is created). Each track creates one or more notebooks + skill files following the existing course style. All content uses public datasets only — no real research data, no researcher names in notebooks.

**Tech Stack:** cooler · cooltools · pymc · bambi · arviz · statsmodels · pysam · pybedtools · squidpy · scanpy · transformers · trl · peft · bitsandbytes · torch · patsy · scipy · matplotlib · seaborn · numpy · pandas · palmerpenguins

---

## Pre-flight: Read These Files Before Starting Any Task

- `Skills/ngs-variant-calling.md` — existing skill format reference
- `Course/Tier_2_Core_Bioinformatics/01_Biological_Databases/01_biological_databases.ipynb` — existing notebook style
- `docs/superpowers/specs/2026-04-10-course-extension-design.md` — full spec

---

## PHASE 1 — Parallel (Tasks 1–6, run simultaneously)

### Task 1: Hi-C Analysis Module (Tier 2, Module 14)

**Files:**
- Create: `Course/Tier_2_Core_Bioinformatics/14_Hi-C_Analysis/14_hic_analysis.ipynb`
- Create: `Course/Tier_2_Core_Bioinformatics/14_Hi-C_Analysis/README.md`
- Create: `Skills/hic-analysis.md`

**Steps:**

- [ ] Step 1: Create the directory structure

```bash
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_2_Core_Bioinformatics/14_Hi-C_Analysis
```

- [ ] Step 2: Create `Skills/hic-analysis.md`

The file content:

```markdown
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
# matrix is a numpy ndarray; NaN where weight is missing
import numpy as np
matrix_log = np.log1p(matrix)  # log-transform for visualization
```

**Pattern 3: Compute contact decay (expected) curve**
```python
import cooltools

view_df = cooler.util.make_chromarms(
    clr.chromsizes, mid_point_flag=False
)
view_df.columns = ["chrom", "start", "end", "name"]

expected = cooltools.expected_cis(
    clr,
    view_df=view_df,
    nproc=4,
    chunksize=1_000_000,
)
# expected has columns: region1, region2, dist, n_valid, balanced.avg
```

**Pattern 4: A/B compartment eigenvectors**
```python
import bioframe

gc_cov = bioframe.load_fasta_sequence(  # or provide GC track separately
    "genome.fa", clr.bins()[["chrom","start","end"]]
)

eigenvalues, eigenvectors = cooltools.eigs_cis(
    clr,
    phasing_track=gc_cov,   # used to orient E1 sign
    view_df=view_df,
    nproc=4,
)
# eigenvectors["E1"] > 0 → A compartment, < 0 → B compartment
```

**Pattern 5: Insulation score and TAD boundaries**
```python
insulation_df = cooltools.insulation(
    clr,
    window_bp=[100_000, 200_000],  # list of window sizes to test
    view_df=view_df,
    nproc=4,
)
# Boundary bins have is_boundary_<window> == True
boundaries = insulation_df[
    insulation_df["is_boundary_100000"] == True
]
```

## Code Templates

**Template 1: Full cooler load + region matrix + heatmap**
```python
import cooler
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def plot_contact_matrix(cool_path: str, region: str, vmax: float = 0.05) -> None:
    clr = cooler.Cooler(cool_path)
    matrix = clr.matrix(balance=True).fetch(region)
    matrix = np.nan_to_num(matrix, nan=0.0)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.matshow(
        np.log1p(matrix),
        cmap="fall",
        norm=mcolors.PowerNorm(gamma=0.5, vmin=0, vmax=np.log1p(vmax)),
    )
    ax.set_title(f"Contact matrix — {region}", pad=12)
    ax.set_xlabel("Genomic bins")
    ax.set_ylabel("Genomic bins")
    plt.tight_layout()
    plt.show()
```

**Template 2: Contact decay log-log plot**
```python
import matplotlib.pyplot as plt

def plot_expected_decay(expected_df, chrom: str = "chr1") -> None:
    chrom_expected = expected_df[expected_df["region1"] == chrom].copy()
    chrom_expected = chrom_expected[chrom_expected["balanced.avg"] > 0]

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.loglog(
        chrom_expected["dist"] + 1,
        chrom_expected["balanced.avg"],
        lw=1.5,
        color="#2166ac",
    )
    ax.set_xlabel("Genomic distance (bins)")
    ax.set_ylabel("Mean balanced contact frequency")
    ax.set_title(f"Contact decay — {chrom}")
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.show()
```

**Template 3: Insulation score track plot with boundary calls**
```python
import matplotlib.pyplot as plt

def plot_insulation_track(
    insulation_df,
    chrom: str,
    start: int,
    end: int,
    window: int = 100_000,
) -> None:
    col_score = f"log2_insulation_score_{window}"
    col_boundary = f"is_boundary_{window}"

    region_df = insulation_df[
        (insulation_df["chrom"] == chrom)
        & (insulation_df["start"] >= start)
        & (insulation_df["end"] <= end)
    ].copy()

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(
        region_df["start"],
        region_df[col_score],
        lw=1.2,
        color="#333333",
    )
    boundaries = region_df[region_df[col_boundary] == True]
    ax.vlines(
        boundaries["start"],
        ymin=region_df[col_score].min(),
        ymax=region_df[col_score].max(),
        color="#d73027",
        lw=0.8,
        alpha=0.6,
        label="TAD boundary",
    )
    ax.set_xlabel("Genomic position (bp)")
    ax.set_ylabel("log2 insulation score")
    ax.set_title(f"Insulation score — {chrom}:{start}-{end}")
    ax.legend()
    plt.tight_layout()
    plt.show()
```

**Template 4: Synthetic cooler for testing without real data**
```python
import numpy as np
import pandas as pd
import cooler

def create_synthetic_cooler(path: str, n_bins: int = 200, resolution: int = 50_000) -> cooler.Cooler:
    bins = pd.DataFrame({
        "chrom": ["chr1"] * n_bins,
        "start": np.arange(n_bins) * resolution,
        "end": (np.arange(n_bins) + 1) * resolution,
    })
    rng = np.random.default_rng(42)
    # Simulate distance-dependent decay
    i_idx, j_idx = np.triu_indices(n_bins, k=0)
    distances = j_idx - i_idx
    values = rng.poisson(lam=np.maximum(0, 50 / (distances + 1)))
    pixels = pd.DataFrame({"bin1_id": i_idx, "bin2_id": j_idx, "count": values})
    pixels = pixels[pixels["count"] > 0]

    cooler.create_cooler(path, bins=bins, pixels=pixels, ordered=True)
    return cooler.Cooler(path)
```

## Common Pitfalls

1. **Unbalanced matrices**: always check `clr.bins()["weight"].isna().sum()` before using `balance=True`. If many bins have no weight, run `cooler balance` first or use `balance=False` with raw counts.

2. **Resolution mismatch with `.mcool`**: multi-resolution files require the `::resolutionN` suffix, e.g., `cooler.Cooler("file.mcool::10000")`. Omitting it raises a confusing error about groups.

3. **Chromosome name format**: `cooler` preserves the name format from the file. UCSC uses `chr1`; Ensembl uses `1`. Always check `clr.chromnames[0]` and prepend `"chr"` if needed before passing regions.

4. **Eigenvector sign ambiguity**: `cooltools.eigs_cis()` returns eigenvectors with arbitrary sign. Always orient E1 using a phasing track (GC content or gene density) — positive E1 should correspond to A compartments (high GC).

5. **NaN propagation in balanced matrices**: `clr.matrix(balance=True).fetch()` returns `NaN` for bins with missing weights. Convert with `np.nan_to_num(matrix, nan=0.0)` before any arithmetic or visualization.
```

- [ ] Step 3: Create `14_hic_analysis.ipynb`

The notebook JSON:

```json
{
 "nbformat": 4,
 "nbformat_minor": 5,
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "name": "python",
   "version": "3.11.0"
  }
 },
 "cells": [
  {
   "cell_type": "markdown",
   "id": "header-cell",
   "metadata": {},
   "source": [
    "# Module 14: Hi-C Analysis\n",
    "\n",
    "Workshop materials inspired by MPIB Hi-C Analysis Workshop\n",
    "\n",
    "---\n",
    "\n",
    "## Learning Objectives\n",
    "\n",
    "By the end of this notebook you will be able to:\n",
    "\n",
    "1. Load and inspect `.cool` contact matrix files using `cooler`\n",
    "2. Fetch and visualize submatrices for genomic regions of interest\n",
    "3. Compute and interpret the contact decay (expected) curve\n",
    "4. Detect A/B compartments using eigenvector decomposition\n",
    "5. Calculate insulation scores and identify TAD boundaries\n",
    "6. Create pileup (aggregate) plots around genomic features\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "imports-cell",
   "metadata": {},
   "outputs": [],
   "source": [
    "import warnings\n",
    "warnings.filterwarnings('ignore')\n",
    "\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt\n",
    "import matplotlib.colors as mcolors\n",
    "\n",
    "import cooler\n",
    "import cooltools\n",
    "import cooltools.lib.plotting  # registers 'fall' colormap\n",
    "\n",
    "print(f'cooler version:    {cooler.__version__}')\n",
    "print(f'cooltools version: {cooltools.__version__}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-cooler",
   "metadata": {},
   "source": [
    "## 1. The Cooler Format\n",
    "\n",
    "Hi-C contact matrices are stored as sparse upper-triangular tables in the **Cooler** format (`.cool` for a single resolution, `.mcool` for multi-resolution). The file is an HDF5 container with three core tables:\n",
    "\n",
    "| Table | Contents |\n",
    "|---|---|\n",
    "| `bins` | Genomic intervals (chrom, start, end) + optional balancing weight |\n",
    "| `pixels` | Non-zero contacts between bin pairs (bin1_id, bin2_id, count) |\n",
    "| `chroms` | Chromosome names and lengths |\n",
    "\n",
    "**Iterative correction (ICE) balancing** removes systematic biases (GC content, mappability, restriction site density) by computing a weight for each bin such that all rows/columns of the balanced matrix sum to the same value.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "create-synthetic-cooler",
   "metadata": {},
   "outputs": [],
   "source": [
    "# We use a synthetic cooler so no external download is required.\n",
    "# The simulation models distance-dependent contact decay (P(s) ~ 1/s)\n",
    "# with three TAD-like blocks of elevated contacts.\n",
    "\n",
    "import tempfile, os\n",
    "\n",
    "RESOLUTION = 50_000   # 50 kb bins\n",
    "N_BINS = 200          # 200 bins = 10 Mb on chr1\n",
    "\n",
    "rng = np.random.default_rng(42)\n",
    "\n",
    "bins = pd.DataFrame({\n",
    "    'chrom': ['chr1'] * N_BINS,\n",
    "    'start': np.arange(N_BINS) * RESOLUTION,\n",
    "    'end':   (np.arange(N_BINS) + 1) * RESOLUTION,\n",
    "})\n",
    "\n",
    "i_idx, j_idx = np.triu_indices(N_BINS, k=0)\n",
    "distances = j_idx - i_idx\n",
    "\n",
    "# Base decay: Poisson-sampled counts proportional to 1/(d+1)\n",
    "lam = 80.0 / (distances + 1)\n",
    "\n",
    "# Add three TAD blocks (bins 0-49, 50-109, 110-199)\n",
    "tad_boundaries = [0, 50, 110, N_BINS]\n",
    "for t in range(len(tad_boundaries) - 1):\n",
    "    lo, hi = tad_boundaries[t], tad_boundaries[t + 1]\n",
    "    mask = (i_idx >= lo) & (j_idx < hi)\n",
    "    lam[mask] *= 3.0   # 3x enrichment within TAD\n",
    "\n",
    "counts = rng.poisson(lam=np.maximum(lam, 0))\n",
    "\n",
    "pixels = pd.DataFrame({'bin1_id': i_idx, 'bin2_id': j_idx, 'count': counts})\n",
    "pixels = pixels[pixels['count'] > 0].reset_index(drop=True)\n",
    "\n",
    "COOL_PATH = os.path.join(tempfile.gettempdir(), 'synthetic_hic.cool')\n",
    "cooler.create_cooler(COOL_PATH, bins=bins, pixels=pixels, ordered=True)\n",
    "\n",
    "print(f'Cooler written to: {COOL_PATH}')\n",
    "print(f'Total non-zero pixels: {len(pixels):,}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-inspect",
   "metadata": {},
   "source": [
    "## 2. Loading and Inspecting a Cooler File\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "inspect-cooler",
   "metadata": {},
   "outputs": [],
   "source": [
    "clr = cooler.Cooler(COOL_PATH)\n",
    "\n",
    "print('=== Cooler info ===')\n",
    "for key, val in clr.info.items():\n",
    "    print(f'  {key}: {val}')\n",
    "\n",
    "print(f'\\nChromosomes: {clr.chromnames}')\n",
    "print(f'Resolution:  {clr.binsize:,} bp')\n",
    "print(f'Matrix shape: {clr.shape}')\n",
    "print(f'\\nBin table (first 5 rows):')\n",
    "display(clr.bins()[:5].to_pandas())"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-matrix-fetch",
   "metadata": {},
   "source": [
    "## 3. Fetching and Visualising a Contact Matrix\n",
    "\n",
    "`clr.matrix(balance=True).fetch(region)` returns a dense numpy array for the specified region. `balance=True` applies ICE weights so each bin's value is divided by its weight.\n",
    "\n",
    "The **'fall'** colormap (registered by `cooltools.lib.plotting`) is a standard Hi-C colormap ranging from white (low) through orange to dark red (high).\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "fetch-matrix",
   "metadata": {},
   "outputs": [],
   "source": [
    "region = 'chr1:0-10000000'   # first 10 Mb = all 200 bins at 50 kb\n",
    "\n",
    "# balance=False because our synthetic cooler has no ICE weights yet\n",
    "matrix = clr.matrix(balance=False).fetch(region).astype(float)\n",
    "\n",
    "fig, axes = plt.subplots(1, 2, figsize=(13, 5))\n",
    "\n",
    "axes[0].matshow(matrix, cmap='fall', vmin=0, vmax=np.percentile(matrix, 99))\n",
    "axes[0].set_title('Raw contact matrix', pad=10)\n",
    "axes[0].set_xlabel('Bin index')\n",
    "axes[0].set_ylabel('Bin index')\n",
    "\n",
    "axes[1].matshow(np.log1p(matrix), cmap='fall')\n",
    "axes[1].set_title('log1p-transformed contact matrix', pad=10)\n",
    "axes[1].set_xlabel('Bin index')\n",
    "axes[1].set_ylabel('Bin index')\n",
    "\n",
    "plt.suptitle(f'Region: {region}', y=1.01)\n",
    "plt.tight_layout()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-expected",
   "metadata": {},
   "source": [
    "## 4. Contact Decay Curve — P(s)\n",
    "\n",
    "The **contact decay** (or P(s) curve) describes how contact frequency decreases with increasing genomic distance s. In a log-log plot:\n",
    "\n",
    "- A straight line with slope ≈ -1 indicates a polymer-like (fractal globule) organisation.\n",
    "- Deviations from this slope at specific distances reveal TAD/compartment boundaries.\n",
    "\n",
    "`cooltools.expected_cis()` computes the mean balanced contact frequency at each genomic separation, averaging over all valid bin pairs on the same chromosome.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "compute-expected",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Build a view DataFrame specifying the genomic regions to analyse.\n",
    "# For our synthetic cooler the whole chromosome is one region.\n",
    "view_df = pd.DataFrame([\n",
    "    {'chrom': 'chr1', 'start': 0, 'end': N_BINS * RESOLUTION, 'name': 'chr1'}\n",
    "])\n",
    "\n",
    "expected = cooltools.expected_cis(\n",
    "    clr,\n",
    "    view_df=view_df,\n",
    "    ignore_diags=2,   # skip the main diagonal and nearest neighbour (self-ligation artefacts)\n",
    "    nproc=1,\n",
    ")\n",
    "\n",
    "print('Expected DataFrame columns:', expected.columns.tolist())\n",
    "display(expected.head())"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "plot-expected",
   "metadata": {},
   "outputs": [],
   "source": [
    "chr1_exp = expected[\n",
    "    (expected['region1'] == 'chr1') & (expected['count.avg'] > 0)\n",
    "].copy()\n",
    "\n",
    "# Convert bin distance to base-pair distance\n",
    "chr1_exp['dist_bp'] = chr1_exp['dist'] * RESOLUTION\n",
    "\n",
    "fig, ax = plt.subplots(figsize=(6, 5))\n",
    "ax.loglog(\n",
    "    chr1_exp['dist_bp'] + 1,\n",
    "    chr1_exp['count.avg'],\n",
    "    lw=1.8,\n",
    "    color='#2166ac',\n",
    "    label='P(s)',\n",
    ")\n",
    "ax.set_xlabel('Genomic distance s (bp)', fontsize=12)\n",
    "ax.set_ylabel('Mean contact frequency', fontsize=12)\n",
    "ax.set_title('Contact decay curve — chr1', fontsize=13)\n",
    "ax.grid(True, which='both', alpha=0.3)\n",
    "ax.legend()\n",
    "plt.tight_layout()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-compartments",
   "metadata": {},
   "source": [
    "## 5. A/B Compartments via Eigenvector Decomposition\n",
    "\n",
    "Chromatin is spatially segregated into:\n",
    "- **A compartments** — open, transcriptionally active, high GC content\n",
    "- **B compartments** — closed, repressed, low GC content\n",
    "\n",
    "The first eigenvector (E1) of the **observed/expected** contact matrix separates these two populations. The sign of E1 is arbitrary and must be oriented using a phasing track such as GC content (A compartment bins have high GC → E1 should be positive there).\n",
    "\n",
    "`cooltools.eigs_cis()` returns eigenvalues and an eigenvector table per chromosome.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "compute-eigenvectors",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Simulate a GC-content phasing track.\n",
    "# In real analysis this comes from reference genome sequence.\n",
    "# A compartment (first and last TAD blocks) → high GC; B compartment (middle) → low GC.\n",
    "gc_track = clr.bins()[['chrom', 'start', 'end']].to_pandas().copy()\n",
    "gc_values = np.full(N_BINS, 0.38)   # B compartment baseline\n",
    "gc_values[:50]    = 0.55            # first TAD block → A compartment\n",
    "gc_values[110:]   = 0.52            # last TAD block → A compartment\n",
    "gc_track['GC'] = gc_values\n",
    "\n",
    "eigenvalues, eigenvectors = cooltools.eigs_cis(\n",
    "    clr,\n",
    "    phasing_track=gc_track,          # used to orient E1 sign\n",
    "    view_df=view_df,\n",
    "    n_eigs=3,\n",
    "    ignore_diags=2,\n",
    ")\n",
    "\n",
    "print('Eigenvalue table:')\n",
    "display(eigenvalues)\n",
    "print('\\nEigenvector table (first 5 rows):')\n",
    "display(eigenvectors.head())"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "plot-eigenvectors",
   "metadata": {},
   "outputs": [],
   "source": [
    "e1 = eigenvectors['E1'].values\n",
    "bin_centers = (eigenvectors['start'].values + eigenvectors['end'].values) / 2\n",
    "\n",
    "fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True)\n",
    "\n",
    "# Panel 1: GC content used for phasing\n",
    "axes[0].fill_between(bin_centers, gc_values, 0.38, alpha=0.7, color='#4dac26')\n",
    "axes[0].axhline(0.38, color='grey', lw=0.8, ls='--')\n",
    "axes[0].set_ylabel('GC content', fontsize=11)\n",
    "axes[0].set_title('GC content phasing track', fontsize=12)\n",
    "\n",
    "# Panel 2: E1 coloured by compartment\n",
    "pos_mask = e1 >= 0\n",
    "axes[1].fill_between(bin_centers, e1, 0, where=pos_mask,  color='#d73027', alpha=0.8, label='A compartment (E1 > 0)')\n",
    "axes[1].fill_between(bin_centers, e1, 0, where=~pos_mask, color='#4575b4', alpha=0.8, label='B compartment (E1 < 0)')\n",
    "axes[1].axhline(0, color='black', lw=0.8)\n",
    "axes[1].set_ylabel('E1', fontsize=11)\n",
    "axes[1].set_xlabel('Genomic position (bp)', fontsize=11)\n",
    "axes[1].set_title('First eigenvector — A/B compartments', fontsize=12)\n",
    "axes[1].legend(loc='upper right', fontsize=9)\n",
    "\n",
    "plt.tight_layout()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-insulation",
   "metadata": {},
   "source": [
    "## 6. Insulation Score and TAD Boundaries\n",
    "\n",
    "The **insulation score** measures how many contacts cross each bin by summing contacts in a sliding diamond window centred on each bin:\n",
    "\n",
    "$$\\text{IS}(i) = \\log_2 \\frac{\\text{contacts in window at } i}{\\text{mean contacts across all windows}}$$\n",
    "\n",
    "- **Minima** of the insulation score correspond to TAD boundaries (few contacts cross that bin).\n",
    "- The `window_bp` parameter controls the size of the diamond; typical values are 100–500 kb.\n",
    "\n",
    "`cooltools.insulation()` returns a DataFrame with one row per bin, including the insulation score and a Boolean boundary call column for each window size.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "compute-insulation",
   "metadata": {},
   "outputs": [],
   "source": [
    "insulation_df = cooltools.insulation(\n",
    "    clr,\n",
    "    window_bp=[100_000, 200_000],\n",
    "    view_df=view_df,\n",
    "    ignore_diags=2,\n",
    "    min_dist_bad_bin=0,\n",
    ")\n",
    "\n",
    "print('Insulation score DataFrame columns:')\n",
    "print(insulation_df.columns.tolist())\n",
    "\n",
    "# Count boundary calls at each window size\n",
    "for w in [100_000, 200_000]:\n",
    "    col = f'is_boundary_{w}'\n",
    "    n = insulation_df[col].sum()\n",
    "    print(f'  Boundaries at window={w:,} bp: {int(n)}')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "plot-insulation",
   "metadata": {},
   "outputs": [],
   "source": [
    "WINDOW = 100_000\n",
    "score_col = f'log2_insulation_score_{WINDOW}'\n",
    "boundary_col = f'is_boundary_{WINDOW}'\n",
    "\n",
    "track = insulation_df[insulation_df['chrom'] == 'chr1'].copy()\n",
    "bin_centers = (track['start'].values + track['end'].values) / 2\n",
    "scores = track[score_col].values\n",
    "boundaries = track[track[boundary_col] == True]\n",
    "\n",
    "fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True,\n",
    "                          gridspec_kw={'height_ratios': [2, 1]})\n",
    "\n",
    "# Top panel: contact matrix\n",
    "matrix = clr.matrix(balance=False).fetch('chr1').astype(float)\n",
    "axes[0].matshow(\n",
    "    np.log1p(matrix), cmap='fall', aspect='auto',\n",
    "    extent=[0, N_BINS * RESOLUTION, N_BINS * RESOLUTION, 0],\n",
    ")\n",
    "for _, row in boundaries.iterrows():\n",
    "    axes[0].axvline(row['start'], color='cyan', lw=1.2, alpha=0.8)\n",
    "axes[0].set_title('Contact matrix with TAD boundaries (cyan)', fontsize=12)\n",
    "axes[0].set_ylabel('Genomic position (bp)')\n",
    "\n",
    "# Bottom panel: insulation score\n",
    "axes[1].plot(bin_centers, scores, lw=1.4, color='#333333')\n",
    "axes[1].vlines(\n",
    "    boundaries['start'].values,\n",
    "    ymin=np.nanmin(scores), ymax=np.nanmax(scores),\n",
    "    color='#d73027', lw=1.0, alpha=0.7, label='TAD boundary',\n",
    ")\n",
    "axes[1].axhline(0, color='grey', lw=0.8, ls='--')\n",
    "axes[1].set_ylabel(f'log2 IS ({WINDOW//1000} kb window)', fontsize=10)\n",
    "axes[1].set_xlabel('Genomic position (bp)', fontsize=11)\n",
    "axes[1].legend()\n",
    "\n",
    "plt.tight_layout()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-pileup",
   "metadata": {},
   "source": [
    "## 7. Pileup (Aggregate) Plots\n",
    "\n",
    "A **pileup** aligns and averages contact sub-matrices centred on a set of genomic features (e.g. CTCF loop anchors, TSS pairs). The resulting average matrix reveals subtle enrichment patterns invisible in individual loci.\n",
    "\n",
    "The workflow:\n",
    "1. Define a feature table (pairs of genomic coordinates).\n",
    "2. For each feature, extract a fixed-size matrix snippet.\n",
    "3. Average all snippets.\n",
    "\n",
    "Here we simulate CTCF-like loop anchor pairs and aggregate them manually without requiring external BED files.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "pileup-plot",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Simulate loop anchor pairs at TAD boundaries\n",
    "# In real analysis these come from ENCODE CTCF ChIP-seq peaks\n",
    "anchor_pairs = [\n",
    "    (0, 49),     # TAD 1 anchors (bin indices)\n",
    "    (50, 109),   # TAD 2 anchors\n",
    "    (110, 199),  # TAD 3 anchors\n",
    "]\n",
    "\n",
    "FLANK = 10   # bins on each side of each anchor\n",
    "snippets = []\n",
    "\n",
    "full_matrix = clr.matrix(balance=False).fetch('chr1').astype(float)\n",
    "\n",
    "for anchor1, anchor2 in anchor_pairs:\n",
    "    r1_lo = max(0, anchor1 - FLANK)\n",
    "    r1_hi = min(N_BINS, anchor1 + FLANK + 1)\n",
    "    r2_lo = max(0, anchor2 - FLANK)\n",
    "    r2_hi = min(N_BINS, anchor2 + FLANK + 1)\n",
    "\n",
    "    snippet = full_matrix[r1_lo:r1_hi, r2_lo:r2_hi]\n",
    "    # Resize to (2*FLANK+1, 2*FLANK+1) with padding\n",
    "    padded = np.zeros((2 * FLANK + 1, 2 * FLANK + 1))\n",
    "    sh, sw = snippet.shape\n",
    "    padded[:sh, :sw] = snippet\n",
    "    snippets.append(padded)\n",
    "\n",
    "pileup = np.mean(snippets, axis=0)\n",
    "\n",
    "fig, ax = plt.subplots(figsize=(5, 5))\n",
    "im = ax.matshow(np.log1p(pileup), cmap='fall')\n",
    "plt.colorbar(im, ax=ax, label='log1p(contact count)')\n",
    "ax.set_title('Pileup around TAD anchor pairs', pad=12)\n",
    "ax.set_xlabel('Bins relative to anchor 2')\n",
    "ax.set_ylabel('Bins relative to anchor 1')\n",
    "# Mark the central pixel (the loop point)\n",
    "ax.plot(FLANK, FLANK, 'c+', markersize=12, markeredgewidth=2, label='Loop point')\n",
    "ax.legend(loc='upper right', fontsize=8)\n",
    "plt.tight_layout()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "exercise-header",
   "metadata": {},
   "source": [
    "## Exercise\n",
    "\n",
    "Using the synthetic cooler created above, complete the following tasks:\n",
    "\n",
    "1. **Extend the cooler to include a chr2** with 150 bins (7.5 Mb at 50 kb resolution) and a different TAD structure: two TADs at bins 0–74 and 75–149.\n",
    "2. **Compute the insulation score on chr2** with a 100 kb window.\n",
    "3. **Identify the top 10 boundary positions on chr2** by boundary strength (use the `boundary_strength_100000` column, or rank by the depth of the insulation score minima).\n",
    "4. **Plot the insulation score track for chr2** with the top 10 boundaries highlighted in red.\n",
    "\n",
    "_Hint_: rebuild the cooler with both chromosomes using the same `cooler.create_cooler()` function, then re-run `cooltools.insulation()` with an updated `view_df` that includes both chromosomes.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "exercise-cell",
   "metadata": {},
   "outputs": [],
   "source": [
    "# --- Your code here ---\n",
    "\n",
    "# Step 1: Build a two-chromosome cooler\n",
    "N_BINS_CHR2 = 150\n",
    "COOL_PATH_2 = os.path.join(tempfile.gettempdir(), 'synthetic_hic_two_chroms.cool')\n",
    "\n",
    "# bins for chr1 (reuse) + chr2\n",
    "# ...\n",
    "\n",
    "# Step 2: Compute insulation on chr2\n",
    "# ...\n",
    "\n",
    "# Step 3: Identify top 10 boundaries\n",
    "# ...\n",
    "\n",
    "# Step 4: Plot\n",
    "# ..."
   ]
  }
 ]
}
```

- [ ] Step 4: Create `14_Hi-C_Analysis/README.md`

```markdown
# Module 14: Hi-C Analysis

**Tier 2 — Core Bioinformatics | Module 14**

## Overview

This module covers the analysis of Hi-C chromatin conformation data using the `cooler` and `cooltools` Python libraries. Hi-C measures genome-wide pairwise contact frequencies between all pairs of genomic loci, revealing the three-dimensional organisation of chromatin.

## Topics

1. The Cooler format — bins, pixels, balancing
2. Loading and inspecting `.cool` files
3. Fetching and visualising contact matrices
4. Contact decay (P(s)) curves
5. A/B compartment detection via eigenvector decomposition
6. TAD boundary calling with insulation score
7. Pileup/aggregate plots around genomic features

## Prerequisites

- Module 11 (ChIP-seq) — familiarity with genomic interval formats
- NumPy, pandas, matplotlib

## Installation

```bash
pip install cooler cooltools bioframe
```

## Key Libraries

| Library | Purpose |
|---|---|
| `cooler` | Read/write `.cool` and `.mcool` files |
| `cooltools` | Hi-C analysis algorithms (expected, eigs, insulation, pileup) |
| `bioframe` | Genomic interval operations, chromosome size utilities |

## Learning Objectives

After completing this module you will be able to:

- Load and inspect `.cool` contact matrix files
- Fetch balanced submatrices and visualise them
- Compute and interpret contact decay curves
- Detect A/B compartments using eigenvector decomposition
- Calculate insulation scores and call TAD boundaries
- Build pileup plots around genomic features

## Skill File

See `Skills/hic-analysis.md` for quick-reference patterns and copy-paste code templates.
```

- [ ] Step 5: Commit

```bash
git add Course/Tier_2_Core_Bioinformatics/14_Hi-C_Analysis/14_hic_analysis.ipynb \
        Course/Tier_2_Core_Bioinformatics/14_Hi-C_Analysis/README.md \
        Skills/hic-analysis.md
git commit -m "feat(tier2): add Hi-C analysis module (14) and skill"
```

---

### Task 2: Motif Discovery Module (Tier 2, Module 15)

**Files:**
- Create: `Course/Tier_2_Core_Bioinformatics/15_Motif_Discovery/15_motif_discovery.ipynb`
- Create: `Course/Tier_2_Core_Bioinformatics/15_Motif_Discovery/README.md`
- Create: `Skills/motif-discovery.md`

**Steps:**

- [ ] Step 1: Create the directory structure

```bash
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_2_Core_Bioinformatics/15_Motif_Discovery
```

- [ ] Step 2: Create `Skills/motif-discovery.md`

```markdown
---
name: motif-discovery
description: Transcription factor motif discovery — build PWMs, compute information content, score sequences, test enrichment, query JASPAR REST API
---

## When to Use

Use this skill when:
- Discovering over-represented sequence motifs in ChIP-seq or ATAC-seq peak sequences
- Building or loading Position Weight Matrices (PWMs) for TF binding site prediction
- Scoring genomic sequences against known motifs (log-odds scoring)
- Testing motif enrichment in foreground vs background sequences (Fisher's exact test)
- Querying JASPAR or HOCOMOCO databases for known TF motifs
- Computing information content to rank or filter motifs

## Quick Reference

| Concept / Tool | Formula / URL | Notes |
|---|---|---|
| PPM shape convention | `(L, 4)` — rows = positions, columns = [A,C,G,T] | Never use `(4, L)` — causes silent transposition bugs |
| PWM (log-odds) | `log2(PPM / background_freq)` | Typical background: uniform 0.25 each |
| Information content per position | `IC_i = 2 + sum_b(p_b * log2(p_b))` | Max 2 bits; pseudocount before computing |
| KDIC | `sum(IC_i) / (2 * L)` | Normalised to [0,1]; 1 = perfect specificity |
| JASPAR 2024 API | `https://jaspar.elixir.no/api/v1/matrix/{matrix_id}/` | Returns JSON with pfm field |
| HOCOMOCO v12 | `https://hocomoco12.autosome.org/` | MEME format download |
| MEME Suite | `meme`, `tomtom`, `fimo` | CLI tools for de novo + comparison |
| HOMER | `findMotifsGenome.pl` | ChIP-seq–specific pipeline |
| STREME | `streme --p foreground.fa --n background.fa` | Modern MEME-Suite replacement |

## Key Patterns

**Pattern 1: Build a PPM (Position Probability Matrix) from aligned sequences**
```python
import numpy as np

BASES = ['A', 'C', 'G', 'T']
BASE_INDEX = {b: i for i, b in enumerate(BASES)}

def build_ppm(sequences: list[str], pseudocount: float = 0.5) -> np.ndarray:
    """Return PPM of shape (L, 4) with columns ordered [A,C,G,T]."""
    L = len(sequences[0])
    counts = np.full((L, 4), pseudocount)
    for seq in sequences:
        for pos, base in enumerate(seq.upper()):
            if base in BASE_INDEX:
                counts[pos, BASE_INDEX[base]] += 1
    ppm = counts / counts.sum(axis=1, keepdims=True)
    return ppm   # shape (L, 4)
```

**Pattern 2: Compute information content per position and KDIC**
```python
def information_content(ppm: np.ndarray, bg: np.ndarray | None = None) -> np.ndarray:
    """Return IC array of shape (L,) in bits. bg defaults to uniform 0.25."""
    if bg is None:
        bg = np.full(4, 0.25)
    # Avoid log(0): ppm already has pseudocounts, but clip for safety
    ppm_safe = np.clip(ppm, 1e-9, 1.0)
    ic = np.sum(ppm_safe * np.log2(ppm_safe / bg), axis=1)  # shape (L,)
    return ic

def kdic(ppm: np.ndarray) -> float:
    """Normalised information content in [0, 1]."""
    ic = information_content(ppm)
    L = ppm.shape[0]
    return float(ic.sum() / (2.0 * L))
```

**Pattern 3: Score a sequence with a PWM (log-odds)**
```python
def build_pwm(ppm: np.ndarray, bg: np.ndarray | None = None) -> np.ndarray:
    """Return PWM of shape (L, 4) in log2 log-odds."""
    if bg is None:
        bg = np.full(4, 0.25)
    return np.log2(np.clip(ppm, 1e-9, 1.0) / bg)

def score_sequence(seq: str, pwm: np.ndarray) -> float:
    """Return the log-odds score of seq against pwm."""
    score = 0.0
    for pos, base in enumerate(seq.upper()):
        if base in BASE_INDEX:
            score += pwm[pos, BASE_INDEX[base]]
    return score
```

**Pattern 4: Fisher's exact test for motif enrichment**
```python
from scipy.stats import fisher_exact

def motif_enrichment_fisher(
    n_fg_with: int, n_fg_total: int,
    n_bg_with: int, n_bg_total: int,
) -> tuple[float, float]:
    """Return (odds_ratio, p_value) for motif enrichment."""
    table = [
        [n_fg_with,             n_fg_total - n_fg_with],
        [n_bg_with,             n_bg_total - n_bg_with],
    ]
    odds_ratio, p_value = fisher_exact(table, alternative='greater')
    return odds_ratio, p_value
```

**Pattern 5: Benjamini-Hochberg multiple testing correction**
```python
from scipy.stats import false_discovery_control

def apply_bh_correction(p_values: list[float], alpha: float = 0.05) -> np.ndarray:
    """Return BH-adjusted p-values. Requires scipy >= 1.11."""
    adjusted = false_discovery_control(p_values, method='bh')
    return adjusted
```

**Pattern 6: Query the JASPAR REST API**
```python
import requests

def fetch_jaspar_motif(matrix_id: str) -> dict:
    """Fetch a JASPAR motif by ID and return a dict with ppm (L,4) and name."""
    url = f'https://jaspar.elixir.no/api/v1/matrix/{matrix_id}/'
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    data = response.json()
    # data['pfm'] is a dict {'A': [...], 'C': [...], 'G': [...], 'T': [...]}
    pfm = data['pfm']
    L = len(pfm['A'])
    counts = np.array([pfm['A'], pfm['C'], pfm['G'], pfm['T']], dtype=float).T  # (L,4)
    row_sums = counts.sum(axis=1, keepdims=True)
    ppm = counts / np.where(row_sums == 0, 1, row_sums)
    return {'name': data['name'], 'matrix_id': matrix_id, 'ppm': ppm}
```

## Code Templates

**Template 1: Full motif discovery pipeline from sequences**
```python
import numpy as np
from scipy.stats import fisher_exact, false_discovery_control

BASES = ['A', 'C', 'G', 'T']
BASE_INDEX = {b: i for i, b in enumerate(BASES)}

def build_ppm(sequences: list[str], pseudocount: float = 0.5) -> np.ndarray:
    L = len(sequences[0])
    counts = np.full((L, 4), pseudocount)
    for seq in sequences:
        for pos, base in enumerate(seq.upper()):
            if base in BASE_INDEX:
                counts[pos, BASE_INDEX[base]] += 1
    return counts / counts.sum(axis=1, keepdims=True)

def build_pwm(ppm: np.ndarray, bg: np.ndarray | None = None) -> np.ndarray:
    if bg is None:
        bg = np.full(4, 0.25)
    return np.log2(np.clip(ppm, 1e-9, 1.0) / bg)

def scan_sequences(sequences: list[str], pwm: np.ndarray) -> list[float]:
    """Return max PWM score across all positions for each sequence."""
    L = pwm.shape[0]
    max_scores = []
    for seq in sequences:
        seq = seq.upper()
        scores = [
            sum(pwm[i, BASE_INDEX.get(seq[pos + i], 0)] for i in range(L))
            for pos in range(len(seq) - L + 1)
            if all(seq[pos + i] in BASE_INDEX for i in range(L))
        ]
        max_scores.append(max(scores) if scores else float('-inf'))
    return max_scores

def run_enrichment_pipeline(
    fg_seqs: list[str],
    bg_seqs: list[str],
    pwm: np.ndarray,
    score_threshold: float,
) -> dict:
    fg_scores = scan_sequences(fg_seqs, pwm)
    bg_scores = scan_sequences(bg_seqs, pwm)
    n_fg_with = sum(s >= score_threshold for s in fg_scores)
    n_bg_with = sum(s >= score_threshold for s in bg_scores)
    table = [
        [n_fg_with, len(fg_seqs) - n_fg_with],
        [n_bg_with, len(bg_seqs) - n_bg_with],
    ]
    odds_ratio, p_value = fisher_exact(table, alternative='greater')
    return {
        'n_fg_with': n_fg_with, 'n_fg_total': len(fg_seqs),
        'n_bg_with': n_bg_with, 'n_bg_total': len(bg_seqs),
        'odds_ratio': odds_ratio, 'p_value': p_value,
    }
```

**Template 2: Exact score distribution for short motifs (L ≤ 8)**
```python
import itertools

def exact_score_distribution(pwm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Enumerate all 4^L k-mers and compute their PWM scores.
    Returns (scores array, frequency array) assuming uniform background.
    Only feasible for L <= 8 (max 65536 k-mers).
    """
    L = pwm.shape[0]
    assert L <= 8, f'L={L} is too large for exhaustive enumeration (max 8)'
    all_scores = []
    for kmer in itertools.product(range(4), repeat=L):
        score = sum(pwm[pos, base] for pos, base in enumerate(kmer))
        all_scores.append(score)
    all_scores = np.array(all_scores)
    return all_scores

def p_value_from_exact_distribution(score: float, all_scores: np.ndarray) -> float:
    """Fraction of k-mers scoring >= score under uniform background."""
    return float((all_scores >= score).mean())
```

**Template 3: IUPAC consensus from PPM**
```python
IUPAC = {
    frozenset('A'):    'A', frozenset('C'):    'C',
    frozenset('G'):    'G', frozenset('T'):    'T',
    frozenset('AC'):   'M', frozenset('AG'):   'R',
    frozenset('AT'):   'W', frozenset('CG'):   'S',
    frozenset('CT'):   'Y', frozenset('GT'):   'K',
    frozenset('ACG'):  'V', frozenset('ACT'):  'H',
    frozenset('AGT'):  'D', frozenset('CGT'):  'B',
    frozenset('ACGT'): 'N',
}

def ppm_to_iupac(ppm: np.ndarray, threshold: float = 0.25) -> str:
    """Return IUPAC consensus string. Positions with p >= threshold are included."""
    consensus = []
    for row in ppm:
        present = frozenset(b for b, p in zip(BASES, row) if p >= threshold)
        consensus.append(IUPAC.get(present, 'N'))
    return ''.join(consensus)
```

## Common Pitfalls

1. **PPM shape ambiguity**: always store PPMs as `(L, 4)` — rows are positions, columns are [A,C,G,T] in alphabetical order. Using `(4, L)` causes silent transposition when you index `ppm[pos, base_idx]`. Always verify with `assert ppm.shape[1] == 4`.

2. **Background frequency assumption**: the default uniform background (0.25 each) is wrong for GC-rich genomes. For human genome, use `bg = [0.295, 0.205, 0.205, 0.295]` (ACGT order matching PPM columns).

3. **Multiple testing**: when scanning thousands of sequences or testing many candidate motifs, always apply Benjamini-Hochberg correction. A raw p-value of 0.01 across 1000 tests expects 10 false positives.

4. **Monte Carlo convergence for long motifs**: exhaustive enumeration is only feasible for L ≤ 8. For longer motifs, use `scipy.stats.norm` approximation or dynamic programming (MOODS library).

5. **Pseudocount magnitude**: a pseudocount of 0.5 (Laplace smoothing) works well for motifs built from ≥ 20 sequences. For fewer sequences, use a larger pseudocount or a Dirichlet prior — a zero count position makes the entire sequence score −∞.
```

- [ ] Step 3: Create `15_motif_discovery.ipynb`

```json
{
 "nbformat": 4,
 "nbformat_minor": 5,
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "name": "python",
   "version": "3.11.0"
  }
 },
 "cells": [
  {
   "cell_type": "markdown",
   "id": "header-cell",
   "metadata": {},
   "source": [
    "# Module 15: Motif Discovery\n",
    "\n",
    "Code patterns inspired by TotipotencyLab Motif Discovery Pipeline (MIT License), Max Planck Institute of Biochemistry\n",
    "\n",
    "---\n",
    "\n",
    "## Learning Objectives\n",
    "\n",
    "By the end of this notebook you will be able to:\n",
    "\n",
    "1. Build Position Probability Matrices (PPMs) and Position Weight Matrices (PWMs) from aligned sequences\n",
    "2. Compute per-position information content and the KDIC summary score\n",
    "3. Score arbitrary sequences against a PWM using log-odds\n",
    "4. Test motif enrichment with Fisher's exact test and correct for multiple testing with Benjamini-Hochberg\n",
    "5. Query the JASPAR 2024 REST API to retrieve known TF motifs\n",
    "6. Understand TomTom-style motif comparison using cosine similarity between PPMs\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "imports-cell",
   "metadata": {},
   "outputs": [],
   "source": [
    "import warnings\n",
    "warnings.filterwarnings('ignore')\n",
    "\n",
    "import itertools\n",
    "import requests\n",
    "\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt\n",
    "import matplotlib.patches as mpatches\n",
    "from scipy.stats import fisher_exact, false_discovery_control\n",
    "\n",
    "# Alphabet convention used throughout this notebook\n",
    "BASES = ['A', 'C', 'G', 'T']\n",
    "BASE_INDEX = {b: i for i, b in enumerate(BASES)}\n",
    "\n",
    "print('All imports successful')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-pwm",
   "metadata": {},
   "source": [
    "## 1. From Sequences to Position Probability Matrix\n",
    "\n",
    "A **Position Probability Matrix (PPM)** — also called a Position Frequency Matrix normalised to probabilities — stores, for each position i in the motif, the probability of observing each nucleotide.\n",
    "\n",
    "**Shape convention**: `(L, 4)` where L = motif length and 4 columns correspond to [A, C, G, T] in alphabetical order.\n",
    "\n",
    "> Critical: always store as (L, 4), never (4, L). Transposing silently produces wrong scores when you index `ppm[position, base_index]`.\n",
    "\n",
    "We add a **pseudocount** of 0.5 (Laplace smoothing) to avoid zero probabilities, which would produce −∞ in log-odds scoring.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "build-ppm",
   "metadata": {},
   "outputs": [],
   "source": [
    "def build_ppm(sequences: list[str], pseudocount: float = 0.5) -> np.ndarray:\n",
    "    \"\"\"\n",
    "    Build a PPM of shape (L, 4) from a list of same-length sequences.\n",
    "    Columns are ordered [A, C, G, T] (alphabetical).\n",
    "    \"\"\"\n",
    "    L = len(sequences[0])\n",
    "    assert all(len(s) == L for s in sequences), 'All sequences must have the same length'\n",
    "    counts = np.full((L, 4), pseudocount)\n",
    "    for seq in sequences:\n",
    "        for pos, base in enumerate(seq.upper()):\n",
    "            if base in BASE_INDEX:\n",
    "                counts[pos, BASE_INDEX[base]] += 1\n",
    "    ppm = counts / counts.sum(axis=1, keepdims=True)\n",
    "    assert ppm.shape == (L, 4), f'Expected shape ({L}, 4), got {ppm.shape}'\n",
    "    return ppm\n",
    "\n",
    "\n",
    "# Demonstrate with the AP-1 consensus TGASTCA (S = C or G)\n",
    "ap1_sequences = [\n",
    "    'TGACTCA', 'TGAGTCA', 'TGACTCA', 'TGAGTCA',\n",
    "    'TGACTCA', 'TGAGTCA', 'TGACTCA', 'TGAGTCA',\n",
    "]\n",
    "\n",
    "ppm = build_ppm(ap1_sequences)\n",
    "\n",
    "ppm_df = pd.DataFrame(ppm, columns=BASES)\n",
    "ppm_df.index.name = 'position'\n",
    "print('PPM shape:', ppm.shape, '  (L, 4) ✓')\n",
    "display(ppm_df.round(3))"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-logo",
   "metadata": {},
   "source": [
    "## 2. Sequence Logo — Visualising the PPM\n",
    "\n",
    "A sequence logo represents each position as a stack of letters whose height is proportional to the letter's probability, scaled by the position's information content (total height = IC in bits).\n",
    "\n",
    "We implement a simplified logo using `matplotlib` rectangles without requiring the `logomaker` library.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "sequence-logo",
   "metadata": {},
   "outputs": [],
   "source": [
    "BASE_COLORS = {'A': '#2ca02c', 'C': '#1f77b4', 'G': '#ff7f0e', 'T': '#d62728'}\n",
    "\n",
    "def plot_sequence_logo(ppm: np.ndarray, title: str = 'Sequence Logo') -> None:\n",
    "    \"\"\"Draw a simplified sequence logo using matplotlib bar chart.\"\"\"\n",
    "    L = ppm.shape[0]\n",
    "    bg = np.full(4, 0.25)\n",
    "    ppm_safe = np.clip(ppm, 1e-9, 1.0)\n",
    "    ic_per_pos = np.sum(ppm_safe * np.log2(ppm_safe / bg), axis=1)\n",
    "    heights = ppm * ic_per_pos[:, np.newaxis]  # (L, 4)\n",
    "\n",
    "    fig, ax = plt.subplots(figsize=(max(6, L * 0.8), 3))\n",
    "    for pos in range(L):\n",
    "        order = np.argsort(heights[pos])   # smallest first (bottom of stack)\n",
    "        bottom = 0.0\n",
    "        for base_idx in order:\n",
    "            h = heights[pos, base_idx]\n",
    "            ax.bar(\n",
    "                pos, h, bottom=bottom, color=BASE_COLORS[BASES[base_idx]],\n",
    "                width=0.8, edgecolor='white', linewidth=0.3,\n",
    "            )\n",
    "            if h > 0.08:\n",
    "                ax.text(\n",
    "                    pos, bottom + h / 2, BASES[base_idx],\n",
    "                    ha='center', va='center', fontsize=9, fontweight='bold', color='white',\n",
    "                )\n",
    "            bottom += h\n",
    "\n",
    "    legend_patches = [mpatches.Patch(color=BASE_COLORS[b], label=b) for b in BASES]\n",
    "    ax.legend(handles=legend_patches, loc='upper right', fontsize=8, ncol=4)\n",
    "    ax.set_xlim(-0.5, L - 0.5)\n",
    "    ax.set_ylim(0, 2.1)\n",
    "    ax.set_xticks(range(L))\n",
    "    ax.set_xticklabels([str(i + 1) for i in range(L)])\n",
    "    ax.set_ylabel('Information content (bits)')\n",
    "    ax.set_title(title)\n",
    "    plt.tight_layout()\n",
    "    plt.show()\n",
    "\n",
    "\n",
    "plot_sequence_logo(ppm, title='AP-1 (TGASTCA) motif logo')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-ic",
   "metadata": {},
   "source": [
    "## 3. Information Content and KDIC\n",
    "\n",
    "**Per-position information content** measures how specific the motif is at each position, relative to a background frequency:\n",
    "\n",
    "$$\\text{IC}_i = \\sum_{b \\in \\{A,C,G,T\\}} p_{i,b} \\log_2 \\frac{p_{i,b}}{q_b}$$\n",
    "\n",
    "where $q_b$ is the background frequency for nucleotide $b$. With uniform background ($q_b = 0.25$), a perfectly conserved position has IC = 2 bits; a random position has IC = 0 bits.\n",
    "\n",
    "**KDIC** (KL-divergence information content) normalises the sum of IC values to [0, 1]:\n",
    "\n",
    "$$\\text{KDIC} = \\frac{\\sum_i \\text{IC}_i}{2 \\cdot L}$$\n",
    "\n",
    "A KDIC of 1.0 means every position is perfectly conserved; 0.0 means the matrix is uniform (no information).\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "compute-ic",
   "metadata": {},
   "outputs": [],
   "source": [
    "def information_content(ppm: np.ndarray, bg: np.ndarray | None = None) -> np.ndarray:\n",
    "    \"\"\"Return IC array of shape (L,) in bits. bg defaults to uniform 0.25.\"\"\"\n",
    "    if bg is None:\n",
    "        bg = np.full(4, 0.25)\n",
    "    ppm_safe = np.clip(ppm, 1e-9, 1.0)\n",
    "    return np.sum(ppm_safe * np.log2(ppm_safe / bg), axis=1)\n",
    "\n",
    "\n",
    "def kdic(ppm: np.ndarray) -> float:\n",
    "    \"\"\"Return KDIC normalised to [0, 1].\"\"\"\n",
    "    ic = information_content(ppm)\n",
    "    L = ppm.shape[0]\n",
    "    return float(ic.sum() / (2.0 * L))\n",
    "\n",
    "\n",
    "ic = information_content(ppm)\n",
    "kd = kdic(ppm)\n",
    "\n",
    "print('Per-position IC (bits):')\n",
    "for i, (ic_val, base_probs) in enumerate(zip(ic, ppm)):\n",
    "    dominant = BASES[np.argmax(base_probs)]\n",
    "    print(f'  Position {i+1}: IC = {ic_val:.3f} bits  (dominant: {dominant})')\n",
    "print(f'\\nTotal IC: {ic.sum():.3f} bits')\n",
    "print(f'KDIC:     {kd:.4f}  (max possible = 1.0)')\n",
    "\n",
    "fig, ax = plt.subplots(figsize=(7, 3))\n",
    "ax.bar(range(1, len(ic) + 1), ic, color='#5aae61', edgecolor='white')\n",
    "ax.axhline(2.0, color='red', ls='--', lw=1, label='Max IC (2 bits)')\n",
    "ax.set_xlabel('Motif position')\n",
    "ax.set_ylabel('IC (bits)')\n",
    "ax.set_title(f'Per-position information content  |  KDIC = {kd:.3f}')\n",
    "ax.legend()\n",
    "plt.tight_layout()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-pwm-scoring",
   "metadata": {},
   "source": [
    "## 4. PWM Scoring and IUPAC Consensus\n",
    "\n",
    "A **Position Weight Matrix (PWM)** converts PPM probabilities to log-odds scores:\n",
    "\n",
    "$$\\text{PWM}_{i,b} = \\log_2 \\frac{p_{i,b}}{q_b}$$\n",
    "\n",
    "The total score of a sequence is the sum of per-position log-odds values. A positive score means the sequence is more likely under the motif model than under background.\n",
    "\n",
    "**IUPAC consensus** collapses the PPM to a single string using the standard degenerate nucleotide codes (e.g., S = C or G, R = A or G).\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "pwm-scoring",
   "metadata": {},
   "outputs": [],
   "source": [
    "IUPAC = {\n",
    "    frozenset('A'):    'A', frozenset('C'):    'C',\n",
    "    frozenset('G'):    'G', frozenset('T'):    'T',\n",
    "    frozenset('AC'):   'M', frozenset('AG'):   'R',\n",
    "    frozenset('AT'):   'W', frozenset('CG'):   'S',\n",
    "    frozenset('CT'):   'Y', frozenset('GT'):   'K',\n",
    "    frozenset('ACG'):  'V', frozenset('ACT'):  'H',\n",
    "    frozenset('AGT'):  'D', frozenset('CGT'):  'B',\n",
    "    frozenset('ACGT'): 'N',\n",
    "}\n",
    "\n",
    "def build_pwm(ppm: np.ndarray, bg: np.ndarray | None = None) -> np.ndarray:\n",
    "    \"\"\"Return PWM of shape (L, 4) in log2 log-odds.\"\"\"\n",
    "    if bg is None:\n",
    "        bg = np.full(4, 0.25)\n",
    "    return np.log2(np.clip(ppm, 1e-9, 1.0) / bg)\n",
    "\n",
    "\n",
    "def score_sequence(seq: str, pwm: np.ndarray) -> float:\n",
    "    \"\"\"Score a sequence of length L against a PWM of shape (L, 4).\"\"\"\n",
    "    return sum(\n",
    "        pwm[pos, BASE_INDEX[base]]\n",
    "        for pos, base in enumerate(seq.upper())\n",
    "        if base in BASE_INDEX\n",
    "    )\n",
    "\n",
    "\n",
    "def ppm_to_iupac(ppm: np.ndarray, threshold: float = 0.25) -> str:\n",
    "    \"\"\"Return IUPAC consensus. Include bases with probability >= threshold.\"\"\"\n",
    "    consensus = []\n",
    "    for row in ppm:\n",
    "        present = frozenset(b for b, p in zip(BASES, row) if p >= threshold)\n",
    "        consensus.append(IUPAC.get(present, 'N'))\n",
    "    return ''.join(consensus)\n",
    "\n",
    "\n",
    "pwm = build_pwm(ppm)\n",
    "consensus = ppm_to_iupac(ppm)\n",
    "print(f'IUPAC consensus: {consensus}')\n",
    "\n",
    "# Score a genuine AP-1 site vs a random sequence\n",
    "test_seqs = {\n",
    "    'TGACTCA (genuine)': 'TGACTCA',\n",
    "    'TGAGTCA (genuine)': 'TGAGTCA',\n",
    "    'AAAAAAA (random)':  'AAAAAAA',\n",
    "    'GCGCGCG (random)':  'GCGCGCG',\n",
    "    'TGAAACA (1 mismatch)': 'TGAAACA',\n",
    "}\n",
    "\n",
    "print('\\nPWM scores:')\n",
    "for label, seq in test_seqs.items():\n",
    "    print(f'  {label}: {score_sequence(seq, pwm):+.3f}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-exact-dist",
   "metadata": {},
   "source": [
    "## 5. Exact Score Distribution for Short Motifs\n",
    "\n",
    "For motifs of length L ≤ 8, we can enumerate all 4^L possible k-mers and compute their PWM scores exactly. This gives the full null distribution under a uniform background, enabling exact p-value calculation without Monte Carlo sampling.\n",
    "\n",
    "For L = 7 (our AP-1 motif): 4^7 = 16 384 k-mers — fast to enumerate.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "exact-distribution",
   "metadata": {},
   "outputs": [],
   "source": [
    "def exact_score_distribution(pwm: np.ndarray) -> np.ndarray:\n",
    "    \"\"\"\n",
    "    Enumerate all 4^L k-mers and return their PWM scores.\n",
    "    Only feasible for L <= 8.\n",
    "    \"\"\"\n",
    "    L = pwm.shape[0]\n",
    "    assert L <= 8, f'L={L} too large for exhaustive enumeration (max 8)'\n",
    "    scores = np.array([\n",
    "        sum(pwm[pos, base] for pos, base in enumerate(kmer))\n",
    "        for kmer in itertools.product(range(4), repeat=L)\n",
    "    ])\n",
    "    return scores\n",
    "\n",
    "\n",
    "def exact_p_value(score: float, all_scores: np.ndarray) -> float:\n",
    "    \"\"\"Fraction of k-mers scoring >= score under uniform background.\"\"\"\n",
    "    return float((all_scores >= score).mean())\n",
    "\n",
    "\n",
    "all_scores = exact_score_distribution(pwm)\n",
    "\n",
    "genuine_score = score_sequence('TGACTCA', pwm)\n",
    "p_val = exact_p_value(genuine_score, all_scores)\n",
    "\n",
    "print(f'Score of TGACTCA:  {genuine_score:.3f}')\n",
    "print(f'Exact p-value:     {p_val:.2e}')\n",
    "print(f'Total k-mers:      {len(all_scores):,}')\n",
    "\n",
    "fig, ax = plt.subplots(figsize=(8, 4))\n",
    "ax.hist(all_scores, bins=60, color='#aec7e8', edgecolor='white', label='All 4^L k-mers')\n",
    "ax.axvline(genuine_score, color='#d62728', lw=2, label=f'TGACTCA score = {genuine_score:.2f}')\n",
    "ax.set_xlabel('PWM score (log2 log-odds)')\n",
    "ax.set_ylabel('Count')\n",
    "ax.set_title(f'Exact null score distribution for L={pwm.shape[0]} motif')\n",
    "ax.legend()\n",
    "plt.tight_layout()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-enrichment",
   "metadata": {},
   "source": [
    "## 6. Motif Enrichment — Fisher's Exact Test\n",
    "\n",
    "To test whether a motif is **enriched in foreground sequences** (e.g. ChIP-seq peaks) relative to background (e.g. GC-matched random regions), we:\n",
    "\n",
    "1. Score every sequence in both sets.\n",
    "2. Call a sequence \"motif-positive\" if its max score ≥ threshold (e.g. 50th percentile of the null distribution).\n",
    "3. Build a 2x2 contingency table and apply Fisher's exact test.\n",
    "\n",
    "With multiple motifs being tested, apply **Benjamini-Hochberg** correction to control the false discovery rate.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "enrichment-test",
   "metadata": {},
   "outputs": [],
   "source": [
    "rng = np.random.default_rng(0)\n",
    "\n",
    "def random_seq(length: int) -> str:\n",
    "    return ''.join(rng.choice(BASES, size=length))\n",
    "\n",
    "def embed_motif(seq: str, motif: str, position: int | None = None) -> str:\n",
    "    \"\"\"Insert motif into seq at a random or specified position.\"\"\"\n",
    "    if position is None:\n",
    "        position = rng.integers(0, len(seq) - len(motif) + 1)\n",
    "    return seq[:position] + motif + seq[position + len(motif):]\n",
    "\n",
    "SEQ_LEN = 200\n",
    "MOTIF = 'TGACTCA'\n",
    "\n",
    "# Foreground: 80% of sequences contain the motif\n",
    "fg_seqs = [\n",
    "    embed_motif(random_seq(SEQ_LEN), MOTIF) if rng.random() < 0.80 else random_seq(SEQ_LEN)\n",
    "    for _ in range(200)\n",
    "]\n",
    "\n",
    "# Background: 5% of sequences contain the motif (background enrichment rate)\n",
    "bg_seqs = [\n",
    "    embed_motif(random_seq(SEQ_LEN), MOTIF) if rng.random() < 0.05 else random_seq(SEQ_LEN)\n",
    "    for _ in range(500)\n",
    "]\n",
    "\n",
    "def scan_max_score(sequences: list[str], pwm: np.ndarray) -> list[float]:\n",
    "    \"\"\"Return max PWM score across all valid L-mer positions for each sequence.\"\"\"\n",
    "    L = pwm.shape[0]\n",
    "    results = []\n",
    "    for seq in sequences:\n",
    "        seq = seq.upper()\n",
    "        window_scores = [\n",
    "            sum(pwm[i, BASE_INDEX.get(seq[pos + i], 0)] for i in range(L))\n",
    "            for pos in range(len(seq) - L + 1)\n",
    "        ]\n",
    "        results.append(max(window_scores) if window_scores else float('-inf'))\n",
    "    return results\n",
    "\n",
    "\n",
    "# Threshold: 50th percentile of null distribution (liberal, for demonstration)\n",
    "threshold = np.percentile(all_scores, 50)\n",
    "\n",
    "fg_scores = scan_max_score(fg_seqs, pwm)\n",
    "bg_scores = scan_max_score(bg_seqs, pwm)\n",
    "\n",
    "n_fg_with = sum(s >= threshold for s in fg_scores)\n",
    "n_bg_with = sum(s >= threshold for s in bg_scores)\n",
    "\n",
    "contingency = [\n",
    "    [n_fg_with,              len(fg_seqs) - n_fg_with],\n",
    "    [n_bg_with,              len(bg_seqs) - n_bg_with],\n",
    "]\n",
    "odds_ratio, p_value = fisher_exact(contingency, alternative='greater')\n",
    "\n",
    "print(f'Threshold score:       {threshold:.3f}')\n",
    "print(f'Foreground positives:  {n_fg_with}/{len(fg_seqs)} ({100*n_fg_with/len(fg_seqs):.1f}%)')\n",
    "print(f'Background positives:  {n_bg_with}/{len(bg_seqs)} ({100*n_bg_with/len(bg_seqs):.1f}%)')\n",
    "print(f'Odds ratio:            {odds_ratio:.2f}')\n",
    "print(f'p-value:               {p_value:.2e}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-bh",
   "metadata": {},
   "source": [
    "## 7. Multiple Testing Correction — Benjamini-Hochberg\n",
    "\n",
    "When testing many motifs simultaneously, raw p-values will produce false positives. **Benjamini-Hochberg (BH)** controls the expected False Discovery Rate (FDR) at level α:\n",
    "\n",
    "1. Sort p-values: $p_{(1)} \\leq p_{(2)} \\leq \\ldots \\leq p_{(m)}$\n",
    "2. Find the largest k such that $p_{(k)} \\leq \\frac{k}{m} \\alpha$\n",
    "3. Reject all hypotheses with rank ≤ k\n",
    "\n",
    "`scipy.stats.false_discovery_control` (scipy ≥ 1.11) implements this directly.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "bh-correction",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Simulate testing 20 motifs; one of them is truly enriched\n",
    "simulated_p_values = list(rng.uniform(0.1, 1.0, size=19)) + [p_value]\n",
    "motif_names = [f'Motif_{i:02d}' for i in range(19)] + ['AP1_TGASTCA']\n",
    "\n",
    "adjusted_p = false_discovery_control(simulated_p_values, method='bh')\n",
    "\n",
    "results_df = pd.DataFrame({\n",
    "    'motif': motif_names,\n",
    "    'p_value': simulated_p_values,\n",
    "    'q_value': adjusted_p,\n",
    "}).sort_values('p_value').reset_index(drop=True)\n",
    "\n",
    "results_df['significant_q05'] = results_df['q_value'] < 0.05\n",
    "\n",
    "print('Top motifs by p-value:')\n",
    "display(results_df.head(10).round(4))\n",
    "\n",
    "# Visualise the correction\n",
    "fig, axes = plt.subplots(1, 2, figsize=(11, 4))\n",
    "\n",
    "axes[0].scatter(\n",
    "    range(len(results_df)), results_df['p_value'],\n",
    "    c=results_df['significant_q05'].map({True: '#d62728', False: '#aec7e8'}),\n",
    "    s=60, zorder=3,\n",
    ")\n",
    "axes[0].axhline(0.05, color='grey', ls='--', lw=1, label='p = 0.05')\n",
    "axes[0].set_xlabel('Motif rank')\n",
    "axes[0].set_ylabel('Raw p-value')\n",
    "axes[0].set_title('Raw p-values')\n",
    "axes[0].legend()\n",
    "\n",
    "axes[1].scatter(\n",
    "    range(len(results_df)), results_df['q_value'],\n",
    "    c=results_df['significant_q05'].map({True: '#d62728', False: '#aec7e8'}),\n",
    "    s=60, zorder=3,\n",
    ")\n",
    "axes[1].axhline(0.05, color='grey', ls='--', lw=1, label='FDR = 0.05')\n",
    "axes[1].set_xlabel('Motif rank')\n",
    "axes[1].set_ylabel('BH-adjusted q-value')\n",
    "axes[1].set_title('BH-corrected q-values')\n",
    "axes[1].legend()\n",
    "\n",
    "red_patch = mpatches.Patch(color='#d62728', label='Significant (q < 0.05)')\n",
    "blue_patch = mpatches.Patch(color='#aec7e8', label='Not significant')\n",
    "fig.legend(handles=[red_patch, blue_patch], loc='upper center', ncol=2, fontsize=9)\n",
    "\n",
    "plt.tight_layout()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-jaspar",
   "metadata": {},
   "source": [
    "## 8. Querying the JASPAR 2024 REST API\n",
    "\n",
    "JASPAR is the reference database of curated TF binding profiles. The REST API at `https://jaspar.elixir.no/api/v1/` allows programmatic access to all matrices.\n",
    "\n",
    "Each matrix has a stable ID such as `MA0139.1` (CTCF) or `MA0099.3` (AP-1/FOS). The `/matrix/{id}/` endpoint returns a JSON object including the raw PFM counts.\n",
    "\n",
    "**Note**: this cell requires an internet connection. If offline, the fallback PPM below is used instead.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "jaspar-api",
   "metadata": {},
   "outputs": [],
   "source": [
    "def fetch_jaspar_ppm(matrix_id: str) -> dict | None:\n",
    "    \"\"\"\n",
    "    Fetch a JASPAR PFM by matrix ID and convert to PPM of shape (L, 4).\n",
    "    Returns None if the request fails (offline, etc.).\n",
    "    \"\"\"\n",
    "    url = f'https://jaspar.elixir.no/api/v1/matrix/{matrix_id}/'\n",
    "    try:\n",
    "        response = requests.get(url, timeout=8)\n",
    "        response.raise_for_status()\n",
    "    except requests.RequestException as exc:\n",
    "        print(f'JASPAR request failed: {exc}')\n",
    "        return None\n",
    "\n",
    "    data = response.json()\n",
    "    pfm = data['pfm']   # {'A': [...], 'C': [...], 'G': [...], 'T': [...]}\n",
    "    L = len(pfm['A'])\n",
    "    counts = np.array([pfm[b] for b in BASES], dtype=float).T  # (L, 4)\n",
    "    row_sums = counts.sum(axis=1, keepdims=True)\n",
    "    ppm_jaspar = counts / np.where(row_sums == 0, 1.0, row_sums)\n",
    "    return {'name': data['name'], 'matrix_id': matrix_id, 'ppm': ppm_jaspar}\n",
    "\n",
    "\n",
    "# Fetch AP-1 (FOS::JUN heterodimer) from JASPAR 2024\n",
    "jaspar_result = fetch_jaspar_ppm('MA0099.3')\n",
    "\n",
    "if jaspar_result is not None:\n",
    "    print(f\"Fetched: {jaspar_result['name']} ({jaspar_result['matrix_id']})\")\n",
    "    print(f\"PPM shape: {jaspar_result['ppm'].shape}\")\n",
    "    plot_sequence_logo(jaspar_result['ppm'], title=f\"JASPAR {jaspar_result['matrix_id']}: {jaspar_result['name']}\")\n",
    "else:\n",
    "    print('Using locally built AP-1 PPM as fallback.')\n",
    "    plot_sequence_logo(ppm, title='Local AP-1 PPM (offline fallback)')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-tomtom",
   "metadata": {},
   "source": [
    "## 9. Motif Comparison — TomTom-style Cosine Similarity\n",
    "\n",
    "**TomTom** (MEME Suite) compares a query motif to a database by computing a similarity score at every alignment offset. We implement a simplified version using **cosine similarity** between flattened PPM vectors:\n",
    "\n",
    "$$\\text{cosine}(\\mathbf{u}, \\mathbf{v}) = \\frac{\\mathbf{u} \\cdot \\mathbf{v}}{\\|\\mathbf{u}\\| \\|\\mathbf{v}\\|}$$\n",
    "\n",
    "For two PPMs of possibly different lengths, we slide the shorter over the longer, compute column-by-column cosine similarity at each offset, and return the maximum.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "tomtom-similarity",
   "metadata": {},
   "outputs": [],
   "source": [
    "def cosine_similarity_ppms(ppm_a: np.ndarray, ppm_b: np.ndarray) -> float:\n",
    "    \"\"\"\n",
    "    Compute maximum cosine similarity between two PPMs at the best alignment offset.\n",
    "    Both PPMs have shape (L, 4). Shorter is slid over longer.\n",
    "    \"\"\"\n",
    "    if ppm_a.shape[0] > ppm_b.shape[0]:\n",
    "        ppm_a, ppm_b = ppm_b, ppm_a   # ensure ppm_a is shorter or equal\n",
    "    La, Lb = ppm_a.shape[0], ppm_b.shape[0]\n",
    "    best = -1.0\n",
    "    for offset in range(Lb - La + 1):\n",
    "        window = ppm_b[offset: offset + La]   # (La, 4)\n",
    "        u = ppm_a.ravel()\n",
    "        v = window.ravel()\n",
    "        norm = np.linalg.norm(u) * np.linalg.norm(v)\n",
    "        if norm > 0:\n",
    "            sim = float(np.dot(u, v) / norm)\n",
    "            best = max(best, sim)\n",
    "    return best\n",
    "\n",
    "\n",
    "# Compare our locally built AP-1 PPM to a few synthetic motifs\n",
    "random_ppm = build_ppm([random_seq(7) for _ in range(50)])\n",
    "similar_ppm = build_ppm(['TGACTCA', 'TGACTCA', 'TGAGTCA', 'TGACTCA'] * 5)\n",
    "\n",
    "print('Cosine similarity — local AP-1 PPM vs:')\n",
    "print(f'  Self:           {cosine_similarity_ppms(ppm, ppm):.4f}  (should be 1.0)')\n",
    "print(f'  Similar PPM:    {cosine_similarity_ppms(ppm, similar_ppm):.4f}')\n",
    "print(f'  Random PPM:     {cosine_similarity_ppms(ppm, random_ppm):.4f}')\n",
    "\n",
    "if jaspar_result is not None:\n",
    "    # JASPAR motif may have different length — handled by offset sliding\n",
    "    sim_jaspar = cosine_similarity_ppms(ppm, jaspar_result['ppm'])\n",
    "    print(f\"  JASPAR {jaspar_result['matrix_id']}: {sim_jaspar:.4f}\")"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-pipeline",
   "metadata": {},
   "source": [
    "## 10. Pipeline Design Pattern — Abstract Base Class and Dataclass Config\n",
    "\n",
    "When wrapping external tools (MEME, HOMER, STREME) in a Python pipeline, a common pattern is:\n",
    "\n",
    "1. A `dataclass` for configuration (avoids long constructor signatures)\n",
    "2. An abstract base class defining the interface (`run`, `parse_output`)\n",
    "3. Concrete subclasses for each tool\n",
    "\n",
    "This lets you swap tools without changing downstream code.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "pipeline-pattern",
   "metadata": {},
   "outputs": [],
   "source": [
    "from abc import ABC, abstractmethod\n",
    "from dataclasses import dataclass, field\n",
    "from pathlib import Path\n",
    "\n",
    "\n",
    "@dataclass\n",
    "class MotifToolConfig:\n",
    "    fasta_fg: Path\n",
    "    fasta_bg: Path\n",
    "    output_dir: Path\n",
    "    n_motifs: int = 5\n",
    "    motif_width: int | None = None   # None = auto-detect\n",
    "    n_threads: int = 4\n",
    "    extra_args: list[str] = field(default_factory=list)\n",
    "\n",
    "\n",
    "@dataclass\n",
    "class DiscoveredMotif:\n",
    "    name: str\n",
    "    ppm: np.ndarray   # shape (L, 4)\n",
    "    e_value: float\n",
    "    n_sites: int\n",
    "    consensus: str\n",
    "\n",
    "\n",
    "class MotifToolBase(ABC):\n",
    "    \"\"\"Abstract base class for motif discovery tool wrappers.\"\"\"\n",
    "\n",
    "    def __init__(self, config: MotifToolConfig) -> None:\n",
    "        self.config = config\n",
    "\n",
    "    @abstractmethod\n",
    "    def build_command(self) -> list[str]:\n",
    "        \"\"\"Return the CLI command as a list of strings.\"\"\"\n",
    "\n",
    "    @abstractmethod\n",
    "    def parse_output(self) -> list[DiscoveredMotif]:\n",
    "        \"\"\"Parse tool output into a list of DiscoveredMotif objects.\"\"\"\n",
    "\n",
    "    def run(self) -> list[DiscoveredMotif]:\n",
    "        \"\"\"Execute the tool and return parsed motifs.\"\"\"\n",
    "        import subprocess\n",
    "        cmd = self.build_command()\n",
    "        print(f'Running: {\" \".join(str(c) for c in cmd)}')\n",
    "        result = subprocess.run(cmd, capture_output=True, text=True, check=True)\n",
    "        return self.parse_output()\n",
    "\n",
    "\n",
    "class StremeWrapper(MotifToolBase):\n",
    "    \"\"\"Wrapper for the STREME motif discovery tool (MEME Suite >= 5.5).\"\"\"\n",
    "\n",
    "    def build_command(self) -> list[str]:\n",
    "        cmd = [\n",
    "            'streme',\n",
    "            '--p', str(self.config.fasta_fg),\n",
    "            '--n', str(self.config.fasta_bg),\n",
    "            '--oc', str(self.config.output_dir),\n",
    "            '--nmotifs', str(self.config.n_motifs),\n",
    "            '--totallength', '4000000',\n",
    "        ]\n",
    "        if self.config.motif_width is not None:\n",
    "            cmd += ['--w', str(self.config.motif_width)]\n",
    "        cmd += self.config.extra_args\n",
    "        return cmd\n",
    "\n",
    "    def parse_output(self) -> list[DiscoveredMotif]:\n",
    "        # In a real implementation: parse streme.txt or streme.xml\n",
    "        # Here we return an empty list as a structural placeholder\n",
    "        streme_txt = self.config.output_dir / 'streme.txt'\n",
    "        if not streme_txt.exists():\n",
    "            return []\n",
    "        # ... XML/text parsing logic here ...\n",
    "        return []\n",
    "\n",
    "\n",
    "print('Pipeline classes defined successfully.')\n",
    "print(f'MotifToolConfig fields: {[f.name for f in MotifToolConfig.__dataclass_fields__.values()]}')\n",
    "print(f'DiscoveredMotif fields: {[f.name for f in DiscoveredMotif.__dataclass_fields__.values()]}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "exercise-header",
   "metadata": {},
   "source": [
    "## Exercise\n",
    "\n",
    "Twenty DNA sequences of length 30 bp are provided below. Exactly 15 of them contain the motif **TGASTCA** (S = C or G; the exact variant is randomly chosen) embedded at a random position. The remaining 5 are pure random sequence.\n",
    "\n",
    "Your task:\n",
    "\n",
    "1. **Build a PPM** from the embedded motif instances. You can extract the 7-mer at its known position for this toy example, or try to recover it de novo by scanning.\n",
    "2. **Compute the KDIC** of your recovered PPM. A correct recovery should give KDIC ≥ 0.85.\n",
    "3. **Identify the top-scoring 7-mer position** in each of the 20 sequences using the PWM you built.\n",
    "4. **Check**: in sequences that contain the motif, does the top-scoring position coincide with the embedded location?\n",
    "\n",
    "_Hint_: build the PPM only from sequences 0–14 (those you know contain the motif). The motif positions are provided in `motif_positions`.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "exercise-cell",
   "metadata": {},
   "outputs": [],
   "source": [
    "exercise_rng = np.random.default_rng(99)\n",
    "\n",
    "SEQ_LEN = 30\n",
    "MOTIF_LEN = 7\n",
    "N_WITH_MOTIF = 15\n",
    "\n",
    "exercise_seqs = []\n",
    "motif_positions = []\n",
    "embedded_motifs = []\n",
    "\n",
    "for i in range(20):\n",
    "    seq = ''.join(exercise_rng.choice(BASES, size=SEQ_LEN))\n",
    "    if i < N_WITH_MOTIF:\n",
    "        pos = int(exercise_rng.integers(0, SEQ_LEN - MOTIF_LEN + 1))\n",
    "        # S = C or G at position 3 (0-indexed)\n",
    "        s_base = exercise_rng.choice(['C', 'G'])\n",
    "        motif = f'TGA{s_base}TCA'\n",
    "        seq = seq[:pos] + motif + seq[pos + MOTIF_LEN:]\n",
    "        motif_positions.append(pos)\n",
    "        embedded_motifs.append(motif)\n",
    "    else:\n",
    "        motif_positions.append(None)\n",
    "        embedded_motifs.append(None)\n",
    "    exercise_seqs.append(seq)\n",
    "\n",
    "print('Exercise sequences (first 5):')\n",
    "for i in range(5):\n",
    "    print(f'  seq[{i}]: {exercise_seqs[i]}  motif at pos {motif_positions[i]} ({embedded_motifs[i]})')\n",
    "\n",
    "print('\\n--- Your solution below ---')\n",
    "\n",
    "# Step 1: Extract the embedded 7-mers from sequences 0-14 and build a PPM\n",
    "# ...\n",
    "\n",
    "# Step 2: Compute KDIC\n",
    "# ...\n",
    "\n",
    "# Step 3: Find top-scoring position in each sequence\n",
    "# ...\n",
    "\n",
    "# Step 4: Check recovery accuracy\n",
    "# ..."
   ]
  }
 ]
}
```

- [ ] Step 4: Create `15_Motif_Discovery/README.md`

```markdown
# Module 15: Motif Discovery

**Tier 2 — Core Bioinformatics | Module 15**

## Overview

This module covers the theory and practice of transcription factor motif discovery from DNA sequences. Starting from first principles (Position Probability Matrices, information content), it builds up to enrichment testing, multiple-testing correction, and database queries against JASPAR 2024.

## Topics

1. Building PPMs and PWMs from aligned sequences
2. Per-position information content and KDIC
3. Sequence logo visualisation
4. Log-odds scoring of arbitrary sequences
5. Exact null score distribution for short motifs (L ≤ 8)
6. Motif enrichment with Fisher's exact test
7. Benjamini-Hochberg FDR correction for multiple motifs
8. Querying the JASPAR 2024 REST API
9. Motif comparison with cosine similarity (TomTom concept)
10. Pipeline design pattern: abstract base class + dataclass config

## Prerequisites

- Modules 10–12 (ChIP-seq, ATAC-seq, peak calling)
- NumPy, scipy, requests

## Installation

```bash
pip install numpy scipy requests matplotlib
# optional: pip install logomaker  # for publication-quality logos
```

## Public Data Sources

| Resource | URL |
|---|---|
| JASPAR 2024 API | https://jaspar.elixir.no/api/v1/ |
| HOCOMOCO v12 | https://hocomoco12.autosome.org/ |
| MEME Suite (STREME/TomTom) | https://meme-suite.org/ |

## Learning Objectives

After completing this module you will be able to:

- Build PPMs with the correct (L, 4) shape convention
- Compute IC and KDIC to quantify motif specificity
- Score sequences with log-odds PWMs
- Test enrichment using Fisher's exact test with BH correction
- Retrieve motifs from JASPAR via REST API
- Design a tool-wrapper pipeline using abstract base classes

## Skill File

See `Skills/motif-discovery.md` for quick-reference patterns and copy-paste code templates.
```

- [ ] Step 5: Commit

```bash
git add Course/Tier_2_Core_Bioinformatics/15_Motif_Discovery/15_motif_discovery.ipynb \
        Course/Tier_2_Core_Bioinformatics/15_Motif_Discovery/README.md \
        Skills/motif-discovery.md
git commit -m "feat(tier2): add motif discovery module (15) and skill"
```

### Task 3: GWAS + Spatial Transcriptomics + Copy Number Analysis (Tier 3, Modules 19–21)

**Files:**
- Create: `Course/Tier_3_Applied_Bioinformatics/19_GWAS/19_gwas.ipynb`
- Create: `Course/Tier_3_Applied_Bioinformatics/19_GWAS/README.md`
- Create: `Course/Tier_3_Applied_Bioinformatics/20_Spatial_Transcriptomics/20_spatial_transcriptomics.ipynb`
- Create: `Course/Tier_3_Applied_Bioinformatics/20_Spatial_Transcriptomics/README.md`
- Create: `Course/Tier_3_Applied_Bioinformatics/21_Copy_Number_Analysis/21_copy_number_analysis.ipynb`
- Create: `Course/Tier_3_Applied_Bioinformatics/21_Copy_Number_Analysis/README.md`
- Create: `Skills/gwas-population-genetics.md`
- Create: `Skills/spatial-transcriptomics.md`

**Steps:**

- [ ] Step 1: Create directory structure

```bash
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_3_Applied_Bioinformatics/19_GWAS
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_3_Applied_Bioinformatics/20_Spatial_Transcriptomics
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_3_Applied_Bioinformatics/21_Copy_Number_Analysis
```

- [ ] Step 2: Create `Skills/gwas-population-genetics.md`

```markdown
---
name: gwas-population-genetics
description: Genome-wide association studies — study design, genotype QC, population stratification PCA, per-SNP logistic regression, Manhattan and QQ plots, LD clumping, and fine-mapping concepts
---

## When to Use

Use this skill when:
- Designing or analyzing a GWAS (case-control or quantitative trait)
- Performing genotype QC (MAF, call rate, HWE filters)
- Detecting and correcting for population stratification
- Running per-SNP association tests with covariate adjustment
- Creating Manhattan and Q-Q plots
- Interpreting GWAS hits with LD clumping

## Quick Reference

| Step | Tool / Method | Notes |
|---|---|---|
| Genotype matrix | `numpy` ndarray (samples × SNPs) | 0/1/2 = ref/het/alt dosage |
| MAF filter | `maf >= 0.01` | Remove rare variants |
| HWE filter | `chi2 test, p > 1e-6` | Apply to controls only |
| PCA | `sklearn.decomposition.PCA` | Standardize genotypes first |
| Association | `statsmodels.Logit` / `OLS` | Add PC1-10 as covariates |
| GW significance | `p < 5e-8` | Bonferroni for ~1M tests |
| Manhattan plot | `-log10(p)` vs genomic position | Color alternating chromosomes |
| QQ plot | observed vs expected `-log10(p)` | Slope = genomic inflation λ |
| LD clumping | group nearby SNPs by R² | Keep index SNP per block |

## Key Patterns

**Pattern 1: Simulate genotype data**
```python
import numpy as np

rng = np.random.default_rng(42)
N_SAMPLES, N_SNPS = 500, 2000
mafs = rng.uniform(0.05, 0.45, N_SNPS)

genotypes = np.column_stack([
    rng.choice([0, 1, 2], N_SAMPLES,
               p=[(1-f)**2, 2*f*(1-f), f**2])
    for f in mafs
]).astype(np.int8)

causal_idx = rng.choice(N_SNPS, 5, replace=False)
logit = genotypes[:, causal_idx].dot(rng.normal(0.6, 0.1, 5)) - 1.5
y = rng.binomial(1, 1/(1+np.exp(-logit))).astype(float)
```

**Pattern 2: SNP QC (MAF + HWE)**
```python
from scipy import stats

maf_obs = np.minimum(genotypes.mean(0)/2, 1 - genotypes.mean(0)/2)
maf_pass = maf_obs >= 0.01

def hwe_pvalue(col):
    counts = np.bincount(col.astype(int), minlength=3)
    n = counts.sum()
    p = (2*counts[2] + counts[1]) / (2*n)
    q = 1 - p
    exp = np.array([n*q**2, 2*n*p*q, n*p**2])
    chi2 = ((counts - exp)**2 / (exp + 1e-10)).sum()
    return float(stats.chi2.sf(chi2, df=1))

ctrl = y == 0
hwe_p = np.array([hwe_pvalue(genotypes[ctrl, j]) for j in range(N_SNPS)])
keep = maf_pass & (hwe_p > 1e-6)
geno_qc = genotypes[:, keep]
```

**Pattern 3: PCA for stratification**
```python
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

pcs = PCA(10).fit_transform(StandardScaler().fit_transform(geno_qc))
```

**Pattern 4: Per-SNP logistic regression**
```python
import statsmodels.api as sm

def gwas_logit(geno_col, pheno, covs):
    X = sm.add_constant(np.column_stack([geno_col, covs]))
    try:
        fit = sm.Logit(pheno, X, missing='drop').fit(
            disp=False, maxiter=100, method='bfgs')
        return fit.params[1], fit.pvalues[1]
    except Exception:
        return np.nan, 1.0

covariates = pcs[:, :5]
results = [gwas_logit(geno_qc[:, j], y, covariates)
           for j in range(geno_qc.shape[1])]
betas, pvals = map(np.array, zip(*results))
```

**Pattern 5: Manhattan and QQ plots**
```python
import matplotlib.pyplot as plt

neg_log10_p = -np.log10(np.maximum(pvals, 1e-300))

fig, axes = plt.subplots(1, 2, figsize=(16, 4))

# Manhattan
colors = ['#2166ac', '#b2182b']
x_off = 0
xticks, xlabs = [], []
for chrom in range(1, 23):
    mask = chr_qc == chrom
    if not mask.any(): continue
    xs = pos_qc[mask] + x_off
    axes[0].scatter(xs, neg_log10_p[mask], c=colors[chrom%2], s=4, alpha=0.8, linewidths=0)
    xticks.append(xs.mean()); xlabs.append(str(chrom))
    x_off += pos_qc[mask].max() + 5_000_000
axes[0].axhline(-np.log10(5e-8), color='red', lw=1, ls='--')
axes[0].set_xticks(xticks); axes[0].set_xticklabels(xlabs, fontsize=7)
axes[0].set_xlabel('Chromosome'); axes[0].set_ylabel('-log₁₀(p)')

# QQ
n = len(pvals)
obs = np.sort(neg_log10_p)[::-1]
exp = -np.log10((np.arange(n)+0.5)/n)
chi2_obs = stats.chi2.ppf(1-pvals, df=1)
lam = np.nanmedian(chi2_obs) / stats.chi2.ppf(0.5, df=1)
axes[1].scatter(np.sort(exp), np.sort(obs), s=4, alpha=0.7, c='#333333')
axes[1].plot([0, exp.max()], [0, exp.max()], 'r--', lw=1)
axes[1].set_title(f'QQ plot  λ = {lam:.3f}')
plt.tight_layout(); plt.show()
```

## Code Templates

**Template 1: Complete GWAS pipeline (simulated data)**
```python
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats

rng = np.random.default_rng(0)
N_SAMPLES, N_SNPS = 500, 2000
mafs = rng.uniform(0.05, 0.4, N_SNPS)
genotypes = np.column_stack([
    rng.choice([0,1,2], N_SAMPLES, p=[(1-f)**2, 2*f*(1-f), f**2])
    for f in mafs
])
causal = rng.choice(N_SNPS, 3, replace=False)
y = rng.binomial(1, 1/(1+np.exp(-genotypes[:,causal].sum(1)*0.6+1.5))).astype(float)

maf_obs = np.minimum(genotypes.mean(0)/2, 1-genotypes.mean(0)/2)
keep = maf_obs >= 0.01
geno = genotypes[:, keep]
pcs = PCA(5).fit_transform(StandardScaler().fit_transform(geno))

results = []
for j in range(geno.shape[1]):
    X = sm.add_constant(np.column_stack([geno[:,j], pcs]))
    try:
        fit = sm.Logit(y, X).fit(disp=False, maxiter=50)
        results.append((fit.params[1], fit.pvalues[1]))
    except Exception:
        results.append((np.nan, 1.0))
betas, pvals = map(np.array, zip(*results))
print(f"Hits at p < 5e-8: {(pvals < 5e-8).sum()}")
```

## Common Pitfalls

1. **Stratification without PCs**: Always include PC1–10 as covariates. Omitting them inflates λ and produces false positives from ancestry differences.

2. **HWE filter on cases**: Apply HWE filtering only in controls. Cases may genuinely deviate at disease-associated loci.

3. **Bonferroni vs 5×10⁻⁸**: The 5×10⁻⁸ threshold is Bonferroni for ~1M independent LD blocks genome-wide. For a small simulated set, use `0.05 / N_SNPS` instead.

4. **LD inflation of hit count**: Multiple significant SNPs in a LD block are not independent signals. Always clump before counting.

5. **Log-odds vs odds ratio**: `statsmodels` returns log-odds coefficients. Report `np.exp(beta)` as the odds ratio. Effect sizes < 1 mean the alt allele is protective.
```

- [ ] Step 3: Create `Skills/spatial-transcriptomics.md`

```markdown
---
name: spatial-transcriptomics
description: Spatial transcriptomics analysis — AnnData with spatial coordinates, QC, normalization, spatial neighborhood graphs, spatially variable genes via Moran's I, and visualization using squidpy and scanpy
---

## When to Use

Use this skill when:
- Analyzing 10x Visium, Xenium, MERFISH, or Slide-seq data
- Loading count matrices that include spatial (x, y) coordinates
- Building spatial neighborhood graphs for downstream analyses
- Identifying spatially variable genes (SVGs) via Moran's I
- Visualizing gene expression overlaid on tissue images
- Understanding cell-type deconvolution concepts (RCTD, cell2location)

## Quick Reference

| Tool | Method | Purpose |
|---|---|---|
| `squidpy` | `sq.datasets.visium_hne_adata_crop()` | Built-in demo Visium dataset |
| `scanpy` | `sc.pp.filter_cells/genes` | Standard QC |
| `scanpy` | `sc.pp.normalize_total`, `log1p` | Size-factor normalization |
| `squidpy` | `sq.gr.spatial_neighbors()` | Build spatial graph |
| `squidpy` | `sq.gr.spatial_autocorr()` | Moran's I for SVGs |
| `squidpy` | `sq.pl.spatial_scatter()` | Tissue overlay plot |
| AnnData | `.obsm["spatial"]` | Coordinates array (n_obs × 2) |
| AnnData | `.uns["spatial"]` | Image data + scale factors |

## Key Patterns

**Pattern 1: Load built-in Visium dataset**
```python
import squidpy as sq
import scanpy as sc

adata = sq.datasets.visium_hne_adata_crop()
print(adata)
print(adata.obsm["spatial"][:3])   # (row, col) pixel coordinates
```

**Pattern 2: Spatial QC**
```python
sc.pp.filter_cells(adata, min_counts=200)
sc.pp.filter_genes(adata, min_cells=5)
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
```

**Pattern 3: Normalize, reduce, cluster**
```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=0.4)
```

**Pattern 4: Spatial neighbor graph + Moran's I SVGs**
```python
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
sq.gr.spatial_autocorr(adata, mode="moran", n_perms=100, n_jobs=1)
svgs = adata.uns["moranI"].sort_values("I", ascending=False).head(20)
print(svgs[["I", "pval_norm"]])
```

**Pattern 5: Spatial scatter visualization**
```python
sq.pl.spatial_scatter(adata, color=["leiden", "Ttr"], wspace=0.4)
```

## Code Templates

**Template 1: Full spatial pipeline**
```python
import squidpy as sq
import scanpy as sc

adata = sq.datasets.visium_hne_adata_crop()

sc.pp.filter_cells(adata, min_counts=200)
sc.pp.filter_genes(adata, min_cells=5)
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata)

sq.gr.spatial_neighbors(adata, coord_type="grid")
sq.gr.spatial_autocorr(adata, mode="moran", n_perms=50, n_jobs=1)

sq.pl.spatial_scatter(adata, color="leiden")
```

## Common Pitfalls

1. **`coord_type` mismatch**: Use `"grid"` for Visium (hexagonal spots); use `"generic"` with `radius` or `n_neighs` for Xenium/MERFISH (irregular layouts).

2. **Raw vs normalized for Moran's I**: Run `sq.gr.spatial_autocorr` on log-normalized data, not raw counts. Size-factor variation dominates raw count autocorrelation.

3. **Coordinate units**: `.obsm["spatial"]` stores pixel coordinates; `.uns["spatial"][library_id]["scalefactors"]` has the µm/pixel conversion. Use consistent units for distance-based analyses.

4. **Deconvolution reference**: Cell-type deconvolution (RCTD, cell2location) requires a matching scRNA-seq reference from the same tissue. Cross-species or cross-tissue references give unreliable proportions.

5. **Leiden resolution for spatial data**: Spatial clusters tend to be large (anatomical regions). Start with `resolution=0.3` and visualize spatially before tuning up.
```

- [ ] Step 4: Create `19_gwas.ipynb`

```json
{
 "nbformat": 4,
 "nbformat_minor": 5,
 "metadata": {
  "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
  "language_info": {"name": "python", "version": "3.11.0"}
 },
 "cells": [
  {
   "cell_type": "markdown",
   "id": "header-gwas",
   "metadata": {},
   "source": [
    "# Module 19: Genome-Wide Association Studies (GWAS)\n",
    "\n",
    "Patterns inspired by NGSchool 2023 practical\n",
    "\n",
    "---\n",
    "\n",
    "## Learning Objectives\n",
    "\n",
    "By the end of this notebook you will be able to:\n",
    "\n",
    "1. Explain GWAS study design: cohort, phenotype, SNP density, covariates\n",
    "2. Apply standard genotype QC filters (MAF, call rate, HWE)\n",
    "3. Detect population stratification with PCA and include PCs as covariates\n",
    "4. Run per-SNP logistic regression and interpret odds ratios\n",
    "5. Produce Manhattan and Q-Q plots and compute the genomic inflation factor\n",
    "6. Describe LD clumping and fine-mapping at a conceptual level\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "imports-gwas",
   "metadata": {},
   "outputs": [],
   "source": [
    "import warnings\n",
    "warnings.filterwarnings('ignore')\n",
    "\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt\n",
    "import statsmodels.api as sm\n",
    "from sklearn.decomposition import PCA\n",
    "from sklearn.preprocessing import StandardScaler\n",
    "from scipy import stats\n",
    "\n",
    "rng = np.random.default_rng(42)\n",
    "print('Libraries loaded.')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-design",
   "metadata": {},
   "source": [
    "## 1. Study Design\n",
    "\n",
    "A GWAS scans the genome for common variants (SNPs) associated with a trait. The statistical test at each SNP is simple: 'Is the minor allele more frequent in cases than controls?' Genome-wide the test is repeated for ~500K–10M SNPs, so multiple testing is the central challenge.\n",
    "\n",
    "| Component | Typical value | Note |\n",
    "|---|---|---|\n",
    "| Cohort size | 1,000–1,000,000+ | Larger = more power |\n",
    "| Phenotype | Binary or continuous | Logistic vs linear regression |\n",
    "| Significance threshold | p < 5×10⁻⁸ | Bonferroni for ~1M LD blocks |\n",
    "| Covariates | Age, sex, PC1–10, batch | Reduce confounding |\n",
    "\n",
    "We simulate a **500-sample, 2000-SNP** case-control dataset with 5 causal variants.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "simulate-data",
   "metadata": {},
   "outputs": [],
   "source": [
    "N_SAMPLES, N_SNPS = 500, 2000\n",
    "N_CAUSAL = 5\n",
    "\n",
    "mafs = rng.uniform(0.05, 0.45, N_SNPS)\n",
    "genotypes = np.column_stack([\n",
    "    rng.choice([0, 1, 2], N_SAMPLES,\n",
    "               p=[(1-f)**2, 2*f*(1-f), f**2])\n",
    "    for f in mafs\n",
    "]).astype(np.int8)\n",
    "\n",
    "chromosomes = np.sort(rng.choice(np.arange(1, 23), N_SNPS))\n",
    "positions   = rng.integers(1, 248_000_000, N_SNPS)\n",
    "\n",
    "causal_idx = rng.choice(N_SNPS, N_CAUSAL, replace=False)\n",
    "effects    = rng.normal(0.6, 0.1, N_CAUSAL)\n",
    "logit      = genotypes[:, causal_idx].dot(effects) - 1.5\n",
    "y          = rng.binomial(1, 1/(1+np.exp(-logit))).astype(float)\n",
    "\n",
    "print(f'Cases: {y.sum():.0f}  Controls: {(1-y).sum():.0f}')\n",
    "print(f'Causal SNP indices: {sorted(causal_idx.tolist())}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-qc",
   "metadata": {},
   "source": [
    "## 2. Quality Control\n",
    "\n",
    "Before testing, low-quality SNPs must be removed:\n",
    "\n",
    "| Filter | Threshold | Reason |\n",
    "|---|---|---|\n",
    "| MAF | ≥ 1% | Low-MAF variants underpowered and often artifacts |\n",
    "| SNP call rate | ≥ 95% | High missingness → genotyping failure |\n",
    "| HWE (controls only) | p > 1×10⁻⁶ | Violation suggests batch error |\n",
    "| Sample call rate | ≥ 98% | Remove poorly-genotyped individuals |\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "qc-code",
   "metadata": {},
   "outputs": [],
   "source": [
    "def hwe_pvalue(col: np.ndarray) -> float:\n",
    "    counts = np.bincount(col.astype(int), minlength=3)\n",
    "    n = counts.sum()\n",
    "    p = (2*counts[2] + counts[1]) / (2*n)\n",
    "    q = 1 - p\n",
    "    exp = np.array([n*q**2, 2*n*p*q, n*p**2])\n",
    "    chi2 = ((counts - exp)**2 / (exp + 1e-10)).sum()\n",
    "    return float(stats.chi2.sf(chi2, df=1))\n",
    "\n",
    "maf_obs  = np.minimum(genotypes.mean(0)/2, 1 - genotypes.mean(0)/2)\n",
    "maf_pass = maf_obs >= 0.01\n",
    "\n",
    "ctrl_mask = y == 0\n",
    "hwe_p     = np.array([hwe_pvalue(genotypes[ctrl_mask, j]) for j in range(N_SNPS)])\n",
    "hwe_pass  = hwe_p > 1e-6\n",
    "\n",
    "keep_snps  = maf_pass & hwe_pass\n",
    "geno_qc    = genotypes[:, keep_snps]\n",
    "chr_qc     = chromosomes[keep_snps]\n",
    "pos_qc     = positions[keep_snps]\n",
    "\n",
    "print(f'SNPs before QC: {N_SNPS}')\n",
    "print(f'SNPs after  QC: {keep_snps.sum()}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-pca",
   "metadata": {},
   "source": [
    "## 3. Population Stratification\n",
    "\n",
    "If cases and controls differ in ancestry, allele frequencies differ genome-wide — not due to disease, but due to ancestry. This inflates test statistics.\n",
    "\n",
    "**Fix:** PCA on the genotype matrix captures ancestry axes. Including PC1–10 as covariates in the regression absorbs stratification effects.\n",
    "\n",
    "> In practice: prune SNPs to an LD-independent set (r² < 0.1, window 50 SNPs) before PCA to avoid PCs reflecting local LD rather than ancestry.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "pca-code",
   "metadata": {},
   "outputs": [],
   "source": [
    "geno_std = StandardScaler().fit_transform(geno_qc)\n",
    "pca      = PCA(n_components=10)\n",
    "pcs      = pca.fit_transform(geno_std)\n",
    "\n",
    "fig, axes = plt.subplots(1, 2, figsize=(11, 4))\n",
    "for label, name, c in [(0, 'Control', '#2166ac'), (1, 'Case', '#d73027')]:\n",
    "    m = y == label\n",
    "    axes[0].scatter(pcs[m,0], pcs[m,1], s=12, alpha=0.6, c=c, label=name)\n",
    "axes[0].set_xlabel('PC1'); axes[0].set_ylabel('PC2')\n",
    "axes[0].set_title('PCA — stratification check'); axes[0].legend()\n",
    "\n",
    "axes[1].bar(range(1,11), pca.explained_variance_ratio_*100, color='#4393c3')\n",
    "axes[1].set_xlabel('PC'); axes[1].set_ylabel('Variance explained (%)')\n",
    "axes[1].set_title('Scree plot')\n",
    "plt.tight_layout(); plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-assoc",
   "metadata": {},
   "source": [
    "## 4. Association Testing\n",
    "\n",
    "For a binary phenotype, each SNP is tested with **logistic regression**:\n",
    "\n",
    "$$\\log\\frac{P(\\text{case})}{1-P(\\text{case})} = \\beta_0 + \\beta_1 G_j + \\sum_{k=1}^{10} \\beta_{k+1} \\text{PC}_k$$\n",
    "\n",
    "Test: $H_0: \\beta_1 = 0$ (no association). The effect estimate is the **log odds ratio** for the alt allele; `np.exp(beta)` gives the odds ratio.\n",
    "\n",
    "For quantitative traits, use `statsmodels.OLS` instead.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "assoc-code",
   "metadata": {},
   "outputs": [],
   "source": [
    "def gwas_logit(geno_col, pheno, covs):\n",
    "    X = sm.add_constant(np.column_stack([geno_col, covs]))\n",
    "    try:\n",
    "        fit = sm.Logit(pheno, X, missing='drop').fit(\n",
    "            disp=False, maxiter=100, method='bfgs')\n",
    "        return fit.params[1], fit.pvalues[1]\n",
    "    except Exception:\n",
    "        return np.nan, 1.0\n",
    "\n",
    "covariates = pcs[:, :5]\n",
    "results    = [gwas_logit(geno_qc[:, j], y, covariates)\n",
    "              for j in range(geno_qc.shape[1])]\n",
    "betas, pvals = map(np.array, zip(*results))\n",
    "\n",
    "gw_hits = (pvals < 5e-8).sum()\n",
    "print(f'Genome-wide significant hits (p < 5×10⁻⁸): {gw_hits}')\n",
    "top5 = pd.DataFrame({'beta': betas, 'OR': np.exp(betas), 'p': pvals}).nsmallest(5, 'p')\n",
    "print(top5.to_string())"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "manhattan-qq",
   "metadata": {},
   "outputs": [],
   "source": [
    "neg_log10_p = -np.log10(np.maximum(pvals, 1e-300))\n",
    "\n",
    "fig, axes = plt.subplots(1, 2, figsize=(16, 4))\n",
    "\n",
    "# Manhattan\n",
    "cols = ['#2166ac', '#b2182b']\n",
    "x_off = 0; xticks, xlabs = [], []\n",
    "for chrom in range(1, 23):\n",
    "    mask = chr_qc == chrom\n",
    "    if not mask.any(): continue\n",
    "    xs = pos_qc[mask] + x_off\n",
    "    axes[0].scatter(xs, neg_log10_p[mask], c=cols[chrom%2], s=5, alpha=0.8, linewidths=0)\n",
    "    xticks.append(xs.mean()); xlabs.append(str(chrom))\n",
    "    x_off += pos_qc[mask].max() + 5_000_000\n",
    "axes[0].axhline(-np.log10(5e-8), color='red', lw=1, ls='--', label='5×10⁻⁸')\n",
    "axes[0].set_xticks(xticks); axes[0].set_xticklabels(xlabs, fontsize=7)\n",
    "axes[0].set_xlabel('Chromosome'); axes[0].set_ylabel('-log₁₀(p)')\n",
    "axes[0].set_title('Manhattan plot'); axes[0].legend()\n",
    "\n",
    "# QQ\n",
    "n   = len(pvals)\n",
    "obs = np.sort(neg_log10_p)[::-1]\n",
    "exp = -np.log10((np.arange(n)+0.5)/n)\n",
    "chi2_obs = stats.chi2.ppf(1-np.clip(pvals, 1e-300, 1), df=1)\n",
    "lam = np.nanmedian(chi2_obs) / stats.chi2.ppf(0.5, df=1)\n",
    "axes[1].scatter(np.sort(exp), np.sort(obs), s=5, alpha=0.7, c='#333333')\n",
    "axes[1].plot([0, exp.max()], [0, exp.max()], 'r--', lw=1)\n",
    "axes[1].set_xlabel('Expected -log₁₀(p)'); axes[1].set_ylabel('Observed -log₁₀(p)')\n",
    "axes[1].set_title(f'Q-Q plot  |  λ = {lam:.3f}')\n",
    "plt.tight_layout(); plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-ld",
   "metadata": {},
   "source": [
    "## 5. LD, Clumping, and Fine-mapping\n",
    "\n",
    "**Linkage disequilibrium (LD):** Nearby SNPs are correlated because they are rarely separated by recombination. One causal SNP can make hundreds of correlated neighbours appear significant.\n",
    "\n",
    "**LD clumping (PLINK `--clump`):**\n",
    "1. Sort significant SNPs by p-value\n",
    "2. The top SNP becomes the *index SNP* for that locus\n",
    "3. Remove all SNPs within 250 kb with R² > 0.1 (its LD block)\n",
    "4. Repeat — each remaining SNP is an independent signal\n",
    "\n",
    "**Fine-mapping:** Within the LD block, compute a **95% credible set** — the minimal set of SNPs that contains the true causal SNP with 95% probability. Tools: SuSiE, FINEMAP.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "exercise-gwas",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Exercise 19\n",
    "# ───────────\n",
    "# 1. Re-run association WITHOUT PC covariates.\n",
    "#    How does λ change? Why?\n",
    "#\n",
    "# 2. Our simulation has 2000 SNPs. What is the Bonferroni-corrected p-value\n",
    "#    for this dataset? How does it compare to 5×10⁻⁸?\n",
    "#    When would you use 5×10⁻⁸ vs a dataset-specific threshold?\n",
    "#\n",
    "# 3. (Challenge) Implement simple LD clumping:\n",
    "#    - Find all SNPs with p < 0.01/N_SNPS_QC\n",
    "#    - Sort by p-value; greedily remove SNPs within ±30 indices of each index SNP\n",
    "#    - How many independent hits remain?\n",
    "\n",
    "# Your code here\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "summary-gwas",
   "metadata": {},
   "source": [
    "## Summary\n",
    "\n",
    "| Step | Key tool | Output |\n",
    "|---|---|---|\n",
    "| Simulate / load genotypes | numpy / plink | (samples × SNPs) dosage matrix |\n",
    "| SNP QC | MAF, call rate, HWE | Filtered genotype matrix |\n",
    "| Stratification | PCA on genotypes | PC covariates |\n",
    "| Association | statsmodels.Logit | p-values, odds ratios |\n",
    "| Visualisation | Manhattan, QQ | Genome-wide scan summary |\n",
    "| Localisation | LD clumping | Independent loci |\n",
    "| Fine-mapping | SuSiE / FINEMAP | 95% credible sets |\n",
    "\n",
    "See `Skills/gwas-population-genetics.md` for copy-paste code patterns.\n"
   ]
  }
 ]
}
```

- [ ] Step 5: Create `19_GWAS/README.md`

```markdown
# Module 19: Genome-Wide Association Studies (GWAS)

Patterns inspired by NGSchool 2023 practical

## Overview

This module covers the full GWAS workflow: study design, genotype quality control, population stratification correction via PCA, per-SNP association testing, and result interpretation with Manhattan and Q-Q plots.

## Prerequisites

- Module 02 (Python Essentials), Module 15 (Population Genetics), Module 02 (Variant Calling)

## Learning Objectives

After completing this module you will be able to:

- Apply genotype QC filters (MAF, call rate, HWE)
- Detect population stratification and include PC covariates
- Run per-SNP logistic regression and interpret odds ratios
- Produce Manhattan and Q-Q plots and interpret the genomic inflation factor
- Describe LD clumping and fine-mapping concepts

## Skill File

See `Skills/gwas-population-genetics.md` for quick-reference patterns.
```

- [ ] Step 6: Create `20_spatial_transcriptomics.ipynb`

```json
{
 "nbformat": 4,
 "nbformat_minor": 5,
 "metadata": {
  "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
  "language_info": {"name": "python", "version": "3.11.0"}
 },
 "cells": [
  {
   "cell_type": "markdown",
   "id": "header-spatial",
   "metadata": {},
   "source": [
    "# Module 20: Spatial Transcriptomics\n",
    "\n",
    "Patterns inspired by NGSchool 2023 practical\n",
    "\n",
    "---\n",
    "\n",
    "## Learning Objectives\n",
    "\n",
    "By the end of this notebook you will be able to:\n",
    "\n",
    "1. Explain spatial transcriptomics data formats (AnnData + spatial coordinates)\n",
    "2. Perform QC, normalization, and dimensionality reduction for spatial data\n",
    "3. Build a spatial neighborhood graph with squidpy\n",
    "4. Identify spatially variable genes using Moran's I\n",
    "5. Visualize gene expression and clusters overlaid on tissue coordinates\n",
    "6. Describe cell-type deconvolution approaches (RCTD, cell2location)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "imports-spatial",
   "metadata": {},
   "outputs": [],
   "source": [
    "import warnings\n",
    "warnings.filterwarnings('ignore')\n",
    "\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt\n",
    "import scanpy as sc\n",
    "import squidpy as sq\n",
    "\n",
    "sc.settings.verbosity = 0\n",
    "print(f'scanpy  {sc.__version__}')\n",
    "print(f'squidpy {sq.__version__}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-spatial-formats",
   "metadata": {},
   "source": [
    "## 1. Spatial Data Formats\n",
    "\n",
    "Spatial transcriptomics extends the single-cell AnnData object with two spatial additions:\n",
    "\n",
    "| Field | Location | Contents |\n",
    "|---|---|---|\n",
    "| Spatial coordinates | `adata.obsm[\"spatial\"]` | (n_spots × 2) array of (x, y) positions |\n",
    "| Tissue image | `adata.uns[\"spatial\"][lib_id]` | H&E or fluorescence image + scale factors |\n",
    "\n",
    "**10x Visium** spots are arranged in a hexagonal grid covering a tissue section (~5,000 spots, 55 µm diameter each). Each spot captures mRNA from multiple cells, so downstream analysis includes **deconvolution** to infer cell-type proportions.\n",
    "\n",
    "We load a built-in squidpy demo dataset (mouse cortex, H&E staining).\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "load-data-spatial",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Built-in cropped Visium H&E dataset — no external download needed\n",
    "adata = sq.datasets.visium_hne_adata_crop()\n",
    "print(adata)\n",
    "print(f'\\nSpatial coordinates shape: {adata.obsm[\"spatial\"].shape}')\n",
    "print(f'First 3 spots:\\n{adata.obsm[\"spatial\"][:3]}')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "qc-spatial",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Quality control\n",
    "sc.pp.filter_cells(adata, min_counts=200)\n",
    "sc.pp.filter_genes(adata, min_cells=5)\n",
    "\n",
    "adata.var['mt'] = adata.var_names.str.startswith('mt-')\n",
    "sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)\n",
    "\n",
    "# Remove high-MT spots (likely damaged)\n",
    "adata = adata[adata.obs['pct_counts_mt'] < 20].copy()\n",
    "\n",
    "fig, axes = plt.subplots(1, 3, figsize=(12, 3))\n",
    "for ax, key in zip(axes, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']):\n",
    "    ax.hist(adata.obs[key], bins=40, color='#4393c3', edgecolor='none')\n",
    "    ax.set_title(key); ax.set_xlabel('')\n",
    "plt.tight_layout(); plt.show()\n",
    "print(f'Spots after QC: {adata.n_obs}  Genes: {adata.n_vars}')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "normalize-spatial",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Normalize and reduce dimensions (same as standard scRNA-seq pipeline)\n",
    "sc.pp.normalize_total(adata, target_sum=1e4)\n",
    "sc.pp.log1p(adata)\n",
    "sc.pp.highly_variable_genes(adata, n_top_genes=2000)\n",
    "adata.raw = adata\n",
    "\n",
    "sc.pp.pca(adata, n_comps=30)\n",
    "sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)\n",
    "sc.tl.umap(adata)\n",
    "sc.tl.leiden(adata, resolution=0.4)\n",
    "\n",
    "print(f'Leiden clusters: {adata.obs[\"leiden\"].nunique()}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-neighbors",
   "metadata": {},
   "source": [
    "## 2. Spatial Neighborhood Graph\n",
    "\n",
    "A **spatial neighborhood graph** connects each spot to its physical neighbours (not transcriptome neighbours). This is used to:\n",
    "\n",
    "- Test spatial autocorrelation (Moran's I)\n",
    "- Smooth expression values spatially\n",
    "- Find spatially co-expressed gene modules\n",
    "\n",
    "`squidpy` builds this from `.obsm[\"spatial\"]` using:\n",
    "- `coord_type=\"grid\"` → fixed 6-neighbour Visium grid\n",
    "- `coord_type=\"generic\"` → radius- or kNN-based for irregular layouts (Xenium, MERFISH)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "spatial-graph",
   "metadata": {},
   "outputs": [],
   "source": [
    "sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6)\n",
    "print('Spatial graph adjacency matrix shape:',\n",
    "      adata.obsp['spatial_connectivities'].shape)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-svgs",
   "metadata": {},
   "source": [
    "## 3. Spatially Variable Genes\n",
    "\n",
    "**Moran's I** measures spatial autocorrelation: do nearby spots tend to have similar expression?\n",
    "\n",
    "$$I = \\frac{n}{W} \\cdot \\frac{\\sum_i \\sum_j w_{ij}(x_i - \\bar{x})(x_j - \\bar{x})}{\\sum_i (x_i - \\bar{x})^2}$$\n",
    "\n",
    "where $w_{ij}$ are spatial weights (1 if neighbours, 0 otherwise), $W = \\sum w_{ij}$.\n",
    "\n",
    "- **I ≈ +1**: Strong positive spatial autocorrelation (gene expressed in contiguous regions)\n",
    "- **I ≈ 0**: No spatial pattern (random)\n",
    "- **I ≈ −1**: Checkerboard pattern (very rare in transcriptomics)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "morans-i",
   "metadata": {},
   "outputs": [],
   "source": [
    "sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=1)\n",
    "svgs = adata.uns['moranI'].sort_values('I', ascending=False)\n",
    "print('Top 10 spatially variable genes:')\n",
    "print(svgs[['I', 'pval_norm']].head(10).to_string())"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "visualize-spatial",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Spatial scatter: cluster labels + top SVG\n",
    "top_svg = svgs.index[0]\n",
    "sq.pl.spatial_scatter(\n",
    "    adata,\n",
    "    color=['leiden', top_svg],\n",
    "    wspace=0.4,\n",
    "    title=['Leiden clusters', f'Top SVG: {top_svg}'],\n",
    ")\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-deconv",
   "metadata": {},
   "source": [
    "## 4. Cell-Type Deconvolution (Concept)\n",
    "\n",
    "Each Visium spot captures mRNA from multiple cells (average ~5–10 for mouse brain). **Deconvolution** infers the cell-type composition of each spot.\n",
    "\n",
    "**Popular methods:**\n",
    "\n",
    "| Method | Approach | Notes |\n",
    "|---|---|---|\n",
    "| RCTD (`spacexr`) | Likelihood model with Poisson noise | Fast; requires scRNA-seq reference |\n",
    "| cell2location | Negative-binomial regression | Handles multiple cell types per spot |\n",
    "| CARD | Conditional autoregressive model | Spatially aware |\n",
    "\n",
    "**Requirements:** A scRNA-seq reference dataset from the same tissue + species, with cell-type labels. The deconvolution algorithm learns a signature matrix from the reference and fits proportions to each spot.\n",
    "\n",
    "**Output:** `adata.obs[\"cell_type_X_proportion\"]` — fraction of spot signal attributed to each type.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "exercise-spatial",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Exercise 20\n",
    "# ───────────\n",
    "# 1. Change the Leiden resolution from 0.4 to 0.8. How does the spatial pattern\n",
    "#    of clusters change? Do finer clusters correspond to distinct tissue layers?\n",
    "#\n",
    "# 2. Plot the 3rd and 5th spatially variable genes on the tissue.\n",
    "#    Are they co-localised with specific Leiden clusters?\n",
    "#\n",
    "# 3. Compute Moran's I using coord_type='generic' (requires specifying n_neighs=8).\n",
    "#    Do the top SVGs change? Why might the neighborhood definition matter?\n",
    "\n",
    "# Your code here\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "summary-spatial",
   "metadata": {},
   "source": [
    "## Summary\n",
    "\n",
    "| Step | Tool | Key output |\n",
    "|---|---|---|\n",
    "| Load data | `sq.datasets.*` | AnnData with `.obsm['spatial']` |\n",
    "| QC + normalize | scanpy | Filtered, log-normalized AnnData |\n",
    "| Cluster | scanpy Leiden | Transcriptomic clusters |\n",
    "| Spatial graph | `sq.gr.spatial_neighbors` | `.obsp['spatial_connectivities']` |\n",
    "| SVGs | `sq.gr.spatial_autocorr` | Moran's I per gene |\n",
    "| Visualize | `sq.pl.spatial_scatter` | Tissue overlay plots |\n",
    "\n",
    "See `Skills/spatial-transcriptomics.md` for copy-paste patterns.\n"
   ]
  }
 ]
}
```

- [ ] Step 7: Create `20_Spatial_Transcriptomics/README.md`

```markdown
# Module 20: Spatial Transcriptomics

Patterns inspired by NGSchool 2023 practical

## Overview

This module covers spatial transcriptomics data analysis using squidpy and scanpy: loading AnnData with spatial coordinates, quality control, normalization, spatial neighborhood graphs, spatially variable gene detection via Moran's I, and tissue-overlay visualization.

## Prerequisites

- Module 03 (RNA-seq Analysis), Module 12 (Single-Cell fundamentals)

## Learning Objectives

After completing this module you will be able to:

- Load and inspect spatial AnnData objects (coordinates, images, scale factors)
- Apply QC and normalization for spatial data
- Build spatial neighborhood graphs and detect spatially variable genes
- Visualize gene expression and clusters on tissue coordinates
- Explain cell-type deconvolution approaches

## Skill File

See `Skills/spatial-transcriptomics.md` for quick-reference patterns.
```

- [ ] Step 8: Create `21_copy_number_analysis.ipynb`

```json
{
 "nbformat": 4,
 "nbformat_minor": 5,
 "metadata": {
  "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
  "language_info": {"name": "python", "version": "3.11.0"}
 },
 "cells": [
  {
   "cell_type": "markdown",
   "id": "header-cnv",
   "metadata": {},
   "source": [
    "# Module 21: DNA Copy Number Analysis\n",
    "\n",
    "Patterns inspired by NGSchool 2023 practical\n",
    "\n",
    "---\n",
    "\n",
    "## Learning Objectives\n",
    "\n",
    "By the end of this notebook you will be able to:\n",
    "\n",
    "1. Define copy number variation (CNV): gains, losses, LOH\n",
    "2. Normalize read-depth data to remove GC content and mappability biases\n",
    "3. Apply circular binary segmentation (CBS) to detect breakpoints\n",
    "4. Call integer copy number states from segment log-ratios\n",
    "5. Visualize genome-wide copy number profiles\n",
    "6. Annotate segments with gene-level copy number calls\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "imports-cnv",
   "metadata": {},
   "outputs": [],
   "source": [
    "import warnings\n",
    "warnings.filterwarnings('ignore')\n",
    "\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt\n",
    "from scipy import stats\n",
    "\n",
    "rng = np.random.default_rng(42)\n",
    "print('Libraries loaded.')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-cnv",
   "metadata": {},
   "source": [
    "## 1. Copy Number Variation\n",
    "\n",
    "In cancer genomes, chromosomal regions are frequently **gained** (more than 2 copies) or **lost** (fewer than 2 copies). **Loss of heterozygosity (LOH)** occurs when one allele of a diploid region is lost, sometimes exposing a deleterious recessive mutation.\n",
    "\n",
    "**Detection strategy (WGS/WES):**\n",
    "1. Compute read depth in non-overlapping bins across the genome\n",
    "2. Normalize for GC content, mappability, and sample coverage\n",
    "3. Compute log₂(tumour depth / normal depth) → the **log-ratio** track\n",
    "4. Segment the log-ratio track into piecewise-constant regions\n",
    "5. Call copy number states: log-ratio ≈ 0 → diploid; ≈ +1 → one extra copy; ≈ −∞ → homozygous deletion\n",
    "\n",
    "We simulate a tumour vs. normal log-ratio track with planted copy number events.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "simulate-cnv",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Simulate a genome: 22 chromosomes, 100 bins each (10 kb bins → 1 Mb per chrom)\n",
    "N_CHROMS, N_BINS = 22, 100\n",
    "BIN_SIZE = 10_000\n",
    "\n",
    "bins = pd.DataFrame({\n",
    "    'chrom': np.repeat(np.arange(1, N_CHROMS+1), N_BINS),\n",
    "    'start': np.tile(np.arange(N_BINS) * BIN_SIZE, N_CHROMS),\n",
    "    'end':   np.tile((np.arange(N_BINS)+1) * BIN_SIZE, N_CHROMS),\n",
    "})\n",
    "n_total = len(bins)\n",
    "\n",
    "# Diploid log-ratio baseline + noise\n",
    "log_ratio = rng.normal(0.0, 0.15, n_total)\n",
    "\n",
    "# Plant copy number events\n",
    "events = [\n",
    "    (1,  10, 50, +1.0),   # chr1 bins 10-50: gain (+1 copy)\n",
    "    (3,  30, 70, -1.0),   # chr3 bins 30-70: loss (1 copy)\n",
    "    (7,  60, 90, +2.0),   # chr7 bins 60-90: amplification\n",
    "    (12, 20, 45, -8.0),   # chr12 bins 20-45: homozygous deletion\n",
    "    (18, 50, 80, +0.6),   # chr18 bins 50-80: subclonal gain\n",
    "]\n",
    "\n",
    "for chrom, b_start, b_end, delta in events:\n",
    "    idx = (bins['chrom'] == chrom) & \\\n",
    "          (bins['start'] >= b_start*BIN_SIZE) & \\\n",
    "          (bins['end']   <= b_end*BIN_SIZE)\n",
    "    log_ratio[idx] += delta\n",
    "\n",
    "bins['log_ratio'] = log_ratio\n",
    "print(f'Total bins: {n_total}')\n",
    "print(bins.head())"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-segmentation",
   "metadata": {},
   "source": [
    "## 2. Segmentation — Circular Binary Segmentation\n",
    "\n",
    "The **CBS algorithm** (Olshen et al. 2004) recursively splits the log-ratio track at positions where the mean changes significantly. It uses a permutation test to assess each candidate breakpoint.\n",
    "\n",
    "The classic implementation is in the R `DNAcopy` package. In Python we can approximate CBS with the `ruptures` package (change-point detection) or with a simple recursive mean-shift algorithm.\n",
    "\n",
    "Here we implement a lightweight segmentation using `ruptures` with the L2 cost function (equivalent to minimizing within-segment variance).\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "segment-cnv",
   "metadata": {},
   "outputs": [],
   "source": [
    "try:\n",
    "    import ruptures as rpt\n",
    "    HAS_RUPTURES = True\n",
    "except ImportError:\n",
    "    HAS_RUPTURES = False\n",
    "    print('ruptures not installed — using a simple rolling-window segmentation fallback')\n",
    "\n",
    "def segment_chromosome(lr: np.ndarray, n_bkps_max: int = 10) -> list[int]:\n",
    "    \"\"\"Return breakpoint indices using ruptures Pelt or a simple fallback.\"\"\"\n",
    "    if HAS_RUPTURES:\n",
    "        model = rpt.Pelt(model='l2', min_size=5, jump=1)\n",
    "        model.fit(lr.reshape(-1, 1))\n",
    "        bkps = model.predict(pen=0.5)\n",
    "        return bkps\n",
    "    else:\n",
    "        # Fallback: detect jumps > 0.5 in 10-bin rolling mean\n",
    "        rolling = pd.Series(lr).rolling(10, center=True).mean().fillna(method='bfill').fillna(method='ffill').values\n",
    "        diff = np.abs(np.diff(rolling))\n",
    "        bkps = list(np.where(diff > 0.4)[0]) + [len(lr)]\n",
    "        return sorted(set(bkps))\n",
    "\n",
    "# Segment per chromosome\n",
    "segment_records = []\n",
    "for chrom in range(1, N_CHROMS+1):\n",
    "    mask = bins['chrom'] == chrom\n",
    "    lr   = bins.loc[mask, 'log_ratio'].values\n",
    "    bkps = segment_chromosome(lr)\n",
    "    prev = 0\n",
    "    for bkp in bkps:\n",
    "        seg_lr = lr[prev:bkp]\n",
    "        if len(seg_lr) == 0: continue\n",
    "        segment_records.append({\n",
    "            'chrom': chrom,\n",
    "            'start': bins.loc[mask].iloc[prev]['start'],\n",
    "            'end':   bins.loc[mask].iloc[min(bkp-1, len(lr)-1)]['end'],\n",
    "            'log_ratio_mean': seg_lr.mean(),\n",
    "            'n_bins': len(seg_lr),\n",
    "        })\n",
    "        prev = bkp\n",
    "\n",
    "segments = pd.DataFrame(segment_records)\n",
    "print(f'Segments detected: {len(segments)}')\n",
    "print(segments.head(10).to_string(index=False))"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-cn-calling",
   "metadata": {},
   "source": [
    "## 3. Copy Number State Calling\n",
    "\n",
    "Convert segment log-ratios to integer copy number states:\n",
    "\n",
    "| Log₂ ratio | Interpretation | CN state |\n",
    "|---|---|---|\n",
    "| ≤ −3 | Homozygous deletion | 0 |\n",
    "| −3 to −0.7 | Heterozygous deletion | 1 |\n",
    "| −0.3 to +0.3 | Diploid (neutral) | 2 |\n",
    "| +0.4 to +0.8 | Single copy gain | 3 |\n",
    "| > +0.8 | Amplification | ≥ 4 |\n",
    "\n",
    "These thresholds assume a pure tumour (no normal contamination). In practice, tumour purity and ploidy corrections are applied first (ASCAT, PURPLE).\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "cn-calling",
   "metadata": {},
   "outputs": [],
   "source": [
    "def call_cn(log2r: float) -> int:\n",
    "    if   log2r <= -3.0: return 0\n",
    "    elif log2r <= -0.7: return 1\n",
    "    elif log2r <=  0.3: return 2\n",
    "    elif log2r <=  0.8: return 3\n",
    "    else:               return 4\n",
    "\n",
    "segments['cn'] = segments['log_ratio_mean'].apply(call_cn)\n",
    "print('CN state distribution:')\n",
    "print(segments['cn'].value_counts().sort_index())"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "plot-cnv",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Genome-wide CN profile\n",
    "cn_colors = {0: '#d73027', 1: '#fc8d59', 2: '#ffffbf', 3: '#91cf60', 4: '#1a9850'}\n",
    "\n",
    "fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6), sharex=False)\n",
    "\n",
    "x_off = 0; xticks, xlabs = [], []\n",
    "for chrom in range(1, N_CHROMS+1):\n",
    "    b_mask   = bins['chrom'] == chrom\n",
    "    s_mask   = segments['chrom'] == chrom\n",
    "    xs       = bins.loc[b_mask, 'start'].values + x_off\n",
    "    lr       = bins.loc[b_mask, 'log_ratio'].values\n",
    "\n",
    "    # Raw log-ratio points\n",
    "    ax1.scatter(xs, lr, s=1, alpha=0.4, c='#555555', linewidths=0)\n",
    "\n",
    "    # Segment means\n",
    "    for _, row in segments.loc[s_mask].iterrows():\n",
    "        seg_xs = [row['start'] + x_off, row['end'] + x_off]\n",
    "        ax1.plot(seg_xs, [row['log_ratio_mean']]*2,\n",
    "                 color=cn_colors[row['cn']], lw=2.5)\n",
    "        ax2.barh(chrom, row['end']-row['start'],\n",
    "                 left=row['start']+x_off, color=cn_colors[row['cn']],\n",
    "                 height=0.8, align='center')\n",
    "\n",
    "    xticks.append(xs.mean()); xlabs.append(str(chrom))\n",
    "    x_off += bins.loc[b_mask, 'end'].values.max() + BIN_SIZE*5\n",
    "\n",
    "ax1.axhline(0, color='black', lw=0.5, ls='--')\n",
    "ax1.set_ylabel('log₂ ratio'); ax1.set_title('Genome-wide copy number profile')\n",
    "ax1.set_xticks(xticks); ax1.set_xticklabels(xlabs, fontsize=7)\n",
    "\n",
    "ax2.set_yticks(range(1, N_CHROMS+1))\n",
    "ax2.set_yticklabels([str(c) for c in range(1, N_CHROMS+1)], fontsize=7)\n",
    "ax2.set_ylabel('Chromosome'); ax2.set_title('CN state map')\n",
    "ax2.set_xticks([])\n",
    "\n",
    "from matplotlib.patches import Patch\n",
    "legend_elements = [Patch(facecolor=cn_colors[k], label=f'CN={k}') for k in sorted(cn_colors)]\n",
    "ax2.legend(handles=legend_elements, loc='lower right', fontsize=8)\n",
    "plt.tight_layout(); plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "gene-annotation",
   "metadata": {},
   "source": [
    "## 4. Gene-Level Copy Number Annotation\n",
    "\n",
    "In practice, gene coordinates from Ensembl or UCSC are intersected with segments to assign each gene its copy number call. Genes in deleted regions (CN=0) may be tumour suppressor genes; genes in amplified regions (CN≥4) may be oncogenes.\n",
    "\n",
    "**Pattern:**\n",
    "```python\n",
    "import pyranges as pr  # or pybedtools\n",
    "\n",
    "genes = pr.read_gtf('Homo_sapiens.GRCh38.gtf')\n",
    "segs  = pr.from_dict({'Chromosome': segments['chrom'].astype(str),\n",
    "                       'Start': segments['start'],\n",
    "                       'End':   segments['end'],\n",
    "                       'CN':    segments['cn']})\n",
    "annotated = genes.join(segs)\n",
    "```\n",
    "\n",
    "Below we simulate a simple gene annotation step.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "gene-level-cnv",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Simulate gene coordinates (100 genes, random positions)\n",
    "fake_genes = pd.DataFrame({\n",
    "    'gene': [f'GENE{i:04d}' for i in range(100)],\n",
    "    'chrom': rng.integers(1, 23, 100),\n",
    "    'start': rng.integers(0, N_BINS-5, 100) * BIN_SIZE,\n",
    "    'end':   None,\n",
    "})\n",
    "fake_genes['end'] = fake_genes['start'] + 5*BIN_SIZE\n",
    "\n",
    "def assign_cn_to_gene(row, segs):\n",
    "    overlapping = segs[\n",
    "        (segs['chrom'] == row['chrom']) &\n",
    "        (segs['start'] <= row['end']) &\n",
    "        (segs['end']   >= row['start'])\n",
    "    ]\n",
    "    if overlapping.empty:\n",
    "        return 2\n",
    "    return int(overlapping.iloc[0]['cn'])\n",
    "\n",
    "fake_genes['cn'] = fake_genes.apply(\n",
    "    lambda r: assign_cn_to_gene(r, segments), axis=1\n",
    ")\n",
    "print('Genes by CN state:')\n",
    "print(fake_genes['cn'].value_counts().sort_index())\n",
    "print('\\nDeleted genes (CN=0):')\n",
    "print(fake_genes[fake_genes['cn'] == 0][['gene', 'chrom', 'start']].to_string(index=False))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "exercise-cnv",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Exercise 21\n",
    "# ───────────\n",
    "# 1. Add noise level rng.normal(0, 0.30) instead of 0.15.\n",
    "#    How does this affect the number of detected segments?\n",
    "#    What does this teach about the importance of tumour purity and sequencing depth?\n",
    "#\n",
    "# 2. Plot a single chromosome (e.g. chr3) showing both the raw log-ratio\n",
    "#    scatter and the CBS segment means, with the true breakpoints marked as\n",
    "#    vertical dashed lines.\n",
    "#\n",
    "# 3. (Challenge) Compute the fraction of the genome altered (FGA):\n",
    "#    sum of lengths of segments with CN ≠ 2 / total genome length.\n",
    "#    High FGA is associated with genomic instability in cancer.\n",
    "\n",
    "# Your code here\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "summary-cnv",
   "metadata": {},
   "source": [
    "## Summary\n",
    "\n",
    "| Step | Tool | Output |\n",
    "|---|---|---|\n",
    "| Read depth in bins | samtools / mosdepth | Per-bin depth |\n",
    "| Normalize log-ratio | GC / mappability correction | Centred log₂ track |\n",
    "| Segmentation | ruptures / DNAcopy | Segment table |\n",
    "| CN calling | Threshold rules | Integer CN per segment |\n",
    "| Visualisation | matplotlib | Genome-wide CN profile |\n",
    "| Gene annotation | pyranges / pybedtools | Per-gene CN |\n",
    "\n",
    "Copy number analysis extends the variant calling concepts in Module 02. For somatic SV detection, see Module 17 (Genome Assembly and Advanced NGS).\n"
   ]
  }
 ]
}
```

- [ ] Step 9: Create `21_Copy_Number_Analysis/README.md`

```markdown
# Module 21: DNA Copy Number Analysis

Patterns inspired by NGSchool 2023 practical

## Overview

This module covers the complete workflow for detecting and visualizing somatic copy number alterations from WGS/WES data: log-ratio computation, segmentation via circular binary segmentation, CN state calling, genome-wide visualization, and gene-level annotation.

## Prerequisites

- Module 02 (Variant Calling and SNP Analysis), Module 01 (NGS Fundamentals)

## Learning Objectives

After completing this module you will be able to:

- Explain CNV types: gains, losses, LOH, amplifications, homozygous deletions
- Normalize read-depth tracks and compute log-ratios
- Apply CBS segmentation to detect breakpoints
- Call integer copy number states from segment log-ratios
- Visualize genome-wide CN profiles and annotate at gene level

## Skill File

Copy number analysis extends `Skills/ngs-variant-calling.md`. No separate skill file is provided for this module.
```

- [ ] Step 10: Commit

```bash
git add Course/Tier_3_Applied_Bioinformatics/19_GWAS/ \
        Course/Tier_3_Applied_Bioinformatics/20_Spatial_Transcriptomics/ \
        Course/Tier_3_Applied_Bioinformatics/21_Copy_Number_Analysis/ \
        Skills/gwas-population-genetics.md \
        Skills/spatial-transcriptomics.md
git commit -m "feat(tier3): add GWAS (19), Spatial (20), CNV (21) modules and skills"
```

### Task 4: Bayesian Statistics in Python (Tier 3, Module 22)

**Files:**
- Create: `Course/Tier_3_Applied_Bioinformatics/22_Bayesian_Statistics_Python/22_bayesian_statistics_python.ipynb`
- Create: `Course/Tier_3_Applied_Bioinformatics/22_Bayesian_Statistics_Python/README.md`
- Create: `Skills/bayesian-python.md`

**Steps:**

- [ ] Step 1: Create directory

```bash
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_3_Applied_Bioinformatics/22_Bayesian_Statistics_Python
```

- [ ] Step 2: Create `Skills/bayesian-python.md`

```markdown
---
name: bayesian-python
description: Bayesian statistics in Python — linear models with pymc and statsmodels, prior specification, model comparison with WAIC/LOO-CV via arviz, mixed-effects models with bambi, GLMs (Poisson, NegBin, Bernoulli), and advanced topics (GLMM, zero-inflation, GAMs)
---

## When to Use

Use this skill when:
- Fitting Bayesian linear or generalized linear models (pymc, bambi)
- Specifying informative or weakly informative priors and running prior predictive checks
- Comparing models with information criteria: WAIC, LOO-CV (arviz)
- Fitting linear mixed-effects models (random intercepts/slopes) with bambi
- Fitting GLMs: Poisson, Negative Binomial, Bernoulli, Binomial (statsmodels or pymc)
- Handling zero-inflated count data or fitting GAM/GAMM models

## Quick Reference

| Task | Tool | Key function |
|---|---|---|\n| Bayesian linear model | `pymc` | `pm.Model()`, `pm.sample()` |
| Posterior summary | `arviz` | `az.summary()`, `az.plot_posterior()` |
| Prior predictive | `arviz` / pymc | `pm.sample_prior_predictive()` |
| Model comparison | `arviz` | `az.compare()` (WAIC or LOO) |
| Mixed-effects model | `bambi` | `bmb.Model(formula, data)` |
| GLM (freq.) | `statsmodels` | `smf.glm(formula, family=sm.families.*)` |
| GLM (Bayesian) | `pymc` | `pm.Poisson`, `pm.NegativeBinomial` |
| Zero-inflated | `pymc` | `pm.ZeroInflatedPoisson` |
| Collinearity | `statsmodels` | `variance_inflation_factor()` |

## Key Patterns

**Pattern 1: Bayesian linear model (pymc)**
```python
import pymc as pm
import arviz as az
import numpy as np

with pm.Model() as linear_model:
    # Weakly informative priors
    alpha = pm.Normal('alpha', mu=0, sigma=10)
    beta  = pm.Normal('beta',  mu=0, sigma=10)
    sigma = pm.HalfNormal('sigma', sigma=1)

    mu = alpha + beta * X
    y_obs = pm.Normal('y_obs', mu=mu, sigma=sigma, observed=y)

    idata = pm.sample(1000, tune=1000, target_accept=0.9, return_inferencedata=True)

az.summary(idata, var_names=['alpha', 'beta', 'sigma'])
```

**Pattern 2: Prior predictive check**
```python
with linear_model:
    prior_pred = pm.sample_prior_predictive(samples=200)

az.plot_ppc(az.from_pymc3(prior=prior_pred, model=linear_model),
            group='prior')
```

**Pattern 3: Model comparison with LOO-CV**
```python
with model_1:
    idata_1 = pm.sample(1000, tune=1000, return_inferencedata=True)
    pm.compute_log_likelihood(idata_1)

with model_2:
    idata_2 = pm.sample(1000, tune=1000, return_inferencedata=True)
    pm.compute_log_likelihood(idata_2)

comparison = az.compare({'model_1': idata_1, 'model_2': idata_2}, ic='loo')
az.plot_compare(comparison)
```

**Pattern 4: Mixed-effects model (bambi)**
```python
import bambi as bmb

model = bmb.Model('bill_length_mm ~ body_mass_g + (1|species)',
                  data=penguins.dropna())
idata = model.fit(draws=1000, tune=1000)
model.plot_posteriors()
```

**Pattern 5: Poisson GLM (statsmodels)**
```python
import statsmodels.formula.api as smf
import statsmodels as sm

glm = smf.glm('count ~ x1 + x2',
              data=df,
              family=sm.families.Poisson()).fit()
print(glm.summary())
print(np.exp(glm.params))  # incident rate ratios
```

**Pattern 6: Zero-inflated Poisson (pymc)**
```python
with pm.Model() as zip_model:
    psi   = pm.Beta('psi', alpha=1, beta=1)       # zero-inflation mixing weight
    mu    = pm.Gamma('mu', alpha=2, beta=1)        # Poisson rate
    y_obs = pm.ZeroInflatedPoisson('y_obs',
                                   psi=psi, mu=mu,
                                   observed=counts)
    idata = pm.sample(1000, tune=1000, return_inferencedata=True)
```

## Code Templates

**Template 1: Bayesian vs frequentist linear model comparison**
```python
import pymc as pm, arviz as az, numpy as np
import statsmodels.api as sm

# Frequentist
X_const = sm.add_constant(X)
ols = sm.OLS(y, X_const).fit()
print(ols.summary())

# Bayesian
with pm.Model() as bayes_lm:
    a = pm.Normal('a', 0, 10)
    b = pm.Normal('b', 0, 10)
    s = pm.HalfNormal('s', 1)
    pm.Normal('obs', a + b*X, s, observed=y)
    idata = pm.sample(1000, tune=500, progressbar=False)

print(az.summary(idata))
az.plot_posterior(idata)
```

## Common Pitfalls

1. **NUTS divergences**: Divergences indicate the sampler is exploring regions of high curvature. Increase `target_accept` to 0.95 or reparameterise (e.g., non-centered parameterisation for hierarchical models).

2. **R-hat > 1.01**: Indicates chains have not converged. Run more tuning steps (`tune=2000`) or increase `chains=4`.

3. **Prior too wide**: A flat or very wide prior allows implausible parameter values (e.g., negative standard deviations). Use `pm.HalfNormal` or `pm.Exponential` for scale parameters, never `pm.Uniform(0, 1000)`.

4. **LOO warning (high Pareto k)**: k > 0.7 for some observations means LOO is unreliable. Use `az.loo` with `pointwise=True` to identify influential observations; consider robust regression.

5. **bambi formula syntax**: bambi uses R-style formulae. Fixed effect: `y ~ x`. Random intercept: `y ~ x + (1|group)`. Random slope: `y ~ x + (x|group)`. Interaction: `y ~ x1:x2` or `y ~ x1*x2`.
```

- [ ] Step 3: Create `22_bayesian_statistics_python.ipynb`

```json
{
 "nbformat": 4,
 "nbformat_minor": 5,
 "metadata": {
  "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
  "language_info": {"name": "python", "version": "3.11.0"}
 },
 "cells": [
  {
   "cell_type": "markdown",
   "id": "header-bayes",
   "metadata": {},
   "source": [
    "# Module 22: Bayesian Statistics in Python\n",
    "\n",
    "Statistical concepts based on Fränzi Korner-Nievergelt's Applied Statistics course; Python implementation by course authors.\n",
    "\n",
    "---\n",
    "\n",
    "## Learning Objectives\n",
    "\n",
    "By the end of this notebook you will be able to:\n",
    "\n",
    "1. Contrast frequentist and Bayesian inference: posterior = likelihood × prior\n",
    "2. Specify informative and weakly informative priors; run prior predictive checks\n",
    "3. Diagnose convergence: R-hat, ESS, trace plots\n",
    "4. Fit multiple regression models and diagnose collinearity with VIF\n",
    "5. Compare models using WAIC and LOO-CV via arviz\n",
    "6. Fit linear mixed-effects models with random intercepts/slopes using bambi\n",
    "7. Fit Poisson, Negative Binomial, and Zero-Inflated GLMs with pymc\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "imports-bayes",
   "metadata": {},
   "outputs": [],
   "source": [
    "import warnings\n",
    "warnings.filterwarnings('ignore')\n",
    "\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt\n",
    "import statsmodels.api as sm\n",
    "import statsmodels.formula.api as smf\n",
    "import pymc as pm\n",
    "import arviz as az\n",
    "import bambi as bmb\n",
    "from palmerpenguins import load_penguins\n",
    "\n",
    "rng = np.random.default_rng(42)\n",
    "penguins = load_penguins().dropna()\n",
    "print(f'pymc {pm.__version__}  |  arviz {az.__version__}  |  bambi {bmb.__version__}')\n",
    "print(f'Penguins dataset: {len(penguins)} rows')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-bayes-freq",
   "metadata": {},
   "source": [
    "## 1. Frequentist vs Bayesian Framing\n",
    "\n",
    "**Frequentist:** Parameters are fixed, unknown constants. Data are random. A 95% confidence interval means 'if we repeated the experiment many times, 95% of the intervals would contain the true parameter.' The interval itself either does or does not contain the parameter — it is not a probability statement about the parameter.\n",
    "\n",
    "**Bayesian:** Parameters are uncertain. We encode prior knowledge as a *prior distribution* $p(\\theta)$. After observing data $y$, we update via Bayes' theorem:\n",
    "\n",
    "$$p(\\theta | y) = \\frac{p(y | \\theta) \\cdot p(\\theta)}{p(y)} \\propto \\underbrace{p(y|\\theta)}_{\\text{likelihood}} \\cdot \\underbrace{p(\\theta)}_{\\text{prior}}$$\n",
    "\n",
    "The result is a **posterior distribution** — a complete probability distribution over the parameter. A 95% **credible interval** literally means 'there is 95% probability the parameter lies in this range, given our model and data.'\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "bayes-vs-freq",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Predict penguin bill length from body mass\n",
    "X = penguins['body_mass_g'].values.astype(float)\n",
    "y = penguins['bill_length_mm'].values.astype(float)\n",
    "\n",
    "# Standardize\n",
    "X_mean, X_std = X.mean(), X.std()\n",
    "X_s = (X - X_mean) / X_std\n",
    "y_mean, y_std = y.mean(), y.std()\n",
    "y_s = (y - y_mean) / y_std\n",
    "\n",
    "# Frequentist OLS\n",
    "ols = sm.OLS(y_s, sm.add_constant(X_s)).fit()\n",
    "print('OLS coefficient for body_mass (standardized):', ols.params[1].round(4))\n",
    "print('OLS 95% CI:', ols.conf_int().iloc[1].round(4).tolist())\n",
    "\n",
    "# Bayesian linear model\n",
    "with pm.Model() as lm:\n",
    "    alpha = pm.Normal('alpha', mu=0, sigma=2)\n",
    "    beta  = pm.Normal('beta',  mu=0, sigma=2)\n",
    "    sigma = pm.HalfNormal('sigma', sigma=1)\n",
    "    mu    = alpha + beta * X_s\n",
    "    pm.Normal('obs', mu=mu, sigma=sigma, observed=y_s)\n",
    "    idata = pm.sample(1000, tune=1000, chains=2,\n",
    "                      target_accept=0.9, progressbar=False,\n",
    "                      random_seed=42)\n",
    "\n",
    "print('\\nBayesian posterior summary:')\n",
    "print(az.summary(idata, var_names=['alpha', 'beta', 'sigma'])[['mean','sd','hdi_3%','hdi_97%']])"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-priors",
   "metadata": {},
   "source": [
    "## 2. Prior Specification and Prior Predictive Checks\n",
    "\n",
    "**Weakly informative priors** are the sensible default: they rule out impossible values (negative standard deviations, probabilities > 1) without strongly constraining plausible effects. A good strategy:\n",
    "\n",
    "| Parameter type | Suggested prior | Reason |\n",
    "|---|---|---|\n",
    "| Intercept | `Normal(0, 10)` | On standardized scale, large |\n",
    "| Regression slope | `Normal(0, 2)` | On std. scale, slopes rarely > 2 |\n",
    "| Standard deviation | `HalfNormal(1)` or `Exponential(1)` | Positive-only; regularizes |\n",
    "| Probability | `Beta(2, 2)` | Centred on 0.5; allows extremes |\n",
    "| Count rate | `Gamma(2, 1)` | Positive, right-skewed |\n",
    "\n",
    "**Prior predictive check:** Sample from the prior and simulate data. Does the simulated data look plausible (correct scale, sign, range)? If not, tighten or shift your priors.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "prior-predictive",
   "metadata": {},
   "outputs": [],
   "source": [
    "with lm:\n",
    "    prior_pred = pm.sample_prior_predictive(samples=200, random_seed=0)\n",
    "\n",
    "fig, axes = plt.subplots(1, 2, figsize=(12, 4))\n",
    "\n",
    "# Prior predictive lines\n",
    "x_grid = np.linspace(X_s.min(), X_s.max(), 50)\n",
    "prior_a = prior_pred.prior['alpha'].values.flatten()[:200]\n",
    "prior_b = prior_pred.prior['beta'].values.flatten()[:200]\n",
    "for a, b in zip(prior_a[:50], prior_b[:50]):\n",
    "    axes[0].plot(x_grid, a + b*x_grid, alpha=0.1, color='#4393c3', lw=0.8)\n",
    "axes[0].scatter(X_s, y_s, s=8, alpha=0.5, c='black', zorder=5)\n",
    "axes[0].set_title('Prior predictive lines'); axes[0].set_xlabel('body_mass (std)')\n",
    "axes[0].set_ylabel('bill_length (std)')\n",
    "\n",
    "# Posterior trace\n",
    "az.plot_trace(idata, var_names=['beta'], axes=axes[1:], compact=True)\n",
    "plt.tight_layout(); plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-multiple-regression",
   "metadata": {},
   "source": [
    "## 3. Multiple Regression and Collinearity\n",
    "\n",
    "When predictors are correlated, their individual regression coefficients become unstable — the **collinearity** problem. The **variance inflation factor (VIF)** quantifies this: VIF > 5 warrants attention; VIF > 10 indicates severe collinearity.\n",
    "\n",
    "$$\\text{VIF}_j = \\frac{1}{1 - R^2_j}$$\n",
    "\n",
    "where $R^2_j$ is the R² from regressing predictor $j$ on all other predictors.\n",
    "\n",
    "**Solutions:** Center predictors, drop one of the correlated pair, use principal components regression, or use regularisation (Lasso/Ridge).\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "collinearity",
   "metadata": {},
   "outputs": [],
   "source": [
    "from statsmodels.stats.outliers_influence import variance_inflation_factor\n",
    "\n",
    "feature_cols = ['body_mass_g', 'flipper_length_mm', 'bill_depth_mm']\n",
    "X_multi = penguins[feature_cols].values.astype(float)\n",
    "X_multi_std = (X_multi - X_multi.mean(0)) / X_multi.std(0)\n",
    "X_multi_const = sm.add_constant(X_multi_std)\n",
    "\n",
    "# Frequentist multiple regression\n",
    "ols_multi = sm.OLS(y_s, X_multi_const).fit()\n",
    "print(ols_multi.summary().tables[1])\n",
    "\n",
    "# VIF\n",
    "vifs = [variance_inflation_factor(X_multi_std, i) for i in range(X_multi_std.shape[1])]\n",
    "print('\\nVIF per predictor:')\n",
    "for name, vif in zip(feature_cols, vifs):\n",
    "    print(f'  {name:25s}: VIF = {vif:.2f}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-model-comparison",
   "metadata": {},
   "source": [
    "## 4. Model Comparison: WAIC and LOO-CV\n",
    "\n",
    "**WAIC (Widely Applicable Information Criterion)** and **LOO-CV (leave-one-out cross-validation)** are fully Bayesian alternatives to AIC/BIC. They estimate out-of-sample predictive accuracy using the posterior.\n",
    "\n",
    "- LOO-CV is generally preferred over WAIC (more stable for small samples)\n",
    "- Lower is better (ELPD is the expected log predictive density; displayed as negative ELPD so lower = worse; `az.compare` sorts by best model first)\n",
    "- A difference > 2 SE between models is considered meaningful\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "model-comparison",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Model A: body_mass only\n",
    "with pm.Model() as model_a:\n",
    "    a = pm.Normal('a', 0, 2)\n",
    "    b = pm.Normal('b', 0, 2)\n",
    "    s = pm.HalfNormal('s', 1)\n",
    "    pm.Normal('obs', a + b*X_s, s, observed=y_s)\n",
    "    idata_a = pm.sample(1000, tune=500, chains=2, progressbar=False, random_seed=1)\n",
    "    pm.compute_log_likelihood(idata_a)\n",
    "\n",
    "# Model B: body_mass + flipper_length\n",
    "X2 = (penguins['flipper_length_mm'].values - penguins['flipper_length_mm'].mean()) / penguins['flipper_length_mm'].std()\n",
    "with pm.Model() as model_b:\n",
    "    a  = pm.Normal('a',  0, 2)\n",
    "    b1 = pm.Normal('b1', 0, 2)\n",
    "    b2 = pm.Normal('b2', 0, 2)\n",
    "    s  = pm.HalfNormal('s', 1)\n",
    "    pm.Normal('obs', a + b1*X_s + b2*X2, s, observed=y_s)\n",
    "    idata_b = pm.sample(1000, tune=500, chains=2, progressbar=False, random_seed=2)\n",
    "    pm.compute_log_likelihood(idata_b)\n",
    "\n",
    "comparison = az.compare({'body_mass_only': idata_a, 'body_mass+flipper': idata_b}, ic='loo')\n",
    "print(comparison)\n",
    "az.plot_compare(comparison, insample_dev=False)\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-mixed",
   "metadata": {},
   "source": [
    "## 5. Linear Mixed-Effects Models (bambi)\n",
    "\n",
    "When observations are grouped (penguins from 3 species, patients from multiple hospitals), treating groups as fixed effects wastes degrees of freedom and fails to generalise. **Random effects** model the group-level variation as draws from a population distribution.\n",
    "\n",
    "**Random intercept:** Each group gets its own intercept, but they share a common distribution — partial pooling between 'no pooling' and 'complete pooling.'\n",
    "\n",
    "**Random slope:** Each group can also have its own slope for a predictor.\n",
    "\n",
    "`bambi` provides an R-formula interface over pymc:\n",
    "- `y ~ x + (1|group)` → random intercepts per group\n",
    "- `y ~ x + (x|group)` → random intercepts and slopes\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "mixed-effects",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Random intercepts per species\n",
    "model_me = bmb.Model(\n",
    "    'bill_length_mm ~ body_mass_g + (1|species)',\n",
    "    data=penguins,\n",
    "    dropna=True,\n",
    ")\n",
    "idata_me = model_me.fit(draws=1000, tune=1000, chains=2,\n",
    "                         target_accept=0.9, progressbar=False,\n",
    "                         random_seed=42)\n",
    "\n",
    "# Species-specific intercepts\n",
    "species_intercepts = az.summary(\n",
    "    idata_me,\n",
    "    var_names=['1|species']\n",
    ")[['mean', 'sd', 'hdi_3%', 'hdi_97%']]\n",
    "print('Species random intercepts:')\n",
    "print(species_intercepts)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-glm",
   "metadata": {},
   "source": [
    "## 6. Generalized Linear Models\n",
    "\n",
    "GLMs extend linear regression to non-normal outcomes via a **link function** and **error distribution**:\n",
    "\n",
    "| Outcome | Distribution | Link | Use case |\n",
    "|---|---|---|---|\n",
    "| 0/1 binary | Bernoulli | logit | Classification, disease presence |\n",
    "| Count (no overdispersion) | Poisson | log | Gene counts, event counts |\n",
    "| Count (overdispersed) | Negative Binomial | log | RNA-seq counts, species counts |\n",
    "| Proportion | Binomial | logit | Fraction of successes |\n",
    "\n",
    "**Bayesian approach:** Place priors on the linear predictor coefficients. The link function is deterministic; the likelihood is the chosen distribution.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "glm-code",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Simulate overdispersed count data (e.g., species abundance per plot)\n",
    "n = 200\n",
    "x_pred  = rng.normal(0, 1, n)\n",
    "true_mu = np.exp(1.5 + 0.8 * x_pred)\n",
    "alpha_nb = 3.0  # overdispersion (lower = more overdispersed)\n",
    "counts  = rng.negative_binomial(alpha_nb, alpha_nb / (alpha_nb + true_mu))\n",
    "\n",
    "# Frequentist Negative Binomial\n",
    "glm_nb = smf.glm(\n",
    "    'counts ~ x_pred',\n",
    "    data=pd.DataFrame({'counts': counts, 'x_pred': x_pred}),\n",
    "    family=sm.families.NegativeBinomial()\n",
    ").fit()\n",
    "print('NegBin GLM coefficients (log scale):')\n",
    "print(glm_nb.params.round(3))\n",
    "print('IRR (incidence rate ratios):')\n",
    "print(np.exp(glm_nb.params).round(3))\n",
    "\n",
    "# Bayesian Negative Binomial (pymc)\n",
    "with pm.Model() as nb_model:\n",
    "    a     = pm.Normal('a', 0, 2)\n",
    "    b     = pm.Normal('b', 0, 2)\n",
    "    alpha = pm.Exponential('alpha', 1)\n",
    "    mu    = pm.math.exp(a + b * x_pred)\n",
    "    pm.NegativeBinomial('obs', mu=mu, alpha=alpha, observed=counts)\n",
    "    idata_nb = pm.sample(1000, tune=500, chains=2, progressbar=False, random_seed=7)\n",
    "\n",
    "print('\\nBayesian posterior (a, b):')\n",
    "print(az.summary(idata_nb, var_names=['a', 'b'])[['mean', 'sd', 'hdi_3%', 'hdi_97%']])"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "theory-advanced",
   "metadata": {},
   "source": [
    "## 7. Advanced Topics\n",
    "\n",
    "### Zero-Inflated Models\n",
    "\n",
    "When count data has excess zeros (more than Poisson/NegBin predicts), use a zero-inflated model: a mixture of a point mass at zero and a count distribution.\n",
    "\n",
    "```python\n",
    "with pm.Model():\n",
    "    psi = pm.Beta('psi', 1, 1)   # P(structural zero)\n",
    "    mu  = pm.Gamma('mu', 2, 1)\n",
    "    pm.ZeroInflatedPoisson('obs', psi=psi, mu=mu, observed=counts)\n",
    "```\n",
    "\n",
    "### Generalized Additive Models (GAMs)\n",
    "\n",
    "GAMs replace linear predictors with smooth functions:\n",
    "\n",
    "$$g(E[y]) = \\beta_0 + f_1(x_1) + f_2(x_2) + \\ldots$$\n",
    "\n",
    "The smooth functions $f_j$ are typically represented as B-splines with a smoothing penalty. In Python: `pygam` for frequentist GAMs; `pymc` with B-spline basis for Bayesian GAMs.\n",
    "\n",
    "### Meta-analysis\n",
    "\n",
    "A Bayesian hierarchical model where each study contributes an estimate with known uncertainty. The global effect is estimated with partial pooling across studies:\n",
    "\n",
    "```python\n",
    "with pm.Model():\n",
    "    mu_global = pm.Normal('mu_global', 0, 1)\n",
    "    tau       = pm.HalfNormal('tau', 1)\n",
    "    theta     = pm.Normal('theta', mu_global, tau, shape=n_studies)\n",
    "    pm.Normal('obs', theta, sigma=se, observed=effect_sizes)\n",
    "```\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "zero-inflated",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Zero-inflated Poisson example\n",
    "# Simulate: 30% structural zeros + Poisson(4)\n",
    "n_zip = 300\n",
    "structural_zero = rng.binomial(1, 0.3, n_zip).astype(bool)\n",
    "zip_counts = np.where(structural_zero, 0, rng.poisson(4, n_zip))\n",
    "\n",
    "print(f'Observed zeros:  {(zip_counts == 0).sum()} / {n_zip}')\n",
    "print(f'Expected from Poisson(4): {int(np.exp(-4)*n_zip)} structural zeros')\n",
    "\n",
    "with pm.Model() as zip_model:\n",
    "    psi  = pm.Beta('psi', 1, 5)\n",
    "    mu   = pm.Gamma('mu', alpha=4, beta=1)\n",
    "    pm.ZeroInflatedPoisson('obs', psi=psi, mu=mu, observed=zip_counts)\n",
    "    idata_zip = pm.sample(1000, tune=500, chains=2, progressbar=False, random_seed=3)\n",
    "\n",
    "print('\\nZero-inflation parameter ψ:')\n",
    "print(az.summary(idata_zip, var_names=['psi', 'mu'])[['mean', 'sd', 'hdi_3%', 'hdi_97%']])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "exercise-bayes",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Exercise 22\n",
    "# ───────────\n",
    "# 1. In the linear model (Section 1), change the prior on beta from Normal(0, 2)\n",
    "#    to Normal(0, 0.1). Run the prior predictive check. How does this change the\n",
    "#    implied range of regression lines? Is this prior appropriate?\n",
    "#\n",
    "# 2. Fit a random-intercepts model with bambi using island as the grouping factor\n",
    "#    instead of species. Compare the group-level intercepts. Which factor explains\n",
    "#    more variation in bill length?\n",
    "#\n",
    "# 3. (Challenge) Fit a Bayesian Poisson model to the simulated count data.\n",
    "#    Then fit a Negative Binomial model. Compare them with LOO-CV.\n",
    "#    When does the NB model win by more than 2 SE over Poisson?\n",
    "\n",
    "# Your code here\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "summary-bayes",
   "metadata": {},
   "source": [
    "## Summary\n",
    "\n",
    "| Topic | Key tool | Key concept |\n",
    "|---|---|---|\n",
    "| Bayesian linear model | pymc | Posterior = likelihood × prior |\n",
    "| Prior predictive check | pm.sample_prior_predictive | Validate prior plausibility |\n",
    "| Convergence | arviz R-hat, ESS | R-hat < 1.01, ESS > 400 |\n",
    "| Multiple regression | statsmodels OLS | VIF to diagnose collinearity |\n",
    "| Model comparison | az.compare LOO | ELPD; diff > 2 SE = meaningful |\n",
    "| Mixed effects | bambi | Random intercepts/slopes |\n",
    "| GLMs | statsmodels + pymc | NegBin for overdispersed counts |\n",
    "| Zero-inflated | pymc ZeroInflatedPoisson | Excess zeros |\n",
    "\n",
    "See `Skills/bayesian-python.md` for copy-paste code patterns.\n"
   ]
  }
 ]
}
```

- [ ] Step 4: Create `22_Bayesian_Statistics_Python/README.md`

```markdown
# Module 22: Bayesian Statistics in Python

Statistical concepts based on Fränzi Korner-Nievergelt's Applied Statistics course; Python implementation by course authors.

## Overview

This module converts key concepts from applied Bayesian statistics to Python. It covers frequentist vs. Bayesian framing, prior specification, multiple regression and collinearity diagnostics, model comparison with LOO-CV/WAIC, linear mixed-effects models with bambi, and GLMs (Poisson, Negative Binomial, Zero-Inflated).

## Prerequisites

- Module 06 (Statistics for Bioinformatics), Module 07 (Machine Learning for Biology)

## Learning Objectives

After completing this module you will be able to:

- Fit Bayesian linear models with pymc and interpret posterior summaries
- Run prior predictive checks and diagnose convergence
- Compare models using LOO-CV via arviz
- Fit random-intercept and random-slope models with bambi
- Fit Poisson, NegBin, and Zero-Inflated GLMs

## Skill File

See `Skills/bayesian-python.md` for quick-reference patterns.
```

- [ ] Step 5: Commit

```bash
git add Course/Tier_3_Applied_Bioinformatics/22_Bayesian_Statistics_Python/ \
        Skills/bayesian-python.md
git commit -m "feat(tier3): add Bayesian statistics Python module (22) and skill"
```

### Task 5: TF Footprinting & Chromatin Accessibility (Tier 3, Module 23)

**Files:**
- Create: `Course/Tier_3_Applied_Bioinformatics/23_TF_Footprinting/23_tf_footprinting.ipynb`
- Create: `Course/Tier_3_Applied_Bioinformatics/23_TF_Footprinting/README.md`

**Steps:**

- [ ] Step 1: Create directory

```bash
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_3_Applied_Bioinformatics/23_TF_Footprinting
```

- [ ] Step 2: Create `23_tf_footprinting.ipynb`

Notebook cells (write as standard `.ipynb` JSON with nbformat=4, kernelspec Python 3):

**Cell 1 — Markdown header:**
```
# Module 23: TF Footprinting & Chromatin Accessibility
Code patterns inspired by TotipotencyLab chromatin analysis workflows, Max Planck Institute of Biochemistry
Learning Objectives: (1) ATAC-seq fragment size distribution + NFR; (2) TF footprinting concept; (3) cut-site profiles around motifs; (4) footprint score; (5) pybedtools interval arithmetic; (6) accumulation plots.
```

**Cell 2 — Imports:**
```python
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
rng = np.random.default_rng(42)
```

**Cell 3 — Markdown: ATAC-seq Fragment Size Distribution (theory)**
Explain: NFR < 150 bp, mono-nucleosomal ~200 bp, di ~400 bp, tri ~600 bp. Nucleosomal ladder = QC pass.

**Cell 4 — Code: Simulate and plot fragment size distribution**
```python
N_FRAGS = 100_000
components = [(0.60, 80, 30), (0.25, 200, 25), (0.10, 400, 30), (0.05, 600, 35)]
frags = []
for frac, mu, sd in components:
    n = int(frac * N_FRAGS)
    sizes = rng.normal(mu, sd, n).astype(int)
    frags.extend(sizes[sizes > 30].tolist())
fragment_sizes = np.array(frags)
fig, ax = plt.subplots(figsize=(9, 4))
ax.hist(fragment_sizes, bins=np.arange(30, 800, 5), color='#4393c3', edgecolor='none', alpha=0.8)
ax.axvline(150, color='red', ls='--', lw=1, label='NFR cutoff (150 bp)')
ax.set_xlabel('Fragment size (bp)'); ax.set_ylabel('Count')
ax.set_title('ATAC-seq fragment size distribution'); ax.legend()
plt.tight_layout(); plt.show()
print(f"NFR fraction: {(fragment_sizes < 150).mean():.1%}")
```

**Cell 5 — Markdown: TF Footprinting (theory)**
Explain: TF bound → Tn5 excluded from motif site → depletion within motif + elevation in flanks. Profile = average Tn5 insertions around all motif occurrences.

**Cell 6 — Code: Simulate footprint profiles**
```python
WINDOW = 200; N_MOTIFS = 5000; MOTIF_LEN = 12
positions = np.arange(-WINDOW, WINDOW + 1)

def make_profile(n_motifs, depth=0.35):
    profile = np.exp(-0.5*(positions/80)**2)*2.5 + 0.3
    motif_mask = np.abs(positions) <= MOTIF_LEN//2
    profile[motif_mask] *= (1 - depth)
    profile += 0.15*np.cos(2*np.pi*positions/200)
    counts = rng.poisson(lam=np.maximum(profile, 0.01), size=(n_motifs, len(positions)))
    return counts.mean(0)

observed = make_profile(N_MOTIFS, depth=0.35)
control  = make_profile(N_MOTIFS, depth=0.0)

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(positions, gaussian_filter1d(observed, 3), label='Observed', c='#d73027', lw=2)
ax.plot(positions, gaussian_filter1d(control,  3), label='Control (shuffled)', c='#2166ac', lw=1.5, ls='--')
ax.axvspan(-MOTIF_LEN//2, MOTIF_LEN//2, alpha=0.12, color='gray', label='Motif')
ax.set_xlabel('Position relative to motif centre (bp)')
ax.set_ylabel('Mean Tn5 insertions')
ax.set_title('TF footprint profile'); ax.legend()
plt.tight_layout(); plt.show()
```

**Cell 7 — Code: Footprint score**
```python
def footprint_score(profile, positions, motif_half=6, flank_min=20, flank_max=100):
    centre = profile[np.abs(positions) <= motif_half].mean()
    flank  = profile[(np.abs(positions) >= flank_min) & (np.abs(positions) <= flank_max)].mean()
    return flank / (centre + 1e-9)

tf_depths = [('TF_A strong', 0.55), ('TF_B moderate', 0.35), ('TF_C weak', 0.15), ('TF_D none', 0.0)]
for name, d in tf_depths:
    p = gaussian_filter1d(make_profile(3000, depth=d), 2)
    print(f'{name:20s}: score = {footprint_score(p, positions):.3f}')
```

**Cell 8 — Markdown: pybedtools interval arithmetic (theory + syntax)**
Show key operations: `.intersect()`, `.slop()`, `.subtract()`, `.merge()`, `.coverage()`. Note: requires BEDTools installed; use pyranges as alternative.

**Cell 9 — Code: pybedtools patterns (show as string)**
```python
PATTERNS = """
import pybedtools
peaks  = pybedtools.BedTool('atac_peaks.bed')
motifs = pybedtools.BedTool('motif_hits.bed')
black  = pybedtools.BedTool('blacklist_hg38.bed')

# Remove blacklist, intersect with motifs, extend ±200 bp
clean = peaks.subtract(black)
peaks_with_motif = clean.intersect(motifs, u=True)
windows = motifs.slop(b=200, g=pybedtools.genome_registry.hg38).merge()
coverage = windows.coverage('atac.bam', d=True)
"""
print(PATTERNS)
```

**Cell 10 — Code: Accumulation plot**
```python
HALF_WIN = 300; xs = np.arange(-HALF_WIN, HALF_WIN+1)

def sim_accum(n, amp):
    signal = amp * np.exp(-0.5*(xs/50)**2)
    counts = rng.poisson(np.maximum(signal, 0.05), (n, len(xs)))
    return counts.mean(0), counts.std(0)/np.sqrt(n)

mc, sec = sim_accum(2000, 2.0)
mt, set_ = sim_accum(2000, 3.5)

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(xs, gaussian_filter1d(mt, 4), label='Treated', c='#d73027', lw=2)
ax.fill_between(xs, gaussian_filter1d(mt-set_, 4), gaussian_filter1d(mt+set_, 4), alpha=0.2, color='#d73027')
ax.plot(xs, gaussian_filter1d(mc, 4), label='Control', c='#2166ac', lw=2)
ax.fill_between(xs, gaussian_filter1d(mc-sec, 4), gaussian_filter1d(mc+sec, 4), alpha=0.2, color='#2166ac')
ax.axvline(0, color='black', lw=0.8, ls='--')
ax.set_xlabel('Position relative to feature centre (bp)'); ax.set_ylabel('Mean Tn5 insertions')
ax.set_title('Accumulation plot'); ax.legend()
plt.tight_layout(); plt.show()
```

**Cell 11 — Markdown: Exercise 23**
1. Change motif half-length to 18 bp — does footprint score change at same depth?
2. Compute footprint scores for treated vs control accumulation profiles.
3. Challenge: simulate 20% occupancy TF (80% unbound). What is the average footprint score?

- [ ] Step 3: Create `23_TF_Footprinting/README.md`

```markdown
# Module 23: TF Footprinting & Chromatin Accessibility

Code patterns inspired by TotipotencyLab chromatin analysis workflows, Max Planck Institute of Biochemistry

## Overview

ATAC-seq footprinting workflow: fragment size QC, Tn5 cut-site profiles around motif sites, footprint scoring, pybedtools interval arithmetic, and accumulation plots.

## Prerequisites

Module 10 (Sequence Motifs), Module 15 (Motif Discovery), Module 01 (NGS Fundamentals)

## Skill File

Extends `Skills/motif-discovery.md`. No separate skill file.
```

- [ ] Step 4: Commit

```bash
git add Course/Tier_3_Applied_Bioinformatics/23_TF_Footprinting/
git commit -m "feat(tier3): add TF footprinting module (23)"
```

---

### Task 6: Tier 5 — Modern AI for Science (LLM Fine-tuning, Vision RAG, Diffusion Models)

**Files:**
- Create: `Course/Tier_5_Modern_AI_for_Science/01_LLM_Finetuning/01_LLM_Finetuning.ipynb`
- Create: `Course/Tier_5_Modern_AI_for_Science/02_Vision_RAG/02_Vision_RAG.ipynb`
- Create: `Course/Tier_5_Modern_AI_for_Science/03_Diffusion_Generative_Models/03_Diffusion_Generative_Models.ipynb`
- Create: `Skills/llm-finetuning.md`
- Create: `Skills/vision-rag.md`
- Create: `Skills/diffusion-generative.md`

**Steps:**

- [ ] Step 1: Create directories

```bash
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_5_Modern_AI_for_Science/01_LLM_Finetuning
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_5_Modern_AI_for_Science/02_Vision_RAG
mkdir -p /Users/pavel/Documents/Python-scripts-for-everything/Course/Tier_5_Modern_AI_for_Science/03_Diffusion_Generative_Models
```

- [ ] Step 2: Create `Skills/llm-finetuning.md`

```markdown
---
name: llm-finetuning
description: LLM fine-tuning — LoRA adapter math, 4-bit NF4 quantization with bitsandbytes, chat template formatting, SFTTrainer workflow with trl, synthetic instruction data generation, training tips
---

## When to Use

Use this skill when:
- Fine-tuning a pre-trained LLM on a custom instruction/chat dataset
- Configuring LoRA (rank, alpha, target_modules) for parameter-efficient fine-tuning
- Applying 4-bit or 8-bit quantization to fit large models on consumer GPUs
- Formatting data with chat templates (system/user/assistant structure)
- Running SFTTrainer from `trl` for supervised fine-tuning
- Generating synthetic instruction data for domain adaptation

## Quick Reference

| Concept | Tool / Class | Key parameter |
|---|---|---|
| Quantization | `BitsAndBytesConfig` | `load_in_4bit=True`, `bnb_4bit_quant_type="nf4"` |
| LoRA config | `LoraConfig` (peft) | `r=16`, `lora_alpha=32`, `target_modules` |
| Apply LoRA | `get_peft_model(model, config)` | Adds trainable adapters |
| Training args | `SFTConfig` | `num_train_epochs`, `learning_rate`, `fp16` |
| SFT trainer | `SFTTrainer` (trl) | `dataset_text_field`, `max_seq_length` |
| Chat template | `tokenizer.apply_chat_template` | `tokenize=False` for preview |
| Save adapter | `model.save_pretrained(path)` | LoRA weights only (~10 MB) |
| Merge adapter | `model.merge_and_unload()` | Merge into base weights |

## Key Patterns

**Pattern 1: Load model with 4-bit quantization**
```python
from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig
import torch

bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_quant_type="nf4",
    bnb_4bit_compute_dtype=torch.float16,
    bnb_4bit_use_double_quant=True,
)
model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Meta-Llama-3-8B",
    quantization_config=bnb_config,
    device_map="auto",
)
tokenizer = AutoTokenizer.from_pretrained("meta-llama/Meta-Llama-3-8B")
tokenizer.pad_token = tokenizer.eos_token
```

**Pattern 2: Configure and apply LoRA**
```python
from peft import LoraConfig, get_peft_model, TaskType

lora_config = LoraConfig(
    r=16,
    lora_alpha=32,
    target_modules=["q_proj", "v_proj", "k_proj", "o_proj"],
    lora_dropout=0.05,
    bias="none",
    task_type=TaskType.CAUSAL_LM,
)
model.enable_input_require_grads()  # required for 4-bit + gradient checkpointing
model = get_peft_model(model, lora_config)
model.print_trainable_parameters()
```

**Pattern 3: Chat template formatting**
```python
messages = [
    {"role": "system",    "content": "You are a helpful bioinformatics assistant."},
    {"role": "user",      "content": "What is a p-value?"},
    {"role": "assistant", "content": "A p-value is the probability of observing..."},
]
formatted = tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=False)
```

**Pattern 4: SFTTrainer training loop**
```python
from trl import SFTTrainer, SFTConfig

trainer = SFTTrainer(
    model=model,
    args=SFTConfig(
        output_dir="./lora-checkpoint",
        num_train_epochs=3,
        per_device_train_batch_size=4,
        gradient_accumulation_steps=4,
        learning_rate=2e-4,
        fp16=True,
        max_seq_length=2048,
        report_to="none",
    ),
    train_dataset=dataset,
    tokenizer=tokenizer,
    dataset_text_field="text",
)
trainer.train()
```

**Pattern 5: Synthetic instruction data**
```python
from datasets import Dataset

def make_example(topic, question, answer):
    messages = [
        {"role": "system",    "content": f"You are an expert in {topic}."},
        {"role": "user",      "content": question},
        {"role": "assistant", "content": answer},
    ]
    return {"text": tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=False)}

rows = [make_example("genomics", "What is RNA-seq?", "RNA-seq measures transcript abundance...")]
dataset = Dataset.from_list(rows)
```

**Pattern 6: Inference after fine-tuning**
```python
# Switch model to inference mode (no gradient tracking needed)
model.train(False)
prompt = tokenizer.apply_chat_template(test_messages, tokenize=False, add_generation_prompt=True)
inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
with torch.no_grad():
    output = model.generate(**inputs, max_new_tokens=200, do_sample=False)
response = tokenizer.decode(output[0][inputs.input_ids.shape[1]:], skip_special_tokens=True)
```

## Code Templates

**Template 1: Minimal LoRA fine-tuning scaffold (Colab-ready)**
```python
# pip install transformers peft trl bitsandbytes accelerate datasets
from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig
from peft import LoraConfig, get_peft_model, TaskType
from trl import SFTTrainer, SFTConfig
from datasets import Dataset
import torch

MODEL_ID = "TinyLlama/TinyLlama-1.1B-Chat-v1.0"
bnb = BitsAndBytesConfig(load_in_4bit=True, bnb_4bit_compute_dtype=torch.float16)
model = AutoModelForCausalLM.from_pretrained(MODEL_ID, quantization_config=bnb, device_map="auto")
tokenizer = AutoTokenizer.from_pretrained(MODEL_ID)
tokenizer.pad_token = tokenizer.eos_token

lora = LoraConfig(r=8, lora_alpha=16, target_modules=["q_proj","v_proj"],
                  lora_dropout=0.05, task_type=TaskType.CAUSAL_LM)
model.enable_input_require_grads()
model = get_peft_model(model, lora)
model.print_trainable_parameters()

data = [{"text": "<|user|>What is RNA?<|assistant|>RNA carries genetic info..."}] * 50
ds = Dataset.from_list(data)

trainer = SFTTrainer(
    model=model, tokenizer=tokenizer, train_dataset=ds,
    args=SFTConfig(output_dir="/tmp/lora", num_train_epochs=1, fp16=True,
                   max_seq_length=256, report_to="none"),
    dataset_text_field="text",
)
trainer.train()
```

## Common Pitfalls

1. **No pad token**: LLaMA tokenizers lack a pad token. Always set `tokenizer.pad_token = tokenizer.eos_token` before training.

2. **LoRA + 4-bit gradients**: Call `model.enable_input_require_grads()` BEFORE `get_peft_model()` when using 4-bit quantization.

3. **Chat template mismatch**: Each model family has its own template. Never hardcode `<|user|>` for a different model family — always use the model's own tokenizer.

4. **Train on completions only**: Use `DataCollatorForCompletionOnlyLM` to set user/system token labels to `-100` so loss is computed only on assistant turns.

5. **LoRA alpha scaling**: The effective adapter magnitude is `lora_alpha / r`. A common heuristic: set `lora_alpha = 2*r`. Doubling alpha without changing r increases the learning rate for the adapter.
```

- [ ] Step 3: Create `Skills/vision-rag.md`

```markdown
---
name: vision-rag
description: Vision RAG pipelines — VLM architecture, ColPali late-interaction (MaxSim) document retrieval, RAG pipeline patterns (retrieve → inject → generate), Qwen2-VL inference, recall@k evaluation
---

## When to Use

Use this skill when:
- Building document Q&A systems that retrieve over page images (PDFs, scanned docs)
- Using vision-language models (VLMs) for multi-modal generation
- Implementing ColPali late interaction retrieval (MaxSim scoring)
- Constructing a RAG pipeline: embed pages → query → retrieve → generate
- Evaluating retrieval (recall@k) and generation faithfulness

## Quick Reference

| Component | Tool | Notes |
|---|---|---|
| Page embeddings | ColPali (`vidore/colpali-v1.2`) | Late interaction; patch token matrices |
| VLM inference | Qwen2-VL or InternVL | Open weights; T4 Colab viable |
| PDF → images | `pdf2image` / `pymupdf` | ≥ 150 DPI for ColPali |
| MaxSim score | `einsum('qd,pd->qp').max(1).sum()` | Per-query-token max over patches |
| Vector store | `faiss` or `usearch` | Store page embedding matrices |

## Key Patterns

**Pattern 1: ColPali embed a page**
```python
from colpali_engine.models import ColPali, ColPaliProcessor
import torch

model = ColPali.from_pretrained("vidore/colpali-v1.2", torch_dtype=torch.bfloat16, device_map="auto")
proc  = ColPaliProcessor.from_pretrained("vidore/colpali-v1.2")

inputs = proc.process_images([page_image]).to(model.device)
with torch.no_grad():
    emb = model(**inputs)   # (1, n_patches, dim)
```

**Pattern 2: MaxSim late-interaction score**
```python
def maxsim(q: torch.Tensor, d: torch.Tensor) -> float:
    # q: (q_tokens, dim)  d: (n_patches, dim)  — both L2-normalised
    q_n = q / (q.norm(dim=1, keepdim=True) + 1e-8)
    d_n = d / (d.norm(dim=1, keepdim=True) + 1e-8)
    return torch.einsum('qd,pd->qp', q_n, d_n).max(dim=1).values.sum().item()
```

**Pattern 3: Retrieve top-k pages**
```python
scores  = [maxsim(query_emb, d) for d in page_embeddings]
top_idx = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)[:3]
```

**Pattern 4: Generate with Qwen2-VL**
```python
from transformers import Qwen2VLForConditionalGeneration, AutoProcessor

vlm  = Qwen2VLForConditionalGeneration.from_pretrained("Qwen/Qwen2-VL-2B-Instruct",
                                                        torch_dtype=torch.float16, device_map="cuda")
proc = AutoProcessor.from_pretrained("Qwen/Qwen2-VL-2B-Instruct")

messages = [{"role": "user", "content":
    [{"type": "image"} for _ in top_images] + [{"type": "text", "text": query}]
}]
prompt = proc.apply_chat_template(messages, add_generation_prompt=True)
inputs = proc(text=prompt, images=top_images, return_tensors="pt").to(vlm.device)
output = vlm.generate(**inputs, max_new_tokens=512)
```

**Pattern 5: Recall@k evaluation**
```python
def recall_at_k(ranked: list[int], relevant: set[int], k: int) -> float:
    return len(set(ranked[:k]) & relevant) / len(relevant)
```

## Code Templates

**Template 1: CPU mock RAG pipeline**
```python
import torch, numpy as np

rng = np.random.default_rng(0)
N_PAGES, N_PATCHES, DIM = 20, 64, 128

page_embeddings = [torch.tensor(rng.normal(0,1,(N_PATCHES,DIM)), dtype=torch.float32) for _ in range(N_PAGES)]
query_emb = torch.tensor(rng.normal(0,1,(16,DIM)), dtype=torch.float32)

# Make page 7 relevant
page_embeddings[7][:16] = query_emb + torch.randn(16, DIM)*0.1

def maxsim(q, d):
    q_n = q / (q.norm(dim=1, keepdim=True) + 1e-8)
    d_n = d / (d.norm(dim=1, keepdim=True) + 1e-8)
    return torch.einsum('qd,pd->qp', q_n, d_n).max(dim=1).values.sum().item()

scores = [maxsim(query_emb, d) for d in page_embeddings]
ranked = sorted(range(N_PAGES), key=lambda i: scores[i], reverse=True)
print(f"Top-3: {ranked[:3]}  (page 7 should be rank 1)")
print(f"Recall@1: {recall_at_k(ranked, {7}, 1):.0%}  Recall@3: {recall_at_k(ranked, {7}, 3):.0%}")
```

## Common Pitfalls

1. **ColPali is NOT a bi-encoder**: Do not mean-pool patch embeddings and use cosine similarity — this loses the spatial information that makes ColPali work.

2. **PDF resolution**: Render at ≥ 150 DPI. Low-res pages lose text tokens and degrade retrieval significantly.

3. **VLM context overflow**: Passing > 5 high-res pages may exceed context. Use `max_pixels` in Qwen2-VL processor to cap resolution.

4. **`add_generation_prompt=True`**: Always set this during inference with `apply_chat_template`, otherwise the model repeats the user prompt.

5. **Faithfulness vs. fluency**: A fluent answer can still be hallucinated. Evaluate whether the answer is grounded in the retrieved pages separately from reading quality.
```

- [ ] Step 4: Create `Skills/diffusion-generative.md`

```markdown
---
name: diffusion-generative
description: Diffusion and score-based generative models — DDIM forward/reverse process, linear and cosine noise schedules, score matching intuition, SVD-based degradation operators for imaging inverse problems, score field visualization
---

## When to Use

Use this skill when:
- Implementing DDIM/DDPM forward and reverse processes from scratch
- Understanding or comparing linear vs. cosine noise schedules
- Formulating imaging inverse problems (denoising, inpainting, super-resolution) as posterior sampling
- Implementing SVD-based degradation operators (H = UΣVᵀ)
- Visualizing score fields for toy distributions
- Scientific applications: cryo-EM denoising, MRI reconstruction, fluorescence deconvolution

## Quick Reference

| Concept | Formula / Code | Notes |
|---|---|---|
| Forward (single step) | `xt = sqrt(ab)*x0 + sqrt(1-ab)*eps` | `ab = alpha_bars[t]` |
| DDIM reverse step | `x0_pred = (xt - sqrt(1-ab)*eps) / sqrt(ab)` | Clip x0_pred to [-1,1] |
| Linear schedule | `betas = linspace(1e-4, 0.02, T)` | |
| Cosine schedule | `ab = cos²(π/2·(t/T+s)/(1+s))` | s=0.008; smoother near t=0 |
| SNR | `10*log10(ab / (1-ab))` | Useful for schedule comparison |
| Degradation | `y = H @ x + noise` | H = UΣVᵀ via np.linalg.svd |
| Score function | `∇ₓ log p(x) ≈ -eps_theta(xt, t) / sqrt(1-ab)` | |

## Key Patterns

**Pattern 1: Linear noise schedule**
```python
import numpy as np

def linear_schedule(T, beta_start=1e-4, beta_end=0.02):
    betas = np.linspace(beta_start, beta_end, T)
    alphas = 1.0 - betas
    alpha_bars = np.cumprod(alphas)
    return {"betas": betas, "alphas": alphas, "alpha_bars": alpha_bars, "T": T}
```

**Pattern 2: Cosine schedule**
```python
def cosine_schedule(T, s=0.008):
    t  = np.linspace(0, T, T+1)
    f  = np.cos(((t/T + s) / (1+s)) * np.pi/2)**2
    ab = (f/f[0])[1:]
    betas = np.clip(1 - ab / np.concatenate([[1.0], ab[:-1]]), 0, 0.999)
    return {"betas": betas, "alphas": 1-betas, "alpha_bars": ab, "T": T}
```

**Pattern 3: Forward diffusion**
```python
def q_sample(x0, t, sched, rng):
    ab  = sched["alpha_bars"][t]
    eps = rng.standard_normal(x0.shape)
    return np.sqrt(ab)*x0 + np.sqrt(1-ab)*eps, eps
```

**Pattern 4: DDIM deterministic reverse step**
```python
def ddim_step(xt, eps_pred, t, t_prev, sched):
    ab_t    = sched["alpha_bars"][t]
    ab_prev = sched["alpha_bars"][t_prev] if t_prev >= 0 else 1.0
    x0_pred = np.clip((xt - np.sqrt(1-ab_t)*eps_pred) / np.sqrt(ab_t), -1, 1)
    return np.sqrt(ab_prev)*x0_pred + np.sqrt(1-ab_prev)*eps_pred
```

**Pattern 5: SVD degradation operator**
```python
def make_degradation(image_size, frac=0.25):
    n = image_size**2
    H = np.eye(n)[:int(n*frac)]   # observe frac of pixels
    U, s, Vt = np.linalg.svd(H, full_matrices=False)
    return {"H": H, "U": U, "s": s, "Vt": Vt}

def degrade(x_flat, op, sigma=0.1, rng=None):
    if rng is None: rng = np.random.default_rng(0)
    return op["H"] @ x_flat + rng.normal(0, sigma, op["H"].shape[0])
```

**Pattern 6: Score field visualization**
```python
import numpy as np, matplotlib.pyplot as plt

def gaussian_score(xy, mu, sigma=0.8):
    return -(xy - mu) / sigma**2

grid = 18
xg, yg = np.linspace(-4,4,grid), np.linspace(-3,3,grid)
X, Y = np.meshgrid(xg, yg)
pts  = np.stack([X.ravel(), Y.ravel()], axis=1)
sc   = gaussian_score(pts, np.array([1.0, 0.5]))

plt.figure(figsize=(6,5))
plt.quiver(X, Y, sc[:,0].reshape(grid,grid), sc[:,1].reshape(grid,grid), alpha=0.7)
plt.scatter(1.0, 0.5, s=100, c='red'); plt.title('Score field ∇ₓ log p(x)')
plt.show()
```

## Code Templates

**Template 1: Full DDIM forward+reverse demo (CPU)**
```python
import numpy as np, matplotlib.pyplot as plt

T   = 1000
rng = np.random.default_rng(42)
sched = linear_schedule(T)

# 16x16 checkerboard
size = 16
x0 = np.array([[1 if (i+j)%2==0 else -1 for j in range(size)] for i in range(size)], float).ravel()
xT, _ = q_sample(x0, T-1, sched, rng)

xt = xT.copy()
ts = list(range(T-1, -1, -T//50))
for i in range(len(ts)-1):
    _, eps = q_sample(x0, ts[i], sched, rng)   # oracle noise
    xt = ddim_step(xt, eps, ts[i], ts[i+1], sched)

fig, axes = plt.subplots(1,3, figsize=(9,3))
for ax, d, t in zip(axes, [x0, xT, xt], ['Original','Noisy','Reconstructed']):
    ax.imshow(d.reshape(size,size), cmap='gray', vmin=-1, vmax=1); ax.set_title(t); ax.axis('off')
plt.tight_layout(); plt.show()
```

## Common Pitfalls

1. **Off-by-one in schedule indexing**: `alpha_bars` must be indexed with `t ∈ {0,...,T-1}`. Check `schedule["alpha_bars"].shape == (T,)`.

2. **No clipping of x₀ prediction**: Without `np.clip(x0_pred, -1, 1)` the predicted image drifts in early steps, causing compounding errors.

3. **DDIM eta parameter**: Standard DDIM is deterministic (η=0). η>0 adds stochasticity; η=1 recovers DDPM. Match η to training assumptions.

4. **Inverse problem vs. unconditional**: Unconditional sampling ignores the measurement `y`. You need a data-consistency term (DDRM, DPS) to actually solve the inverse problem.

5. **SVD memory for large images**: Explicit H matrix is O(n²). For n > 64×64, use implicit operators (FFT-based convolution, masking) instead of dense SVD.
```

- [ ] Step 5: Create `01_LLM_Finetuning.ipynb`

Notebook cells (write as standard `.ipynb` JSON with nbformat=4, kernelspec Python 3, colab GPU tag):

**Header cell (markdown):**
```
# Module 01: LLM Fine-tuning
Inspired by Unsloth AI and Manuel Faysse's fine-tuning tutorials
GPU required: Colab T4. Enable: Runtime → Change runtime type → T4 GPU.
Learning Objectives: (1) base vs instruction vs chat models; (2) LoRA math;
(3) 4-bit NF4 quantization; (4) chat templates; (5) SFTTrainer; (6) synthetic data.
```

**Install cell (code, commented):**
```python
# !pip install -q transformers peft trl bitsandbytes accelerate datasets
```

**Imports cell (code):**
```python
import warnings; warnings.filterwarnings('ignore')
import torch, numpy as np
from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig
from peft import LoraConfig, get_peft_model, TaskType
from trl import SFTTrainer, SFTConfig
from datasets import Dataset
print(f'PyTorch {torch.__version__} | GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else "not available"}')
```

**Theory cell: Base vs Instruction vs Chat (markdown)**
Table: Base (causal LM pre-training → text continuation), Instruction (+ SFT → follows instructions), Chat (+ RLHF/DPO → multi-turn, aligned). Fine-tune from instruction model for domain adaptation.

**Theory cell: LoRA math (markdown)**
`W' = W + ΔW = W + BA` where B∈ℝ^(d×r), A∈ℝ^(r×k), r≪min(d,k). Parameters: r (rank), lora_alpha (=2r typical), target_modules, lora_dropout. Table: r=16 on 8B model → ~8M trainable (0.1%).

**LoRA parameter count visualisation cell (code):**
```python
import matplotlib.pyplot as plt
ranks = [4, 8, 16, 32, 64, 128]
params = [2 * 4096 * r * 32 * 4 for r in ranks]  # approx LLaMA-3-8B
fig, ax = plt.subplots(figsize=(7,4))
ax.bar(ranks, [p/1e6 for p in params], color='#4393c3', width=5)
ax.axhline(8000, color='red', ls='--', lw=1, label='8B base params')
ax.set_xlabel('LoRA rank r'); ax.set_ylabel('Trainable params (M)')
ax.set_title('LoRA parameter count vs rank'); ax.legend()
plt.tight_layout(); plt.show()
for r, p in zip(ranks, params):
    print(f'r={r:3d}: {p/1e6:.1f}M ({100*p/8e9:.3f}% of 8B)')
```

**Theory cell: Quantization (markdown)**
Table: fp32 4B/param 32GB, bf16 2B 16GB, int8 1B 8GB, nf4 0.5B 5GB. NF4 = NormalFloat4: optimal 4-bit values for normal-distribution weights. Double quantization reduces constants overhead to < 0.5 bit/param.

**Load model cell (code):**
```python
MODEL_ID = "TinyLlama/TinyLlama-1.1B-Chat-v1.0"
bnb_config = BitsAndBytesConfig(load_in_4bit=True, bnb_4bit_quant_type="nf4",
                                 bnb_4bit_compute_dtype=torch.float16, bnb_4bit_use_double_quant=True)
model = AutoModelForCausalLM.from_pretrained(MODEL_ID, quantization_config=bnb_config, device_map="auto")
tokenizer = AutoTokenizer.from_pretrained(MODEL_ID)
tokenizer.pad_token = tokenizer.eos_token
model.config.use_cache = False
print(f'Params: {sum(p.numel() for p in model.parameters())/1e9:.2f}B')
```

**Apply LoRA cell (code):**
```python
lora_config = LoraConfig(r=16, lora_alpha=32,
                          target_modules=["q_proj","v_proj","k_proj","o_proj"],
                          lora_dropout=0.05, bias="none", task_type=TaskType.CAUSAL_LM)
model.enable_input_require_grads()
model = get_peft_model(model, lora_config)
model.print_trainable_parameters()
```

**Theory cell: Chat Templates (markdown)**
Explain: each model family has its own template. Use tokenizer.apply_chat_template — never hardcode tags.

**Chat template demo cell (code):**
```python
sample = [
    {"role": "system",    "content": "You are a bioinformatics tutor."},
    {"role": "user",      "content": "What is RNA-seq?"},
    {"role": "assistant", "content": "RNA-seq (RNA sequencing) quantifies gene expression..."},
]
print(tokenizer.apply_chat_template(sample, tokenize=False, add_generation_prompt=False))
```

**Make dataset cell (code):**
```python
qa_pairs = [
    ("What is a p-value?", "A p-value is the probability of results at least as extreme under H₀."),
    ("What does RNA-seq measure?", "RNA-seq measures transcript abundance via read counts per gene."),
    ("What is a SNP?", "A SNP is a single base-pair variation at a specific genomic position."),
    ("Explain k-means.", "k-means assigns each point to one of k clusters minimising within-cluster variance."),
    ("What is CRISPR?", "CRISPR-Cas9 uses a guide RNA to direct Cas9 to a target DNA site for editing."),
] * 4

def fmt(q, a):
    msgs = [{"role":"system","content":"You are a bioinformatics tutor."},
            {"role":"user","content":q}, {"role":"assistant","content":a}]
    return {"text": tokenizer.apply_chat_template(msgs, tokenize=False, add_generation_prompt=False)}

dataset = Dataset.from_list([fmt(q,a) for q,a in qa_pairs])
print(f'Dataset: {len(dataset)} examples')
```

**Train cell (code):**
```python
trainer = SFTTrainer(
    model=model,
    args=SFTConfig(output_dir="/tmp/lora-bio", num_train_epochs=2,
                   per_device_train_batch_size=2, gradient_accumulation_steps=4,
                   learning_rate=2e-4, fp16=torch.cuda.is_available(),
                   logging_steps=5, save_strategy="no", max_seq_length=512, report_to="none"),
    train_dataset=dataset, tokenizer=tokenizer, dataset_text_field="text",
)
trainer.train()
print('Training complete.')
```

**Inference cell (code):**
```python
# Switch to inference mode
model.train(False)
test = [{"role":"system","content":"You are a bioinformatics tutor."},
        {"role":"user",  "content":"What is RNA-seq?"}]
prompt = tokenizer.apply_chat_template(test, tokenize=False, add_generation_prompt=True)
inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
with torch.no_grad():
    out = model.generate(**inputs, max_new_tokens=100, do_sample=False)
print(tokenizer.decode(out[0][inputs.input_ids.shape[1]:], skip_special_tokens=True))
```

**Exercise cell (code, commented):**
```python
# Exercise T5-01
# 1. Change LoRA r from 16 to 4. How does trainable param count change?
# 2. Add gate_proj, up_proj to target_modules. How many extra params?
# 3. Challenge: add 10 more Q&A pairs on a new topic; retrain and compare responses.
```

- [ ] Step 6: Create `02_Vision_RAG.ipynb`

Notebook cells (nbformat=4, Colab T4 tag):

**Header cell (markdown):** Module 02 Vision RAG. Attribution: Unsloth AI, Manuel Faysse. GPU recommended (T4). Learning Objectives: (1) VLM architecture; (2) page retrieval vs token retrieval; (3) ColPali MaxSim; (4) RAG pipeline; (5) recall@k evaluation.

**Imports cell (code):**
```python
import warnings; warnings.filterwarnings('ignore')
import numpy as np, torch, matplotlib.pyplot as plt
from PIL import Image
rng = np.random.default_rng(42)
print(f'PyTorch {torch.__version__} | GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else "not available"}')
```

**Theory cell: VLM architecture (markdown):** Table: visual encoder (CLIP ViT / SigLIP → patch tokens), projector (MLP → LLM hidden dim), LLM decoder (Qwen2-VL, InternVL). Patches: 14×14px each → ~500–1000 tokens per page. Full-corpus retrieval needs a retrieval step before generation.

**Theory cell: ColPali MaxSim (markdown):** Standard dense retrieval = single vector per doc (loses spatial info). ColPali: document page = matrix (n_patches × dim); query = matrix (q_tokens × dim). MaxSim: for each query token find max-sim patch, sum across tokens. Formula: score(q,d) = Σᵢ maxⱼ qᵢᵀdⱼ.

**ColPali CPU demo cell (code):**
```python
N_PAGES, N_PATCHES, DIM = 20, 64, 128
page_embeddings = [torch.tensor(rng.normal(0,1,(N_PATCHES,DIM)), dtype=torch.float32) for _ in range(N_PAGES)]
query_emb = torch.tensor(rng.normal(0,1,(16,DIM)), dtype=torch.float32)
# Make page 7 relevant
page_embeddings[7][:16] = query_emb + torch.randn(16, DIM)*0.1

def maxsim(q, d):
    q_n = q / (q.norm(dim=1,keepdim=True)+1e-8)
    d_n = d / (d.norm(dim=1,keepdim=True)+1e-8)
    return torch.einsum('qd,pd->qp',q_n,d_n).max(dim=1).values.sum().item()

def recall_at_k(ranked, relevant, k):
    return len(set(ranked[:k]) & set(relevant)) / len(relevant)

scores = [maxsim(query_emb, d) for d in page_embeddings]
ranked = sorted(range(N_PAGES), key=lambda i: scores[i], reverse=True)
print(f'Top-3: {ranked[:3]}')
for k in [1,3,5]:
    print(f'Recall@{k}: {recall_at_k(ranked, [7], k):.0%}')
```

**Visualise retrieval scores cell (code):**
```python
fig, ax = plt.subplots(figsize=(9,4))
ax.bar(range(N_PAGES), scores, color=['#d73027' if i==7 else '#4393c3' for i in range(N_PAGES)])
ax.set_xlabel('Page index'); ax.set_ylabel('MaxSim score')
ax.set_title('ColPali MaxSim scores (red = relevant page)')
plt.tight_layout(); plt.show()
```

**RAG scaffold cell (code):**
```python
class MockVLM:
    def generate(self, images, query, page_ids):
        return f"[Mock VLM] From {len(images)} pages — answer to '{query}': ..."

class DocumentRAG:
    def __init__(self, page_embs, vlm):
        self.embs, self.vlm = page_embs, vlm

    def retrieve(self, q_emb, k=3):
        sc = [maxsim(q_emb, d) for d in self.embs]
        return sorted(range(len(sc)), key=lambda i: sc[i], reverse=True)[:k]

    def answer(self, query, q_emb, page_images, k=3):
        idx = self.retrieve(q_emb, k)
        return self.vlm.generate([page_images[i] for i in idx], query, idx), idx

rag = DocumentRAG(page_embeddings, MockVLM())
mock_pages = [Image.new('RGB', (64,64), (i*10,50,100)) for i in range(N_PAGES)]
ans, retrieved = rag.answer("What is on page 7?", query_emb, mock_pages)
print(f'Retrieved: {retrieved}\nAnswer: {ans}')
```

**Real ColPali code cell (code, as string):**
```python
REAL_CODE = """
# pip install colpali-engine qwen-vl-utils transformers
# Run on Colab T4 GPU

from colpali_engine.models import ColPali, ColPaliProcessor
from transformers import Qwen2VLForConditionalGeneration, AutoProcessor
import torch

retriever = ColPali.from_pretrained("vidore/colpali-v1.2", torch_dtype=torch.bfloat16, device_map="cuda")
rproc = ColPaliProcessor.from_pretrained("vidore/colpali-v1.2")

# Embed pages
page_embs = []
for img in page_images:
    inp = rproc.process_images([img]).to(retriever.device)
    with torch.no_grad():
        page_embs.append(retriever(**inp)[0])

# Embed query
query = "What are the main findings in Figure 3?"
qinp = rproc.process_queries([query]).to(retriever.device)
with torch.no_grad():
    q_emb = retriever(**qinp)[0]

# Retrieve and generate
top3 = sorted(range(len(page_embs)), key=lambda i: maxsim(q_emb.cpu(), page_embs[i].cpu()), reverse=True)[:3]
vlm  = Qwen2VLForConditionalGeneration.from_pretrained("Qwen/Qwen2-VL-2B-Instruct",
                                                        torch_dtype=torch.float16, device_map="cuda")
vproc = AutoProcessor.from_pretrained("Qwen/Qwen2-VL-2B-Instruct")
msgs = [{"role":"user","content":[{"type":"image"} for _ in top3]+[{"type":"text","text":query}]}]
prompt = vproc.apply_chat_template(msgs, add_generation_prompt=True)
inputs = vproc(text=prompt, images=[page_images[i] for i in top3], return_tensors="pt").to(vlm.device)
out = vlm.generate(**inputs, max_new_tokens=512)
print(vproc.decode(out[0], skip_special_tokens=True))
"""
print(REAL_CODE)
```

**Exercise cell (code, commented):**
```python
# Exercise T5-02
# 1. Make pages 2 AND 11 relevant. Recompute Recall@1, @3, @5.
#    Why does Recall@1 drop vs single-relevant-page case?
# 2. Remove L2 normalisation from maxsim. How do scores change?
# 3. Challenge: implement re-ranker: retrieve top-10 with MaxSim,
#    re-rank by mean-pooled cosine similarity. Does Recall@3 improve?
```

- [ ] Step 7: Create `03_Diffusion_Generative_Models.ipynb`

Notebook cells (nbformat=4, GPU optional):

**Header cell (markdown):** Module 03 Diffusion & Generative Models. Attribution: adapted from DDRM (Kawar et al. 2022), bahjat-kawar/ddrm. GPU optional — all demos run on CPU. Learning Objectives: (1) score function + score matching; (2) DDIM forward/reverse; (3) linear vs cosine schedules; (4) inverse problems formulation; (5) SVD degradation operators; (6) score field visualization.

**Imports cell (code):**
```python
import warnings; warnings.filterwarnings('ignore')
import numpy as np, matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
rng = np.random.default_rng(42)
print('All demos run on CPU.')
```

**Theory cell: Score-based generative models (markdown):** Score = ∇ₓ log p(x): gradient of log-density, points toward higher probability regions. Score matching: train εθ(xt,t) to predict noise ε added to x0 → equivalent to estimating score. Table: DDPM (stochastic, ~1000 steps) vs DDIM (deterministic ODE, ~50 steps, same model). Key insight: same trained denoiser works for both.

**Noise schedules cell (code):**
```python
def linear_schedule(T, b0=1e-4, bT=0.02):
    betas = np.linspace(b0, bT, T)
    ab = np.cumprod(1-betas)
    return {"betas": betas, "alpha_bars": ab, "T": T}

def cosine_schedule(T, s=0.008):
    t  = np.linspace(0, T, T+1)
    f  = np.cos(((t/T + s)/(1+s)) * np.pi/2)**2
    ab = (f/f[0])[1:]
    betas = np.clip(1 - ab/np.concatenate([[1.0], ab[:-1]]), 0, 0.999)
    return {"betas": betas, "alpha_bars": ab, "T": T}

T = 1000
lin, cos = linear_schedule(T), cosine_schedule(T)

fig, axes = plt.subplots(1, 3, figsize=(14,4))
ts = np.arange(T)
axes[0].plot(ts, lin['betas'], label='Linear'); axes[0].plot(ts, cos['betas'], label='Cosine')
axes[0].set_title('β schedule'); axes[0].legend()
axes[1].plot(ts, lin['alpha_bars']); axes[1].plot(ts, cos['alpha_bars'])
axes[1].set_title('ᾱ (cumul. product)')
snr = lambda sched: 10*np.log10(sched['alpha_bars']/(1-sched['alpha_bars']+1e-8))
axes[2].plot(ts, snr(lin)); axes[2].plot(ts, snr(cos))
axes[2].set_title('SNR (dB)')
plt.tight_layout(); plt.show()
```

**Forward process cell (code):**
```python
SIZE = 16
cx, cy = SIZE//2, SIZE//2
ii, jj = np.mgrid[:SIZE, :SIZE]
x0 = np.cos(np.sqrt((ii-cx)**2+(jj-cy)**2)*0.8).ravel()

def q_sample(x0, t, sched):
    ab = sched['alpha_bars'][t]
    eps = rng.standard_normal(x0.shape)
    return np.sqrt(ab)*x0 + np.sqrt(1-ab)*eps, eps

show_t = [0, 100, 250, 500, 750, 999]
fig, axes = plt.subplots(2, len(show_t), figsize=(14,5))
for col, t in enumerate(show_t):
    for row, sched in enumerate([lin, cos]):
        xt, _ = q_sample(x0, t, sched)
        axes[row,col].imshow(xt.reshape(SIZE,SIZE), cmap='gray', vmin=-2, vmax=2)
        axes[row,col].set_title(f't={t}', fontsize=9); axes[row,col].axis('off')
axes[0,0].set_ylabel('Linear', fontsize=9); axes[1,0].set_ylabel('Cosine', fontsize=9)
plt.suptitle('Forward diffusion noise progression'); plt.tight_layout(); plt.show()
```

**DDIM reverse cell (code):**
```python
def ddim_step(xt, eps_pred, t, t_prev, sched):
    ab_t    = sched['alpha_bars'][t]
    ab_prev = sched['alpha_bars'][t_prev] if t_prev >= 0 else 1.0
    x0_pred = np.clip((xt - np.sqrt(1-ab_t)*eps_pred)/np.sqrt(ab_t), -1, 1)
    return np.sqrt(ab_prev)*x0_pred + np.sqrt(1-ab_prev)*eps_pred

xT, _ = q_sample(x0, T-1, lin)
xt = xT.copy()
ts = list(range(T-1, -1, -T//50))
for i in range(len(ts)-1):
    _, eps = q_sample(x0, ts[i], lin)
    xt = ddim_step(xt, eps, ts[i], ts[i+1], lin)

fig, axes = plt.subplots(1,3, figsize=(9,3))
for ax, d, title in zip(axes, [x0, xT, xt], ['Original','Noisy xT','Reconstructed']):
    ax.imshow(d.reshape(SIZE,SIZE), cmap='gray', vmin=-1, vmax=1); ax.set_title(title); ax.axis('off')
plt.tight_layout(); plt.show()
```

**Theory cell: Inverse problems (markdown):** y = Hx + n. Table: H=I denoising, mask inpainting, downsample super-res, PSF conv deblurring. DDRM combines diffusion prior with data-consistency projection at each reverse step via SVD of H. Scientific apps: cryo-EM, MRI, fluorescence microscopy.

**SVD degradation cell (code):**
```python
N = SIZE**2
H = np.eye(N)[:N//4]   # observe 25% of pixels
U, s, Vt = np.linalg.svd(H, full_matrices=False)
y = H @ x0 + rng.normal(0, 0.1, N//4)
x_pinv = Vt.T @ np.diag(1/(s+1e-8)) @ U.T @ y

fig, axes = plt.subplots(1,3, figsize=(9,3))
for ax, d, title in zip(axes, [x0, rng.normal(0,0.4,N)+x0, x_pinv],
                         ['Original','Noisy','Pseudo-inverse (25%)']):
    ax.imshow(d.reshape(SIZE,SIZE), cmap='gray', vmin=-1, vmax=1); ax.set_title(title); ax.axis('off')
plt.tight_layout(); plt.show()
```

**Score field visualization cell (code):**
```python
def bimodal_score(xy, sigma=0.8):
    mu1, mu2 = np.array([-2.,0.]), np.array([2.,0.])
    s1 = -(xy-mu1)/sigma**2; s2 = -(xy-mu2)/sigma**2
    w1 = np.exp(-0.5*((xy-mu1)**2/sigma**2).sum(1,keepdims=True))
    w2 = np.exp(-0.5*((xy-mu2)**2/sigma**2).sum(1,keepdims=True))
    return (w1*s1 + w2*s2)/(w1+w2+1e-8)

g = 18; xg, yg = np.linspace(-4,4,g), np.linspace(-3,3,g)
X, Y = np.meshgrid(xg, yg)
pts  = np.stack([X.ravel(), Y.ravel()], axis=1)
sc   = bimodal_score(pts)

fig, ax = plt.subplots(figsize=(8,5))
ax.quiver(X, Y, sc[:,0].reshape(g,g), sc[:,1].reshape(g,g), alpha=0.8, color='#4393c3')
for mu in [[-2,0],[2,0]]:
    ax.scatter(*mu, s=120, c='red', zorder=5)
ax.set_title('Score field ∇ₓ log p(x) — bimodal Gaussian')
ax.set_xlabel('x₁'); ax.set_ylabel('x₂')
plt.tight_layout(); plt.show()
```

**Exercise cell (code, commented):**
```python
# Exercise T5-03
# 1. At what t does ᾱₜ drop below 0.01 for linear vs cosine schedule?
# 2. Run DDIM with cosine schedule + oracle noise. Compare reconstruction quality.
# 3. Challenge: implement DDRM data-consistency for denoising (H=I):
#    x0_proj = x0_pred + (sigma_obs² / (sigma_obs² + (1-ab)/ab)) * (y_obs - x0_pred)
#    Use x0_proj in ddim_step instead of x0_pred. Compare results.
```

- [ ] Step 8: Commit

```bash
git add Course/Tier_5_Modern_AI_for_Science/ \
        Skills/llm-finetuning.md \
        Skills/vision-rag.md \
        Skills/diffusion-generative.md
git commit -m "feat(tier5): add LLM fine-tuning (01), Vision RAG (02), Diffusion (03) and skills"
```

---

## PHASE 2 — Sequential (Task 7, runs after Tasks 1–6)

### Task 7: README Updates (runs after Tasks 1–6 complete)

**Files:**
- Create: `Course/Tier_5_Modern_AI_for_Science/README.md`
- Update: `README.md` (root)
- Update: `Course/README.md`
- Update: `Skills/README.md`

**Prerequisite:** All notebooks, module READMEs, and skill files from Tasks 1–6 must exist. Verify with:
```bash
ls Course/Tier_2_Core_Bioinformatics/14_Hi-C_Analysis/
ls Course/Tier_2_Core_Bioinformatics/15_Motif_Discovery/
ls Course/Tier_3_Applied_Bioinformatics/19_GWAS/
ls Course/Tier_3_Applied_Bioinformatics/20_Spatial_Transcriptomics/
ls Course/Tier_3_Applied_Bioinformatics/21_Copy_Number_Analysis/
ls Course/Tier_3_Applied_Bioinformatics/22_Bayesian_Statistics_Python/
ls Course/Tier_3_Applied_Bioinformatics/23_TF_Footprinting/
ls Course/Tier_5_Modern_AI_for_Science/01_LLM_Finetuning/
ls Course/Tier_5_Modern_AI_for_Science/02_Vision_RAG/
ls Course/Tier_5_Modern_AI_for_Science/03_Diffusion_Generative_Models/
```

**Steps:**

- [ ] Step 1: Create `Course/Tier_5_Modern_AI_for_Science/README.md`

```markdown
# Tier 5: Modern AI for Science

GPU-optional modules covering contemporary AI methods applied to scientific research. Every notebook is designed to run on free-tier Google Colab. Theory and code-pattern cells run without a GPU; hands-on training cells require a T4 or A100 runtime.

## Prerequisites

- Tier 3, Module 10 (Deep Learning for Biology) — transformer architectures and PyTorch basics
- Comfort with numpy, pandas, and matplotlib

## Modules

| Module | Topic | GPU Required |
|--------|-------|-------------|
| [01 LLM Fine-tuning](01_LLM_Finetuning/01_LLM_Finetuning.ipynb) | LoRA, quantization, SFTTrainer, instruction datasets | Yes (T4 or better) |
| [02 Vision RAG](02_Vision_RAG/02_Vision_RAG.ipynb) | VLMs, ColPali, document retrieval, RAG pipeline | Optional (CPU feasible for inference) |
| [03 Diffusion & Generative Models](03_Diffusion_Generative_Models/03_Diffusion_Generative_Models.ipynb) | Score matching, DDIM, inverse problems, noise schedules | Optional (CPU feasible for small examples) |

## Running on Colab

Each notebook contains a setup cell that installs all required packages. Recommended runtime: **T4 GPU** (free tier). Switch to CPU runtime for theory-only reading.

```python
# Standard Colab setup (at top of each notebook)
!pip install -q unsloth trl peft bitsandbytes transformers accelerate
```

## Attribution

| Module | Source | Attribution |
|--------|--------|-------------|
| LLM Fine-tuning | FinetuningSmallLanguageModels workshop patterns | Inspired by Unsloth AI & Manuel Faysse |
| Vision RAG | VisionRag workshop patterns | Inspired by Unsloth AI & Manuel Faysse |
| Diffusion & Generative Models | Lundi workshop patterns | Adapted from DDRM (Kawar et al., 2022), github.com/bahjat-kawar/ddrm |

All content uses public datasets (HuggingFace Hub, torchvision). No proprietary research data.
```

- [ ] Step 2: Update `README.md` (root) — read it first with the Read tool, then apply these exact edits

**Edit 2a — Stats line:** find and replace:
```
**96 notebooks** · **5 tiers** · **30 interactive visualizations** · **108 glossary terms** · **[28 Claude Code skills](#claude-code-skills)**
```
→
```
**106 notebooks** · **6 tiers** · **30 interactive visualizations** · **108 glossary terms** · **[36 Claude Code skills](#claude-code-skills)**
```

**Edit 2b — Tier 2 line in Structure block:** find and replace:
```
Tier 2  Core Bioinformatics            15 notebooks
        Databases · BioPython · Alignment · BLAST · MSA ·
        Phylogenetics · Protein Structure · Nucleic Acids ·
        Chromatograms · Motifs · GO/Pathways · Comparative Genomics ·
        Computational Genetics
```
→
```
Tier 2  Core Bioinformatics            17 notebooks
        Databases · BioPython · Alignment · BLAST · MSA ·
        Phylogenetics · Protein Structure · Nucleic Acids ·
        Chromatograms · Motifs · GO/Pathways · Comparative Genomics ·
        Computational Genetics · Hi-C Analysis · Motif Discovery
```

**Edit 2c — Tier 3 lines in Structure block:** find and replace:
```
Tier 3  Applied Bioinformatics         21 notebooks
        NGS · Variant Calling · RNA-seq · Microbial Diversity ·
        Promoters · Statistics · Machine Learning · Deep Learning ·
        Molecular Modeling · Clinical Genomics · Capstone Project ·
        Biochemistry & Enzyme Kinetics · Genetic Engineering ·
        Population Genetics · Numerical Methods ·
        Genome Assembly · Proteomics & Structural Methods
```
→
```
Tier 3  Applied Bioinformatics         26 notebooks
        NGS · Variant Calling · RNA-seq · Microbial Diversity ·
        Promoters · Statistics · Machine Learning · Deep Learning ·
        Molecular Modeling · Clinical Genomics · Capstone Project ·
        Biochemistry & Enzyme Kinetics · Genetic Engineering ·
        Population Genetics · Numerical Methods ·
        Genome Assembly · Proteomics & Structural Methods ·
        GWAS · Spatial Transcriptomics · Copy Number Analysis ·
        Bayesian Statistics · TF Footprinting
```

**Edit 2d — Add Tier 5 to Structure block.** Find the closing triple-backtick that ends the structure block (after the Tier 4 lines). The Tier 4 block ends with:
```
        KMP · Rabin-Karp · Tries · Suffix Trees · Graphs · DP
```
Replace with:
```
        KMP · Rabin-Karp · Tries · Suffix Trees · Graphs · DP

Tier 5  Modern AI for Science          3 notebooks
        LLM Fine-tuning · Vision RAG · Diffusion & Generative Models
```

**Edit 2e — Update the "See the full table" line:** find and replace:
```
See the full table of contents in [Course/README.md](Course/README.md) and [Tier 4 README](Course/Tier_4_Algorithms_and_Data_Structures/README.md).
```
→
```
See the full table of contents in [Course/README.md](Course/README.md), [Tier 4 README](Course/Tier_4_Algorithms_and_Data_Structures/README.md), and [Tier 5 README](Course/Tier_5_Modern_AI_for_Science/README.md).
```

**Edit 2f — Update skills intro sentence:** find and replace:
```
The entire course is compressed into **28 modular skill files** for [Claude Code](https://claude.com/claude-code) — maximum knowledge density, minimum tokens. Each skill provides quick-reference tables, copy-paste code templates, and common pitfalls for a focused topic.
```
→
```
The entire course is compressed into **36 modular skill files** for [Claude Code](https://claude.com/claude-code) — maximum knowledge density, minimum tokens. Each skill provides quick-reference tables, copy-paste code templates, and common pitfalls for a focused topic.
```

**Edit 2g — Add Modern AI row to Claude Code Skills table.** The table currently ends with:
```
| **Biology & Computation** | `probability-statistics-python` · `genetics-computational` · `biochemistry-enzymology` · `genetic-engineering-insilico` · `population-genetics-evolution` · `numerical-methods-bio` |
```
Replace with:
```
| **Biology & Computation** | `probability-statistics-python` · `genetics-computational` · `biochemistry-enzymology` · `genetic-engineering-insilico` · `population-genetics-evolution` · `numerical-methods-bio` |
| **Tier 2 Depth** | `hic-analysis` · `motif-discovery` |
| **Applied Bio Depth** | `gwas-population-genetics` · `spatial-transcriptomics` · `bayesian-python` |
| **Modern AI** | `llm-finetuning` · `vision-rag` · `diffusion-generative` |
```

- [ ] Step 3: Update `Course/README.md` — read it first with the Read tool (it is 1364 lines; read in two passes: lines 1–100 for header, then the tier sections you need), then apply these edits

**Edit 3a — Stats line at top:** find and replace:
```
`96 notebooks` | `5 tiers` | `108 glossary terms` | `12 sample data files` | `30 interactive visualizations`
```
→
```
`106 notebooks` | `6 tiers` | `108 glossary terms` | `12 sample data files` | `30 interactive visualizations`
```

**Edit 3b — Visual map Tier 2 count:** find and replace:
```
│  TIER 2: CORE BIOINFORMATICS                                   15 notebooks  │
```
→
```
│  TIER 2: CORE BIOINFORMATICS                                   17 notebooks  │
```

**Edit 3c — Visual map Tier 3 count:** find and replace:
```
│  TIER 3: APPLIED BIOINFORMATICS                                21 notebooks  │
```
→
```
│  TIER 3: APPLIED BIOINFORMATICS                                26 notebooks  │
```

**Edit 3d — Add Tier 5 box to visual map.** The map ends with Tier 4's closing box:
```
└──────────────────────────────────────────────────────────────────────────────┘

  ENTRY POINTS ──────────────────────────────────────────────────────────────
```
Replace with:
```
└──────────────────────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  TIER 5: MODERN AI FOR SCIENCE                                  3 notebooks  │
│  ──────────────────────────────────────────────────────────────────────────  │
│  LLM Fine-tuning │ Vision RAG │ Diffusion & Generative Models                │
│                                                                              │
│  Entry: Tiers 1–3 complete. GPU-optional; runs on free-tier Colab.          │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘

  ENTRY POINTS ──────────────────────────────────────────────────────────────
```

**Edit 3e — "Who This Course Is For" table:** find and replace:
```
| **Complete beginner** -- no programming experience, new to biology or coming from a wet-lab background | Tier 0, Module 01 (Linux Fundamentals) | All five tiers (0 through 4) | 150--200 hours |
```
→
```
| **Complete beginner** -- no programming experience, new to biology or coming from a wet-lab background | Tier 0, Module 01 (Linux Fundamentals) | All six tiers (0 through 5) | 160--215 hours |
```

**Edit 3f — Tier 2 section header:** find and replace:
```
### Tier 2: Core Bioinformatics -- 15 notebooks
```
→
```
### Tier 2: Core Bioinformatics -- 17 notebooks
```

**Edit 3g — Add modules 2.14 and 2.15 after module 2.13.** The 2.13 section ends with:
```
`genetic code` `codon usage` `CAI` `restriction enzymes` `ORF finding` `genetic mapping` `Ts/Tv` `Hardy-Weinberg`

---

### Tier 3: Applied Bioinformatics -- 21 notebooks
```
Replace with:
```
`genetic code` `codon usage` `CAI` `restriction enzymes` `ORF finding` `genetic mapping` `Ts/Tv` `Hardy-Weinberg`

---

#### 2.14 Hi-C Analysis

[14_hic_analysis.ipynb](Tier_2_Core_Bioinformatics/14_Hi-C_Analysis/14_hic_analysis.ipynb)

3D genome organization from Hi-C experiments using the cooler and cooltools Python stack. cooler file format: loading, inspecting metadata, and slicing contact matrices. Contact decay (expected) curves to normalize distance effects. Eigenvector decomposition for A/B compartment identification. Insulation score calculation and TAD boundary detection. Saddle plots for compartment strength visualization. Pileup (aggregate) analysis around genomic features such as CTCF sites and loop anchors. Uses public Hi-C data from the 4DN Data Portal or ENCODE.

`cooler` `cooltools` `contact matrices` `A/B compartments` `TADs` `insulation score` `saddle plots` `pileup`

---

#### 2.15 Motif Discovery

[15_motif_discovery.ipynb](Tier_2_Core_Bioinformatics/15_Motif_Discovery/15_motif_discovery.ipynb)

Quantitative motif analysis from position frequency matrices to enrichment testing. PPM/PWM construction, normalization, and pseudocount handling. Information content per position and total IC calculation. KDIC score (mean IC normalized to [0,1]). IUPAC consensus sequence generation. Score distributions: exact enumeration for short motifs (length ≤ 10), Monte Carlo with confidence intervals for longer. Motif enrichment using Fisher's exact test with Benjamini-Hochberg correction. TomTom matching concept against JASPAR 2024 and HOCOMOCO public databases. Pipeline design patterns using abstract interfaces and dataclasses; BED/FASTA I/O patterns.

`PWM` `PPM` `information content` `KDIC` `IUPAC` `Fisher enrichment` `Benjamini-Hochberg` `TomTom` `JASPAR`

---

### Tier 3: Applied Bioinformatics -- 26 notebooks
```

**Edit 3h — Update Tier 3 section description paragraph:** find and replace:
```
Advanced topics and real-world analysis pipelines. Each notebook covers a complete workflow from raw data to biological conclusions. Includes a capstone project integrating skills from every tier, plus specialized modules on molecular modeling, deep learning, clinical genomics, and modern bioinformatics workflows (single-cell analysis, pipeline engines, testing/CI-CD).
```
→
```
Advanced topics and real-world analysis pipelines. Each notebook covers a complete workflow from raw data to biological conclusions. Includes a capstone project integrating skills from every tier, plus specialized modules on molecular modeling, deep learning, clinical genomics, modern bioinformatics workflows, GWAS, spatial transcriptomics, copy number analysis, Bayesian statistics, and TF footprinting.
```

**Edit 3i — Add modules 3.19–3.23 after module 3.18.** The 3.18 section ends with:
```
`mass spectrometry` `proteomics` `MS/MS` `peptide identification` `FDR` `TMT` `protein engineering` `X-ray crystallography` `cryo-EM`

---

### Tier 4: Algorithms & Data Structures -- 30 notebooks, 927 cells
```
Replace with:
```
`mass spectrometry` `proteomics` `MS/MS` `peptide identification` `FDR` `TMT` `protein engineering` `X-ray crystallography` `cryo-EM`

---

#### 3.19 Genome-Wide Association Studies (GWAS)

[19_gwas.ipynb](Tier_3_Applied_Bioinformatics/19_GWAS/19_gwas.ipynb)

GWAS study design from first principles. Case/control phenotype definition and confounder identification. Quality control: SNP and sample filtering, Hardy-Weinberg equilibrium testing, MAF thresholds. Population stratification detection via PCA on genotype data. Association testing using logistic and linear regression per SNP. Multiple testing correction with the genome-wide significance threshold (5×10⁻⁸). Manhattan and QQ plot generation from scratch. Linkage disequilibrium and clumping concepts. Downstream: fine-mapping principles and GWAS catalog lookup. Uses 1000 Genomes public data or simulated genotype matrices.

`GWAS` `population stratification` `PCA` `Manhattan plot` `QQ plot` `LD` `clumping` `fine-mapping` `5e-8`

---

#### 3.20 Spatial Transcriptomics

[20_spatial_transcriptomics.ipynb](Tier_3_Applied_Bioinformatics/20_Spatial_Transcriptomics/20_spatial_transcriptomics.ipynb)

Spatial gene expression analysis with Squidpy and Scanpy. AnnData structure with spatial coordinates; Visium and Xenium layout conventions. Quality control for spatial data: mitochondrial fraction and spot-level filtering. Normalization and dimensionality reduction in spatial context. Spatial neighborhood graph construction. Spatially variable gene detection. Cell-type deconvolution concepts (RCTD, cell2location patterns). Visualization: spatial scatter plots and expression overlays on tissue sections. Uses the public 10x Visium mouse brain dataset available through Squidpy.

`spatial transcriptomics` `Squidpy` `AnnData` `Visium` `spatially variable genes` `deconvolution` `neighborhood graph` `scanpy`

---

#### 3.21 DNA Copy Number Analysis

[21_copy_number_analysis.ipynb](Tier_3_Applied_Bioinformatics/21_Copy_Number_Analysis/21_copy_number_analysis.ipynb)

Copy number variation analysis from sequencing data. CNV concepts: gains, losses, and loss of heterozygosity. Read depth normalization approaches across genomic windows. Segmentation using the Circular Binary Segmentation (CBS) algorithm concept. Copy number state calling from segments. Genome-wide CN profile visualization. Gene-level annotation of CN events. Extends the variant calling pipeline concepts from Module 3.02.

`copy number` `CNV` `CBS segmentation` `read depth normalization` `LOH` `genome-wide profile` `somatic variants`

---

#### 3.22 Bayesian Statistics in Python

[22_bayesian_statistics_python.ipynb](Tier_3_Applied_Bioinformatics/22_Bayesian_Statistics_Python/22_bayesian_statistics_python.ipynb)

Bayesian statistical modeling in Python across seven sub-sections. Frequentist vs. Bayesian framing: posterior = likelihood × prior, credible intervals. Prior specification: informative vs. weakly informative priors, prior predictive checks. Multiple regression: collinearity diagnostics and variance inflation factor. Model comparison with WAIC and LOO-CV using ArviZ. Linear mixed-effects models with random intercepts and slopes using Bambi. GLMs in a Bayesian framework: Bernoulli, Binomial, Poisson, Negative Binomial with PyMC. Advanced: GLMM, zero-inflated models, GAM concepts, and Bayesian meta-analysis. Converted from Fränzi Korner-Nievergelt's applied statistics R course; uses the public palmerpenguins dataset.

`pymc` `bambi` `arviz` `LOO-CV` `WAIC` `credible intervals` `GLM` `mixed-effects` `prior predictive` `zero-inflated`

---

#### 3.23 TF Footprinting & Chromatin Accessibility

[23_tf_footprinting.ipynb](Tier_3_Applied_Bioinformatics/23_TF_Footprinting/23_tf_footprinting.ipynb)

Transcription factor footprinting from ATAC-seq data. ATAC-seq recap: fragment size distribution and nucleosome-free region identification. TF footprinting concept: Tn5 insertion bias around motif binding sites. Expected vs. observed cut-site profiles around motif centers. Footprint score calculation and interpretation across conditions. Genomic interval arithmetic with pybedtools: intersection, subtraction, closest-feature queries. Accumulation plots: average signal enrichment around genomic features. Extends the ngs-variant-calling and motif-discovery skills; uses public ENCODE ATAC-seq data.

`ATAC-seq` `TF footprinting` `Tn5 bias` `cut-site profiles` `pybedtools` `genomic intervals` `accumulation plots` `chromatin accessibility`

---

### Tier 4: Algorithms & Data Structures -- 30 notebooks, 927 cells
```

**Edit 3j — Add new Tier 5 section after Tier 4 (before `## Skills Check Guide`).** Find:
```
## Skills Check Guide
```
Replace with:
```
### Tier 5: Modern AI for Science -- 3 notebooks

GPU-optional modules covering contemporary AI methods for scientific research. Each notebook is designed to run on free-tier Google Colab. Theory cells run without GPU; hands-on training cells require a T4 or better. See the [Tier 5 README](Tier_5_Modern_AI_for_Science/README.md) for GPU setup and Colab instructions.

---

#### 5.01 LLM Fine-tuning

[01_LLM_Finetuning.ipynb](Tier_5_Modern_AI_for_Science/01_LLM_Finetuning/01_LLM_Finetuning.ipynb)

Fine-tuning large language models for domain-specific instruction following. Base vs. instruction/chat models: what changes during fine-tuning and why. LoRA: low-rank adapter mathematics, rank selection, and target module identification. Quantization: 4-bit NF4 with bitsandbytes and trade-offs with output quality. Chat template formatting: system/user/assistant structure. SFTTrainer workflow: dataset preparation, training loop configuration, and evaluation. Synthetic data generation for instruction tuning. Practical tips: gradient checkpointing, batch size scheduling, learning rate warm-up. Inspired by Unsloth AI and Manuel Faysse fine-tuning patterns; uses public instruction datasets from HuggingFace Hub.

`LoRA` `quantization` `NF4` `SFTTrainer` `trl` `peft` `bitsandbytes` `chat templates` `instruction tuning` `synthetic data`

---

#### 5.02 Vision RAG

[02_Vision_RAG.ipynb](Tier_5_Modern_AI_for_Science/02_Vision_RAG/02_Vision_RAG.ipynb)

Retrieval-augmented generation with vision-language models for document understanding. VLM architecture overview: visual encoder + LLM decoder. Document understanding: page-level vs. token-level retrieval approaches. ColPali: late-interaction document retrieval concept and API pattern. RAG pipeline: retrieval → context injection → generation. Qwen2-VL inference pattern for multi-page document Q&A. Evaluation: retrieval recall and generation faithfulness metrics. Uses public PDF documents (arXiv papers, open-access reports). Inspired by Unsloth AI and Manuel Faysse patterns.

`VLM` `ColPali` `late interaction` `RAG` `document retrieval` `Qwen2-VL` `retrieval recall` `generation faithfulness`

---

#### 5.03 Diffusion & Generative Models

[03_Diffusion_Generative_Models.ipynb](Tier_5_Modern_AI_for_Science/03_Diffusion_Generative_Models/03_Diffusion_Generative_Models.ipynb)

Score-based generative models and diffusion for scientific imaging applications. Score matching intuition: score field as ∇ₓ log p(x), learned denoising. DDIM: deterministic sampling, linear and cosine noise schedules, reverse process equations. Inverse problems in imaging: denoising, inpainting, and colorization as special cases of DDRM. SVD-based degradation operators and the pseudoinverse projection step. Linear and cosine scheduler implementation patterns with explicit tensor shapes. Score field visualization with quiver plots. Scientific applications: cryo-EM denoising and medical image restoration concepts. Adapted from DDRM (Kawar et al., 2022), github.com/bahjat-kawar/ddrm; uses public MNIST/CIFAR-10 data.

`diffusion models` `DDIM` `score matching` `noise schedule` `inverse problems` `DDRM` `SVD` `cryo-EM` `image restoration`

---

## Skills Check Guide
```

- [ ] Step 4: Update `Skills/README.md` — read it first, then apply these edits

**Edit 4a — Intro line:** find and replace:
```
28 skill files compressing a 96-notebook, 5-tier bioinformatics course into actionable Claude Code references. Each skill provides key patterns, code templates, complexity tables, and common pitfalls for a focused topic area.
```
→
```
36 skill files compressing a 106-notebook, 6-tier bioinformatics course into actionable Claude Code references. Each skill provides key patterns, code templates, complexity tables, and common pitfalls for a focused topic area.
```

**Edit 4b — Course Origin section:** find and replace:
```
Based on the **Bioinformatics with Python** course — 96 notebooks across 5 tiers, built from materials by:
```
→
```
Based on the **Bioinformatics with Python** course — 106 notebooks across 6 tiers, built from materials by:
```

**Edit 4c — Add three new skill sections after the last existing section.** The `### Biology & Computation (Skills 22–28)` section ends with the `numerical-methods-bio` row, then `---`. Find:
```
| [`numerical-methods-bio`](numerical-methods-bio.md) | Interpolation, curve fitting, optimization, FFT for bioinformatics data |

---

## Skill File Format
```
Replace with:
```
| [`numerical-methods-bio`](numerical-methods-bio.md) | Interpolation, curve fitting, optimization, FFT for bioinformatics data |

### Tier 2 Depth (Skills 29–30)

| Skill | Use When... |
|-------|-------------|
| [`hic-analysis`](hic-analysis.md) | Analyzing 3D genome organization: loading cooler/mcool files, computing A/B compartments, TAD boundaries, pileup plots with cooltools |
| [`motif-discovery`](motif-discovery.md) | Building PWMs from aligned binding sites, IC/KDIC scoring, Fisher enrichment testing, TomTom database matching |

### Applied Bioinformatics Depth (Skills 31–33)

| Skill | Use When... |
|-------|-------------|
| [`gwas-population-genetics`](gwas-population-genetics.md) | GWAS study design, QC pipelines, PCA for stratification, Manhattan/QQ plots, LD clumping, fine-mapping |
| [`spatial-transcriptomics`](spatial-transcriptomics.md) | AnnData with spatial coordinates, spatially variable gene detection, Squidpy neighborhood graphs, cell-type deconvolution |
| [`bayesian-python`](bayesian-python.md) | Bayesian linear and generalized linear models with pymc/bambi, model comparison via arviz LOO-CV/WAIC, prior specification |

### Modern AI for Science (Skills 34–36)

| Skill | Use When... |
|-------|-------------|
| [`llm-finetuning`](llm-finetuning.md) | Fine-tuning LLMs with LoRA and quantization, configuring SFTTrainer, designing chat templates, building instruction datasets |
| [`vision-rag`](vision-rag.md) | Vision-language models, ColPali document retrieval, RAG pipeline construction, Qwen2-VL inference patterns |
| [`diffusion-generative`](diffusion-generative.md) | DDIM sampling, linear/cosine noise schedules, score matching intuition, inverse imaging problems, scientific image restoration |

---

## Skill File Format
```

- [ ] Step 5: Commit

```bash
cd /Users/pavel/Documents/Python-scripts-for-everything
git add README.md \
        Course/README.md \
        Skills/README.md \
        Course/Tier_5_Modern_AI_for_Science/README.md
git commit -m "docs: update READMEs for 10 new modules, 8 skills, and Tier 5 (task 7)"
```

---
