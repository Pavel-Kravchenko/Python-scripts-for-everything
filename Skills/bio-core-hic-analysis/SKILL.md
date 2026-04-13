---
name: bio-core-hic-analysis
description: "*Prerequisites: Tier 2 Modules 1–13. Familiarity with numpy arrays and matplotlib.*"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/14_Hi-C_Analysis/14_hic_analysis.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Hi-C Analysis: 3D Genome Organization

*Source: Course notebook `Tier_2_Core_Bioinformatics/14_Hi-C_Analysis/14_hic_analysis.ipynb`*

# Hi-C Analysis: 3D Genome Organization

**Tier 2 — Core Bioinformatics | Module 14**

*Prerequisites: Tier 2 Modules 1–13. Familiarity with numpy arrays and matplotlib.*

---

This module covers 3D genome organization analysis using Hi-C sequencing data. You will learn to:

1. Load and inspect Hi-C contact matrices in cooler format
2. Visualize contact matrices and identify TADs visually
3. Compute the P(s) contact decay curve for quality assessment
4. Detect A/B compartments by eigenvector decomposition
5. Calculate insulation scores for TAD boundary calling
6. Build aggregate pileup plots for weak signal detection

**Tools covered:** `cooler`, `cooltools`, `numpy`, `pandas`, `matplotlib`

## Why this notebook matters

Hi-C (chromosome conformation capture followed by high-throughput sequencing) reveals the 3D organization of the genome inside the nucleus. The genome is not a linear string of DNA — it folds into a hierarchy of structures (compartments, TADs, loops) that directly control gene expression. Enhancers can act over megabases of linear distance by looping to their target promoters. TAD disruption is a mechanism of oncogene activation in cancer. Understanding how to load, visualize, normalize, and analyze Hi-C contact matrices is increasingly essential for understanding regulatory genomics and chromatin biology.

## Complicated moments explained

- **Why cooler format?**: Raw Hi-C produces read pairs that are mapped to the genome and converted into a contact matrix (each cell = number of ligations between two genomic bins). The cooler format (.cool, .mcool) stores these sparse contact matrices efficiently with metadata. Use `mcool` files for multi-resolution analysis — you can access specific resolutions with the `::resolutions/25000` suffix.
- **Balancing (ICE normalization)**: Raw contact counts are biased by GC content, mappability, restriction fragment density, and bin size. ICE (Iterative Correction and Eigenvector decomposition) balancing removes these biases by assuming that all genomic loci should have equal total contacts. The balanced matrix is stored in the `weight` column of cooler files. Use `balance=True` for compartment and insulation analysis; `balance=False` for visualization.
- **Diagonal artifacts**: The first few diagonals of a contact matrix (very short genomic distances) are dominated by unligated fragments and self-ligation products — they are not true contacts. Always `ignore_diags=2` (or more) in cooltools analyses to exclude these.
- **A/B compartment sign flip**: The first eigenvector (E1) from PCA of the O/E matrix captures the A/B pattern, but its sign is arbitrary. Flip so that A compartment (gene-dense, active) has positive values. Use GC content or gene density as a reference to determine the correct orientation.
- **TAD calling depends on resolution**: TADs appear at 25-40 kb resolution. At 5-10 kb you see sub-TADs; at 100 kb you see compartment-scale domains. The insulation score window size controls the scale of boundaries detected. Always report the resolution and window size used.

## Background: 3D Genome Organization

DNA is not a linear string — it folds into a hierarchy of three-dimensional structures inside the nucleus:

| Level | Scale | Feature | Detection |
|---|---|---|---|
| Compartments | ~1–10 Mb | A (active) / B (inactive) chromatin | Eigenvector decomposition |
| TADs | ~100 kb–3 Mb | Topologically Associating Domains | Insulation score |
| Loops | ~10–300 kb | CTCF-anchored contacts | Pileup analysis |

**Hi-C protocol:** Cells are cross-linked, DNA is digested with a restriction enzyme, proximity-ligated, and sequenced. Read pairs reflect spatial proximity: reads mapping far apart in sequence but close in 3D space show up as off-diagonal contacts.

**Contact matrix:** After mapping, read pairs are binned into a genome-wide matrix where entry (i, j) counts contacts between bins i and j. Balancing (iterative correction / KR normalization) removes biases from mappability and GC content.

**cooler format:** A compact HDF5-based container for contact matrices at one or more resolutions. The Python `cooler` library provides random access to any genomic region.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import warnings
warnings.filterwarnings("ignore")

# Hi-C specific libraries
try:
    import cooler
    import cooltools
    print(f"cooler {cooler.__version__}, cooltools {cooltools.__version__}")
except ImportError:
    print("Install with: pip install cooler cooltools")
    print("On Colab: !pip install cooler cooltools")

# For download helper
import urllib.request, os

plt.rcParams.update({"figure.dpi": 120, "axes.spines.top": False, "axes.spines.right": False})
```

## 1. Loading Hi-C Data with cooler

We use a small public Hi-C dataset from the 4DN Data Portal. The cell below downloads a `.cool` file at 25 kb resolution (chromosome 1 only, ~10 MB) for demonstration.

**Note:** In a real analysis you would use `.mcool` files containing multiple resolutions. Access a specific resolution with `cooler.Cooler("file.mcool::resolutions/25000")`.

```python
# Download a small public Hi-C cool file (chr1 only, 25 kb resolution)
# Source: 4DN Data Portal — GM12878 in situ Hi-C (Rao et al. 2014)
# This is a publicly available subset for educational use

DATA_URL = "https://dl.4dnucleome.org/files-processed/4DNFIYECESRC/"  # Ready
COOL_FILE = "demo_hic_25kb.cool"

# Simulate a small contact matrix if real download is unavailable
def make_demo_cooler(path, n_bins=200, binsize=25_000):
    """Create a minimal demo .cool file with synthetic contacts."""
    import cooler
    import scipy.sparse as sp
    import tempfile

    chroms = pd.DataFrame({"name": ["chr1"], "length": [n_bins * binsize]})
    bins = pd.DataFrame({
        "chrom": ["chr1"] * n_bins,
        "start": np.arange(n_bins) * binsize,
        "end":   np.arange(1, n_bins + 1) * binsize,
    })

    # Synthetic contacts: distance-decay pattern
    rng = np.random.default_rng(42)
    rows, cols, vals = [], [], []
    for i in range(n_bins):
        for j in range(i, min(i + 80, n_bins)):
            w = np.exp(-0.1 * (j - i)) * rng.poisson(20)
            if w > 0:
                rows.append(i); cols.append(j); vals.append(int(w))
    pixels = pd.DataFrame({"bin1_id": rows, "bin2_id": cols, "count": vals})

    cooler.create_cooler(path, bins=bins, pixels=pixels, dtypes={"count": np.int32})
    print(f"Demo cooler created: {path}")

if not os.path.exists(COOL_FILE):
    print("Creating synthetic demo cooler (real 4DN download requires authentication)...")
    make_demo_cooler(COOL_FILE)

clr = cooler.Cooler(COOL_FILE)
print(f"Resolution: {clr.binsize:,} bp")
print(f"Chromosomes: {clr.chromnames}")
print(f"Matrix shape: {clr.shape}")
print(f"\nBin table (first 5 rows):")
print(clr.bins()[:5])
```

## 2. Visualizing the Contact Matrix

A contact matrix is a symmetric 2D array. Entries near the diagonal = short-range contacts (always high); off-diagonal blocks = long-range contacts, which reveal compartments and TADs.

**Color scale:** Log-transformed counts. Strong diagonal = expected from polymer physics. Square blocks off-diagonal = TADs.

```python
region = "chr1:0-5000000"  # first 5 Mb

mat = clr.matrix(balance=False).fetch(region).astype(float)
mat[mat == 0] = np.nan

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Raw counts
im0 = axes[0].matshow(mat, cmap="YlOrRd", norm=mcolors.LogNorm(vmin=1, vmax=np.nanpercentile(mat, 99)))
plt.colorbar(im0, ax=axes[0], fraction=0.046, label="contacts")
axes[0].set_title(f"Raw contact matrix\n{region}")
axes[0].set_xlabel("Genomic bin"); axes[0].set_ylabel("Genomic bin")

# Log1p transform
im1 = axes[1].matshow(np.log1p(np.nan_to_num(mat)), cmap="RdPu")
plt.colorbar(im1, ax=axes[1], fraction=0.046, label="log1p(contacts)")
axes[1].set_title(f"Log1p transformed\n{region}")
axes[1].set_xlabel("Genomic bin")

plt.tight_layout()
plt.show()
print(f"Matrix shape: {mat.shape} | Max: {np.nanmax(mat):.0f} | NaN bins: {np.isnan(mat).sum()}")
```

## 3. Contact Decay Curve (P(s) Curve)

The *contact probability as a function of genomic distance* (P(s) curve) is a fundamental Hi-C quality metric and normalization anchor. Contacts decay roughly as a power law with distance: P(s) ~ s^α.

- **Active (A) compartment:** steeper decay, more short-range contacts
- **Inactive (B) compartment:** flatter decay

We compute *expected contacts* by averaging over all pairs of bins at the same genomic distance. This creates a distance-specific baseline for normalization.

```python
# Build view DataFrame (required by cooltools)
view_df = pd.DataFrame({
    "chrom": clr.chromnames,
    "start": [0] * len(clr.chromnames),
    "end":   list(clr.chromsizes.values),
    "name":  clr.chromnames,
})

expected = cooltools.expected_cis(clr, view_df=view_df, ignore_diags=2)
print(expected.columns.tolist())
print(expected.head())

# Plot P(s)
fig, ax = plt.subplots(figsize=(7, 4))
# Use raw count average (balanced not available for synthetic data)
count_col = "count.avg" if "count.avg" in expected.columns else expected.filter(like="avg").columns[0]
dist_bp = expected["dist"] * clr.binsize

mask = dist_bp > 0
ax.loglog(dist_bp[mask], expected[count_col][mask], "o-", color="steelblue", ms=4, lw=1.5)
ax.set_xlabel("Genomic distance (bp)")
ax.set_ylabel("Mean contact frequency")
ax.set_title("Contact decay curve (P(s))")
ax.grid(True, alpha=0.3, which="both")
plt.tight_layout()
plt.show()
```

## 4. A/B Compartment Detection

Chromatin segregates into two compartments:
- **A compartment** — transcriptionally active, gene-dense, open chromatin (positive E1)
- **B compartment** — silent, gene-poor, heterochromatic (negative E1)

The method: compute the *observed/expected* (O/E) matrix (divide each contact by the expected at that distance), then apply eigenvector decomposition (PCA). The first eigenvector (E1) captures the A/B pattern.

**Sign convention:** E1 sign is arbitrary. Flip so positive = high GC content = A compartment.

```python
try:
    eigvals, eigvecs = cooltools.eigs_cis(
        clr,
        view_df=view_df,
        n_eigs=3,
        ignore_diags=2,
    )
    print("Eigenvector DataFrame columns:", eigvecs.columns.tolist())
    print(eigvecs[["chrom", "start", "end", "E1"]].head(10))

    # Plot E1 along chromosome
    ev = eigvecs[eigvecs["chrom"] == "chr1"].copy()
    fig, ax = plt.subplots(figsize=(12, 3))
    pos = (ev["start"] + ev["end"]) / 2 / 1e6
    ax.fill_between(pos, ev["E1"], where=ev["E1"] > 0, color="#e74c3c", alpha=0.7, label="A (E1 > 0)")
    ax.fill_between(pos, ev["E1"], where=ev["E1"] < 0, color="#3498db", alpha=0.7, label="B (E1 < 0)")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Position (Mb)"); ax.set_ylabel("E1")
    ax.set_title("A/B compartments — chr1 (E1)")
    ax.legend(frameon=False)
    plt.tight_layout(); plt.show()

except Exception as e:
    print(f"Eigenvector computation skipped for synthetic data: {e}")
    print("In real data: cooltools.eigs_cis(clr, view_df=view_df, n_eigs=3)")
```

## 5. Insulation Score and TAD Boundary Detection

Topologically Associating Domains (TADs) are genomic intervals within which contacts are enriched. TAD boundaries are detected as local minima of the *insulation score*.

**Insulation score at bin i:** average contact frequency in a square window centred on the diagonal at bin i. Minima = regions where contacts across the boundary are depleted = TAD boundaries.

**Window size:** typically 100–500 kb. Use multiple windows to capture hierarchy.

```python
windows = [5 * clr.binsize, 10 * clr.binsize]  # 5 and 10 bins in bp

try:
    insulation = cooltools.insulation(clr, window_bp=windows, view_df=view_df, ignore_diags=2)
    print(insulation.columns.tolist())

    w_col = f"log2_insulation_score_{windows[0]}"
    b_col = f"is_boundary_{windows[0]}"

    if w_col in insulation.columns:
        ins = insulation[insulation["chrom"] == "chr1"].copy()
        pos = (ins["start"] + ins["end"]) / 2 / 1e6

        fig, ax = plt.subplots(figsize=(12, 3))
        ax.plot(pos, ins[w_col], color="black", lw=1.2)
        if b_col in ins.columns:
            bnd = ins[ins[b_col]]
            ax.scatter((bnd["start"] + bnd["end"]) / 2 / 1e6,
                       bnd[w_col], color="red", zorder=5, label="TAD boundary", s=30)
        ax.set_xlabel("Position (Mb)"); ax.set_ylabel("Insulation score (log2)")
        ax.set_title(f"Insulation score — window {windows[0]//1000} kb")
        ax.legend(frameon=False)
        plt.tight_layout(); plt.show()

except Exception as e:
    print(f"Insulation score skipped for synthetic data: {e}")
    print("In real data: cooltools.insulation(clr, window_bp=[200_000], view_df=view_df)")
```
