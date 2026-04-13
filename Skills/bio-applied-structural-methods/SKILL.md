---
name: bio-applied-structural-methods
description: "This notebook covers mass spectrometry-based proteomics, from raw data to protein identification and quantification, together with computational protein engineering and the structural determination me"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/18_Proteomics_and_Structural_Methods/02_structural_methods.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Proteomics and Structural Methods

*Source: Course notebook `Tier_3_Applied_Bioinformatics/18_Proteomics_and_Structural_Methods/02_structural_methods.ipynb`*


**Tier 3 -- Applied Bioinformatics**

This notebook covers mass spectrometry-based proteomics, from raw data to protein identification and quantification, together with computational protein engineering and the structural determination methods that underpin it all. These topics are genuinely complementary: you need structural data to engineer proteins rationally, and you need proteomics to validate your engineered variants at scale.

The material is grounded in the Физико-химические методы and Белковая инженерия courses taught at the Faculty of Bioengineering and Bioinformatics (ФББ), Moscow State University.

**Prerequisites:** Tier 2.07 (PDB format, Bio.PDB, Ramachandran, DSSP), Tier 3.09 (force fields, MD concepts, docking)  
**Libraries:** `numpy`, `matplotlib`, `pandas`, `collections`  
**Topics NOT repeated from earlier modules:** PDB parsing, RMSD, DSSP, py3Dmol, force fields, energy minimization, AutoDock Vina

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from collections import Counter, defaultdict

%matplotlib inline
plt.rcParams['figure.figsize'] = (12, 5)
plt.rcParams['font.size'] = 12
np.random.seed(42)
```python

---
## Part 1: Mass Spectrometry for Proteomics

### 1.1 Ionization Methods

Mass spectrometry measures the **mass-to-charge ratio (m/z)** of gas-phase ions. To get proteins and peptides into the gas phase without destroying them, two "soft" ionization methods are used:

| Feature | MALDI | ESI |
|---------|-------|-----|
| Full name | Matrix-Assisted Laser Desorption/Ionization | Electrospray Ionization |
| Physical principle | UV laser pulse vaporizes analyte embedded in crystalline matrix | High-voltage needle sprays droplets; solvent evaporates leaving multiply charged ions |
| Charge state | Mostly singly charged (+1) | Multiple charges; large proteins appear at lower m/z |
| Sample state | Solid (co-crystallized with matrix) | Liquid (coupled directly to LC) |
| Best for | Peptide mass fingerprinting, MALDI-TOF | LC-MS/MS proteomics, intact protein analysis |
| Throughput | Medium; requires sample preparation | High; online LC coupling |

**Key formula for ESI multiply-charged ions:**

$$m/z = \frac{M + z \cdot m_H}{z}$$

where $M$ is the neutral molecular mass, $z$ is the charge state, and $m_H = 1.00728$ Da (proton mass). A 50 kDa protein ionized to $z=20$ appears at $m/z \approx 2502$.

### 1.2 Mass Analyzers

| Analyzer | Principle | Mass accuracy | Resolution | Typical use |
|----------|-----------|--------------|------------|-------------|
| **TOF** (Time-of-Flight) | Flight time ∝ √(m/z) through field-free tube | ~5–50 ppm | 10,000–40,000 | MALDI-TOF, fast survey scans |
| **Orbitrap** | Ions orbit around central electrode; frequency ∝ 1/√(m/z) | <2 ppm | 100,000–500,000 | High-resolution proteomics (Thermo Fusion, Exploris) |
| **Q-TOF** | Quadrupole selects precursor → TOF measures fragments | ~5 ppm | 20,000–50,000 | MS/MS, metabolomics |
| **Ion trap** (IT) | Ions trapped by oscillating RF field; sequential mass-selective ejection | ~100–200 ppm | ~4,000 | Fast MS/MS, MSⁿ fragmentation trees |
| **Triple quadrupole** (QQQ) | Q1 selects precursor, Q2 fragments, Q3 scans products | ~100 ppm | ~4,000 | Targeted quantification (SRM/MRM) |

Modern instruments like the Orbitrap Fusion combine multiple analyzers: the Orbitrap provides high-resolution survey scans, while the ion trap simultaneously acquires MS/MS spectra for identified precursors — a strategy called **data-dependent acquisition (DDA)**.

### 1.3 MS/MS Fragmentation: b and y Ion Series

In tandem mass spectrometry (MS/MS), a selected peptide precursor ion is fragmented by **CID** (collision-induced dissociation) or **HCD** (higher-energy collisional dissociation). Fragmentation primarily breaks the peptide backbone at amide bonds, generating two complementary series:

- **b ions** (N-terminal fragments): contain the N-terminus up to position k.
  Singly charged m/z: $m_{b_k}^{+1} = \sum_{i=1}^{k} m_{\text{res},i} + m_H$
  (sum of residue masses + one proton; b ions retain the acylium structure — no hydroxyl group)

- **y ions** (C-terminal fragments): contain the C-terminus from position k+1.
  Singly charged m/z: $m_{y_{n-k}}^{+1} = \sum_{i=k+1}^{n} m_{\text{res},i} + m_{H_2O} + m_H$
  (sum of residue masses + water + one proton; y ions carry the C-terminal OH and N-terminal H)

Where $m_{\text{res},i}$ is the **residue mass** of amino acid $i$ (monoisotopic), $m_H = 1.00728$ Da (proton mass), and $m_{H_2O} = 18.01056$ Da.

**Note:** A common source of confusion — the full neutral peptide mass = $\sum m_{\text{res}} + H_2O$ (for the two termini), while b ions do NOT add water (they lose it at the C-terminus). This means $m_{b_k} + m_{y_{n-k}} = M_{\text{peptide}} + 2m_H$ (both carry a proton).

```python
  b1   b2   b3   b4
  ↓    ↓    ↓    ↓
H-AA1-AA2-AA3-AA4-AA5-OH
          ↑    ↑    ↑
          y3   y2   y1
```python

The mass difference between adjacent b (or y) ions directly reveals the **residue mass** of the intervening amino acid, enabling de novo sequencing. For example: $m_{b_3} - m_{b_2} = m_{\text{res},3}$.

**Other fragment ion types:**
- **a ions**: b ions − CO (−27.99 Da); common in ECD/ETD
- **c/z ions**: produced by ETD/EThcD fragmentation; useful for labile modifications (phospho, glyco)

```python
# Monoisotopic residue masses (Da)
AA_MONO_MASS = {
    'A': 71.03711,  'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}

H2O = 18.01056
H   = 1.00728   # proton mass


def peptide_mass(seq: str) -> float:
    """Neutral monoisotopic mass of a peptide."""
    return sum(AA_MONO_MASS[aa] for aa in seq) + H2O


def by_ions(seq: str) -> tuple[list[float], list[float]]:
    """Calculate singly-charged b and y ion masses for a peptide."""
    n = len(seq)
    b_ions, y_ions = [], []
    for i in range(1, n):
        # b ion: N-terminal fragment (no OH, one H added as proton)
        b_mass = sum(AA_MONO_MASS[seq[j]] for j in range(i)) + H
        b_ions.append(b_mass)
        # y ion: C-terminal fragment (water + proton)
        y_mass = sum(AA_MONO_MASS[seq[j]] for j in range(i, n)) + H2O + H
        y_ions.append(y_mass)
    return b_ions, y_ions


# Demonstrate with a tryptic peptide
peptide = "ACDEFGHIK"
b, y = by_ions(peptide)

print(f"Peptide: {peptide}")
print(f"Neutral mass: {peptide_mass(peptide):.4f} Da")
print(f"[M+H]+: {peptide_mass(peptide) + H:.4f} Da")
print()
print(f"{'Ion':<8} {'m/z (b)':<12} {'Ion':<8} {'m/z (y)':<12}")
for i, (bi, yi) in enumerate(zip(b, y), 1):
    print(f"b{i:<7} {bi:<12.4f} y{len(peptide)-i:<7} {yi:<12.4f}")
```python

```python
# Visualize a simulated MS/MS spectrum
fig, ax = plt.subplots(figsize=(13, 5))

b_ions, y_ions = by_ions(peptide)

for i, mass in enumerate(b_ions, 1):
    ax.bar(mass, 100 * (0.4 + 0.6 * np.random.random()), width=0.8,
           color='steelblue', alpha=0.8)
    ax.text(mass, 100 * (0.4 + 0.6 * np.random.random()) + 3,
            f'b{i}', ha='center', fontsize=8, color='steelblue')

for i, mass in enumerate(reversed(y_ions), 1):
    ax.bar(mass, 100 * (0.4 + 0.6 * np.random.random()), width=0.8,
           color='tomato', alpha=0.8)
    ax.text(mass, 100 * (0.4 + 0.6 * np.random.random()) + 3,
            f'y{i}', ha='center', fontsize=8, color='tomato')

# Precursor [M+2H]2+
precursor_mz = (peptide_mass(peptide) + 2 * H) / 2
ax.axvline(precursor_mz, color='gray', linestyle='--', alpha=0.5, label=f'[M+2H]²⁺ = {precursor_mz:.1f}')

ax.set_xlabel('m/z')
ax.set_ylabel('Relative intensity (%)')
ax.set_title(f'Simulated MS/MS spectrum: {peptide}')
b_patch = mpatches.Patch(color='steelblue', label='b ions (N-terminal)')
y_patch = mpatches.Patch(color='tomato', label='y ions (C-terminal)')
ax.legend(handles=[b_patch, y_patch])
ax.set_ylim(0, 130)
plt.tight_layout()
plt.show()
```python

---
## Part 2: Protein Identification and Database Searching

### 2.1 Bottom-Up Proteomics Workflow

The dominant proteomics strategy is **bottom-up (shotgun) proteomics**:

```python
Protein mixture
    ↓  (1) Reduction & alkylation (DTT + iodoacetamide)
Denatured protein
    ↓  (2) Enzymatic digestion (trypsin: cleaves after K, R)
Peptide mixture
    ↓  (3) LC separation (reversed-phase nano-LC, 60–120 min gradient)
Eluting peptides
    ↓  (4) ESI-MS survey scan → select precursor ions (DDA)
MS1 spectrum
    ↓  (5) CID/HCD fragmentation
MS/MS spectra
    ↓  (6) Database search (Mascot, X!Tandem, MaxQuant)
Peptide-spectrum matches (PSMs)
    ↓  (7) FDR filtering → peptide / protein lists
Quantitative protein abundance table
```python

**Trypsin** is the enzyme of choice: it cleaves specifically on the C-terminal side of lysine (K) and arginine (R), **except when followed by proline (P)**. This generates peptides of 6–25 residues — the optimal size range for LC-MS/MS.

```python
def trypsin_digest(sequence: str, missed_cleavages: int = 0) -> list[str]:
    """
    In silico trypsin digestion.
    Cleaves after K or R, unless followed by P.
    missed_cleavages: number of allowed missed cleavages (joined fragments).
    """
    seq = sequence.upper()
    # Find cleavage sites: after K/R not followed by P
    sites = [0]
    for i in range(len(seq) - 1):
        if seq[i] in ('K', 'R') and seq[i + 1] != 'P':
            sites.append(i + 1)
    sites.append(len(seq))

    # Generate peptide fragments (zero missed cleavages)
    fragments = [seq[sites[i]:sites[i + 1]] for i in range(len(sites) - 1)]
    fragments = [f for f in fragments if f]  # remove empty strings

    if missed_cleavages == 0:
        return fragments

    # Add missed cleavage peptides
    result = list(fragments)
    for mc in range(1, missed_cleavages + 1):
        for i in range(len(fragments) - mc):
            joined = ''.join(fragments[i:i + mc + 1])
            result.append(joined)
    return sorted(set(result), key=lambda x: sequence.index(x))


# Example: human ubiquitin (76 aa, well-characterized)
ubiquitin = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)

peptides_0mc = trypsin_digest(ubiquitin, missed_cleavages=0)
peptides_1mc = trypsin_digest(ubiquitin, missed_cleavages=1)

print(f"Ubiquitin ({len(ubiquitin)} aa)")
print(f"Tryptic peptides (0 MC): {len(peptides_0mc)}")
print(f"Tryptic peptides (1 MC): {len(peptides_1mc)}")
print()
masses = [(pep, peptide_mass(pep)) for pep in peptides_0mc]
print(f"{'Peptide':<30} {'Mass (Da)':<12} {'Length':<8}")
print("-" * 52)
for pep, mass in sorted(masses, key=lambda x: -x[1]):
    print(f"{pep:<30} {mass:<12.2f} {len(pep):<8}")
```python

### 2.2 Peptide Mass Fingerprinting (PMF)

**PMF** (originally used with MALDI-TOF) identifies a protein by matching the observed set of tryptic peptide masses to theoretical masses from a database. No MS/MS is required — the pattern of masses is the "fingerprint".

Steps:
1. Digest with trypsin → measure peptide masses by MALDI-TOF
2. For each database protein, compute theoretical tryptic masses
3. Count matching masses within tolerance (typically ±50–100 ppm)
4. Score and rank candidates; highest score = best candidate

PMF works well for simple mixtures (single protein from a gel band). For complex mixtures, **LC-MS/MS database searching** (Mascot, X!Tandem, MaxQuant) is required.

### 2.3 Database Search Engines

| Engine | Algorithm | Key feature |
|--------|-----------|-------------|
| **Mascot** (Matrix Science) | Probability-based scoring (Mowse) | Commercial; widely adopted; good statistics |
| **X!Tandem** | Hyperscore (dot product of matched ions) | Free; fast; good for large databases |
| **MaxQuant** | Andromeda search engine built in | Free; integrated LFQ; SILAC quantification |
| **MSFragger** | Fragment index search | Ultra-fast; open modification searches |
| **Comet** | Open-source Sequest-style | Fast; command-line friendly |

All engines compare observed MS/MS spectra against **theoretical spectra** predicted for every peptide in the database at the given precursor mass (±10–20 ppm for Orbitrap).

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
