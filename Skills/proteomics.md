---
name: proteomics
description: Mass spectrometry-based proteomics — peptide mass calculation, b/y ion series, in silico trypsin digestion, target-decoy FDR, LFQ/SILAC/TMT quantification, PMF search, volcano plots
---

# Proteomics & Mass Spectrometry

## When to Use
- Mass spectrometry data processing: peptide mass calculation, b/y ions, in silico digestion
- Protein identification: peptide mass fingerprinting (PMF), FDR estimation with target-decoy
- Quantitative proteomics: LFQ normalization, SILAC/TMT workflows, volcano/MA plots
- Protein engineering computational design: conservation scoring from MSA, ΔΔG-guided mutagenesis

---

## Quick Reference

### Mass Spectrometry Instruments
| Instrument | Analyzer | Resolution | Mass accuracy | Best use |
|-----------|----------|-----------|--------------|----------|
| MALDI-TOF | TOF | 10,000–40,000 | 5–50 ppm | PMF, microbial ID |
| Q-TOF | Quadrupole + TOF | 20,000–50,000 | ~5 ppm | MS/MS, metabolomics |
| Orbitrap Exploris / Fusion | Orbitrap (+ ion trap) | 100,000–500,000 | < 2 ppm | High-res proteomics |
| Triple quadrupole (QQQ) | Q1-Q2-Q3 | ~4,000 | ~100 ppm | Targeted MRM/SRM |

### Proteomics Workflow (Bottom-Up Shotgun)
```
Protein mixture
  → Reduction & alkylation (DTT + iodoacetamide)
  → Trypsin digestion (cleaves after K/R, not before P)
  → LC-MS/MS (nano-LC + Orbitrap)
  → Database search (Mascot / MSFragger / Andromeda)
  → FDR filtering (target-decoy, 1% PSM + protein)
  → Quantification (LFQ / SILAC / TMT)
  → Statistical analysis (Perseus / R)
```

### Proteomics Software
| Tool | Category | Key feature |
|------|----------|-------------|
| MaxQuant | Search + LFQ | Free; Andromeda engine; iBAQ, SILAC, LFQ |
| Mascot | Database search | Commercial; probability-based scoring |
| MSFragger | Database search | Fragment index; ultra-fast; open mod searches |
| Perseus | Statistical analysis | R-like workflows for proteomics matrices |
| Spectronaut | DIA analysis | Spectral library-based; handles chimeric MS/MS |
| DIA-NN | DIA analysis | Neural-network-based; works without library |
| PRIDE | Data repository | European proteomics archive; `PXD` accessions |

### Quantification Methods
| Method | Label | MS level | Ratio accuracy | Note |
|--------|-------|----------|---------------|------|
| LFQ (label-free) | None | MS1 | Good | MaxQuant LFQ intensity |
| SILAC | ¹³C/¹⁵N amino acids | MS1 | Excellent | Mix before digestion |
| TMT / iTRAQ | Isobaric tags | MS2 | Good | Up to 16-plex (TMTpro) |
| DIA (label-free) | None | MS1+MS2 | Good | All-ion fragmentation |

---

## Key Patterns

### Monoisotopic Masses & Peptide Mass
```python
AA_MONO_MASS = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}
H2O = 18.01056   # Da
H   = 1.00728    # proton mass (Da)

def peptide_mass(seq: str) -> float:
    return sum(AA_MONO_MASS[aa] for aa in seq.upper()) + H2O

def mz(mass: float, z: int) -> float:
    """ESI m/z for charge state z."""
    return (mass + z * H) / z
```

### b/y Ion Series
```python
def by_ions(seq: str) -> tuple[list[float], list[float]]:
    """Return (b_ions, y_ions) for singly charged fragment ions."""
    n = len(seq)
    b = [sum(AA_MONO_MASS[seq[j]] for j in range(i)) + H
         for i in range(1, n)]
    y = [sum(AA_MONO_MASS[seq[j]] for j in range(i, n)) + H2O + H
         for i in range(1, n)]
    return b, y
```

### In Silico Trypsin Digestion
```python
def trypsin_digest(sequence: str, missed_cleavages: int = 0) -> list[str]:
    """Cleave after K/R, not before P."""
    seq = sequence.upper()
    sites = [0]
    for i in range(len(seq) - 1):
        if seq[i] in ('K', 'R') and seq[i + 1] != 'P':
            sites.append(i + 1)
    sites.append(len(seq))
    frags = [seq[sites[i]:sites[i + 1]]
             for i in range(len(sites) - 1) if sites[i + 1] > sites[i]]
    if missed_cleavages == 0:
        return frags
    result = list(frags)
    for mc in range(1, missed_cleavages + 1):
        for i in range(len(frags) - mc):
            result.append(''.join(frags[i:i + mc + 1]))
    return sorted(set(result), key=lambda x: sequence.upper().index(x))
```

### Target-Decoy FDR
```python
import numpy as np

def fdr_curve(scores: np.ndarray, is_decoy: np.ndarray) -> np.ndarray:
    """Compute running FDR at each score threshold (descending order)."""
    order = np.argsort(scores)[::-1]
    t = d = 0
    fdrs = np.ones(len(scores))
    for rank, i in enumerate(order):
        if is_decoy[i]:
            d += 1
        else:
            t += 1
        fdrs[rank] = d / t if t > 0 else 1.0
    return fdrs  # indexed in sorted order

def score_at_fdr(scores: np.ndarray, is_decoy: np.ndarray,
                 fdr_cutoff: float = 0.01) -> float:
    """Return lowest score that keeps FDR ≤ fdr_cutoff."""
    order = np.argsort(scores)[::-1]
    fdrs = fdr_curve(scores, is_decoy)
    passing = [scores[order[i]] for i, f in enumerate(fdrs) if f <= fdr_cutoff]
    return min(passing) if passing else float('inf')
```

---

## Code Templates

### Simple PMF Database Search
```python
def pmf_search(
    observed_masses: list[float],
    database: dict[str, str],       # {protein_name: sequence}
    tolerance_ppm: float = 20.0,
    missed_cleavages: int = 1,
) -> list[tuple[str, int, float]]:
    """Return (name, n_matched, coverage) sorted best-first."""
    results = []
    obs = np.array(sorted(observed_masses))
    for name, seq in database.items():
        peptides = trypsin_digest(seq, missed_cleavages)
        theo = np.array(sorted(peptide_mass(p) for p in peptides
                               if all(aa in AA_MONO_MASS for aa in p.upper())))
        matched = sum(
            1 for m in obs
            if np.any(np.abs(theo - m) <= m * tolerance_ppm * 1e-6)
        )
        results.append((name, matched, matched / len(obs)))
    return sorted(results, key=lambda x: -x[1])
```

### LFQ Median Normalization
```python
def median_normalize(log2_matrix: np.ndarray) -> np.ndarray:
    """
    Column-wise median centering of log2 intensity matrix.
    Rows = proteins, columns = samples.
    """
    return log2_matrix - np.nanmedian(log2_matrix, axis=0, keepdims=True)
```

### Conservation Scoring from MSA
```python
from collections import Counter

def shannon_entropy_col(column: list[str]) -> float:
    """Per-column Shannon entropy (bits) from MSA column."""
    col = [aa for aa in column if aa not in ('-', 'X', '.')]
    if not col:
        return 0.0
    counts = Counter(col)
    n = len(col)
    return -sum((c / n) * np.log2(c / n) for c in counts.values())

def mutability_scores(msa: list[str]) -> np.ndarray:
    """Per-column mutability in [0, 1]. 0 = invariant, 1 = fully random."""
    aln_len = len(msa[0])
    H_max = np.log2(20)
    return np.array([
        shannon_entropy_col([s[i] for s in msa]) / H_max
        for i in range(aln_len)
    ])
```

### Volcano Plot for Proteomics
```python
import matplotlib.pyplot as plt

def volcano_plot(log2fc: np.ndarray, pvalues: np.ndarray,
                 fc_thresh: float = 1.0, p_thresh: float = 0.05) -> None:
    neg_log10p = -np.log10(np.clip(pvalues, 1e-300, 1.0))
    sig = (np.abs(log2fc) > fc_thresh) & (pvalues < p_thresh)
    colors = np.where(sig, 'tomato', 'steelblue')
    plt.figure(figsize=(8, 6))
    plt.scatter(log2fc, neg_log10p, c=colors, alpha=0.6, s=15)
    plt.axvline( fc_thresh, color='gray', ls='--', lw=0.8)
    plt.axvline(-fc_thresh, color='gray', ls='--', lw=0.8)
    plt.axhline(-np.log10(p_thresh), color='gray', ls='--', lw=0.8)
    plt.xlabel('log₂ fold change')
    plt.ylabel('-log₁₀(p-value)')
    plt.title(f'Volcano — {sig.sum()} significant proteins')
    plt.tight_layout(); plt.show()
```

---

## Common Pitfalls

- **FDR at PSM level ≠ protein level**: after roll-up, protein-level FDR is always higher. Apply FDR filtering at each level independently.
- **Normalization pitfall**: never normalize between samples without checking that total protein content is comparable. Use median or quantile normalization; avoid sum normalization with abundant contaminants.
- **Missing values in LFQ**: do not impute with zeros — use minimum-value imputation or left-censored imputation (Perseus: width = 0.3, shift = 1.8 in log2 space).
- **SILAC ratio compression**: high-abundance proteins show compressed heavy/light ratios in DDA due to co-isolated interfering ions; verify with MS3.
- **Isobaric amino acids**: Leucine and Isoleucine are mass-identical (113.084 Da). MS/MS cannot distinguish them. Annotate as "I/L" unless ETD data is available.
- **Trypsin missed cleavages**: always allow at least 1 missed cleavage; consecutive K/R residues are common and 0 MC loses real peptides.
- **Variable modifications inflate search space**: limit to ≤ 3 variable mods per PSM in standard searches; too many mods exponentially increase search time and FDR.
- **TMT ratio compression in DDA**: isobaric interference from co-fragmented peptides compresses fold changes; use SPS-MS3 (Fusion Lumos) or DIA-based acquisition for accurate quantification.

---

## Related Skills
- `genome-assembly` — de novo assembly metrics (companion Module 17 skill)
- `structural-bioinformatics` — PDB parsing, RMSD, DSSP, Ramachandran plots
- `biochemistry-enzymology` — enzyme kinetics (Km, Vmax), thermodynamic stability
- `sequence-alignment` — MSA used for conservation scoring and protein engineering
