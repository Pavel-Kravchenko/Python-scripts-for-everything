---
name: genome-assembly-proteomics
description: Genome assembly metrics, de Bruijn graphs, mass spectrometry data analysis, and protein identification workflows
---

# Genome Assembly & Proteomics

## When to Use
- De novo genome assembly or evaluating existing assemblies (N50, L50, BUSCO, assembly graph)
- Mass spectrometry data processing: peptide mass calculation, b/y ions, in silico digestion
- Protein identification: peptide mass fingerprinting, FDR estimation with target-decoy
- Quantitative proteomics: LFQ normalization, SILAC/TMT workflows, volcano/MA plots
- Protein engineering computational design: conservation scoring from MSA, ΔΔG-guided mutagenesis

---

## Quick Reference

### Assembly Quality Metrics
| Metric | Definition | Good value (typical eukaryote) |
|--------|-----------|-------------------------------|
| **N50** | Length L such that 50% of assembled bases are in contigs ≥ L | > 1 Mb (chromosome-scale) |
| **L50** | Number of contigs whose combined length = N50 | Lower is better |
| **NG50** | N50 normalized to genome size (uses estimated genome size) | Comparable across species |
| **BUSCO** | Fraction of conserved single-copy ortholog genes present | > 95% complete |
| **QV score** | Phred-like per-base quality of assembly | > 40 (Merqury) |
| **k-mer completeness** | Fraction of k-mers from reads present in assembly | > 95% |

### Mass Spectrometry Instruments
| Instrument | Analyzer | Resolution | Mass accuracy | Best use |
|-----------|----------|-----------|--------------|----------|
| MALDI-TOF | TOF | 10,000–40,000 | 5–50 ppm | PMF, microbial ID |
| Q-TOF | Quadrupole + TOF | 20,000–50,000 | ~5 ppm | MS/MS, metabolomics |
| Orbitrap Exploris / Fusion | Orbitrap (+ ion trap) | 100,000–500,000 | < 2 ppm | High-res proteomics |
| Triple quadrupole (QQQ) | Q1-Q2-Q3 | ~4,000 | ~100 ppm | Targeted MRM/SRM |

### Proteomics Tools
| Tool | Category | Key feature |
|------|----------|-------------|
| MaxQuant | Search + LFQ quantification | Free; integrated Andromeda; iBAQ, SILAC, LFQ |
| Mascot | Database search | Commercial; probability-based scoring; industry standard |
| MSFragger | Database search | Fragment index; ultra-fast; open mod searches |
| Perseus | Statistical analysis | R-like workflows for proteomics matrices |
| Spectronaut | DIA analysis | Spectral library-based; handles chimeric MS/MS |
| DIA-NN | DIA analysis | Fast neural-network-based; works without library |
| PRIDE | Data repository | European proteomics archive; `PXD` accession numbers |

---

## Key Patterns

### N50 Calculation
```python
def n50(contig_lengths: list[int], genome_size: int | None = None) -> tuple[int, int]:
    """Return (N50, L50). If genome_size given, computes NG50."""
    lengths = sorted(contig_lengths, reverse=True)
    target = (genome_size or sum(lengths)) / 2
    cumsum = 0
    for i, length in enumerate(lengths, 1):
        cumsum += length
        if cumsum >= target:
            return length, i  # (N50, L50)
    return lengths[-1], len(lengths)
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
    fragments = [seq[sites[i]:sites[i+1]] for i in range(len(sites)-1) if sites[i+1] > sites[i]]
    if missed_cleavages == 0:
        return fragments
    result = list(fragments)
    for mc in range(1, missed_cleavages + 1):
        for i in range(len(fragments) - mc):
            result.append(''.join(fragments[i:i+mc+1]))
    return sorted(set(result), key=lambda x: sequence.index(x))
```

### Peptide Mass and b/y Ions
```python
AA_MONO_MASS = {
    'A':71.03711,'R':156.10111,'N':114.04293,'D':115.02694,'C':103.00919,
    'E':129.04259,'Q':128.05858,'G':57.02146, 'H':137.05891,'I':113.08406,
    'L':113.08406,'K':128.09496,'M':131.04049,'F':147.06841,'P':97.05276,
    'S':87.03203, 'T':101.04768,'W':186.07931,'Y':163.06333,'V':99.06841,
}
H2O, H = 18.01056, 1.00728  # H = proton mass

def peptide_mass(seq: str) -> float:
    return sum(AA_MONO_MASS[aa] for aa in seq) + H2O

def mz(mass: float, z: int) -> float:
    """ESI m/z for charge state z."""
    return (mass + z * H) / z

def by_ions(seq: str) -> tuple[list[float], list[float]]:
    n = len(seq)
    b = [sum(AA_MONO_MASS[seq[j]] for j in range(i)) + H for i in range(1, n)]
    y = [sum(AA_MONO_MASS[seq[j]] for j in range(i, n)) + H2O + H for i in range(1, n)]
    return b, y
```

### Target-Decoy FDR
```python
def estimate_fdr(target_hits: int, decoy_hits: int) -> float:
    """FDR = decoy_hits / target_hits at a given score threshold."""
    return decoy_hits / target_hits if target_hits > 0 else 1.0

# Apply at decreasing score thresholds to build FDR curve:
# psms sorted by score (descending); decoy flag for each PSM
def fdr_curve(scores, is_decoy):
    import numpy as np
    order = np.argsort(scores)[::-1]
    t = d = 0
    fdrs = []
    for i in order:
        if is_decoy[i]: d += 1
        else: t += 1
        fdrs.append(d / t if t > 0 else 1.0)
    # 1% FDR threshold
    cutoff_score = scores[order[next(i for i, f in enumerate(fdrs) if f > 0.01) - 1]]
    return cutoff_score
```

### Conservation Scoring from MSA
```python
from collections import Counter
import numpy as np

def shannon_entropy(column: list[str]) -> float:
    col = [aa for aa in column if aa not in ('-', 'X', '.')]
    if not col: return 0.0
    counts = Counter(col)
    n = len(col)
    return -sum((c/n) * np.log2(c/n) for c in counts.values())

def mutability_scores(msa: list[str]) -> np.ndarray:
    """Per-column mutability in [0, 1]. 0 = invariant, 1 = fully random."""
    aln_len = len(msa[0])
    H_max = np.log2(20)
    return np.array([shannon_entropy([s[i] for s in msa]) / H_max
                     for i in range(aln_len)])
```

### LFQ Normalization
```python
def median_normalize(log2_matrix: np.ndarray) -> np.ndarray:
    """
    Column-wise median centering of log2 intensity matrix.
    Rows = proteins, columns = samples.
    """
    return log2_matrix - np.nanmedian(log2_matrix, axis=0, keepdims=True)
```

---

## Code Templates

### Assembly Statistics Report
```python
def assembly_stats(contig_lengths: list[int], genome_size: int | None = None) -> dict:
    lengths = sorted(contig_lengths, reverse=True)
    total = sum(lengths)
    n50_val, l50_val = n50(lengths, genome_size)
    return {
        'n_contigs':     len(lengths),
        'total_bases':   total,
        'largest':       lengths[0],
        'N50':           n50_val,
        'L50':           l50_val,
        'N90':           n50(lengths, genome_size=int(total * 0.9 / 0.5))[0],  # approx
        'mean_length':   int(np.mean(lengths)),
        'median_length': int(np.median(lengths)),
    }
```

### Simple PMF Search
```python
def pmf_search(
    observed_masses: list[float],
    database: dict[str, str],        # {protein_name: sequence}
    tolerance_ppm: float = 20.0,
    missed_cleavages: int = 1,
) -> list[tuple[str, int, float]]:
    """Returns (name, n_matched, coverage) sorted best-first."""
    import numpy as np
    results = []
    obs = np.array(sorted(observed_masses))
    for name, seq in database.items():
        peptides = trypsin_digest(seq, missed_cleavages)
        theo = np.array(sorted(peptide_mass(p) for p in peptides
                               if all(aa in AA_MONO_MASS for aa in p)))
        matched = sum(1 for m in obs if np.any(np.abs(theo - m) <= m * tolerance_ppm * 1e-6))
        results.append((name, matched, matched / len(obs)))
    return sorted(results, key=lambda x: -x[1])
```

### Volcano Plot
```python
def volcano_plot(log2fc: np.ndarray, pvalues: np.ndarray,
                 fc_thresh: float = 1.0, p_thresh: float = 0.05) -> None:
    import matplotlib.pyplot as plt
    import numpy as np
    neg_log10p = -np.log10(np.clip(pvalues, 1e-300, 1.0))
    sig = (np.abs(log2fc) > fc_thresh) & (pvalues < p_thresh)
    colors = np.where(sig, 'tomato', 'steelblue')
    plt.figure(figsize=(8, 6))
    plt.scatter(log2fc, neg_log10p, c=colors, alpha=0.6, s=15)
    plt.axvline( fc_thresh, color='gray', linestyle='--', lw=0.8)
    plt.axvline(-fc_thresh, color='gray', linestyle='--', lw=0.8)
    plt.axhline(-np.log10(p_thresh), color='gray', linestyle='--', lw=0.8)
    plt.xlabel('log₂ fold change'); plt.ylabel('-log₁₀(p-value)')
    plt.title(f'Volcano plot — {sig.sum()} significant proteins')
    plt.tight_layout(); plt.show()
```

---

## Common Pitfalls

### Assembly
- **k-mer size**: too small → repetitive regions collapse; too large → sparse coverage misses low-cov regions. For Illumina 150 bp reads, k=51–77 is typical; for long reads use different assemblers (Flye, Hifiasm).
- **Coverage requirement**: minimum ~30× for reliable de Bruijn graph; < 10× means many k-mers appear only once and are discarded as errors.
- **N50 alone is insufficient**: a 10-contig assembly with N50=1 Mb but 5% BUSCO is worthless. Always pair N50 with BUSCO completeness.
- **Heterozygosity**: diploid organisms produce bubble structures in the assembly graph. Use `--haplotype` modes (hifiasm) or purge haplotigs (purge_dups).
- **Contamination**: always check for unexpected organisms with Kraken2 or BlobTools before submitting.

### Proteomics
- **FDR 1% at PSM level ≠ 1% at protein level**: after roll-up, protein-level FDR is always higher. Apply FDR filtering at each level independently.
- **Normalization**: never normalize between samples without checking that total protein content is comparable. Use median or quantile normalization; avoid sum normalization with highly abundant contaminants.
- **Missing values in LFQ**: do not impute with zeros — use minimum-value imputation or left-censored imputation (Perseus: width=0.3, shift=1.8 in log2 space).
- **SILAC ratio compression**: high-abundance proteins show compressed heavy/light ratios in DDA due to co-isolated interfering ions. Verify with orthogonal method or use high-resolution MS3.
- **Isobaric amino acids**: Leucine and Isoleucine are mass-identical (113.084 Da). MS/MS cannot distinguish them. Annotate as "I/L" unless ETD data available.
- **Trypsin missed cleavages**: always allow at least 1 missed cleavage in the database search. Consecutive K/R residues are common; 0 MC loses real peptides.
- **Variable modifications inflate search space**: limit to ≤ 3 variable mods per PSM in standard searches; too many mods exponentially increase search time and false positives.

---

## Related Skills
- `ngs-variant-calling` — NGS fundamentals, short-read alignment, variant calling (prerequisite for assembly context)
- `structural-bioinformatics` — PDB parsing, RMSD, DSSP, Ramachandran plots, protein structure analysis
- `biochemistry-enzymology` — enzyme kinetics (Km, Vmax, inhibition), thermodynamic stability, linked to protein engineering
