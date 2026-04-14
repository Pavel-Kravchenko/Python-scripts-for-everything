---
name: bio-applied-structural-methods
description: "Proteomics and structural methods: mass spectrometry ionization, MS/MS b/y ion calculation, trypsin digestion, database search engines, peptide mass fingerprinting"
tool_type: python
primary_tool: Matplotlib
---

# Proteomics and Structural Methods

## Mass Spectrometry Ionization

| Feature | MALDI | ESI |
|---------|-------|-----|
| Charge state | Mostly +1 | Multiple charges |
| Sample state | Solid (matrix) | Liquid (LC-coupled) |
| Best for | Peptide mass fingerprinting | LC-MS/MS proteomics |

ESI m/z formula: `m/z = (M + z * 1.00728) / z`

## Mass Analyzers

| Analyzer | Mass accuracy | Resolution | Use |
|----------|--------------|------------|-----|
| TOF | 5-50 ppm | 10k-40k | MALDI-TOF, fast scans |
| Orbitrap | <2 ppm | 100k-500k | High-res proteomics |
| Q-TOF | ~5 ppm | 20k-50k | MS/MS, metabolomics |
| Ion trap | ~100 ppm | ~4k | Fast MS/MS, MSn |
| Triple quad (QQQ) | ~100 ppm | ~4k | Targeted SRM/MRM |

## MS/MS: b and y Ion Series

```python
AA_MONO_MASS = {
    'A': 71.03711,  'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}
H2O = 18.01056; H = 1.00728

def peptide_mass(seq):
    return sum(AA_MONO_MASS[aa] for aa in seq) + H2O

def by_ions(seq):
    n = len(seq)
    b_ions, y_ions = [], []
    for i in range(1, n):
        b_ions.append(sum(AA_MONO_MASS[seq[j]] for j in range(i)) + H)
        y_ions.append(sum(AA_MONO_MASS[seq[j]] for j in range(i, n)) + H2O + H)
    return b_ions, y_ions
```

Key: `b_k + y_(n-k) = M_peptide + 2*H`. Mass difference between adjacent b (or y) ions = residue mass of the intervening amino acid.

## Trypsin Digestion

Trypsin cleaves after K/R, **except when followed by P**. Optimal peptide range: 6-25 residues.

```python
def trypsin_digest(sequence, missed_cleavages=0):
    seq = sequence.upper()
    sites = [0]
    for i in range(len(seq) - 1):
        if seq[i] in ('K', 'R') and seq[i + 1] != 'P':
            sites.append(i + 1)
    sites.append(len(seq))
    fragments = [seq[sites[i]:sites[i+1]] for i in range(len(sites)-1)]
    if missed_cleavages == 0:
        return [f for f in fragments if f]
    result = list(fragments)
    for mc in range(1, missed_cleavages + 1):
        for i in range(len(fragments) - mc):
            result.append(''.join(fragments[i:i+mc+1]))
    return sorted(set(result), key=lambda x: sequence.index(x))
```

## Database Search Engines

| Engine | Key feature |
|--------|-------------|
| Mascot | Probability-based (Mowse), commercial |
| MaxQuant | Integrated LFQ/SILAC, Andromeda engine |
| MSFragger | Ultra-fast, open modification searches |
| X!Tandem | Free, hyperscore dot product |
| Comet | Open-source Sequest-style |

## Pitfalls

- **b ions have no water**: full peptide mass = sum(residues) + H2O, but b ions lose the C-terminal OH
- **FDR in proteomics**: target-decoy approach; 1% FDR at PSM level standard
- **Trypsin missed cleavages**: allow 1-2 in database search; KP/RP sites are not cleaved
