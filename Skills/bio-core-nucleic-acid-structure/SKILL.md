---
name: bio-core-nucleic-acid-structure
description: DNA helix forms (A/B/Z), groove geometry, nearest-neighbor thermodynamics, RNA secondary structure elements and dot-bracket notation
tool_type: python
primary_tool: NumPy
---

# Nucleic Acid Structure

## Pitfalls

- **Strand polarity**: DNA synthesized 5'->3'; duplex strands are antiparallel. Always search both strands for motifs.
- **RNA MFE is only one possibility**: Same RNA can fold multiple ways. Use ViennaRNA/RNAfold for MFE; consider suboptimal structures for regulatory RNAs.
- **Dot-bracket notation**: `(` pairs with `)`, `.` = unpaired. Pseudoknots need `[]{}` brackets; many tools cannot handle them.
- **CpG islands**: CpG O/E > 0.6, GC > 55%, length > 200 bp. Most CpGs in mammals are methylated and depleted by 5mC->T deamination.
- **Coordinate systems**: BED = 0-based half-open; VCF/GFF = 1-based inclusive.

## DNA Helix Forms

| Parameter | A-DNA | B-DNA | Z-DNA |
|-----------|-------|-------|-------|
| Handedness | Right | Right | **Left** |
| bp/turn | 11 | 10.5 | 12 |
| Rise/bp (A) | 2.6 | 3.4 | 3.7 |
| Diameter (A) | 23 | 20 | 18 |
| Major groove | Narrow, deep | Wide, deep | Flat |
| Minor groove | Wide, shallow | Narrow, deep | Narrow, deep |
| Sugar pucker | C3'-endo | C2'-endo | Alternating |
| Conditions | Dehydrated, RNA-DNA hybrids | Physiological | High salt, alt pur-pyr |

B-DNA is the standard form. A-form is adopted by RNA duplexes (2'-OH forces C3'-endo). Z-DNA forms in alternating purine-pyrimidine under high salt.

## Groove Recognition

All four base pairs are distinguishable in the **major groove** (unique H-bond donor/acceptor pattern) but A-T vs T-A and G-C vs C-G are hard to distinguish in the minor groove. Most TFs bind the major groove.

## Nearest-Neighbor Thermodynamics

Stacking interactions contribute more to duplex stability than H-bonding. Energy depends on adjacent base pair identity.

```python
import numpy as np

# SantaLucia (1998) stacking free energies (kcal/mol, 37C, 1M NaCl)
NN_ENERGIES = {
    'AA': -1.0, 'AT': -0.88, 'AG': -1.28, 'AC': -1.44,
    'TA': -0.58, 'TT': -1.0,  'TG': -1.45, 'TC': -1.30,
    'GA': -1.30, 'GT': -1.44, 'GG': -1.84, 'GC': -2.24,
    'CA': -1.45, 'CT': -1.28, 'CG': -2.17, 'CC': -1.84,
}

def nearest_neighbor_dG(seq):
    """Estimate duplex free energy. Initiation: +1.96 kcal/mol."""
    seq = seq.upper()
    dG = 1.96
    for i in range(len(seq) - 1):
        dG += NN_ENERGIES.get(seq[i:i+2], -1.0)
    return dG
```

## Helix Dimension Calculator

```python
HELIX_PARAMS = {
    'A-DNA': {'bp_per_turn': 11,   'rise': 2.6, 'diameter': 23},
    'B-DNA': {'bp_per_turn': 10.5, 'rise': 3.4, 'diameter': 20},
    'Z-DNA': {'bp_per_turn': 12,   'rise': 3.7, 'diameter': 18},
}

def helix_dimensions(num_bp, form='B-DNA'):
    p = HELIX_PARAMS[form]
    length_A = num_bp * p['rise']
    return {'length_nm': length_A / 10.0, 'turns': num_bp / p['bp_per_turn']}
# Human genome (3.2 Gbp B-DNA): ~2 meters packed into a ~6 um nucleus
```

## RNA Secondary Structure Elements

1. **Stem** -- paired region
2. **Hairpin loop** -- stem capped by unpaired loop (most common)
3. **Bulge** -- unpaired bases on one strand
4. **Internal loop** -- mismatches on both strands
5. **Multi-stem junction** -- 3+ stems meet
6. **Pseudoknot** -- base pairs cross; requires extended dot-bracket `[]`

Dot-bracket example: `((((....))))` = 4-bp stem + 4-nt loop
