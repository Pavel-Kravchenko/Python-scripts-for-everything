---
name: bio-applied-genetic-engineering-in-silico
description: "This notebook implements core molecular cloning and genome-editing workflows entirely in Python. All algorithms are built from first principles — no external bioinformatics libraries required."
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/14_Genetic_Engineering_In_Silico/01_genetic_engineering_in_silico.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Genetic Engineering In Silico

*Source: Course notebook `Tier_3_Applied_Bioinformatics/14_Genetic_Engineering_In_Silico/01_genetic_engineering_in_silico.ipynb`*


## Tier 3 - Applied Bioinformatics

This notebook implements core molecular cloning and genome-editing workflows entirely in Python.
All algorithms are built from first principles — no external bioinformatics libraries required.

## Learning Objectives

By the end of this notebook you will be able to:

1. Simulate restriction digestion and predict gel electrophoresis banding patterns
2. Design PCR primers using Tm models ranging from the 4+2 rule to nearest-neighbor thermodynamics
3. Evaluate primer quality: GC content, 3' runs, dimer formation, and hairpin potential
4. Model traditional RE-based cloning and verify reading-frame continuity after ligation
5. Identify CRISPR guide RNA candidates for SpCas9, SaCas9, and Cas12a and rank them by on-target score
6. Optimize a coding sequence for heterologous expression using organism-specific codon usage tables
7. Design Gibson Assembly overlaps for multi-fragment constructs

## 1. In Silico Restriction Digestion

### 1.1 Restriction Enzymes: Recognition Sites and Cut Geometry

Restriction endonucleases recognise specific short DNA sequences (typically 4–8 bp) and cleave both strands.
The recognition sequences of Type II enzymes — the workhorses of cloning — are **palindromic**: the sequence reads the same 5'→3' on both strands.

Cut geometry determines the overhang type:

| Overhang | Example | Description |
|----------|---------|-------------|
| 5' overhang (sticky) | EcoRI | Top strand cut upstream of bottom strand |
| 3' overhang (sticky) | KpnI  | Top strand cut downstream of bottom strand |
| Blunt | SmaI | Both strands cut at the same position |

We represent each enzyme as a recognition sequence plus two cut positions: the position on the **top strand** and the position on the **bottom strand** (both counted from the start of the recognition sequence, 0-indexed).

```python
import re
import math
import itertools
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Restriction enzyme database
# Each entry: (recognition_seq, cut_top, cut_bottom)
# cut positions are 0-indexed from start of recognition site on top strand
# cut_bottom is on the complementary strand read 5'->3' (i.e. counted from end of site)
RESTRICTION_ENZYMES = {
    'EcoRI':  ('GAATTC', 1, 5),   # 5' AATT overhang
    'BamHI':  ('GGATCC', 1, 5),   # 5' GATC overhang
    'HindIII':('AAGCTT', 1, 5),   # 5' AGCT overhang
    'SalI':   ('GTCGAC', 1, 5),   # 5' TCGA overhang
    'XhoI':   ('CTCGAG', 1, 5),   # 5' TCGA overhang (compatible with SalI)
    'NcoI':   ('CCATGG', 1, 5),   # 5' CATG overhang
    'NdeI':   ('CATATG', 2, 4),   # 5' AT overhang
    'SmaI':   ('CCCGGG', 3, 3),   # blunt
    'XmaI':   ('CCCGGG', 1, 5),   # 5' CCGG overhang (isoschizomer of SmaI, different cut)
    'KpnI':   ('GGTACC', 5, 1),   # 3' GTAC overhang
    'SacI':   ('GAGCTC', 5, 1),   # 3' AGCT overhang
    'NotI':   ('GCGGCCGC', 2, 6), # 5' GGCC overhang
    'XbaI':   ('TCTAGA', 1, 5),   # 5' CTAG overhang
    'SpeI':   ('ACTAGT', 1, 5),   # 5' CTAG overhang (compatible with XbaI, scar site)
    'PstI':   ('CTGCAG', 5, 1),   # 3' TGCA overhang
    'ClaI':   ('ATCGAT', 2, 4),   # 5' CG overhang
    'AvaI':   ('CYCGRG', 1, 5),   # degenerate; C[CT]CG[AG]G
}

def complement(base: str) -> str:
    return {'A':'T','T':'A','G':'C','C':'G','N':'N',
            'R':'Y','Y':'R','W':'W','S':'S','M':'K','K':'M'}.get(base, 'N')

def reverse_complement(seq: str) -> str:
    return ''.join(complement(b) for b in reversed(seq.upper()))

def iupac_to_regex(seq: str) -> str:
    iupac = {'A':'A','T':'T','G':'G','C':'C','N':'[ATGC]',
              'R':'[AG]','Y':'[CT]','W':'[AT]','S':'[GC]',
              'M':'[AC]','K':'[GT]','B':'[CGT]','D':'[AGT]',
              'H':'[ACT]','V':'[ACG]'}
    return ''.join(iupac.get(b, b) for b in seq.upper())

print('Restriction enzyme database loaded.')
print(f'Enzymes available: {len(RESTRICTION_ENZYMES)}')
for name, (site, ct, cb) in list(RESTRICTION_ENZYMES.items())[:5]:
    overhang_len = abs(ct - cb)
    overhang_type = '5-prime' if ct < cb else ('3-prime' if ct > cb else 'blunt')
    print(f'  {name:10s}  site={site}  cut_top={ct}  cut_bot={cb}  '
          f'overhang={overhang_len}bp {overhang_type}')
```python

### 1.2 Single-Enzyme Digest: Finding Cut Sites and Computing Fragments

```python
def find_cut_positions(dna: str, enzyme_name: str) -> list[int]:
    """Return top-strand cut positions for all recognition sites of enzyme in dna."""
    site, cut_top, _ = RESTRICTION_ENZYMES[enzyme_name]
    pattern = iupac_to_regex(site)
    positions = []
    for m in re.finditer(f'(?={pattern})', dna.upper()):
        positions.append(m.start() + cut_top)
    return sorted(positions)

def digest(dna: str, enzyme_name: str, circular: bool = False) -> list[str]:
    """Return list of fragment sequences after digestion."""
    cuts = find_cut_positions(dna, enzyme_name)
    if not cuts:
        return [dna]
    if circular:
        # rotate so first cut is at position 0
        first = cuts[0]
        rotated = dna[first:] + dna[:first]
        adjusted = [c - first for c in cuts[1:]] + [len(dna)]
        fragments = []
        prev = 0
        for c in adjusted:
            fragments.append(rotated[prev:c])
            prev = c
        return fragments
    else:
        boundaries = [0] + cuts + [len(dna)]
        return [dna[boundaries[i]:boundaries[i+1]] for i in range(len(boundaries)-1)]

# Demo: digest a synthetic plasmid-like sequence
plasmid_seq = (
    'AAGCTTGCATGCCTGCAGGTCGACGGATCCCCGGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGAC'
    'CTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATT'
    'CCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAAT'
    'TGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCG'
    'CGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGG'
    'CTGCGGCGAGCGGTATCAGGCTACGGCGTAGCGCACGACCGGGCGAGCTCCTCAAGCTTGTGAATTCGGATCC'
)

for enzyme in ['EcoRI', 'BamHI', 'HindIII']:
    cuts = find_cut_positions(plasmid_seq, enzyme)
    frags = digest(plasmid_seq, enzyme)
    print(f'{enzyme}: {len(cuts)} cut(s) at positions {cuts}')
    print(f'  Fragments: {[len(f) for f in frags]} bp')
```python

### 1.3 Multi-Enzyme Double Digest

```python
def double_digest(dna: str, enzyme1: str, enzyme2: str, circular: bool = False) -> list[str]:
    """Digest dna with two enzymes simultaneously."""
    cuts1 = find_cut_positions(dna, enzyme1)
    cuts2 = find_cut_positions(dna, enzyme2)
    all_cuts = sorted(set(cuts1 + cuts2))
    if not all_cuts:
        return [dna]
    if circular:
        first = all_cuts[0]
        rotated = dna[first:] + dna[:first]
        adjusted = [c - first for c in all_cuts[1:]] + [len(dna)]
        frags, prev = [], 0
        for c in adjusted:
            frags.append(rotated[prev:c])
            prev = c
        return frags
    boundaries = [0] + all_cuts + [len(dna)]
    return [dna[boundaries[i]:boundaries[i+1]] for i in range(len(boundaries)-1)]

frags_single = digest(plasmid_seq, 'EcoRI')
frags_double = double_digest(plasmid_seq, 'EcoRI', 'BamHI')

print('Single digest (EcoRI):')
print(f'  {sorted([len(f) for f in frags_single], reverse=True)} bp')
print('Double digest (EcoRI + BamHI):')
print(f'  {sorted([len(f) for f in frags_double], reverse=True)} bp')
```python

### 1.4 Simulating Gel Electrophoresis

In agarose gel electrophoresis, DNA fragments migrate based on size. Migration distance is approximately proportional to the **log** of fragment length. We simulate a lane by plotting bands at log-scale positions.

```python
def plot_gel(lanes: dict[str, list[str]], ladder_sizes: list[int] | None = None):
    """Simulate an agarose gel image. lanes maps lane label to list of fragment sequences."""
    if ladder_sizes is None:
        ladder_sizes = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250, 100]

    all_lanes = {'Ladder': [str(s) for s in ladder_sizes]}
    all_lanes.update(lanes)
    n_lanes = len(all_lanes)

    fig, ax = plt.subplots(figsize=(2.5 * n_lanes, 5))
    ax.set_facecolor('#f5f0e8')   # agarose colour
    ax.set_ylim(1.8, 4.1)
    ax.set_xlim(0, n_lanes)
    ax.set_yticks([])
    ax.set_xticks([])

    for lane_idx, (label, frags) in enumerate(all_lanes.items()):
        x_center = lane_idx + 0.5
        sizes = [len(f) if not f.isdigit() else int(f) for f in frags]
        for size in sizes:
            if size < 50:
                continue
            y = math.log10(size)
            lw = 3 if label == 'Ladder' else 4
            color = '#444444' if label == 'Ladder' else '#e05c1a'
            ax.plot([x_center - 0.35, x_center + 0.35], [y, y],
                    color=color, linewidth=lw, solid_capstyle='butt')
            if label == 'Ladder':
                ax.text(x_center - 0.38, y, f'{size}', ha='right', va='center', fontsize=7)
            else:
                ax.text(x_center, y + 0.04, f'{size}bp', ha='center', va='bottom', fontsize=7)
        ax.text(x_center, 1.85, label, ha='center', va='top', fontsize=9, fontweight='bold')

    ax.set_title('Simulated Agarose Gel', fontsize=12)
    plt.tight_layout()
    plt.show()

lanes = {
    'EcoRI':        digest(plasmid_seq, 'EcoRI'),
    'BamHI':        digest(plasmid_seq, 'BamHI'),
    'EcoRI+BamHI':  double_digest(plasmid_seq, 'EcoRI', 'BamHI'),
}
plot_gel(lanes)
```python

### 1.5 Compatible Ends: Which Enzymes Can Ligate Together?

Two restriction fragments can be ligated if their overhangs are complementary.
For 5' overhangs the single-stranded tail of one end must equal the reverse complement of the other.
Some enzymes produce **identical** overhangs even though their recognition sequences differ (e.g. BamHI/BglII both give 5'-GATC).

```python
def compute_overhang(enzyme_name: str) -> tuple[str, str]:
    """Return (overhang_sequence, overhang_type) for an enzyme.
    overhang_type is '5prime', '3prime', or 'blunt'.
    Overhang sequence is on the top strand for 5' overhangs,
    bottom strand (5'->3') for 3' overhangs.
    """
    site, cut_top, cut_bot = RESTRICTION_ENZYMES[enzyme_name]
    if cut_top == cut_bot:
        return ('', 'blunt')
    if cut_top < cut_bot:
        # 5' overhang: single-stranded portion is site[cut_top:cut_bot]
        return (site[cut_top:cut_bot], '5prime')
    else:
        # 3' overhang: single-stranded on bottom strand
        overhang = reverse_complement(site)[cut_bot: cut_top]
        return (overhang, '3prime')

def are_compatible(enzyme1: str, enzyme2: str) -> bool:
    """Return True if the ends produced by enzyme1 and enzyme2 can be ligated."""
    oh1, type1 = compute_overhang(enzyme1)
    oh2, type2 = compute_overhang(enzyme2)
    if type1 != type2:
        return False
    if type1 == 'blunt':
        return True  # all blunt ends are compatible
    return oh1 == oh2  # same overhang sequence = compatible

print('Overhang summary:')
for name in RESTRICTION_ENZYMES:
    oh, otype = compute_overhang(name)
    print(f'  {name:10s}  {otype:8s}  5\'-{oh}-3\'')

print('\nCompatibility pairs:')
enzymes = list(RESTRICTION_ENZYMES.keys())
for e1, e2 in itertools.combinations(enzymes, 2):
    if are_compatible(e1, e2) and e1 != e2:
        oh, _ = compute_overhang(e1)
        print(f'  {e1} + {e2}  (shared overhang: {oh})')
```python

## 2. Primer Design

### 2.1 PCR Fundamentals

A PCR reaction uses two primers:
- **Forward primer**: matches the top strand, pointing right
- **Reverse primer**: matches the bottom strand, pointing left (i.e. its sequence is the reverse complement of the bottom strand)

Primer quality depends on several properties:

| Property | Recommended range |
|----------|-------------------|
| Length   | 18–25 bp |
| GC content | 40–60% |
| Tm       | 55–65 °C (matched between pair ±2 °C) |
| 3' end   | No runs of ≥4 identical bases |
| Self-complementarity | Low |

### 2.2 Tm Calculation: Three Models

```python
def gc_content(seq: str) -> float:
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq)

def tm_basic(primer: str) -> float:
    """4+2 rule: Tm = 4*(G+C) + 2*(A+T). Only reliable for <20bp."""
    s = primer.upper()
    return 4 * (s.count('G') + s.count('C')) + 2 * (s.count('A') + s.count('T'))

def tm_salt_adjusted(primer: str, salt_mm: float = 50.0) -> float:
    """Salt-adjusted Tm (Wallace rule variant).
    Tm = 81.5 + 16.6*log10([Na+]) + 41*(GC) - 675/N
    salt_mm: [Na+] in mM.
    """
    n = len(primer)
    gc = gc_content(primer)
    return 81.5 + 16.6 * math.log10(salt_mm / 1000) + 41 * gc - 675 / n

# Nearest-neighbor thermodynamic parameters (SantaLucia 1998)
# dH in kcal/mol, dS in cal/mol/K
NN_PARAMS: dict[str, tuple[float, float]] = {
    'AA': (-7.9, -22.2), 'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3), 'CA': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'CT': (-7.8, -21.0), 'GA': (-8.2, -22.2), 'CG': (-10.6, -27.2),
    'GC': (-9.8, -24.4), 'GG': (-8.0, -19.9), 'AC': (-7.8, -21.0), 'TC': (-8.2, -22.2),
    'AG': (-7.8, -21.0), 'TG': (-8.5, -22.7), 'TT': (-7.9, -22.2), 'CC': (-8.0, -19.9),
}
# Initiation parameters
NN_INIT_GC = (0.1, -2.8)   # dH, dS for GC terminal
NN_INIT_AT = (2.3, 4.1)    # dH, dS for AT terminal

def tm_nearest_neighbor(primer: str, dna_conc_nm: float = 250.0, salt_mm: float = 50.0) -> float:
    """Nearest-neighbor Tm using SantaLucia 1998 parameters.
    dna_conc_nm: total primer concentration in nM (assumes non-self-complementary).
    salt_mm: [Na+] in mM.
    """
    seq = primer.upper()
    R = 1.987  # cal/mol/K
    dH, dS = 0.0, 0.0
    # Initiation
    for end_base in (seq[0], seq[-1]):
        if end_base in 'GC':
            dH += NN_INIT_GC[0]; dS += NN_INIT_GC[1]
        else:
            dH += NN_INIT_AT[0]; dS += NN_INIT_AT[1]
    # Propagation
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        h, s = NN_PARAMS.get(dinuc, (-8.0, -20.0))
        dH += h; dS += s
    # Salt correction (Owczarzy 2004 approximation)
    dS += 0.368 * (len(seq) - 1) * math.log(salt_mm / 1000)
    ct = dna_conc_nm * 1e-9  # mol/L
    tm_kelvin = (dH * 1000) / (dS + R * math.log(ct / 4)) - 273.15
    return tm_kelvin

# Compare all three models on example primers
test_primers = [
    ('gfp_fwd',  'ATGGTGAGCAAGGGCGAGGAG'),
    ('gfp_rev',  'TTACTTGTACAGCTCGTCCATGCC'),
    ('short',    'GCATGCAAGCTT'),
    ('long_at',  'AATAATAATAATAATAATAATAAT'),
]
print(f'{"Name":<12} {"Seq":<26} {"Tm_basic":>8} {"Tm_salt":>8} {"Tm_NN":>8}')
for name, seq in test_primers:
    print(f'{name:<12} {seq:<26} {tm_basic(seq):>8.1f} '
          f'{tm_salt_adjusted(seq):>8.1f} {tm_nearest_neighbor(seq):>8.1f}')
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
