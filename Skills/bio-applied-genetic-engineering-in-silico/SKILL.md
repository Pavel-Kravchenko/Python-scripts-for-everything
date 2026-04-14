---
name: bio-applied-genetic-engineering-in-silico
description: In silico restriction digestion, compatible end detection, primer design (Tm models), and gel simulation
tool_type: python
primary_tool: Matplotlib
---

# Genetic Engineering In Silico

## Restriction Enzyme Reference

| Overhang type | Example | Note |
|---------------|---------|------|
| 5' sticky | EcoRI (GAATTC, cut 1/5) | Most common for cloning |
| 3' sticky | KpnI (GGTACC, cut 5/1) | Less common |
| Blunt | SmaI (CCCGGG, cut 3/3) | Ligation-inefficient |

Compatible end rule: same overhang sequence = ligatable (e.g. BamHI + BglII both give 5'-GATC).

## Core Functions

```python
import re

RESTRICTION_ENZYMES = {
    # name: (recognition_seq, cut_top, cut_bottom)  # 0-indexed from site start
    'EcoRI':   ('GAATTC', 1, 5),   # 5'-AATT overhang
    'BamHI':   ('GGATCC', 1, 5),   # 5'-GATC
    'HindIII': ('AAGCTT', 1, 5),   # 5'-AGCT
    'SalI':    ('GTCGAC', 1, 5),   # 5'-TCGA
    'XhoI':    ('CTCGAG', 1, 5),   # 5'-TCGA (compatible with SalI → scar!)
    'NotI':    ('GCGGCCGC', 2, 6), # 5'-GGCC (8-cutter, rare)
    'SmaI':    ('CCCGGG', 3, 3),   # blunt
    'KpnI':    ('GGTACC', 5, 1),   # 3' overhang
}

IUPAC = {'A':'A','T':'T','G':'G','C':'C','N':'[ATGC]',
         'R':'[AG]','Y':'[CT]','W':'[AT]','S':'[GC]',
         'M':'[AC]','K':'[GT]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]'}

def iupac_to_regex(seq):
    return ''.join(IUPAC.get(b, b) for b in seq.upper())

def reverse_complement(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N',
            'R':'Y','Y':'R','W':'W','S':'S','M':'K','K':'M'}
    return ''.join(comp.get(b,'N') for b in reversed(seq.upper()))

def find_cut_positions(dna, enzyme_name):
    """Return top-strand cut positions (0-indexed) for all recognition sites."""
    site, cut_top, _ = RESTRICTION_ENZYMES[enzyme_name]
    pattern = iupac_to_regex(site)
    return sorted(m.start() + cut_top for m in re.finditer(f'(?={pattern})', dna.upper()))

def digest(dna, enzyme_name, circular=False):
    """Return list of fragment sequences."""
    cuts = find_cut_positions(dna, enzyme_name)
    if not cuts:
        return [dna]
    if circular:
        first = cuts[0]
        rotated = dna[first:] + dna[:first]
        adjusted = [c - first for c in cuts[1:]] + [len(dna)]
        frags, prev = [], 0
        for c in adjusted:
            frags.append(rotated[prev:c]); prev = c
        return frags
    boundaries = [0] + cuts + [len(dna)]
    return [dna[boundaries[i]:boundaries[i+1]] for i in range(len(boundaries)-1)]

def double_digest(dna, e1, e2, circular=False):
    all_cuts = sorted(set(find_cut_positions(dna, e1) + find_cut_positions(dna, e2)))
    if not all_cuts:
        return [dna]
    boundaries = [0] + all_cuts + [len(dna)]
    return [dna[boundaries[i]:boundaries[i+1]] for i in range(len(boundaries)-1)]

def compute_overhang(enzyme_name):
    """Return (overhang_seq, type) where type is '5prime', '3prime', or 'blunt'."""
    site, ct, cb = RESTRICTION_ENZYMES[enzyme_name]
    if ct == cb:
        return ('', 'blunt')
    if ct < cb:
        return (site[ct:cb], '5prime')
    return (reverse_complement(site)[cb:ct], '3prime')

def are_compatible(e1, e2):
    oh1, t1 = compute_overhang(e1)
    oh2, t2 = compute_overhang(e2)
    if t1 != t2:
        return False
    return True if t1 == 'blunt' else oh1 == oh2
```

## Primer Design

### Recommended Primer Properties

| Property | Target |
|----------|--------|
| Length | 18–25 bp |
| GC content | 40–60% |
| Tm | 55–65°C (pair within ±2°C) |
| 3' end | Avoid ≥4 identical bases |
| Self-complementarity | Minimal |

### Tm Models

```python
import math

def tm_basic(primer):
    """4+2 rule — reliable only for <20 bp."""
    s = primer.upper()
    return 4*(s.count('G')+s.count('C')) + 2*(s.count('A')+s.count('T'))

def tm_salt_adjusted(primer, salt_mm=50.0):
    """Wallace rule variant. salt_mm: [Na+] in mM."""
    n, gc = len(primer), (primer.upper().count('G')+primer.upper().count('C'))/len(primer)
    return 81.5 + 16.6*math.log10(salt_mm/1000) + 41*gc - 675/n

# SantaLucia 1998 nearest-neighbor parameters (dH kcal/mol, dS cal/mol/K)
NN_PARAMS = {
    'AA':(-7.9,-22.2),'AT':(-7.2,-20.4),'TA':(-7.2,-21.3),'CA':(-8.5,-22.7),
    'GT':(-8.4,-22.4),'CT':(-7.8,-21.0),'GA':(-8.2,-22.2),'CG':(-10.6,-27.2),
    'GC':(-9.8,-24.4),'GG':(-8.0,-19.9),'AC':(-7.8,-21.0),'TC':(-8.2,-22.2),
    'AG':(-7.8,-21.0),'TG':(-8.5,-22.7),'TT':(-7.9,-22.2),'CC':(-8.0,-19.9),
}

def tm_nearest_neighbor(primer, dna_conc_nm=250.0, salt_mm=50.0):
    """SantaLucia 1998 Tm. Most accurate for 18-30 bp primers."""
    seq, R = primer.upper(), 1.987
    dH = dS = 0.0
    # Initiation
    for end in (seq[0], seq[-1]):
        if end in 'GC': dH += 0.1;  dS += -2.8
        else:           dH += 2.3;  dS += 4.1
    # Propagation
    for i in range(len(seq)-1):
        h, s = NN_PARAMS.get(seq[i:i+2], (-8.0,-20.0))
        dH += h; dS += s
    # Salt correction
    dS += 0.368*(len(seq)-1)*math.log(salt_mm/1000)
    return (dH*1000)/(dS + R*math.log(dna_conc_nm*1e-9/4)) - 273.15
```

## Gel Electrophoresis Simulation

Migration distance ∝ log10(fragment size). Plot bands at `y = log10(size)`.

## Pitfalls

- **Compatible ends create scars**: SalI/XhoI are compatible but produce a hybrid site — check the ligation scar sequence if it falls in a coding region
- **Isoschizomers vs neoschizomers**: SmaI and XmaI recognize the same site but cut differently (blunt vs 5'-CCGG overhang) — verify cut position, not just recognition sequence
- **Nearest-neighbor Tm is not absolute**: add 3–5°C for actual PCR annealing temp; use gradient PCR to optimize
- **IUPAC degenerate bases**: enzyme databases use IUPAC codes (e.g., AvaI = CYCGRG) — translate to regex before pattern search
- **Circular vs linear digest**: for plasmids, n cuts → n fragments; for linear DNA, n cuts → n+1 fragments
- **Coordinate systems**: BED is 0-based; VCF/GFF are 1-based — off-by-one errors when mixing
