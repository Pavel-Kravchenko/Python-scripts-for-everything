---
name: bio-core-nucleic-acid-structure
description: "- Describe the chemical components of nucleotides (base, sugar, phosphate) - Distinguish DNA helical forms (A, B, Z) and their biological significance - Explain major and minor groove geometry and its"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/08_Nucleic_Acid_Structure/01_nucleic_acid_structure.ipynb"
---

# Nucleic Acid Structure

*Source: Course notebook `Tier_2_Core_Bioinformatics/08_Nucleic_Acid_Structure/01_nucleic_acid_structure.ipynb`*

# Nucleic Acid Structure

---

## Learning Objectives

By the end of this notebook, you will be able to:

- Describe the chemical components of nucleotides (base, sugar, phosphate)
- Distinguish DNA helical forms (A, B, Z) and their biological significance
- Explain major and minor groove geometry and its role in protein-DNA recognition
- Understand RNA secondary structure elements (stems, loops, bulges, junctions)
- Conceptually understand RNA secondary structure prediction via minimum free energy
- Analyze DNA-protein interactions: transcription factors, nucleosomes
- Parse nucleic acid structures from PDB files with BioPython
- Use tools for nucleic acid visualization and analysis

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **DNA strand polarity matters everywhere**: DNA is always synthesized 5'→3'. The two strands in a duplex run antiparallel — if one strand runs 5'→3' left to right, the complementary strand runs 3'→5' left to right (or equivalently, 5'→3' right to left). When searching for motifs, always search both strands.
- **B-DNA vs. A-DNA vs. Z-DNA**: B-DNA is the canonical right-handed helix under physiological conditions (~10.5 bp/turn, 3.4 Å rise/bp). A-DNA forms in dehydrated conditions (~11 bp/turn, 2.6 Å rise/bp, wider and shorter). Z-DNA is left-handed (~12 bp/turn) and forms in alternating purine-pyrimidine sequences under high salt or torsional stress. Most genomic DNA is B-form.
- **RNA secondary structure complexity**: Unlike DNA, RNA is predominantly single-stranded and folds into complex secondary and tertiary structures via intramolecular base pairing. The same RNA sequence can fold multiple ways; minimum free energy (MFE) structure is only one possibility. Use ViennaRNA or RNAfold for MFE prediction; consider suboptimal structures for regulatory RNAs.
- **Dot-bracket notation**: In dot-bracket notation, `(` and `)` mark paired bases (the `(` at position i pairs with the `)` at position j), and `.` marks unpaired bases. Pseudoknots require `[`, `]`, `{`, `}` brackets and many tools cannot represent them.
- **CpG islands**: Regions with CpG O/E ratio > 0.6, GC content > 55%, and length > 200 bp are often promoters of housekeeping genes. Most CpG dinucleotides in mammalian genomes are methylated and have been depleted by deamination (5mC → T mutation).

## Environment check (run this first)

```python
# Environment check
from Bio.Seq import Seq
from Bio.PDB import PDBParser, PDBList
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

print("Imports ready.")

# Quick demo: base pairing and complement
dna = Seq("ATGCGAATTCGATCG")
print(f"Sequence (5'->3'):     {dna}")
print(f"Complement (3'->5'):   {dna.complement()}")
print(f"Rev. complement (5'->3'): {dna.reverse_complement()}")
print()

# Helix parameters at a glance
helix_params = {
    'B-DNA': {'bp_per_turn': 10.5, 'rise_A': 3.4, 'handedness': 'Right', 'diameter_A': 20},
    'A-DNA': {'bp_per_turn': 11.0, 'rise_A': 2.6, 'handedness': 'Right', 'diameter_A': 23},
    'Z-DNA': {'bp_per_turn': 12.0, 'rise_A': 3.7, 'handedness': 'Left',  'diameter_A': 18},
}
print(f"{'Form':<8} {'bp/turn':<10} {'Rise(Å)':<10} {'Diameter(Å)':<14} {'Handedness'}")
print("-" * 55)
for form, p in helix_params.items():
    print(f"{form:<8} {p['bp_per_turn']:<10} {p['rise_A']:<10} {p['diameter_A']:<14} {p['handedness']}")
print("\nProceed to Section 1.")
```

---

## 3. DNA Double Helix: A, B, and Z Forms

DNA can adopt different helical conformations depending on sequence, hydration, and salt concentration:

```
A-DNA                   B-DNA                    Z-DNA
(dehydrated)            (physiological)          (high salt, alt. pur-pyr)

  /\                      /\                       /\
 /  \                    /  \                     /  \
| ** |  wider,          |    |  the standard     | zz |  zig-zag
 \ **/  shorter         |    |  Watson-Crick      \zz/   backbone
  |  |                  |    |  helix              |  |
  | **|                 |    |                     |zz|
  \** /                  \  /                      \  /
   \/                    \/                        \/

Right-handed            Right-handed              Left-handed
```

| Parameter | A-DNA | B-DNA | Z-DNA |
|-----------|-------|-------|-------|
| Handedness | Right | Right | **Left** |
| Base pairs / turn | 11 | 10.5 | 12 |
| Rise / base pair | 2.6 A | 3.4 A | 3.7 A |
| Helix diameter | 23 A | 20 A | 18 A |
| Major groove | Narrow, deep | Wide, deep | Flat |
| Minor groove | Wide, shallow | Narrow, deep | Narrow, deep |
| Sugar pucker | C3'-endo | C2'-endo | Alternating |
| Conditions | Dehydrated, RNA-DNA hybrids | Physiological | High salt, alternating purine-pyrimidine |

**B-DNA** is the standard physiological form. **A-form** is adopted by RNA duplexes and DNA-RNA hybrids (the 2'-OH of RNA forces C3'-endo sugar pucker). **Z-DNA** forms in regions with alternating purine-pyrimidine sequences under high salt conditions; it may play roles in gene regulation.

```python
# DNA helix parameters and dimension calculations
HELIX_PARAMS = {
    'A-DNA': {'bp_per_turn': 11,   'rise': 2.6, 'diameter': 23, 'handedness': 'Right'},
    'B-DNA': {'bp_per_turn': 10.5, 'rise': 3.4, 'diameter': 20, 'handedness': 'Right'},
    'Z-DNA': {'bp_per_turn': 12,   'rise': 3.7, 'diameter': 18, 'handedness': 'Left'},
}


def helix_dimensions(num_bp, form='B-DNA'):
    """
    Calculate helix dimensions for a given number of base pairs.
    
    Returns dict with length (A and nm), number of turns, diameter.
    """
    p = HELIX_PARAMS[form]
    length_A = num_bp * p['rise']
    return {
        'form': form,
        'base_pairs': num_bp,
        'length_A': length_A,
        'length_nm': length_A / 10.0,
        'turns': num_bp / p['bp_per_turn'],
        'diameter_A': p['diameter']
    }


# Compare 100 bp in different forms
print("100 bp DNA in different helical forms:\n")
print(f"{'Form':<8s} {'Length (nm)':>12s} {'Turns':>8s} {'Diameter (A)':>13s}")
print("-" * 45)
for form in HELIX_PARAMS:
    d = helix_dimensions(100, form)
    print(f"{d['form']:<8s} {d['length_nm']:12.1f} {d['turns']:8.1f} {d['diameter_A']:13.0f}")

# Human genome scale
print("\n--- Human Genome Scale ---")
genome_bp = 3.2e9  # ~3.2 billion base pairs
d = helix_dimensions(genome_bp, 'B-DNA')
print(f"Human genome ({genome_bp/1e9:.1f} Gbp) as B-DNA:")
print(f"  Length: {d['length_nm']/1e9:.2f} m = {d['length_nm']/1e9*100:.0f} cm")
print(f"  Turns:  {d['turns']/1e6:.0f} million")
print(f"  (That is ~2 meters of DNA packed into a ~6 um nucleus!)")
```

---

## 4. Major and Minor Grooves

The double helix has two grooves of unequal size, created because the glycosidic bonds of a base pair are not diametrically opposite:

```
                  MAJOR GROOVE (~12 A wide in B-DNA)
                /                          \
   5' --------/  Rich in chemical groups    \-------- 3'
      \      |   for protein recognition    |      /
       \     |                              |     /
        \----+--------- base pairs ---------+----/
       /     |                              |     \
      /      |   Fewer distinguishing       |      \
   3' --------\  chemical features          /-------- 5'
                \                          /
                  MINOR GROOVE (~6 A wide in B-DNA)


Groove accessibility in different forms:

         Major groove    Minor groove
A-DNA    Narrow, deep    Wide, shallow
B-DNA    Wide, deep      Narrow, deep
Z-DNA    Flat            Narrow, deep
```

### Why Grooves Matter

The **major groove** exposes a unique pattern of hydrogen bond donors and acceptors for each base pair, allowing proteins to "read" the DNA sequence without unwinding the helix:

| Base Pair | Major Groove Pattern | Minor Groove Pattern |
|-----------|---------------------|---------------------|
| A-T | H-bond acceptor, no group, H-bond acceptor, H, CH3 | H-bond acceptor, no group |
| T-A | CH3, H, H-bond acceptor, no group, H-bond acceptor | no group, H-bond acceptor |
| G-C | H-bond acceptor, H-bond acceptor, H-bond donor, no group | H-bond acceptor, H-bond donor |
| C-G | no group, H-bond donor, H-bond acceptor, H-bond acceptor | H-bond donor, H-bond acceptor |

All four base pairs are distinguishable in the **major groove** but A-T vs T-A (and G-C vs C-G) are hard to tell apart in the minor groove. This is why most transcription factors bind in the major groove.

### Base Stacking

In addition to hydrogen bonding between bases, **stacking interactions** between adjacent base pairs are a major stabilizing force. These arise from:

- Van der Waals forces between the flat aromatic rings
- Dipole-dipole interactions
- Hydrophobic effect (excluding water from the helix interior)

Stacking contributes more to duplex stability than hydrogen bonding. The stacking energy depends on the identity of adjacent base pairs (nearest-neighbor model).

```python
# Nearest-neighbor stacking free energies (kcal/mol) at 37C, 1M NaCl
# From SantaLucia (1998)
NN_ENERGIES = {
    'AA': -1.0, 'AT': -0.88, 'AG': -1.28, 'AC': -1.44,
    'TA': -0.58, 'TT': -1.0,  'TG': -1.45, 'TC': -1.30,
    'GA': -1.30, 'GT': -1.44, 'GG': -1.84, 'GC': -2.24,
    'CA': -1.45, 'CT': -1.28, 'CG': -2.17, 'CC': -1.84,
}


def nearest_neighbor_dG(seq):
    """
    Estimate duplex free energy using the nearest-neighbor model.
    
    More accurate than simple base-pair counting.
    Initiation penalty: +1.96 kcal/mol (for both ends)
    """
    seq = seq.upper()
    dG = 1.96  # Initiation penalty
    
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        dG += NN_ENERGIES.get(dinuc, -1.0)
    
    return dG


# Compare stacking energies for different sequences
test_seqs = [
    'GCGCGCGC',   # Alternating GC
    'GGGGCCCC',   # GC block
    'AAAATTTT',   # AT block
    'ATATATAT',   # Alternating AT
    'GATCGATC',   # Mixed
]

print(f"{'Sequence':<12s} {'GC%':>5s} {'dG (kcal/mol)':>14s} {'Tm (C)':>7s}")
print("-" * 42)
for seq in test_seqs:
    dG = nearest_neighbor_dG(seq)
    tm = melting_temperature(seq)
    print(f"{seq:<12s} {gc_content(seq):5.0f} {dG:14.2f} {tm:7.0f}")
```

---

## 5. RNA Structure

Unlike DNA, RNA is typically single-stranded and folds back on itself to form complex secondary and tertiary structures.

### RNA Secondary Structure Elements

```
1. STEM (paired region)        2. HAIRPIN LOOP (stem-loop)

   G - C                              U  C  A
   A - U                             U      G
   C - G                            G--------C
   U - A                            |        |
   G - C                            C--------G
                                    |        |
                                    A--------U
                                    5'      3'

3. BULGE (unpaired on one side) 4. INTERNAL LOOP (unpaired on both sides)

   G - C                            G - C
   A - U                            A   U   <-- mismatches
   C   |                            G   A
   U - A  <-- C is unpaired         C - G
   G - C                            A - U

5. MULTI-STEM JUNCTION           6. PSEUDOKNOT

        G-C                         5'-A-U---G-C--\
       / | \                              |    |   |
      /  |  \                        3'-U-A---C-G--/
   G-C  A-U  C-G                    (pairs with downstream region)
   | |  | |  | |
   C-G  U-A  G-C
```

### Dot-Bracket Notation

RNA secondary structure is commonly represented in dot-bracket format:
- `.` = unpaired nucleotide
- `(` = 5' partner of a base pair
- `)` = 3' partner of a base pair

Example hairpin: `((((....))))` -- 4-bp stem with 4-nt loop
