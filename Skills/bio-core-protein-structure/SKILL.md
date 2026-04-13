---
name: bio-core-protein-structure
description: "- Describe the four levels of protein structure (primary through quaternary) - Parse PDB files and understand ATOM record format - Use BioPython's Bio.PDB module to navigate the Structure/Model/Chain/"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/07_Protein_Structure/01_protein_structure.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Protein Structure Analysis

*Source: Course notebook `Tier_2_Core_Bioinformatics/07_Protein_Structure/01_protein_structure.ipynb`*

# Protein Structure Analysis

---

## Learning Objectives

By the end of this notebook, you will be able to:

- Describe the four levels of protein structure (primary through quaternary)
- Parse PDB files and understand ATOM record format
- Use BioPython's Bio.PDB module to navigate the Structure/Model/Chain/Residue/Atom hierarchy
- Calculate interatomic distances, bond angles, and RMSD
- Assign and analyze secondary structure with DSSP
- Perform structural alignment and superposition (Kabsch algorithm)
- Understand molecular surfaces (van der Waals, solvent-accessible, electrostatic)
- Generate and interpret Ramachandran plots
- Use visualization tools: PyMOL commands and nglview in Jupyter
- Appreciate modern structure prediction with AlphaFold

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **PDB format uses fixed-width columns, not whitespace**: The ATOM record format is strictly column-delimited. Name (cols 13-16), residue name (cols 18-20), chain (col 22), residue number (cols 23-26), x/y/z (cols 31-54), occupancy (cols 55-60), B-factor (cols 61-66). Never use `.split()` to parse PDB ATOM records — use column slicing or BioPython.
- **The SMCRA hierarchy**: BioPython organizes PDB structures as Structure → Model → Chain → Residue → Atom. For X-ray crystal structures, there is usually one model (index 0). NMR structures have multiple models (each model is one conformer). Residue IDs are tuples `(hetflag, resseq, icode)`: `(' ', 10, ' ')` for standard residues, `('H_NAG', 1, ' ')` for heteroatoms like ligands.
- **Resolution means lower number is better**: 1.5 Å resolution is better than 3.0 Å. At <2.0 Å, individual atoms and water molecules are clearly resolved. At 3.0-4.0 Å, only the backbone trace is reliable; side chain positions are approximate.
- **B-factor vs. PLDDT**: In experimental structures, high B-factor (>60 Å²) indicates flexible/disordered regions. In AlphaFold2 structures, the B-factor column stores pLDDT (per-residue confidence, 0-100); pLDDT>90 is very confident, 50-70 is low confidence.
- **RMSD is sensitive to outliers and alignment**: When comparing structures, RMSD depends on which atoms you include and how you align them. Always superimpose using Cα atoms only (backbone), and report the number of aligned residues alongside the RMSD value.

## Environment check (run this first)

```python
# Environment check
from Bio.PDB import PDBParser, PDBList, Superimposer
from Bio.PDB.DSSP import DSSP
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import os

print("Bio.PDB modules loaded.")
print("\nPDB coordinate format reminder:")
print("ATOM  seqno  name  resname  chain  resseq  x       y       z      occ   bfac  elem")
print("ATOM    145  CA    ALA      A      20      12.045  15.321  9.862  1.00  15.30  C")
print()
print("BioPython SMCRA hierarchy:")
print("  Structure")
print("    -> Model (0 for X-ray; 0,1,2... for NMR)")
print("      -> Chain (e.g., 'A', 'B')")
print("        -> Residue (id = (hetflag, resseq, icode))")
print("          -> Atom (name = 'CA', 'N', 'O', 'CB', ...)")
print("\nProceed to Section 1.")
```

---

## 3. Parsing PDB with BioPython (Bio.PDB)

BioPython provides a powerful, object-oriented PDB parser built around the **SMCRA hierarchy**:

```
Structure
  |-- Model (0, 1, ...)
        |-- Chain ('A', 'B', ...)
              |-- Residue (('', 10, ''), 'ALA')
                    |-- Atom ('CA', 'N', 'C', 'O', ...)

S = Structure   One PDB entry (may have multiple models for NMR)
M = Model       One set of coordinates (X-ray: 1 model; NMR: ~20 models)
C = Chain       One polypeptide or nucleic acid chain
R = Residue     One amino acid or nucleotide
A = Atom        One atom with (x, y, z) coordinates
```

Each level is iterable, so you can use nested loops or list comprehensions to traverse the hierarchy.

```python
from Bio.PDB import PDBParser, PDBList
import warnings
warnings.filterwarnings('ignore')  # Suppress PDB parsing warnings for cleaner output

# Download a structure from RCSB PDB
# Using crambin (1CRN) -- a small, well-resolved plant protein (46 residues)
pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file('1CRN', pdir='pdb_files', file_format='pdb')
print(f"Downloaded: {pdb_file}")
```

```python
# Parse the downloaded structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure('crambin', pdb_file)

# Navigate the SMCRA hierarchy
print("=== Structure Hierarchy ===")
print(f"Structure ID: {structure.id}")

for model in structure:
    print(f"\n  Model {model.id}")
    for chain in model:
        residues = list(chain.get_residues())
        # Filter standard amino acids (exclude water, ligands)
        aa_residues = [r for r in residues if r.id[0] == ' ']
        het_residues = [r for r in residues if r.id[0] != ' ']
        print(f"    Chain {chain.id}: {len(aa_residues)} amino acids, {len(het_residues)} het groups")
        
        # Show first 5 residues
        print(f"    First 5 residues:")
        for res in aa_residues[:5]:
            atoms = list(res.get_atoms())
            print(f"      {res.get_resname()} {res.id[1]}: {len(atoms)} atoms")
```

```python
# Accessing individual atoms and their properties
model = structure[0]
chain = model['A']

# Get a specific residue (residue 10 in chain A)
# Residue ID is a tuple: (hetflag, resseq, insertion_code)
residue = chain[(' ', 10, ' ')]
# Shorthand (works when hetflag=' ' and icode=' '):
residue = chain[10]

print(f"Residue: {residue.get_resname()} {residue.id[1]}")
print(f"\nAtoms in this residue:")
for atom in residue:
    coord = atom.get_vector()
    print(f"  {atom.get_name():4s}  element={atom.element:2s}  "
          f"coords=({coord[0]:7.3f}, {coord[1]:7.3f}, {coord[2]:7.3f})  "
          f"B-factor={atom.get_bfactor():.2f}")
```

```python
# Extract all CA atoms (backbone trace)
import numpy as np

ca_atoms = []
for residue in chain:
    if residue.id[0] != ' ':  # Skip heteroatoms
        continue
    if 'CA' in residue:
        ca = residue['CA']
        ca_atoms.append({
            'res_name': residue.get_resname(),
            'res_num': residue.id[1],
            'coord': ca.get_vector().get_array()
        })

print(f"Extracted {len(ca_atoms)} CA atoms")
print(f"\nSequence from CA trace:")
# Three-letter to one-letter conversion
AA_MAP = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}
sequence = ''.join(AA_MAP.get(ca['res_name'], 'X') for ca in ca_atoms)
print(sequence)
```

---

## 4. Calculating Distances, Angles, and RMSD

Structural analysis frequently requires measuring geometric properties:

- **Distance**: Euclidean distance between two atoms
- **Angle**: Angle formed by three atoms
- **Dihedral angle**: Torsion angle defined by four atoms
- **RMSD**: Root Mean Square Deviation between two coordinate sets

```
RMSD = sqrt( (1/N) * sum( |r_i - r_i'|^2 ) )

Interpretation:
  0 - 1 A    Nearly identical structures
  1 - 2 A    Very similar (same fold)
  2 - 3 A    Similar fold, some variation
  3 - 5 A    Same topology, different details
  > 5 A      Different structures
```

```python
import numpy as np

def distance(coord1, coord2):
    """Euclidean distance between two 3D points."""
    return np.linalg.norm(np.array(coord1) - np.array(coord2))


def angle(coord1, coord2, coord3):
    """Angle (in degrees) at coord2 formed by coord1-coord2-coord3."""
    v1 = np.array(coord1) - np.array(coord2)
    v2 = np.array(coord3) - np.array(coord2)
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    return np.degrees(np.arccos(cos_angle))


def dihedral(p1, p2, p3, p4):
    """Dihedral (torsion) angle in degrees defined by four points."""
    p1, p2, p3, p4 = [np.array(p) for p in (p1, p2, p3, p4)]
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2)))


def calculate_rmsd(coords1, coords2):
    """RMSD between two Nx3 coordinate arrays."""
    coords1, coords2 = np.array(coords1), np.array(coords2)
    if coords1.shape != coords2.shape:
        raise ValueError(f"Shape mismatch: {coords1.shape} vs {coords2.shape}")
    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))


# Demonstrate with crambin residues
model = structure[0]
chain = model['A']

# Distance between CA atoms of residues 1 and 2
ca1 = chain[1]['CA'].get_vector().get_array()
ca2 = chain[2]['CA'].get_vector().get_array()
print(f"CA-CA distance (res 1--2): {distance(ca1, ca2):.2f} A")

# N-CA-C bond angle for residue 5
res5 = chain[5]
n_coord = res5['N'].get_vector().get_array()
ca_coord = res5['CA'].get_vector().get_array()
c_coord = res5['C'].get_vector().get_array()
print(f"N-CA-C angle (res 5):     {angle(n_coord, ca_coord, c_coord):.1f} degrees")

# Peptide bond dihedral (omega) between residues 5 and 6
res6 = chain[6]
omega = dihedral(
    res5['CA'].get_vector().get_array(),
    res5['C'].get_vector().get_array(),
    res6['N'].get_vector().get_array(),
    res6['CA'].get_vector().get_array()
)
print(f"Omega angle (res 5--6):   {omega:.1f} degrees (expected ~180 for trans)")
```

```python
# BioPython also provides built-in distance calculation via the minus operator
atom1 = chain[1]['CA']
atom2 = chain[2]['CA']

# The - operator on Atom objects returns distance in Angstroms
d = atom1 - atom2
print(f"BioPython atom distance: {d:.2f} A")

# Contact map: find all CA pairs within 8 Angstroms
ca_list = []
for residue in chain:
    if residue.id[0] == ' ' and 'CA' in residue:
        ca_list.append(residue['CA'])

contacts = []
for i in range(len(ca_list)):
    for j in range(i + 4, len(ca_list)):  # Skip nearby in sequence
        d = ca_list[i] - ca_list[j]
        if d < 8.0:
            contacts.append((i + 1, j + 1, d))

print(f"\nCA-CA contacts (< 8 A, |i-j| >= 4): {len(contacts)}")
print("First 10 contacts:")
for res_i, res_j, dist in contacts[:10]:
    print(f"  Res {res_i:3d} -- Res {res_j:3d}: {dist:.2f} A")
```

---

## 5. Secondary Structure Assignment with DSSP

DSSP (Define Secondary Structure of Proteins) is the gold-standard algorithm for assigning secondary structure from 3D coordinates. It uses hydrogen bond energies calculated from N-H...O=C geometry.

### DSSP Codes

| Code | Structure | Description |
|------|-----------|-------------|
| H | Alpha-helix | 4-turn helix (i -> i+4 H-bonds) |
| B | Beta-bridge | Isolated beta-bridge |
| E | Beta-strand | Extended strand in beta-sheet |
| G | 3_10-helix | 3-turn helix (i -> i+3) |
| I | Pi-helix | 5-turn helix (i -> i+5) |
| T | Turn | H-bonded turn |
| S | Bend | High curvature region |
| - | Coil | No regular structure |

BioPython wraps the DSSP program (must be installed separately via `conda install -c salilab dssp` or `apt install dssp`).

```python
from Bio.PDB.DSSP import DSSP

# Run DSSP on our structure
# Requires the 'dssp' or 'mkdssp' executable to be installed
try:
    dssp = DSSP(structure[0], pdb_file, dssp='mkdssp')
    
    print(f"DSSP assigned {len(dssp)} residues\n")
    
    # DSSP returns: (dssp_index, amino_acid, ss, acc, phi, psi, ...)
    ss_sequence = ''
    for key in dssp.keys():
        residue_dssp = dssp[key]
        ss = residue_dssp[2]  # Secondary structure code
        ss_sequence += ss
    
    print(f"Secondary structure assignment:")
    print(f"  {ss_sequence}")
    
    # Count secondary structure composition
    total = len(ss_sequence)
    helix = ss_sequence.count('H') + ss_sequence.count('G') + ss_sequence.count('I')
    sheet = ss_sequence.count('E') + ss_sequence.count('B')
    coil = total - helix - sheet
    
    print(f"\nSecondary structure composition:")
    print(f"  Helix (H+G+I): {helix:3d} ({100*helix/total:.1f}%)")
    print(f"  Sheet (E+B):   {sheet:3d} ({100*sheet/total:.1f}%)")
    print(f"  Other:         {coil:3d} ({100*coil/total:.1f}%)")

except Exception as e:
    print(f"DSSP not available: {e}")
    print("Install with: conda install -c salilab dssp")
    print("\nFalling back to PDB HELIX/SHEET records...")
    print("(PDB files contain pre-calculated secondary structure in HELIX and SHEET records)")
```
