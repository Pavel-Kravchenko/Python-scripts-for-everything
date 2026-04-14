---
name: bio-core-protein-structure
description: Protein structure analysis with BioPython Bio.PDB — SMCRA hierarchy, distance/RMSD calculations, DSSP secondary structure assignment.
tool_type: python
primary_tool: NumPy
---

# Protein Structure Analysis

## Pitfalls

- **PDB uses fixed-width columns, not whitespace**: ATOM records are column-delimited (name cols 13-16, chain col 22, resseq 23-26, x/y/z 31-54). Never use `.split()` to parse ATOM lines — use column slicing or BioPython.
- **SMCRA hierarchy**: Structure → Model → Chain → Residue → Atom. X-ray = one model (index 0); NMR = multiple models (each is a conformer). Residue IDs are tuples `(hetflag, resseq, icode)`: `(' ', 10, ' ')` for standard residues, `('H_NAG', 1, ' ')` for heteroatoms.
- **Resolution**: lower number = better. <2.0 Å resolves atoms and waters; 3.0–4.0 Å only backbone is reliable; side chains are approximate.
- **B-factor vs pLDDT**: In experimental structures, B-factor >60 Å² = flexible/disordered. In AlphaFold2, the B-factor column stores pLDDT (0–100); pLDDT >90 is very confident, 50–70 is low confidence.
- **RMSD alignment**: RMSD depends on which atoms and alignment method. Always superimpose on Cα only and report the number of aligned residues alongside the value.

## RMSD Interpretation

| RMSD (Å) | Meaning |
|---|---|
| 0–1 | Nearly identical |
| 1–2 | Very similar (same fold) |
| 2–3 | Similar fold, some variation |
| 3–5 | Same topology, different details |
| >5 | Different structures |

## DSSP Secondary Structure Codes

| Code | Structure |
|---|---|
| H | Alpha-helix (i→i+4 H-bonds) |
| E | Beta-strand |
| G | 3₁₀-helix (i→i+3) |
| I | Pi-helix (i→i+5) |
| B | Isolated beta-bridge |
| T | H-bonded turn |
| S | Bend |
| - | Coil |

Install: `conda install -c salilab dssp`

## Key Patterns

### Parse PDB and navigate SMCRA
```python
from Bio.PDB import PDBParser, PDBList
import warnings
warnings.filterwarnings('ignore')

pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file('1CRN', pdir='pdb_files', file_format='pdb')

parser = PDBParser(QUIET=True)
structure = parser.get_structure('crambin', pdb_file)

model = structure[0]
chain = model['A']

# Residue access — shorthand works when hetflag='' and icode=''
residue = chain[10]          # same as chain[(' ', 10, ' ')]

# Filter standard AA residues (skip water/heteroatoms)
aa_residues = [r for r in chain if r.id[0] == ' ']

# Extract CA coords
ca_coords = []
for res in chain:
    if res.id[0] != ' ': continue
    if 'CA' in res:
        ca_coords.append(res['CA'].get_vector().get_array())
```

### Distance, angle, dihedral, RMSD
```python
import numpy as np

def distance(c1, c2):
    return np.linalg.norm(np.array(c1) - np.array(c2))

def angle(c1, c2, c3):
    v1 = np.array(c1) - np.array(c2)
    v2 = np.array(c3) - np.array(c2)
    cos_a = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1)
    return np.degrees(np.arccos(cos_a))

def dihedral(p1, p2, p3, p4):
    p1, p2, p3, p4 = [np.array(p) for p in (p1, p2, p3, p4)]
    b1, b2, b3 = p2-p1, p3-p2, p4-p3
    n1 = np.cross(b1, b2); n2 = np.cross(b2, b3)
    n1 /= np.linalg.norm(n1); n2 /= np.linalg.norm(n2)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2)))

def calculate_rmsd(coords1, coords2):
    coords1, coords2 = np.array(coords1), np.array(coords2)
    if coords1.shape != coords2.shape:
        raise ValueError(f"Shape mismatch: {coords1.shape} vs {coords2.shape}")
    return np.sqrt(np.mean(np.sum((coords1 - coords2)**2, axis=1)))

# BioPython shorthand: atom1 - atom2 returns distance in Å
d = chain[1]['CA'] - chain[2]['CA']

# Contact map: CA pairs within 8 Å
ca_list = [res['CA'] for res in chain if res.id[0] == ' ' and 'CA' in res]
contacts = [(i+1, j+1, ca_list[i] - ca_list[j])
            for i in range(len(ca_list))
            for j in range(i+4, len(ca_list))
            if ca_list[i] - ca_list[j] < 8.0]
```

### DSSP secondary structure
```python
from Bio.PDB.DSSP import DSSP

dssp = DSSP(structure[0], pdb_file, dssp='mkdssp')
ss_sequence = ''.join(dssp[key][2] for key in dssp.keys())

helix = sum(ss_sequence.count(c) for c in 'HGI')
sheet = sum(ss_sequence.count(c) for c in 'EB')
coil  = len(ss_sequence) - helix - sheet
```

### Three-letter to one-letter AA map
```python
AA_MAP = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
    'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
    'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
    'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
}
```
