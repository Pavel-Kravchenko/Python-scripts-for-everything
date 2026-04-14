---
name: structural-bioinformatics
description: Bio.PDB parsing, distance/RMSD/TM-score, DSSP, Ramachandran, PWM/PROSITE, GO enrichment, KEGG REST API.
tool_type: python
primary_tool: NumPy
---

## When to Use
- Parsing PDB files, computing distances/dihedrals/RMSD (Bio.PDB)
- Secondary structure assignment (DSSP), Ramachandran analysis
- Building and scanning PWMs; scanning with PROSITE patterns
- Pfam/InterPro domain search (hmmscan)
- GO/KEGG enrichment from a gene list

## Quick Reference Tables

### Protein Structure Levels
| Level | Stabilized by |
|-------|--------------|
| Primary | Peptide bonds (N→C sequence) |
| Secondary | Backbone H-bonds (α-helix i→i+4, β-sheet) |
| Tertiary | Hydrophobic core, H-bonds, disulfides (single chain 3D) |
| Quaternary | Same forces, multiple subunits |

### Secondary Structure Geometry
| Element | Phi/Psi | Rise/res |
|---------|---------|----------|
| α-helix | -57/-47° | 1.5 Å |
| 3₁₀-helix | -49/-26° | 2.0 Å |
| β-strand (antiparallel) | -139/+135° | 3.4 Å |

### DSSP Codes
`H`=α-helix `G`=3₁₀ `I`=π `E`=β-strand `B`=β-bridge `T`=turn `S`=bend `-`=coil

### DNA Helical Forms
| | A-DNA | B-DNA | Z-DNA |
|---|-------|-------|-------|
| Handedness | Right | Right | Left |
| bp/turn | 11 | 10.5 | 12 |
| Rise/bp | 2.6 Å | 3.4 Å | 3.7 Å |
| Major groove | Narrow/deep | Wide/deep | Flat |
| Conditions | Dehydrated/RNA hybrids | Physiological | High salt, alt pur-pyr |

### PROSITE Syntax
| Syntax | Meaning |
|--------|---------|
| `x` | Any AA |
| `[LIVMF]` | One of listed |
| `{P}` | Not proline |
| `x(3)` | Repeat 3× |
| `x(2,4)` | 2–4 repeats |

N-glycosylation: `N-{P}-[ST]-{P}`

### GO Evidence Hierarchy
- **Experimental** (EXP, IDA, IMP, IGI): highest quality
- **Computational** (ISS, ISO, IBA): medium
- **Automatic** (IEA): lowest — exclude from stringent analyses

## Key Patterns

### Bio.PDB SMCRA Hierarchy
```python
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
structure = parser.get_structure('name', 'file.pdb')

# Navigate: Structure → Model → Chain → Residue → Atom
for residue in structure[0]['A']:
    if residue.id[0] != ' ':   # skip HETATM (ligands/water)
        continue
    if 'CA' in residue:
        coords = residue['CA'].get_vector().get_array()  # np array [x,y,z]
        bfactor = residue['CA'].get_bfactor()
```

Residue ID is a tuple `(hetflag, resseq, icode)` — standard AA has hetflag `' '`.

### Distances, RMSD, TM-score
```python
import numpy as np

def distance(a, b):
    return np.linalg.norm(np.array(a) - np.array(b))

def rmsd(c1, c2):
    diff = np.array(c1) - np.array(c2)
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

def tm_score(c1, c2):
    L = len(c1)
    d0 = max(1.24*(L-15)**(1/3) - 1.8, 0.5)
    d = np.sqrt(np.sum((np.array(c1)-np.array(c2))**2, axis=1))
    return np.sum(1/(1+(d/d0)**2)) / L
# TM > 0.5 → same fold; < 0.3 → different fold
```

### Kabsch Superposition
```python
def kabsch(mobile, ref):
    mc, rc = mobile - mobile.mean(0), ref - ref.mean(0)
    U, S, Vt = np.linalg.svd(mc.T @ rc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1,1,d]) @ U.T
    return (mc @ R) + ref.mean(0)

# Or use Bio.PDB built-in:
from Bio.PDB import Superimposer
sup = Superimposer()
sup.set_atoms(fixed_atoms, moving_atoms)
sup.apply(moving_structure.get_atoms())
print(sup.rms)
```

### DSSP
```python
from Bio.PDB.DSSP import DSSP
dssp = DSSP(structure[0], 'file.pdb', dssp='mkdssp')
for key in dssp.keys():
    chain_id, res_id = key
    ss  = dssp[key][2]   # H/E/G/I/B/T/S/-
    acc = dssp[key][3]   # solvent accessibility Å²
    phi = dssp[key][4]
    psi = dssp[key][5]
```

### Build and Scan PWM
```python
BASES = ['A', 'C', 'G', 'T']

def build_pwm(seqs, pseudocount=0.1):
    pfm = np.zeros((4, len(seqs[0])))
    for seq in seqs:
        for i, b in enumerate(seq.upper()):
            if b in BASES: pfm[BASES.index(b), i] += 1
    ppm = (pfm + pseudocount) / (len(seqs) + 4*pseudocount)
    return np.log2(ppm / 0.25)   # log-odds PWM

def scan_pwm(pwm, sequence, threshold=None):
    L = pwm.shape[1]
    thresh = threshold or 0.6 * np.sum(np.max(pwm, axis=0))
    hits = []
    for i in range(len(sequence)-L+1):
        s = sum(pwm[BASES.index(b), j]
                for j, b in enumerate(sequence[i:i+L].upper()) if b in BASES)
        if s >= thresh:
            hits.append((i, sequence[i:i+L], s))
    return sorted(hits, key=lambda x: -x[2])
```

### PROSITE → Regex
```python
import re

def prosite_to_regex(pattern):
    parts = []
    for elem in pattern.strip('.').split('-'):
        m = re.match(r'^(.+?)\((\d+)(?:,(\d+))?\)$', elem)
        core, low, high = (m.group(1), m.group(2), m.group(3)) if m else (elem, None, None)
        r = ('.' if core=='x' else core if core.startswith('[')
             else f'[^{core[1:-1]}]' if core.startswith('{') else core)
        if low: r += f'{{{low},{high}}}' if high else f'{{{low}}}'
        parts.append(r)
    return ''.join(parts)
```

### GO Hypergeometric Enrichment
```python
from scipy import stats

def go_enrichment(gene_list, term_to_genes, N=20000):
    n = len(set(gene_list))
    results = []
    for term, tgenes in term_to_genes.items():
        K, k = len(tgenes), len(set(gene_list) & tgenes)
        if k == 0: continue
        p = stats.hypergeom.sf(k-1, N, K, n)
        results.append({'term': term, 'k': k, 'K': K, 'p': p})
    results.sort(key=lambda r: r['p'])
    for i, r in enumerate(results):
        r['fdr'] = min(r['p'] * len(results) / (i+1), 1.0)
    return results
```

### KEGG REST
```python
# GET https://rest.kegg.jp/find/pathway/apoptosis
# GET https://rest.kegg.jp/get/hsa04210
# Organism codes: hsa=human, mmu=mouse, sce=yeast, eco=E.coli
import urllib.request
def kegg_get(op, *args):
    url = 'https://rest.kegg.jp/' + '/'.join([op]+list(args))
    with urllib.request.urlopen(url, timeout=15) as r:
        return r.read().decode()
```

## Pitfalls

- **PDB is fixed-width, not whitespace-delimited.** Use `line[30:38]` for X coordinate, not `line.split()`.
- **Residue ID is a tuple, not an int.** `structure[0]['A'][10]` is shorthand for `(' ', 10, ' ')`. Ligands have a non-space hetflag.
- **NMR structures have multiple models.** Use `structure[0]` for first model; iterate over `structure` for ensemble analysis.
- **DSSP requires the external `mkdssp` binary.** Install via `conda install -c salilab dssp`.
- **RMSD without superposition is meaningless.** Always Kabsch-align first.
- **PWM zero probabilities → -inf log-odds.** Always add pseudocount (α=0.1) before log.
- **GO true path rule must be applied before enrichment.** Propagate each annotation to all ancestor terms first.
- **Multiple testing in GO/pathway.** Use BH FDR, not Bonferroni — terms are correlated.
- **IEA annotations are auto-assigned and lower quality.** Filter out for experimental conclusions.
- **PROSITE `{P}` is a negative class, not a quantifier.** Translates to `[^P]` in regex.
- **hmmscan E-value < 0.01 = significant; > 1 = noise.** Use `--domtblout` for domain-level results.

## Related Skills
- `bio-core-gene-ontology` — GO DAG navigation and term propagation
- `bio-pathway-analysis-kegg-pathways` — KEGG enrichment workflows
- `biopython-databases` — BioPython SeqIO, PDB retrieval
