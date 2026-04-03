---
name: structural-bioinformatics
description: PDB parsing, protein/nucleic acid structure, Sanger chromatograms, sequence motifs/PWM, GO enrichment, and KEGG pathways
---

# Structural Bioinformatics, Motifs & Functional Annotation

## When to Use
- Parsing and navigating PDB files with Bio.PDB (Structure→Model→Chain→Residue→Atom)
- Calculating interatomic distances, dihedral angles, RMSD, TM-score
- Assigning secondary structure with DSSP; generating Ramachandran plots
- Analyzing DNA/RNA structure parameters (A/B/Z forms, groove geometry, MFE folding)
- Reading `.ab1` Sanger chromatograms, Phred quality trimming, heterozygous peak detection
- Building PWMs from aligned sites; scanning sequences with PROSITE patterns
- Searching for Pfam/InterPro domains via hmmscan; visualizing domain architectures
- GO enrichment (hypergeometric test, BH FDR) and KEGG pathway enrichment

---

## Quick Reference

### Protein Structure Levels
| Level | Stabilized by | Description |
|-------|--------------|-------------|
| Primary | Peptide bonds | Amino acid sequence (N→C) |
| Secondary | Backbone H-bonds | α-helix (i→i+4), β-sheet, turns |
| Tertiary | Hydrophobic, H-bonds, disulfides | Full 3D fold of one chain |
| Quaternary | Same forces | Multiple subunits |

### Secondary Structure Geometry
| Element | Phi/Psi | Rise/res | H-bond |
|---------|---------|----------|--------|
| α-helix | -57/-47 | 1.5 Å | i→i+4 |
| 3₁₀-helix | -49/-26 | 2.0 Å | i→i+3 |
| β-strand (antiparallel) | -139/+135 | 3.4 Å | between strands |
| β-strand (parallel) | -119/+113 | 3.2 Å | between strands |

### DSSP Secondary Structure Codes
| Code | Meaning |
|------|---------|
| H | α-helix |
| G | 3₁₀-helix |
| I | π-helix |
| E | β-strand (extended) |
| B | Isolated β-bridge |
| T | H-bonded turn |
| S | Bend |
| - | Coil |

### Ramachandran Regions
- α-helix: phi ≈ −57°, psi ≈ −47°
- β-sheet: phi ≈ −120°, psi ≈ +120°
- Left-handed helix: phi ≈ +60°, psi ≈ +60° (Gly only)

### DNA Helical Forms
| Parameter | A-DNA | B-DNA | Z-DNA |
|-----------|-------|-------|-------|
| Handedness | Right | Right | Left |
| bp/turn | 11 | 10.5 | 12 |
| Rise/bp | 2.6 Å | 3.4 Å | 3.7 Å |
| Diameter | 23 Å | 20 Å | 18 Å |
| Major groove | Narrow, deep | Wide, deep | Flat |
| Minor groove | Wide, shallow | Narrow, deep | Narrow, deep |
| Sugar pucker | C3'-endo | C2'-endo | Alternating |
| Conditions | Dehydrated; RNA-DNA hybrids | Physiological | High salt; alternating pur-pyr |

### Phred Quality Scores
| Q | Error prob | Accuracy |
|---|-----------|----------|
| 10 | 1/10 | 90% |
| 20 | 1/100 | 99% |
| 30 | 1/1,000 | 99.9% |
| 40 | 1/10,000 | 99.99% |

Threshold: Q < 20 → unreliable; Q ≥ 30 → high confidence.

### PROSITE Syntax
| Syntax | Meaning |
|--------|---------|
| `C` | Specific amino acid |
| `x` | Any amino acid |
| `[LIVMF]` | One of listed |
| `{P}` | Not proline |
| `x(3)` | Repeat 3× |
| `x(2,4)` | Repeat 2–4× |

Example: N-glycosylation = `N-{P}-[ST]-{P}`

### GO Evidence Code Hierarchy
- **Experimental** (EXP, IDA, IPI, IMP, IGI, IEP): highest quality
- **Computational** (ISS, ISO, IBA): medium
- **Automatic** (IEA): lowest — exclude from stringent analyses

---

## Key Patterns

### Bio.PDB SMCRA Hierarchy
```
Structure[0]         # first Model (X-ray: 1; NMR: ~20)
  ['A']              # Chain by ID
    [10]             # Residue by seq number (shorthand for (' ', 10, ' '))
      ['CA']         # Atom by name
        .get_vector().get_array()   # numpy array [x, y, z]
        .get_bfactor()
        .element
```

Residue ID tuple: `(hetflag, resseq, icode)` — standard AA has hetflag `' '`; skip with `if residue.id[0] != ' '`.

### PDB ATOM record columns (1-indexed, fixed-width)
```
1-6   record type ("ATOM  " / "HETATM")
7-11  serial;  13-16 atom name;  17 alt loc;  18-20 res name
22    chain ID;  23-26 res seq;  31-38 X;  39-46 Y;  47-54 Z
55-60 occupancy;  61-66 B-factor;  77-78 element
```

### PWM Construction (PFM → PPM → log-odds)
```python
pfm[base_idx, pos] += 1                              # count matrix
ppm = (pfm + pseudocount) / (N + 4 * pseudocount)   # add pseudocount α=0.1
pwm = np.log2(ppm / background)                      # log-odds vs 0.25
score = sum(pwm[BASES.index(b), i] for i, b in enumerate(seq))
ic_per_pos = 2.0 - (-sum(p * np.log2(p) for p in ppm[:, pos]))  # bits
```

### GO Hypergeometric Enrichment
```
N = background genes; K = genes annotated to term; n = study list; k = overlap
p = hypergeom.sf(k-1, N, K, n)   # P(X >= k)
```
Apply **true path rule** first: propagate each annotation to all ancestor terms.
Correct for multiple testing with BH FDR (standard: FDR < 0.05).

### KEGG REST API
```
GET https://rest.kegg.jp/find/pathway/apoptosis
GET https://rest.kegg.jp/get/hsa04210
GET https://rest.kegg.jp/link/hsa/pathway:hsa04110
```
Organism codes: `hsa`=human, `mmu`=mouse, `sce`=yeast, `eco`=E.coli.

---

## Code Templates

### Parse and traverse a PDB file
```python
from Bio.PDB import PDBParser, PDBList
import warnings; warnings.filterwarnings('ignore')

pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file('1CRN', pdir='pdb_files', file_format='pdb')
parser = PDBParser(QUIET=True)
structure = parser.get_structure('crambin', pdb_file)

ca_atoms = []
for residue in structure[0]['A']:
    if residue.id[0] == ' ' and 'CA' in residue:
        ca_atoms.append(residue['CA'].get_vector().get_array())
```

### Distances, angles, dihedral, RMSD
```python
def distance(a, b):
    return np.linalg.norm(np.array(a) - np.array(b))

def angle(p1, p2, p3):
    v1, v2 = np.array(p1)-p2, np.array(p3)-p2
    return np.degrees(np.arccos(np.clip(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)), -1, 1)))

def dihedral(p1, p2, p3, p4):
    b1,b2,b3 = p2-p1, p3-p2, p4-p3
    n1 = np.cross(b1,b2); n2 = np.cross(b2,b3)
    n1/=np.linalg.norm(n1); n2/=np.linalg.norm(n2)
    return np.degrees(np.arctan2(np.dot(np.cross(n1,b2/np.linalg.norm(b2)),n2), np.dot(n1,n2)))

def rmsd(c1, c2):
    diff = np.array(c1) - np.array(c2)
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
```

### DSSP secondary structure assignment
```python
from Bio.PDB.DSSP import DSSP
dssp = DSSP(structure[0], pdb_file, dssp='mkdssp')
# dssp[key] = (dssp_idx, aa, ss_code, acc, phi, psi, ...)
for key in dssp.keys():
    chain_id, res_id = key
    ss   = dssp[key][2]   # H/E/G/I/B/T/S/-
    acc  = dssp[key][3]   # solvent accessibility (Å²)
    phi  = dssp[key][4]
    psi  = dssp[key][5]
```

### Ramachandran phi/psi
```python
from Bio.PDB.Polypeptide import PPBuilder
ppb = PPBuilder()
for pp in ppb.build_peptides(structure[0]['A']):
    for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
        if phi and psi:
            print(res.get_resname(), np.degrees(phi), np.degrees(psi))
```

### Kabsch superposition
```python
def kabsch(mobile, ref):
    mc, rc = mobile - mobile.mean(0), ref - ref.mean(0)
    U, S, Vt = np.linalg.svd(mc.T @ rc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1,1,d]) @ U.T
    return (mc @ R) + ref.mean(0)

# Or use built-in:
from Bio.PDB import Superimposer
sup = Superimposer()
sup.set_atoms(fixed_atoms, moving_atoms)   # Atom objects
sup.apply(moving_structure.get_atoms())
print(sup.rms)
```

### TM-score
```python
def tm_score(c1, c2):
    L = len(c1)
    d0 = max(1.24*(L-15)**(1/3) - 1.8, 0.5)
    d = np.sqrt(np.sum((c1-c2)**2, axis=1))
    return np.sum(1/(1+(d/d0)**2)) / L
# > 0.5 = same fold; > 0.3 = possibly same; < 0.3 = different
```

### Read .ab1 chromatogram
```python
from Bio import SeqIO
import numpy as np

record = SeqIO.read('sample.ab1', 'abi')
seq    = str(record.seq)
quals  = record.letter_annotations['phred_quality']

raw = record.annotations['abif_raw']
base_order = raw.get('FWO_1', b'GATC').decode()   # e.g. "GATC"
channels   = {base: np.array(raw[f'DATA{9+i}']) for i, base in enumerate(base_order)}
peak_locs  = list(raw['PLOC1'])                   # scan position per called base
```

### Quality trimming (sliding-window)
```python
def trim_by_quality(quals, min_q=20, window=10):
    n = len(quals)
    start = next((i for i in range(n-window) if np.mean(quals[i:i+window]) >= min_q), 0)
    end   = next((i for i in range(n-1, window, -1) if np.mean(quals[i-window:i]) >= min_q), n)
    return start, end
```

### Build and scan with PWM
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
    max_s = np.sum(np.max(pwm, axis=0))
    thresh = threshold or 0.6 * max_s
    hits = []
    for i in range(len(sequence)-L+1):
        s = sum(pwm[BASES.index(b), j] for j, b in enumerate(sequence[i:i+L].upper()) if b in BASES)
        if s >= thresh:
            hits.append((i, sequence[i:i+L], s))
    return sorted(hits, key=lambda x: -x[2])
```

### PROSITE pattern → regex
```python
import re

def prosite_to_regex(pattern):
    regex_parts = []
    for elem in pattern.strip('.').split('-'):
        m = re.match(r'^(.+?)\((\d+)(?:,(\d+))?\)$', elem)
        core, low, high = (m.group(1), m.group(2), m.group(3)) if m else (elem, None, None)
        r = ('.' if core=='x' else
             core if core.startswith('[') else
             f'[^{core[1:-1]}]' if core.startswith('{') else core)
        if low: r += f'{{{low},{high}}}' if high else f'{{{low}}}'
        regex_parts.append(r)
    return ''.join(regex_parts)

# Example: scan for N-glycosylation sites
regex = prosite_to_regex('N-{P}-[ST]-{P}')
for m in re.finditer(regex, protein_seq):
    print(m.start()+1, m.group())
```

### Parse hmmscan domtblout
```python
def parse_domtblout(text):
    hits = []
    for line in text.split('\n'):
        if line.startswith('#') or not line.strip(): continue
        f = line.split()
        if len(f) >= 23:
            hits.append({'domain': f[0], 'acc': f[1], 'e_value': float(f[6]),
                         'ali_from': int(f[17]), 'ali_to': int(f[18])})
    return sorted(hits, key=lambda h: h['ali_from'])
```

### GO enrichment
```python
from scipy import stats

def go_enrichment(gene_list, term_to_genes, N=20000):
    n = len(set(gene_list))
    results = []
    for term, tgenes in term_to_genes.items():
        K = len(tgenes)
        k = len(set(gene_list) & tgenes)
        if k == 0: continue
        p = stats.hypergeom.sf(k-1, N, K, n)
        results.append({'term': term, 'k': k, 'K': K, 'p': p})
    results.sort(key=lambda r: r['p'])
    m = len(results)
    for i, r in enumerate(results):
        r['fdr'] = min(r['p'] * m / (i+1), 1.0)
    return results
```

### KEGG pathway enrichment
```python
import urllib.request

def kegg_get(operation, *args):
    url = 'https://rest.kegg.jp/' + '/'.join([operation]+list(args))
    with urllib.request.urlopen(url, timeout=15) as r:
        return r.read().decode()

# Pathway enrichment reuses the same hypergeometric logic as GO enrichment
# Replace term_to_genes with {pathway_id: set_of_gene_symbols}
```

### g:Profiler (recommended for production)
```python
from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)
result = gp.profile(organism='hsapiens', query=['TP53','BAX','CASP3'])
# Returns DataFrame with: source, term_id, term_name, p_value, intersection
```

### GOATOOLS (local, offline)
```python
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
obodag = GODag('go-basic.obo')
goe = GOEnrichmentStudy(background_genes, associations, obodag, methods=['fdr_bh'])
results = goe.run_study(study_genes)
```

---

## Common Pitfalls

- **PDB format is fixed-width, not whitespace-delimited.** Use `line[30:38]` not `line.split()` for coordinates.
- **Residue ID is a tuple `(hetflag, resseq, icode)`**, not an integer. Skip ligands/water with `if residue.id[0] != ' '`.
- **NMR structures have multiple models.** Use `structure[0]` for the first model; iterate over all for ensemble analysis.
- **DSSP requires the external `mkdssp` binary.** Install via `conda install -c salilab dssp`.
- **RMSD without superposition is meaningless.** Always Kabsch-align first.
- **PWM zero probabilities cause -inf log-odds.** Always add a pseudocount (α=0.1) before taking log.
- **GO true path rule must be applied before enrichment.** A gene annotated to "G1/S transition" is implicitly annotated to "cell cycle" and all ancestors.
- **Multiple testing in GO/pathway enrichment.** Use BH FDR, not Bonferroni (too conservative for thousands of correlated terms).
- **IEA (electronic) annotations are automatically assigned and lower quality.** Filter to `min_quality >= 3` for experimental/computational evidence.
- **Sanger quality ramps up at 5' end and decays at 3' end.** Always trim before assembly; the first ~20 bp and last ~40 bp are typically low-quality.
- **Double peaks ≠ heterozygosity by default.** Check if the entire trace is noisy (mixed template) vs. single positions (true het).
- **PROSITE `{P}` is a negative class, not a quantifier.** Translates to `[^P]` in Python regex.
- **hmmscan E-value < 0.01 = significant hit; E > 1 = noise.** Use `--domtblout` for domain-level results, not the default output.

---

## Related Skills
- `numpy-pandas-wrangling` — expression matrix operations and PWM scoring
- `rnaseq-metagenomics` — GO/KEGG enrichment is the standard downstream step
- `biopython-databases` — BioPython SeqIO for FASTA/FASTQ preceding structure work
