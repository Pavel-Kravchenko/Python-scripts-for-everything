---
name: alphafold-structure-prediction
description: AlphaFold/ESMFold structure prediction and confidence interpretation.
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


## When to Use

Use this atomic skill for focused work on **alphafold-structure-prediction** without bundling unrelated topics.

## Quick Reference

This skill was split from `alphafold-protein-design.md` to keep topics independent and self-contained.

## Core Patterns

Use the parent material below as the source reference, then keep implementations specific to this topic.

## Source Reference (from merged skill)

---
name: alphafold-protein-design
description: AlphaFold2 structure prediction, pLDDT/PAE interpretation, ColabFold, ESMFold, RFdiffusion backbone design, ProteinMPNN inverse folding, structure quality evaluation
---

# AlphaFold & Protein Design

## When to Use

Use this skill when:
- Predicting protein 3D structure from sequence
- Interpreting AlphaFold2 confidence metrics (pLDDT, PAE)
- Designing novel protein binders or de novo proteins
- Comparing predicted vs experimental structures
- Querying the AlphaFold DB for precomputed structures

## Quick Reference

| Task | Tool | How |
|------|------|-----|
| Structure prediction | ColabFold | Colab notebook or `colabfold_batch` CLI |
| Fast prediction (no MSA) | ESMFold | REST API or `esm.pretrained.esmfold_v1()` |
| Protein complex prediction | AF2-Multimer / ColabFold | `--templates --num-recycle 3` |
| Fetch from AlphaFold DB | REST API | `GET https://alphafold.ebi.ac.uk/files/AF-{UNIPROT}-F1-model_v4.pdb` |
| Structural comparison | TM-align | `TMalign pred.pdb ref.pdb` |
| RMSD | Biopython | `Superimposer().set_atoms()` |
| Backbone design | RFdiffusion | `run_inference.py` + contig maps |
| Sequence design | ProteinMPNN | `protein_mpnn_run.py` |
| Visualize 3D | py3Dmol / nglview | `py3Dmol.view()` |

## Key Patterns

**Pattern 1: Interpret pLDDT and PAE**
```python
import json
import numpy as np
import matplotlib.pyplot as plt

with open('result_model_1.json') as f:
    data = json.load(f)

plddt = data['plddt']     # per-residue confidence (0-100)
pae = np.array(data['predicted_aligned_error'])  # N×N matrix

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
ax1.plot(plddt, color='steelblue')
ax1.axhline(90, color='green', ls='--', label='Very high')
ax1.axhline(70, color='orange', ls='--', label='Confident')
ax1.set(xlabel='Residue', ylabel='pLDDT', title='Per-residue confidence')
ax1.legend()

im = ax2.imshow(pae, cmap='Greens_r', vmin=0, vmax=30)
plt.colorbar(im, ax=ax2, label='Error (Å)')
ax2.set(xlabel='Scored residue', ylabel='Aligned residue', title='PAE')
plt.tight_layout()
```

**Pattern 2: ESMFold API prediction**
```python
import requests

sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY'
response = requests.post(
    'https://api.esmatlas.com/foldSequence/v1/pdb/',
    data=sequence,
    headers={'Content-Type': 'application/x-www-form-urlencoded'}
)
with open('esmfold_pred.pdb', 'w') as f:
    f.write(response.text)
```

**Pattern 3: Download from AlphaFold DB**
```python
import requests

uniprot_id = 'P00533'  # EGFR
url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
pdb_text = requests.get(url).text
with open(f'{uniprot_id}_alphafold.pdb', 'w') as f:
    f.write(pdb_text)
```

**Pattern 4: RMSD with Biopython**
```python
from Bio.PDB import PDBParser, Superimposer

parser = PDBParser(QUIET=True)
s1 = parser.get_structure('pred', 'predicted.pdb')
s2 = parser.get_structure('ref',  'crystal.pdb')

sup = Superimposer()
ca_pred = [a for a in s1.get_atoms() if a.name == 'CA']
ca_ref  = [a for a in s2.get_atoms() if a.name == 'CA']
n = min(len(ca_pred), len(ca_ref))
sup.set_atoms(ca_ref[:n], ca_pred[:n])
print(f'Cα RMSD: {sup.rms:.2f} Å')
```

**Pattern 5: RFdiffusion binder design**
```bash
# Design binders to a target hotspot region
python RFdiffusion/scripts/run_inference.py \
    inference.input_pdb=target.pdb \
    'contigmap.contigs=[A1-100/0 50-100]' \
    'ppi.hotspot_res=[A45,A67,A89]' \
    inference.num_designs=20 \
    inference.output_prefix=designs/binder

# ProteinMPNN: design sequences for the best backbone
python ProteinMPNN/protein_mpnn_run.py \
    --pdb_path designs/binder_0.pdb \
    --out_folder mpnn_seqs/ \
    --num_seq_per_target 8 \
    --sampling_temp 0.1
```

## pLDDT Interpretation

| pLDDT | Confidence | Typical Region |
|-------|-----------|----------------|
| > 90 | Very high | Well-folded structured domain |
| 70–90 | Confident | Structured, minor flexibility |
| 50–70 | Low | Could be disordered or flexible loop |
| < 50 | Very low | Likely intrinsically disordered |

## Common Pitfalls

- **pLDDT ≠ correctness** — high pLDDT only means the model is internally consistent; validate against experimental data
- **PAE for multimers** — PAE between chains indicates confidence in relative domain/chain positioning
- **MSA depth matters** — AF2 performs poorly for orphan sequences with no homologs; ESMFold is better for these
- **Disordered regions** — regions with pLDDT < 50 should not be used for docking or structural analysis
- **ColabFold vs full AF2** — ColabFold uses faster MSA search (MMseqs2); slightly lower accuracy but adequate for most applications
- **RFdiffusion GPU** — requires A100 for reasonable speed; smaller designs (<100 residues) feasible on T4

## Code Templates

### Batch AlphaFold DB Retrieval
```python
import requests
import time
from pathlib import Path

def fetch_alphafold_pdbs(uniprot_ids, out_dir='af_structures', version=4):
    """Download AlphaFold predicted structures for a list of UniProt IDs."""
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    failed = []
    for uid in uniprot_ids:
        url = f'https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{version}.pdb'
        resp = requests.get(url)
        if resp.status_code == 200:
            (out_dir / f'{uid}.pdb').write_text(resp.text)
        else:
            failed.append(uid)
        time.sleep(0.2)  # rate limit
    if failed:
        print(f"Failed: {failed}")
    return out_dir

fetch_alphafold_pdbs(['P00533', 'P04637', 'P01308'])
```

### Extract Disordered Regions (pLDDT < 50)
```python
from Bio.PDB import PDBParser
import numpy as np

def get_disordered_regions(pdb_path, plddt_threshold=50):
    """Return residue ranges with pLDDT < threshold (stored in B-factor column)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('model', pdb_path)
    disordered = []
    run_start = None
    prev_res = None
    for chain in structure.get_chains():
        for residue in chain:
            for atom in residue:
                if atom.name == 'CA':
                    plddt = atom.bfactor  # AF2 stores pLDDT in B-factor
                    res_id = residue.get_id()[1]
                    if plddt < plddt_threshold:
                        if run_start is None:
                            run_start = res_id
                    else:
                        if run_start is not None:
                            disordered.append((run_start, prev_res))
                            run_start = None
                    prev_res = res_id
    if run_start is not None:
        disordered.append((run_start, prev_res))
    return disordered

regions = get_disordered_regions('AF-P00533-F1-model_v4.pdb')
for start, end in regions:
    print(f"Disordered: residues {start}–{end}")
```

### TM-score and RMSD Comparison
```python
from Bio.PDB import PDBParser, Superimposer
import numpy as np

def compare_structures(pred_pdb, ref_pdb):
    """Compute Cα RMSD between predicted and reference structures."""
    parser = PDBParser(QUIET=True)
    pred = parser.get_structure('pred', pred_pdb)
    ref  = parser.get_structure('ref', ref_pdb)

    ca_pred = [a for a in pred.get_atoms() if a.name == 'CA']
    ca_ref  = [a for a in ref.get_atoms()  if a.name == 'CA']
    n = min(len(ca_pred), len(ca_ref))

    sup = Superimposer()
    sup.set_atoms(ca_ref[:n], ca_pred[:n])
    rmsd = sup.rms

    print(f"Aligned residues: {n}")
    print(f"Cα RMSD: {rmsd:.2f} Å")
    print(f"Interpretation: {'Excellent' if rmsd < 1 else 'Good' if rmsd < 2 else 'Moderate' if rmsd < 4 else 'Poor'}")
    return rmsd, sup

rmsd, superimposer = compare_structures('predicted.pdb', 'crystal.pdb')
```

## Related Skills
- `structural-bioinformatics` — PDB format, secondary structure, molecular visualization
- `biopython-databases` — BioPython PDB module, sequence retrieval from UniProt
- `cheminformatics-drug-discovery` — docking predicted structures, binding site analysis
- `ml-deep-learning-bio` — protein language models, ESM embeddings, GNN for proteins


## Related Skills

- `alphafold-structure-prediction` (this file)
- `alphafold-protein-design` (legacy merged skill)
