---
name: alphafold-protein-design
description: AlphaFold2 structure prediction, pLDDT/PAE interpretation, ColabFold, ESMFold, RFdiffusion backbone design, ProteinMPNN inverse folding, structure quality evaluation
---

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
