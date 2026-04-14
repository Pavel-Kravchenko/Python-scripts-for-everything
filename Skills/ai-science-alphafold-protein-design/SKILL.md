---
name: ai-science-alphafold-protein-design
description: "AlphaFold2/3 architecture, confidence metrics (pLDDT/PAE), structure fetching, and protein design ranking workflows"
tool_type: python
primary_tool: NumPy
---

# AlphaFold and Protein Design

## Key Gotchas

- **AF2 does not predict dynamics** — single static conformation, typically most stable
- **High pLDDT ≠ biologically correct** — disordered regions can predict high pLDDT if they resemble structured training data
- **PAE is directional**: PAE[i,j] asks "if residue i is placed correctly, how confident about j?" — not symmetric
- **AF2 vs AF3**: AF2 is for single proteins/homomers. AF3 handles protein-DNA, protein-RNA, protein-ligand, modified residues

## AF2 Architecture

**Inputs**: MSA representation (N_seq × L × 256) + Pair representation (L × L × 128)

**Evoformer** (48 blocks): Row-wise gated self-attention → column-wise attention → outer product mean (MSA→pair) → triangular multiplicative updates (pair→pair) → triangular attention

**Structure Module** (8 iterations): Invariant Point Attention (IPA) — attention using both scalar pair features AND local 3D point attention within each residue's rigid frame. Equivariant to global rotations/translations.

**Outputs**: 3D coordinates (all-atom), per-residue pLDDT, distogram

**AF3 changes**: Replaces Structure Module with a **Diffusion Module** generating all-atom coordinates for proteins, nucleic acids, small molecules, and ions.

## Confidence Metrics

- **pLDDT** (0–100): per-residue local confidence. >90 high, 70–90 useful, <50 often disordered
- **PAE** (angstroms): pairwise relative placement uncertainty. Low inter-domain PAE = reliable domain orientation

```python
import numpy as np
from Bio.PDB import PDBParser, Superimposer

def ca_rmsd(pdb_a: str, pdb_b: str) -> float:
    parser = PDBParser(QUIET=True)
    s1 = parser.get_structure('a', pdb_a)
    s2 = parser.get_structure('b', pdb_b)
    ca1 = [a for a in s1.get_atoms() if a.name == 'CA']
    ca2 = [a for a in s2.get_atoms() if a.name == 'CA']
    n = min(len(ca1), len(ca2))
    if n == 0:
        raise ValueError('No CA atoms found')
    sup = Superimposer()
    sup.set_atoms(ca1[:n], ca2[:n])
    return float(sup.rms)
```

## Fetching from AlphaFold DB

```python
import requests
from pathlib import Path

def fetch_alphafold_pdb(uniprot_id: str, out_dir: str = 'pdb_files'):
    out = Path(out_dir)
    out.mkdir(exist_ok=True)
    url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        raise RuntimeError(f'Fetch failed: {uniprot_id} ({r.status_code})')
    path = out / f'{uniprot_id}.pdb'
    path.write_text(r.text)
    return str(path)
```

## Design Priority Ranking

Pipeline: RFdiffusion (backbones) → ProteinMPNN (sequences) → AF2/AF3 (validation)

```python
def design_priority(mean_plddt: float, interface_pae: float, seq_novelty: float) -> float:
    conf = 0.7 * (mean_plddt / 100.0) + 0.3 * (1.0 - np.clip(interface_pae / 30.0, 0, 1))
    return float(0.8 * conf + 0.2 * seq_novelty)
```

## Sources

- [Jumper et al. 2021, AlphaFold2 (Nature)](https://www.nature.com/articles/s41586-021-03819-2)
- [Abramson et al. 2024, AlphaFold3 (Nature)](https://www.nature.com/articles/s41586-024-07487-w)
- [AlphaFold repo](https://github.com/google-deepmind/alphafold) / [AF3 repo](https://github.com/google-deepmind/alphafold3)
