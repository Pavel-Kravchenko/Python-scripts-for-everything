---
name: ai-science-alphafold-protein-design
description: "*Prerequisites: Tier 2 Module 07 (Protein Structure), Module T5-01 (LLM Fine-tuning)*"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/04_AlphaFold_Protein_Design/04_AlphaFold_Protein_Design.ipynb"
---

# Module T5-04: AlphaFold & Protein Design

*Source: Course notebook `Tier_5_Modern_AI_for_Science/04_AlphaFold_Protein_Design/04_AlphaFold_Protein_Design.ipynb`*

# Module T5-04: AlphaFold & Protein Design

**Tier 5 — Modern AI for Science | Module 04**

*Prerequisites: Tier 2 Module 07 (Protein Structure), Module T5-01 (LLM Fine-tuning)*

---

**By the end of this notebook you will be able to:**
1. Explain AF2 architecture at a practical level
2. Interpret pLDDT and PAE confidence metrics
3. Retrieve and inspect predicted structures from AlphaFold DB
4. Compare structures with simple RMSD and confidence-aware triage
5. Understand RFdiffusion + ProteinMPNN design loops and AF2/AF3 validation

## Why this notebook matters

AlphaFold2 was a paradigm shift for structural biology — it solved the protein folding problem well enough to predict accurate structures for essentially all known proteins. The AlphaFold Protein Structure Database now contains over 200 million predicted structures, covering most of UniProt. Understanding what AF2 actually does (and does not do), how to read its confidence metrics, and how it connects to protein design workflows is essential for modern computational biology.

## How to work through this notebook

1. Read the architecture section (Section 1) carefully — understanding the Evoformer is essential context for interpreting outputs.
2. The pLDDT and PAE sections (Section 2) are the most practically important — most AF2 users should start here.
3. The AlphaFold DB fetching code (Section 3) requires internet access; all other cells are offline.
4. The design workflow (Section 5) is conceptual — no GPU or model download required.

## Common sticking points

- **AF2 does not predict dynamics**: the structure is a single static conformation — typically the most stable one. Disorder, conformational flexibility, and alternate states are not modeled.
- **High pLDDT ≠ biologically correct**: disordered regions can predict with low pLDDT. But some disordered regions also predict high pLDDT if they resemble structured proteins in the training data — always check experimental context.
- **PAE is directional**: PAE[i,j] asks "if residue i is placed correctly, how confident are we about residue j?" It is not symmetric. Low cross-domain PAE means the domain orientation is reliable.
- **AF3 vs AF2**: AF2 is optimized for single-protein monomers (and homomers). AF3 predicts complexes including protein–DNA, protein–RNA, protein–ligand, and modified residues. Use AF3 for any complex-involving question.

## 1. AlphaFold2 Architecture

AlphaFold2 has three main stages: **input preparation**, the **Evoformer trunk**, and the **Structure Module**.

### 1a. Input Representations

AF2 builds two complementary representations of the protein:

**MSA representation** (multiple sequence alignment):
- Search UniRef90 and MGnify databases for homologous sequences
- Stack them into an MSA matrix: rows = sequences, columns = positions
- Shape: (N_seq × L × c_m) where c_m = 256 channels
- Captures evolutionary co-variation: positions that co-mutate are likely in contact

**Pair representation**:
- Tracks pairwise relationships between every pair of residues i, j
- Shape: (L × L × c_z) where c_z = 128 channels
- Initialized from relative positional encodings + MSA outer products

### 1b. Evoformer Block (48 iterations)

The Evoformer is a cross-attention stack that simultaneously updates both representations. Each block contains:

1. **Row-wise gated self-attention with column-wise bias** (MSA → MSA):
   Updates each row (sequence) attending to other rows, biased by the pair representation. This lets evolutionary co-variation signals propagate.

2. **Column-wise gated self-attention** (MSA → MSA):
   Updates each position by attending across the MSA column. This aggregates information from all homologs at the same position.

3. **MSA transition** (feed-forward on MSA representation)

4. **Outer product mean** (MSA → pair update):
   Computes outer products of MSA column embeddings and averages over sequences. This transfers co-evolutionary signal into the pair representation.

5. **Triangular multiplicative updates** (pair → pair):
   Two operations: "outgoing" and "incoming" edges. For residue pair (i,j), the outgoing update gathers information from all paths i→k combined with k→j. This enforces geometric consistency (triangle inequality in 3D space).

6. **Triangular attention** (pair → pair):
   Self-attention on the pair representation using starting/ending node attention.

7. **Pair transition** (feed-forward on pair representation)

The Evoformer runs 48 times. By the end, the pair representation encodes implicit distance and orientation constraints between every residue pair.

### 1c. Structure Module (8 iterations)

The Structure Module converts the pair representation into 3D coordinates. It uses **Invariant Point Attention (IPA)**:

- Each residue has a rigid frame (rotation + translation) representing its local coordinate system
- IPA computes attention using both: global scalar attention on the pair representation, AND local 3D point attention within each residue's frame
- The attention is equivariant to global rotations/translations — moving the whole protein does not change the attention pattern
- Each IPA iteration refines the rigid frame for each residue
- After 8 iterations, backbone frames are decoded into full-atom coordinates

**Key outputs:**
- 3D coordinates (Cα and all-atom)
- Per-residue pLDDT confidence (predicted as a classification head during training)
- Distogram for predicted distance distribution (used during training)

### 1d. AlphaFold3 Changes

AF3 replaces the Structure Module with a **Diffusion Module**:
- Generates all-atom coordinates for proteins, nucleic acids, small molecules, and ions in a single model
- Uses a diffusion process over atom positions conditioned on Evoformer-derived features
- The Evoformer in AF3 is adapted to handle non-protein entities via a token-based representation where each token can be a residue or a small molecule atom group

```python
import numpy as np
import matplotlib.pyplot as plt
import json
from Bio.PDB import PDBParser, Superimposer

np.random.seed(5)
```

## 2. Confidence Metrics: pLDDT and PAE

- **pLDDT** (0..100): per-residue local confidence
  - >90 high, 70–90 useful, <50 often disordered
- **PAE** (angstroms): pairwise relative placement uncertainty
  - low inter-domain PAE means more reliable domain orientation

```python
# synthetic example for visualization workflow
L = 180
plddt = np.clip(np.r_[np.random.normal(88, 4, 120), np.random.normal(58, 8, 60)], 0, 100)
pae = np.random.uniform(2, 18, size=(L, L))
pae = (pae + pae.T) / 2

fig, ax = plt.subplots(1, 2, figsize=(12, 4))
ax[0].plot(plddt, lw=1.2)
ax[0].axhline(90, ls='--', c='g', lw=1)
ax[0].axhline(70, ls='--', c='orange', lw=1)
ax[0].set_title('pLDDT profile')
ax[0].set_xlabel('Residue')
ax[0].set_ylabel('pLDDT')
im = ax[1].imshow(pae, cmap='viridis', vmin=0, vmax=30)
ax[1].set_title('PAE heatmap')
plt.colorbar(im, ax=ax[1], fraction=0.046)
plt.tight_layout()
```

## 3. Fetching Structures from AlphaFold DB

Use UniProt IDs to fetch precomputed structures. This avoids rerunning heavy inference for many exploratory tasks.

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

# Example (uncomment to run online)
# pdb_path = fetch_alphafold_pdb('P04637')  # TP53
# print(pdb_path)
```

## 4. Structure Comparison (RMSD)

When reference structures exist, compare predicted vs reference Cα coordinates.

```python
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

# Example usage:
# print(ca_rmsd('predicted.pdb', 'reference.pdb'))
```

## 5. Protein Design Workflow (Conceptual)

A common loop:
1. **RFdiffusion** proposes backbones
2. **ProteinMPNN** designs sequences for those backbones
3. **AF2/AF3** validates whether designed sequence recovers intended fold/interface

Use confidence filters to reduce false positives before experimental testing.

```python
def design_priority(mean_plddt: float, interface_pae: float, seq_novelty: float) -> float:
    conf = 0.7 * (mean_plddt / 100.0) + 0.3 * (1.0 - np.clip(interface_pae / 30.0, 0, 1))
    return float(0.8 * conf + 0.2 * seq_novelty)

candidates = [
    {'id': 'd1', 'plddt': 86, 'pae': 7.5, 'novelty': 0.40},
    {'id': 'd2', 'plddt': 74, 'pae': 15.0, 'novelty': 0.70},
    {'id': 'd3', 'plddt': 91, 'pae': 5.0, 'novelty': 0.20},
]
for c in candidates:
    c['priority'] = design_priority(c['plddt'], c['pae'], c['novelty'])
print(sorted(candidates, key=lambda x: x['priority'], reverse=True))
```

## Summary

- AF2 outputs are most useful when combined with confidence metrics.
- AlphaFold DB is a fast entry point for structure-driven hypothesis generation.
- Design pipelines should include confidence-aware ranking before experiments.

## Source-backed Context

- AlphaFold2 documentation and paper position AF2 as a high-accuracy monomer predictor with confidence metrics central to interpretation.
- AlphaFold3 introduces broader biomolecular-complex prediction scope including proteins, nucleic acids, ions, and ligands.
- RoseTTAFold remains an official open implementation for structure/interactions with a 3-track architecture.

## Validated Sources

Checked online during content expansion.

- [Jumper et al. 2021, AlphaFold2 (Nature)](https://www.nature.com/articles/s41586-021-03819-2)
- [Abramson et al. 2024, AlphaFold3 (Nature)](https://www.nature.com/articles/s41586-024-07487-w)
- [AlphaFold official repository](https://github.com/google-deepmind/alphafold)
- [AlphaFold3 official repository](https://github.com/google-deepmind/alphafold3)
- [RoseTTAFold official repository](https://github.com/RosettaCommons/RoseTTAFold)
