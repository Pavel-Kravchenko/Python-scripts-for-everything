---
name: bio-applied-virtual-screening
description: "**Tier 3 — Applied Bioinformatics | Module 29 · Notebook 3**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/29_Cheminformatics_Drug_Discovery/03_virtual_screening.ipynb"
primary_tool: Pandas
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Virtual Screening & ADMET Prediction

*Source: Course notebook `Tier_3_Applied_Bioinformatics/29_Cheminformatics_Drug_Discovery/03_virtual_screening.ipynb`*

# Virtual Screening & ADMET Prediction

**Tier 3 — Applied Bioinformatics | Module 29 · Notebook 3**

*Prerequisites: Notebook 2 (QSAR Modeling)*

---

**By the end of this notebook you will be able to:**
1. Explain ligand-based vs structure-based virtual screening strategies
2. Prepare a protein structure for docking (protonation, grid generation)
3. Run AutoDock Vina to dock a small molecule library
4. Filter compounds by ADMET properties using SwissADME and pkCSM
5. Prioritize hits from combined docking + QSAR + ADMET scores



**Key resources:**
- [AutoDock Vina documentation](https://vina.scripps.edu/)
- [SwissADME web tool](http://www.swissadme.ch/)
- [TeachOpenCADD docking talktorial](https://github.com/volkamerlab/teachopencadd)

## 1. Virtual Screening Overview

> Diagram of drug discovery pipeline: target → hit → lead → candidate. Two VS strategies: ligand-based (pharmacophore, similarity, QSAR) vs structure-based (molecular docking). Throughput/accuracy trade-offs.

## 2. Protein Structure Preparation

> Download EGFR crystal structure (PDB: 1IEP). Remove water/ligands. Add hydrogen atoms with Open Babel or PDBFixer. Define docking box around the ATP binding site.

```python
# Example: Download and prepare protein structure
# from Bio.PDB import PDBParser, PDBIO, Select
# import py3Dmol

# !wget -q https://files.rcsb.org/download/1IEP.pdb
# !pdbfixer 1IEP.pdb --output 1IEP_fixed.pdb --add-hydrogens --remove-heterogens

# Example: Visualize structure with py3Dmol
# view = py3Dmol.view()
# view.addModel(open('1IEP_fixed.pdb').read(), 'pdb')
# view.setStyle({'cartoon': {'color': 'spectrum'}})
# view.show()
```

## 3. Molecular Docking with AutoDock Vina

> Prepare ligand SDF to PDBQT with meeko. Set docking box (center_x/y/z, size_x/y/z) from known active site. Run AutoDock Vina. Parse best docking score.

```python
# Example: Prepare and dock ligand
# !mk_prepare_ligand.py -i gefitinib.sdf -o gefitinib.pdbqt

# vina_config = """
# receptor = 1IEP_fixed.pdbqt
# ligand = gefitinib.pdbqt
# center_x = 22.5; center_y = 5.0; center_z = 18.0
# size_x = 20; size_y = 20; size_z = 20
# exhaustiveness = 8
# """
# with open('vina_config.txt', 'w') as f: f.write(vina_config)
# !vina --config vina_config.txt --out gefitinib_docked.pdbqt --log docking.log
# !grep 'REMARK VINA RESULT' gefitinib_docked.pdbqt | head
```

## 4. ADMET Property Prediction

> Compute ADMET properties: Lipinski Ro5, LogS (solubility), Pgp substrate, CYP inhibition, hERG toxicity. Use SwissADME API or DeepChem ADMET models.

```python
# Example: ADMET with DeepChem
# import deepchem as dc
# smiles_list = ['CC1=C(C=C(C=C1)NC2=NC=CC(=N2)N3CCN(CC3)C)NC4=NC=CC(=N4)Cl']  # Imatinib
# featurizer = dc.feat.CircularFingerprint(size=2048)
# X = featurizer.featurize(smiles_list)
# # Example: Use pre-trained Tox21 model for toxicity prediction
```

## 5. Hit Prioritization

> Combine docking score, QSAR pIC50 prediction, and ADMET pass/fail into a composite score. Rank library of 100 compounds. Visualize top 10 structures with their scores.

```python
import pandas as pd

# Example: Composite scoring
# hits = pd.DataFrame({
#     'SMILES': smiles_library,
#     'vina_score': docking_scores,
#     'qsar_pIC50': qsar_predictions,
#     'admet_pass': admet_flags
# })
# hits_filtered = hits[hits['admet_pass']]
# hits_filtered['composite'] = -hits_filtered['vina_score'] + hits_filtered['qsar_pIC50']
# top10 = hits_filtered.nlargest(10, 'composite')
# print(top10[['SMILES', 'vina_score', 'qsar_pIC50', 'composite']])
```

## Summary

> Recap VS pipeline: library preparation → docking → ADMET filter → QSAR scoring → prioritization. Discuss limitations of rigid docking (induced fit, water molecules). Link to GNN notebook for end-to-end learned scoring.
