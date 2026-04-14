---
name: bio-applied-virtual-screening
description: "Virtual screening for drug discovery: pharmacophore modeling, docking score filtering, and ADMET prediction. Use when computationally screening compound libraries."
tool_type: python
primary_tool: Pandas
---

# Virtual Screening  ADMET Prediction

- [AutoDock Vina documentation](https://vina.scripps.edu/)
- [SwissADME web tool](http://www.swissadme.ch/)
- [TeachOpenCADD docking talktorial](https://github.com/volkamerlab/teachopencadd)

## Virtual Screening Overview

> Diagram of drug discovery pipeline: target → hit → lead → candidate. Two VS strategies: ligand-based (pharmacophore, similarity, QSAR) vs structure-based (molecular docking). Throughput/accuracy trade-offs.

## Protein Structure Preparation

> Download EGFR crystal structure (PDB: 1IEP). Remove water/ligands. Add hydrogen atoms with Open Babel or PDBFixer. Define docking box around the ATP binding site.


## Molecular Docking with AutoDock Vina

> Prepare ligand SDF to PDBQT with meeko. Set docking box (center_x/y/z, size_x/y/z) from known active site. Run AutoDock Vina. Parse best docking score.


## ADMET Property Prediction

> Compute ADMET properties: Lipinski Ro5, LogS (solubility), Pgp substrate, CYP inhibition, hERG toxicity. Use SwissADME API or DeepChem ADMET models.


## Hit Prioritization

> Combine docking score, QSAR pIC50 prediction, and ADMET pass/fail into a composite score. Rank library of 100 compounds. Visualize top 10 structures with their scores.

```python
import pandas as pd

# Example: Composite scoring
# hits  pd.DataFrame(
#     'SMILES': smiles_library,
#     'vina_score': docking_scores,
#     'qsar_pIC50': qsar_predictions,
#     'admet_pass': admet_flags
# )
# hits_filtered  hitshits'admet_pass'
# hits_filtered'composite'  -hits_filtered'vina_score' + hits_filtered'qsar_pIC50'
# top10  hits_filtered.nlargest(10, 'composite')
# print(top10'SMILES', 'vina_score', 'qsar_pIC50', 'composite')
```python

## Summary

> Recap VS pipeline: library preparation → docking → ADMET filter → QSAR scoring → prioritization. Discuss limitations of rigid docking (induced fit, water molecules). Link to GNN notebook for end-to-end learned scoring.

## Pitfalls

- **Docking score ≠ binding affinity**: Docking scores correlate weakly with experimental Kd; always validate top hits experimentally
- **Receptor flexibility**: Rigid docking misses induced-fit binding; use ensemble docking for flexible targets
- **PAINS compounds**: Filter pan-assay interference compounds before virtual screening
