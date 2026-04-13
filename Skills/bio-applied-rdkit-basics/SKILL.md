---
name: bio-applied-rdkit-basics
description: "**Tier 3 — Applied Bioinformatics | Module 29 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/29_Cheminformatics_Drug_Discovery/01_rdkit_basics.ipynb"
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, rdkit 2024.03+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Molecular Representations & RDKit

*Source: Course notebook `Tier_3_Applied_Bioinformatics/29_Cheminformatics_Drug_Discovery/01_rdkit_basics.ipynb`*

# Molecular Representations & RDKit

**Tier 3 — Applied Bioinformatics | Module 29 · Notebook 1**

*Prerequisites: Module 07 (Machine Learning), Module 10 (Deep Learning)*

---

**By the end of this notebook you will be able to:**
1. Parse and validate SMILES strings with RDKit
2. Draw and visualize 2D molecular structures
3. Compute molecular descriptors (MW, LogP, HBD/HBA, TPSA)
4. Generate Morgan/ECFP and MACCS fingerprints
5. Calculate Tanimoto similarity and perform chemical similarity search



**Key resources:**
- [RDKit Getting Started in Python](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [TeachOpenCADD — Jupyter drug design course](https://github.com/volkamerlab/teachopencadd)
- [Practical Cheminformatics blog](https://practicalcheminformatics.blogspot.com/)

## 1. Molecular Representations

> Explain SMILES notation (linear text), InChI (International Chemical Identifier), and molecular graph representation. Examples: aspirin, caffeine, ibuprofen.

```python
# Example: Install RDKit
# !pip install rdkit

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem import rdMolDescriptors
import pandas as pd
import numpy as np

# Example: Parse SMILES
# smiles_dict = {
#     'Aspirin':   'CC(=O)Oc1ccccc1C(=O)O',
#     'Caffeine':  'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
#     'Ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
# }
# mols = {name: Chem.MolFromSmiles(smi) for name, smi in smiles_dict.items()}

# Example: Validate (None = invalid SMILES)
# for name, mol in mols.items():
#     print(f'{name}: {"valid" if mol else "INVALID"}')
```

## 2. Structure Visualization

> Draw 2D structures with RDKit Draw. Create a grid image of all three molecules.

```python
from IPython.display import display

# Example: Draw molecule grid
# img = Draw.MolsToGridImage(list(mols.values()), legends=list(mols.keys()),
#                            molsPerRow=3, subImgSize=(300, 250))
# display(img)
```

## 3. Molecular Descriptors & Lipinski's Rule of Five

> Compute MW, LogP, HBD, HBA, TPSA for each molecule. Check Lipinski's Rule of Five for oral bioavailability prediction.

```python
# Example: Compute descriptors
# rows = []
# for name, mol in mols.items():
#     rows.append({
#         'Name':  name,
#         'MW':    Descriptors.MolWt(mol),
#         'LogP':  Descriptors.MolLogP(mol),
#         'HBD':   rdMolDescriptors.CalcNumHBD(mol),
#         'HBA':   rdMolDescriptors.CalcNumHBA(mol),
#         'TPSA':  Descriptors.TPSA(mol),
#     })
# desc_df = pd.DataFrame(rows)
# desc_df['Ro5_pass'] = (desc_df['MW'] <= 500) & (desc_df['LogP'] <= 5) & \
#                       (desc_df['HBD'] <= 5) & (desc_df['HBA'] <= 10)
# print(desc_df)
```

## 4. Molecular Fingerprints

> Compute Morgan (ECFP4) and MACCS fingerprints. Visualize as bit vectors. Explain radius and nBits parameters.

```python
from rdkit.Chem import AllChem
from rdkit import DataStructs

# Example: Morgan fingerprints (ECFP4: radius=2, nBits=2048)
# fps = {name: AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
#        for name, mol in mols.items()}

# Example: MACCS keys
# from rdkit.Chem.MACCSkeys import GenMACCSKeys
# maccs_fps = {name: GenMACCSKeys(mol) for name, mol in mols.items()}
```

## 5. Tanimoto Similarity Search

> Compute pairwise Tanimoto similarity matrix for the three molecules. Query a small library and find the most similar compound.

```python
import seaborn as sns
import matplotlib.pyplot as plt

# Example: Pairwise Tanimoto similarity
# names = list(fps.keys())
# sim_matrix = pd.DataFrame(index=names, columns=names, dtype=float)
# for n1 in names:
#     for n2 in names:
#         sim_matrix.loc[n1, n2] = DataStructs.TanimotoSimilarity(fps[n1], fps[n2])

# sns.heatmap(sim_matrix.astype(float), annot=True, cmap='Blues', vmin=0, vmax=1)
# plt.title('Tanimoto Similarity (ECFP4)')
# plt.show()
```

## Summary

> Recap molecular representations, descriptors, and fingerprints. Discuss use cases: virtual screening, chemical diversity analysis, lead optimization. Link to QSAR modeling notebook.
