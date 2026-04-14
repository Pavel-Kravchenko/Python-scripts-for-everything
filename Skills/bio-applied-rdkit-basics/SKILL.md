---
name: bio-applied-rdkit-basics
description: "RDKit fundamentals: SMILES parsing, molecular properties, substructure search, fingerprints, and chemical similarity. Use when performing cheminformatics operations in Python."
tool_type: python
primary_tool: RDKit
---

# Molecular Representations  RDKit

- [RDKit Getting Started in Python](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [TeachOpenCADD — Jupyter drug design course](https://github.com/volkamerlab/teachopencadd)
- [Practical Cheminformatics blog](https://practicalcheminformatics.blogspot.com/)

## Molecular Representations

> Explain SMILES notation (linear text), InChI (International Chemical Identifier), and molecular graph representation. Examples: aspirin, caffeine, ibuprofen.

```python
# Example: Install RDKit
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem import rdMolDescriptors
import pandas as pd
import numpy as np

# Example: Parse SMILES
# smiles_dict
#     'Aspirin':   'CC(O)Oc1ccccc1C(O)O',
#     'Caffeine':  'Cn1cnc2c1c(O)n(C)c(O)n2C',
#     'Ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(O)O',
#
# mols  name: Chem.MolFromSmiles(smi) for name, smi in smiles_dict.items()

# Example: Validate (None  invalid SMILES)
# for name, mol in mols.items():
#     print(f'name: valid if mol else INVALID')
```python

## Structure Visualization

> Draw 2D structures with RDKit Draw. Create a grid image of all three molecules.

```python
from IPython.display import display

# Example: Draw molecule grid
# img  Draw.MolsToGridImage(list(mols.values()), legendslist(mols.keys()),
#                            molsPerRow3, subImgSize(300, 250))
# display(img)
```python

## Molecular Descriptors & Lipinski's Rule of Five

> Compute MW, LogP, HBD, HBA, TPSA for each molecule. Check Lipinski's Rule of Five for oral bioavailability prediction.


## Molecular Fingerprints

> Compute Morgan (ECFP4) and MACCS fingerprints. Visualize as bit vectors. Explain radius and nBits parameters.

```python
from rdkit.Chem import AllChem
from rdkit import DataStructs

# Example: Morgan fingerprints (ECFP4: radius2, nBits2048)
# fps  name: AllChem.GetMorganFingerprintAsBitVect(mol, radius2, nBits2048)
#        for name, mol in mols.items()

# Example: MACCS keys
# from rdkit.Chem.MACCSkeys import GenMACCSKeys
# maccs_fps  name: GenMACCSKeys(mol) for name, mol in mols.items()
```python

## Tanimoto Similarity Search

> Compute pairwise Tanimoto similarity matrix for the three molecules. Query a small library and find the most similar compound.

```python
import seaborn as sns
import matplotlib.pyplot as plt

# Example: Pairwise Tanimoto similarity
# names  list(fps.keys())
# sim_matrix  pd.DataFrame(indexnames, columnsnames, dtypefloat)
# for n1 in names:
#     for n2 in names:
#         sim_matrix.locn1, n2  DataStructs.TanimotoSimilarity(fpsn1, fpsn2)

# sns.heatmap(sim_matrix.astype(float), annotTrue, cmap'Blues', vmin0, vmax1)
# plt.title('Tanimoto Similarity (ECFP4)')
# plt.show()
```python

## Summary

> Recap molecular representations, descriptors, and fingerprints. Discuss use cases: virtual screening, chemical diversity analysis, lead optimization. Link to QSAR modeling notebook.

## Pitfalls

- **SMILES canonicalization**: Different SMILES can represent the same molecule; always canonicalize before comparison
- **Stereochemistry**: Ignoring chirality in fingerprints can merge enantiomers with different bioactivity
- **Descriptor scaling**: Molecular descriptors span different ranges; always standardize before ML
