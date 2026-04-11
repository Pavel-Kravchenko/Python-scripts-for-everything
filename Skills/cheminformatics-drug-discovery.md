---
name: cheminformatics-drug-discovery
description: Cheminformatics and drug discovery — RDKit molecular representations, fingerprints, Tanimoto similarity, QSAR modeling with ChEMBL data, AutoDock Vina docking, ADMET prediction, graph neural networks for molecules
---

## When to Use

Use this skill when:
- Working with chemical structures (SMILES, SDF, PDB ligands)
- Computing molecular descriptors or fingerprints
- Building QSAR models for bioactivity prediction
- Running virtual screening or molecular docking
- Predicting ADMET properties (absorption, distribution, metabolism, excretion, toxicity)
- Using graph neural networks for molecular property prediction

## Quick Reference

| Task | Tool | Key Method |
|------|------|-----------|
| Parse SMILES | RDKit | `Chem.MolFromSmiles(smi)` |
| Draw structure | RDKit | `Draw.MolToImage(mol)` |
| Molecular weight | RDKit | `Descriptors.MolWt(mol)` |
| LogP | RDKit | `Descriptors.MolLogP(mol)` |
| Morgan FP (ECFP4) | RDKit | `AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)` |
| MACCS keys | RDKit | `MACCSkeys.GenMACCSKeys(mol)` |
| Tanimoto similarity | RDKit | `DataStructs.TanimotoSimilarity(fp1, fp2)` |
| ChEMBL query | chembl_webresource_client | `activity.filter(target_chembl_id=...)` |
| Molecular docking | AutoDock Vina | `vina --config vina.txt` |
| GNN property prediction | DeepChem / PyG | `dc.models.AttentiveFPModel()` |

## Key Patterns

**Pattern 1: RDKit basics**
```python
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem, rdMolDescriptors
from rdkit import DataStructs

# Parse and validate
mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')  # Aspirin
assert mol is not None, "Invalid SMILES"

# Compute Lipinski descriptors
descriptors = {
    'MW':   Descriptors.MolWt(mol),
    'LogP': Descriptors.MolLogP(mol),
    'HBD':  rdMolDescriptors.CalcNumHBD(mol),
    'HBA':  rdMolDescriptors.CalcNumHBA(mol),
    'TPSA': Descriptors.TPSA(mol),
}
ro5_pass = (descriptors['MW'] <= 500 and descriptors['LogP'] <= 5 and
            descriptors['HBD'] <= 5 and descriptors['HBA'] <= 10)
```

**Pattern 2: Morgan fingerprints and similarity**
```python
from rdkit.Chem import AllChem
from rdkit import DataStructs

# ECFP4 fingerprint (radius=2, 2048 bits)
fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

# Tanimoto similarity between two molecules
sim = DataStructs.TanimotoSimilarity(fp1, fp2)

# Similarity search in a library
library_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in library_mols]
similarities = DataStructs.BulkTanimotoSimilarity(query_fp, library_fps)
top_hits = sorted(enumerate(similarities), key=lambda x: -x[1])[:10]
```

**Pattern 3: ChEMBL QSAR dataset**
```python
from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np

activity = new_client.activity
data = activity.filter(
    target_chembl_id='CHEMBL203',  # EGFR
    standard_type='IC50'
).only(['molecule_chembl_id', 'standard_value', 'canonical_smiles'])

df = pd.DataFrame(list(data)).dropna()
df = df[df['standard_value'].astype(float) > 0]
df['pIC50'] = -np.log10(df['standard_value'].astype(float) * 1e-9)
df['active'] = (df['pIC50'] >= 6).astype(int)
```

**Pattern 4: QSAR model**
```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score

# Featurize
X = np.array([list(AllChem.GetMorganFingerprintAsBitVect(
    Chem.MolFromSmiles(smi), 2, 2048))
    for smi in df['canonical_smiles']])
y = df['active'].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)
rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
rf.fit(X_train, y_train)
print(f"AUC: {roc_auc_score(y_test, rf.predict_proba(X_test)[:,1]):.3f}")
```

**Pattern 5: Molecular docking (AutoDock Vina)**
```bash
# Prepare receptor (remove water/heterogens, add H)
pdbfixer receptor.pdb --output receptor_prep.pdb --add-hydrogens --remove-heterogens

# Prepare ligand
mk_prepare_ligand.py -i ligand.sdf -o ligand.pdbqt

# Run docking
vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
    --center_x 22.5 --center_y 5.0 --center_z 18.0 \
    --size_x 20 --size_y 20 --size_z 20 \
    --exhaustiveness 8 --out docked.pdbqt --log docking.log
```

## Lipinski's Rule of Five

A compound is predicted orally bioavailable if it satisfies ≥4 of:
- MW ≤ 500 Da
- LogP ≤ 5
- H-bond donors ≤ 5
- H-bond acceptors ≤ 10
- TPSA ≤ 140 Å²

## Common Pitfalls

- **SMILES validation** — always check `Chem.MolFromSmiles() is not None` before processing
- **Scaffold split** — use scaffold-based train/test splits (not random) to avoid data leakage in QSAR
- **pIC50 vs IC50** — convert IC50 (nM) to pIC50 = -log10(IC50 × 1e-9) for regression targets
- **Docking scoring functions** — Vina score is an approximation; always validate top hits with MD or experimental assay
- **ADMET early filtering** — apply Lipinski/Veber filters before docking to reduce computational cost
