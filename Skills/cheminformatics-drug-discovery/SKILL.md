---
name: cheminformatics-drug-discovery
description: Cheminformatics and drug discovery — RDKit molecular representations, fingerprints, Tanimoto similarity, QSAR modeling with ChEMBL data, AutoDock Vina docking, ADMET prediction, graph neural networks for molecules
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, pandas 2.1+, rdkit 2024.03+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Cheminformatics & Drug Discovery

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
```python

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
```python

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
```python

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
```python

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
```python

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

## Code Templates

### Lipinski Filter for a Compound Library
```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import pandas as pd

def lipinski_filter(smiles_list):
    """Return DataFrame with Ro5 properties and pass/fail."""
    results = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            results.append({'smiles': smi, 'valid': False})
            continue
        props = {
            'smiles': smi,
            'valid': True,
            'MW':   Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'HBD':  rdMolDescriptors.CalcNumHBD(mol),
            'HBA':  rdMolDescriptors.CalcNumHBA(mol),
            'TPSA': Descriptors.TPSA(mol),
        }
        props['ro5_pass'] = (
            props['MW'] <= 500 and props['LogP'] <= 5 and
            props['HBD'] <= 5 and props['HBA'] <= 10
        )
        results.append(props)
    return pd.DataFrame(results)

df = lipinski_filter(smiles_list)
print(f"Ro5 pass rate: {df['ro5_pass'].mean():.1%}")
```python

### Scaffold-Based Train/Test Split
```python
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupShuffleSplit
import numpy as np

def get_scaffold(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ''
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)

scaffolds = [get_scaffold(smi) for smi in df['canonical_smiles']]
# Encode unique scaffolds as group IDs
unique_scaffolds = {s: i for i, s in enumerate(set(scaffolds))}
groups = np.array([unique_scaffolds[s] for s in scaffolds])

gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
train_idx, test_idx = next(gss.split(X, y, groups=groups))
X_train, X_test = X[train_idx], X[test_idx]
y_train, y_test = y[train_idx], y[test_idx]
```python

### Parse SDF File with RDKit
```python
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

def read_sdf_to_df(sdf_path):
    """Load SDF and compute basic descriptors for each molecule."""
    records = []
    for mol in Chem.SDMolSupplier(sdf_path, removeHs=True):
        if mol is None:
            continue
        rec = {'smiles': Chem.MolToSmiles(mol)}
        # Copy SD properties
        rec.update(mol.GetPropsAsDict())
        rec['MW']   = Descriptors.MolWt(mol)
        rec['LogP'] = Descriptors.MolLogP(mol)
        records.append(rec)
    return pd.DataFrame(records)

compounds = read_sdf_to_df('compound_library.sdf')
print(f"Loaded {len(compounds)} molecules")
```python

### Virtual Screening Results Parser
```python
import re
import pandas as pd

def parse_vina_log(log_path):
    """Parse AutoDock Vina log to get docking scores."""
    results = []
    with open(log_path) as f:
        for line in f:
            m = re.match(r'\s+(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)', line)
            if m:
                results.append({
                    'mode': int(m.group(1)),
                    'affinity_kcal_mol': float(m.group(2)),
                    'rmsd_lb': float(m.group(3)),
                    'rmsd_ub': float(m.group(4)),
                })
    return pd.DataFrame(results)

scores = parse_vina_log('docking.log')
print(f"Best pose: {scores.iloc[0]['affinity_kcal_mol']} kcal/mol")
```python

## Related Skills
- `structural-bioinformatics` — protein structure, PDB parsing, binding sites
- `ml-deep-learning-bio` — graph neural networks, molecular property prediction
- `python-core-bio` — data structures, file I/O, sequence handling
- `numerical-methods-bio` — optimization, differential equations for pharmacokinetics
