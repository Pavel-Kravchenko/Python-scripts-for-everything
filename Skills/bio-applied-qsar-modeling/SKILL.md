---
name: bio-applied-qsar-modeling
description: "**Tier 3 — Applied Bioinformatics | Module 29 · Notebook 2**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/29_Cheminformatics_Drug_Discovery/02_qsar_modeling.ipynb"
---

# QSAR Modeling

*Source: Course notebook `Tier_3_Applied_Bioinformatics/29_Cheminformatics_Drug_Discovery/02_qsar_modeling.ipynb`*

# QSAR Modeling

**Tier 3 — Applied Bioinformatics | Module 29 · Notebook 2**

*Prerequisites: Notebook 1 (RDKit Basics), Module 07 (Machine Learning)*

---

**By the end of this notebook you will be able to:**
1. Retrieve bioactivity data from ChEMBL for a target protein
2. Curate activity data (IC50 → pIC50, activity cliffs, structural duplicates)
3. Featurize molecules with Morgan fingerprints and RDKit descriptors
4. Train and evaluate Random Forest and XGBoost QSAR models
5. Define and apply an applicability domain (AD) check



**Key resources:**
- [TeachOpenCADD QSAR notebook](https://github.com/volkamerlab/teachopencadd/tree/main/teachopencadd/talktorials/T022_ligand_based_screening_neural_network)
- [ChEMBL database](https://www.ebi.ac.uk/chembl/)
- [chembl_webresource_client documentation](https://github.com/chembl/chembl_webresource_client)

## 1. Data Retrieval from ChEMBL

> Query ChEMBL for EGFR (CHEMBL203) bioactivity data (IC50 assays, Ki assays). Filter for single-protein assays with standard IC50 values.

```python
# Example: Install ChEMBL client
# !pip install chembl_webresource_client

from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np

# Example: Fetch EGFR bioactivity data
# activity = new_client.activity
# egfr_activities = activity.filter(
#     target_chembl_id='CHEMBL203',
#     standard_type='IC50'
# ).only(['molecule_chembl_id', 'standard_value', 'standard_units', 'canonical_smiles'])
# df = pd.DataFrame(list(egfr_activities))
# print(df.shape)
```

## 2. Data Curation

> Convert IC50 (nM) to pIC50. Remove duplicates (keep median). Filter out molecules with SMILES errors. Split into active (pIC50 > 6) vs inactive.

```python
# Example: Convert to pIC50 and curate
# df = df.dropna(subset=['standard_value', 'canonical_smiles'])
# df = df[df['standard_value'] > 0]
# df['pIC50'] = -np.log10(df['standard_value'].astype(float) * 1e-9)
# df['active'] = (df['pIC50'] >= 6).astype(int)
# print(df['active'].value_counts())
```

## 3. Feature Generation

> Parse SMILES with RDKit. Compute Morgan fingerprints (ECFP4) and 200 RDKit 2D descriptors. Combine into a feature matrix.

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

# Example: Compute Morgan FP + 2D descriptors
# calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors.descList])
# features = []
# for smi in df['canonical_smiles']:
#     mol = Chem.MolFromSmiles(smi)
#     fp = list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048))
#     desc = list(calc.CalcDescriptors(mol))
#     features.append(fp + desc)
# X = np.array(features)
# y = df['active'].values
```

## 4. Model Training and Evaluation

> Scaffold-based train/test split (avoid data leakage). Train Random Forest. Evaluate with AUC-ROC, Matthews Correlation Coefficient, precision-recall curve.

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, matthews_corrcoef, classification_report

# Example: Train/test split and model
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
# rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
# rf.fit(X_train, y_train)
# y_pred = rf.predict(X_test)
# y_prob = rf.predict_proba(X_test)[:, 1]
# print(f'AUC: {roc_auc_score(y_test, y_prob):.3f}')
# print(f'MCC: {matthews_corrcoef(y_test, y_pred):.3f}')
```

## 5. Applicability Domain

> Define applicability domain using k-NN distance threshold. Flag test compounds outside the AD as unreliable predictions.

```python
from sklearn.neighbors import NearestNeighbors

# Example: k-NN applicability domain
# knn = NearestNeighbors(n_neighbors=5)
# knn.fit(X_train)
# distances, _ = knn.kneighbors(X_test)
# avg_distances = distances.mean(axis=1)
# ad_threshold = avg_distances.mean() + 3 * avg_distances.std()
# in_ad = avg_distances <= ad_threshold
# print(f'Compounds in AD: {in_ad.sum()}/{len(in_ad)}')
```

## Summary

> Recap QSAR workflow: data → curation → features → model → AD. Discuss scaffold splits vs random splits for realistic evaluation. Link to virtual screening notebook.
