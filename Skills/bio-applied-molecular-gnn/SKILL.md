---
name: bio-applied-molecular-gnn
description: "**Tier 3 — Applied Bioinformatics | Module 29 · Notebook 4**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/29_Cheminformatics_Drug_Discovery/04_molecular_gnn.ipynb"
---

# Graph Neural Networks for Molecular Property Prediction

*Source: Course notebook `Tier_3_Applied_Bioinformatics/29_Cheminformatics_Drug_Discovery/04_molecular_gnn.ipynb`*

# Graph Neural Networks for Molecular Property Prediction

**Tier 3 — Applied Bioinformatics | Module 29 · Notebook 4**

*Prerequisites: Notebook 2 (QSAR Modeling), Module 10 (Deep Learning)*

---

**By the end of this notebook you will be able to:**
1. Represent a molecule as a graph (atoms = nodes, bonds = edges)
2. Build a Message Passing Neural Network (MPNN) with PyTorch Geometric
3. Train a GNN for molecular property prediction (blood-brain barrier, solubility)
4. Compare GNN performance vs fingerprint-based Random Forest
5. Visualize atom-level GNN attention / feature importance



**Key resources:**
- [PyTorch Geometric documentation](https://pytorch-geometric.readthedocs.io/)
- [DeepChem documentation](https://deepchem.io/)
- [MoleculeNet benchmark](https://moleculenet.org/)

## 1. Molecules as Graphs

> Explain atom features (atomic number, degree, aromaticity, charge) and bond features (bond type, conjugated, ring). Convert SMILES to PyG Data object with node and edge feature matrices.

```python
# Example: Install PyG
# !pip install torch-geometric

from rdkit import Chem
import torch
from torch_geometric.data import Data
import numpy as np

# Example: Convert molecule to PyG graph
# def smiles_to_graph(smiles):
#     mol = Chem.MolFromSmiles(smiles)
#     # Node features: [atomic_num, degree, formal_charge, is_aromatic]
#     node_features = [[a.GetAtomicNum(), a.GetDegree(),
#                        a.GetFormalCharge(), int(a.GetIsAromatic())]
#                       for a in mol.GetAtoms()]
#     # Edge indices
#     edges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
#     edges += [(j, i) for i, j in edges]  # undirected
#     edge_index = torch.tensor(edges, dtype=torch.long).t()
#     x = torch.tensor(node_features, dtype=torch.float)
#     return Data(x=x, edge_index=edge_index)
```

## 2. MPNN Architecture

> Implement a 3-layer GCN (GraphConv) with global mean pooling and a 2-layer MLP head. Use ReLU activations and dropout.

```python
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool

# Example: GNN model
# class MolGNN(nn.Module):
#     def __init__(self, in_channels, hidden_channels, out_channels):
#         super().__init__()
#         self.conv1 = GCNConv(in_channels, hidden_channels)
#         self.conv2 = GCNConv(hidden_channels, hidden_channels)
#         self.conv3 = GCNConv(hidden_channels, hidden_channels)
#         self.lin = nn.Linear(hidden_channels, out_channels)

#     def forward(self, x, edge_index, batch):
#         x = self.conv1(x, edge_index).relu()
#         x = self.conv2(x, edge_index).relu()
#         x = self.conv3(x, edge_index).relu()
#         x = global_mean_pool(x, batch)
#         return self.lin(x)
```

## 3. Training on MoleculeNet BBBP Dataset

> Load BBBP (blood-brain barrier penetration) dataset via DeepChem or direct download. Train/test split by scaffold. Train GNN for 50 epochs. Track train/val loss and AUC.

```python
import torch.optim as optim

# Example: Training loop
# model = MolGNN(in_channels=4, hidden_channels=64, out_channels=1)
# optimizer = optim.Adam(model.parameters(), lr=1e-3)
# criterion = nn.BCEWithLogitsLoss()

# for epoch in range(50):
#     model.train()
#     total_loss = 0
#     for batch in train_loader:
#         optimizer.zero_grad()
#         out = model(batch.x, batch.edge_index, batch.batch).squeeze()
#         loss = criterion(out, batch.y.float())
#         loss.backward()
#         optimizer.step()
#         total_loss += loss.item()
#     if epoch % 10 == 0:
#         print(f'Epoch {epoch:3d}, Loss: {total_loss/len(train_loader):.4f}')
```

## 4. GNN vs Fingerprint Baseline

> Compare GNN AUC vs Random Forest on Morgan fingerprints. Discussion of when GNNs outperform hand-crafted features.

```python
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

# Example: Evaluate GNN and RF
# model.eval()
# gnn_preds = []
# gnn_labels = []
# with torch.no_grad():
#     for batch in test_loader:
#         out = torch.sigmoid(model(batch.x, batch.edge_index, batch.batch).squeeze())
#         gnn_preds.extend(out.numpy())
#         gnn_labels.extend(batch.y.numpy())

# gnn_auc = roc_auc_score(gnn_labels, gnn_preds)
# print(f'GNN AUC: {gnn_auc:.3f}')
# # Example: RF baseline AUC for comparison
```

## Summary

> Recap GNN-based molecular property prediction. Discuss attention mechanisms (GAT) for atom importance. Point to GNN extensions: 3D-aware models (DimeNet, SchNet), generative models for drug design.
