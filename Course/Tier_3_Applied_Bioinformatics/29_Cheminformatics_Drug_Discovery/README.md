# Module 29: Cheminformatics & Drug Discovery

**Tier 3 — Applied Bioinformatics | Module 29**

Molecular representations, fingerprints, QSAR modeling, virtual screening, and deep learning for drug discovery using RDKit and TeachOpenCADD.

## What You'll Learn
- Molecular representations: SMILES, InChI, molecular graphs
- RDKit: loading molecules, drawing structures, computing descriptors
- Molecular fingerprints: Morgan/ECFP, MACCS keys, topological
- Chemical similarity search (Tanimoto coefficient)
- QSAR (Quantitative Structure-Activity Relationship) modeling with scikit-learn
- Virtual screening: docking concepts and AutoDock Vina
- ADMET property prediction (Lipinski's Rule of Five, bioavailability)
- Target-ligand interaction databases: ChEMBL, PubChem
- Deep learning for molecules: graph neural networks (PyTorch Geometric/DGL-LifeSci)

## Prerequisites
- Module 07 (Machine Learning for Biology) — scikit-learn, model evaluation
- Module 10 (Deep Learning for Biology) — neural networks, PyTorch
- Module 02 (Protein Structure) — structural biology basics

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [Molecular Representations & RDKit](01_rdkit_basics.ipynb) | SMILES parsing, descriptors, fingerprints, similarity |
| 2 | [QSAR Modeling](02_qsar_modeling.ipynb) | ChEMBL data retrieval, feature engineering, ML models, applicability domain |
| 3 | [Virtual Screening & ADMET](03_virtual_screening.ipynb) | Docking concepts, AutoDock Vina, ADMET filters |
| 4 | [Graph Neural Networks for Molecules](04_molecular_gnn.ipynb) | Molecular graphs, GNN architectures, property prediction |

## Key Tools
| Tool | Purpose |
|------|---------|
| RDKit | Core cheminformatics toolkit |
| chembl_webresource_client | ChEMBL API access |
| DeepChem | ML/DL for drug discovery |
| AutoDock Vina | Molecular docking |
| scikit-learn | QSAR modeling |
| PyTorch Geometric / DGL | Graph neural networks |
| py3Dmol | 3D molecular visualization |

## Resources
- [RDKit documentation and getting started](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [TeachOpenCADD — Jupyter-based drug design course](https://github.com/volkamerlab/teachopencadd)
- [Practical Cheminformatics blog — Deep Learning for Drug Discovery](https://practicalcheminformatics.blogspot.com/)
- [DeepChem documentation](https://deepchem.io/)
- [ChEMBL database](https://www.ebi.ac.uk/chembl/)
- [PubChem documentation](https://pubchemdocs.ncbi.nlm.nih.gov/)

## Related Skill
`cheminformatics-drug-discovery.md` *(planned)*

---

[← Previous: 3.28 Network Biology](../28_Network_Biology/) | [Course Overview](../../README.md)
