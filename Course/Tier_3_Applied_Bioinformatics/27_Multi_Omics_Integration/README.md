# Module 27: Multi-Omics Integration

**Tier 3 — Applied Bioinformatics | Module 27**

Integrating RNA-seq, proteomics, metabolomics, and other omics layers to uncover systems-level biology.

## What You'll Learn
- Multi-omics data types and their complementarity
- Data harmonization: scaling, missing value imputation across platforms
- Pairwise correlation analysis between omics layers
- Latent factor methods: MOFA2 (Multi-Omics Factor Analysis)
- Partial Least Squares — Discriminant Analysis (PLS-DA) with mixOmics
- Multi-block PLS (DIABLO) for supervised multi-omics classification
- Network-based integration: regulatory networks across omics
- Visualization: loadings plots, STACKED integration plots, circos
- Pathway enrichment across integrated layers

## Prerequisites
- Module 03 (RNA-seq Analysis) — transcriptomics
- Module 07 (Machine Learning for Biology) — dimensionality reduction
- Module 18 (Proteomics & Structural Methods) — proteomics basics

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [Data Harmonization](01_data_harmonization.ipynb) | Scaling, batch correction, missing data strategies |
| 2 | [MOFA2 Factor Analysis](02_mofa2.ipynb) | Latent factor decomposition, variance explained per omic |
| 3 | [mixOmics Integration](03_mixomics.ipynb) | PLS-DA, DIABLO, sPLS, MINT for multi-study |

## Key Tools
| Tool | Purpose |
|------|---------|
| MOFA2 | Multi-omics factor analysis (Python + R) |
| mixOmics | PLS-based multi-omics (R) |
| mofapy2 | Python interface for MOFA2 |
| inmoose | Python port of mixOmics tools |
| pandas / scikit-learn | Data preprocessing |
| seaborn / plotly | Visualization |

## Resources
- [MOFA2 documentation and tutorials](https://biofam.github.io/MOFA2/)
- [mixOmics documentation](http://mixomics.org/)
- [Galaxy Training — Multi-Omics](https://training.galaxyproject.org/)
- [Multi-Omics Factor Analysis (Argelaguet et al., 2018)](https://www.embopress.org/doi/10.15252/msb.20178124)
- [mixOmics: An R package for omics feature selection (Rohart et al., 2017)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005752)

## Related Skill
`multi-omics-integration.md` *(planned)*

---

[← Previous: 3.26 Shotgun Metagenomics](../26_Metagenomics_Shotgun/) | [Course Overview](../../README.md) | [Next: 3.28 Network Biology →](../28_Network_Biology/)
