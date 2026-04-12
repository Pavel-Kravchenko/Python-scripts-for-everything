# Module 30: Single-Cell RNA-seq Analysis

**Tier 3 — Applied Bioinformatics | Module 30**

End-to-end scRNA-seq analysis: quality control, clustering, cell type annotation, and trajectory inference.

## What You'll Learn
- scRNA-seq experimental workflow: 10x Chromium, Smart-seq2, cell barcodes, UMIs
- Count matrix generation with Cell Ranger / STARsolo
- Quality control: mitochondrial content, doublet detection, ambient RNA removal
- Normalization, HVG selection, PCA, UMAP, and Leiden clustering
- Differential expression testing between clusters (Wilcoxon, DESeq2 pseudobulk)
- Manual and automated cell type annotation (SingleR, CellTypist)
- Pseudotime and trajectory analysis (DPT, Monocle 3)
- RNA velocity with scVelo to infer differentiation direction

## Prerequisites
- Module 03 (RNA-seq Analysis) — bulk RNA-seq concepts, normalization, DEseq2
- Module 07 (Machine Learning for Biology) — dimensionality reduction, clustering

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [QC & Preprocessing](01_scrna_preprocessing.ipynb) | Cell Ranger, QC metrics, normalization, HVG selection |
| 2 | [Dimensionality Reduction & Clustering](02_dimensionality_reduction.ipynb) | PCA, UMAP, Leiden, marker genes |
| 3 | [Cell Type Annotation](03_cell_type_annotation.ipynb) | Manual annotation, SingleR, CellTypist |
| 4 | [Trajectory Analysis & RNA Velocity](04_trajectory_analysis.ipynb) | DPT, Monocle 3, scVelo |

## Key Tools
| Tool | Purpose |
|------|---------|
| Scanpy | Python scRNA-seq analysis framework |
| Seurat | R scRNA-seq framework |
| Cell Ranger | 10x read processing → count matrix |
| scVelo | RNA velocity |
| Monocle 3 | Trajectory inference |
| SingleR | Reference-based cell type annotation |
| CellTypist | Python automated annotation |
| DoubletFinder / scDblFinder | Doublet detection |

## Resources
- [Scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html)
- [Single-cell best practices book (Theis lab)](https://www.sc-best-practices.org/)
- [Harvard HBC scRNA-seq training](https://hbctraining.github.io/scRNA-seq_online/)
- [OSCA (Orchestrating Single-Cell Analysis)](https://bioconductor.org/books/release/OSCA/)
- [scVelo documentation](https://scvelo.readthedocs.io/)

## Related Skill
`scrna-seq-analysis.md` *(planned)*

---

[← Previous: 3.29 Cheminformatics & Drug Discovery](../29_Cheminformatics_Drug_Discovery/) | [Course Overview](../../README.md) | [Next: 3.31 Single-Cell Multi-Omics →](../31_Single_Cell_Multi_Omics/)
