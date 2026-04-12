# Module 31: Single-Cell Multi-Omics

**Tier 3 — Applied Bioinformatics | Module 31**

Single-cell chromatin accessibility (scATAC-seq), CITE-seq protein co-profiling, multiome integration, and batch correction.

## What You'll Learn
- scATAC-seq workflow: fragment files, TF-IDF normalization, LSI dimensionality reduction
- Peak calling from single-cell data; co-accessibility with Cicero
- TF motif activity scoring with chromVAR
- CITE-seq: ADT protein data, DSB normalization, WNN joint embedding
- 10x Multiome: simultaneous RNA + ATAC per cell
- Batch effect diagnosis and correction (Harmony, scVI, Seurat CCA)
- Label transfer and reference atlas mapping (Azimuth)
- Differential abundance analysis (Milo)

## Prerequisites
- Module 23 (TF Footprinting & ATAC-seq) — bulk ATAC-seq concepts
- Module 30 (Single-Cell RNA-seq) — scRNA-seq workflow and AnnData

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [scATAC-seq & Chromatin Accessibility](01_scatac_chromatin.ipynb) | Fragment files, LSI, peak-gene links, chromVAR |
| 2 | [CITE-seq & Multiome Integration](02_cite_seq_integration.ipynb) | ADT/DSB normalization, WNN, paired RNA+ATAC |
| 3 | [Batch Correction & Integration](03_sc_integration.ipynb) | Harmony, scVI, label transfer, Azimuth |

## Key Tools
| Tool | Purpose |
|------|---------|
| Signac | scATAC-seq analysis (R/Bioconductor) |
| SnapATAC2 | Python scATAC-seq analysis |
| ArchR | Fast scATAC-seq pipeline |
| muon | Multi-modal AnnData (Python) |
| harmonypy | Batch correction |
| scvi-tools | Deep generative integration |
| Azimuth | Reference atlas mapping |
| Milo | Differential abundance testing |

## Resources
- [Signac documentation](https://stuartlab.org/signac/)
- [SnapATAC2 documentation](https://snapatac2.readthedocs.io/)
- [muon documentation](https://muon.readthedocs.io/)
- [scib benchmarking (Luecken et al. 2022)](https://www.nature.com/articles/s41592-021-01336-8)

## Related Skill
`sc-multiomics.md` *(planned)*

---

[← Previous: 3.30 Single-Cell RNA-seq](../30_Single_Cell_RNA_seq/) | [Course Overview](../../README.md) | [Next: 3.32 DNA Methylation Analysis →](../32_DNA_Methylation_Analysis/)
