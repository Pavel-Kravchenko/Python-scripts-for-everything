# Module 34: Immunogenomics

**Tier 3 — Applied Bioinformatics | Module 34**

Adaptive immune receptor repertoire analysis, V(D)J recombination, HLA typing, and neoantigen prediction.

## What You'll Learn
- TCR and BCR V(D)J recombination biology and CDR3 diversity
- IMGT nomenclature for immunoglobulin and TCR gene segments
- scTCR-seq / scBCR-seq analysis with scirpy (Python)
- Clonotype diversity metrics: Shannon entropy, clonal expansion index
- Bulk immune repertoire from RNA-seq using TRUST4 / MiXCR
- B-cell lineage tree inference from somatic hypermutation
- HLA typing with OptiType from NGS data
- Neoantigen prediction pipeline: VEP → NetMHCpan → pVACseq

## Prerequisites
- Module 13 (Computational Genetics) — genetic variation concepts
- Module 30 (Single-Cell RNA-seq) — for VDJ + transcriptome co-analysis

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [V(D)J Biology & Clonotype Analysis](01_vdj_biology.ipynb) | Recombination, AIRR format, diversity metrics, V-gene usage |
| 2 | [Immune Repertoire Sequencing](02_immune_repertoire.ipynb) | TRUST4, clonal tracking, B-cell phylogenetics |
| 3 | [HLA Typing & Neoantigen Prediction](03_hla_typing.ipynb) | OptiType, NetMHCpan, pVACseq pipeline |

## Key Tools
| Tool | Purpose |
|------|---------|
| scirpy | scTCR/BCR analysis (Python) |
| TRUST4 | BCR/TCR reconstruction from RNA-seq |
| MiXCR | Immune repertoire analysis |
| OptiType | HLA class I typing |
| HLA*LA | HLA class I+II typing |
| NetMHCpan | Peptide-MHC binding prediction |
| pVACseq | Neoantigen prioritization pipeline |
| dowser | B-cell lineage trees (R) |

## Resources
- [IMGT database](https://www.imgt.org/)
- [AIRR Community standards](https://docs.airr-community.org/)
- [scirpy documentation](https://scirpy.scverse.org/)
- [VDJdb antigen database](https://vdjdb.cdr3.net/)
- [pVACtools documentation](https://pvactools.readthedocs.io/)

## Related Skill
`immunogenomics.md` *(planned)*

---

[← Previous: 3.33 CRISPR Screen Analysis](../33_CRISPR_Screen_Analysis/) | [Course Overview](../../README.md) | [Next: 3.35 Small RNA & ncRNA →](../35_Small_RNA_and_ncRNA/)
