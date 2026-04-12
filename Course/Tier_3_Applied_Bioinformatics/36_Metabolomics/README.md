# Module 36: Metabolomics

**Tier 3 — Applied Bioinformatics | Module 36**

Untargeted LC-MS metabolomics: preprocessing, metabolite identification, pathway enrichment, and constraint-based metabolic modeling.

## What You'll Learn
- LC-MS and GC-MS experimental workflows; mzML open data format
- XCMS peak picking, retention time alignment, and feature correspondence
- Normalization strategies: PQN, internal standard, LOESS batch correction
- Metabolite identification levels (MSI 1–4); exact mass matching, MS/MS spectral library
- GNPS molecular networking for unknown metabolite families
- SIRIUS + CSI:FingerID for de novo structure prediction
- Differential abundance analysis and metabolite set enrichment (MSEA)
- Genome-scale metabolic models (GEMs) and flux balance analysis with COBRApy
- 13C metabolic flux analysis (13C-MFA) fundamentals

## Prerequisites
- Module 06 (Statistics for Bioinformatics) — linear models, multiple testing
- Module 18 (Proteomics & Structural Methods) — mass spectrometry concepts

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [LC-MS Preprocessing](01_lc_ms_preprocessing.ipynb) | XCMS, peak picking, normalization, PCA |
| 2 | [Metabolite Identification & Statistics](02_metabolite_identification.ipynb) | Exact mass, MS/MS matching, GNPS, MSEA |
| 3 | [Metabolic Flux Analysis](03_metabolic_flux.ipynb) | COBRApy FBA, gene knockout, transcriptomics integration |

## Key Tools
| Tool | Purpose |
|------|---------|
| XCMS | Peak detection and alignment (R) |
| MZmine 3 | Peak processing (Java, GUI+CLI) |
| pyOpenMS | Python mass spec processing |
| GNPS | Spectral library + molecular networking |
| SIRIUS | Structure elucidation from MS/MS |
| MetaboAnalyst | Statistical analysis and pathways |
| COBRApy | Flux balance analysis (Python) |

## Resources
- [HMDB (Human Metabolome Database)](https://hmdb.ca/)
- [MetaboAnalyst 6.0](https://www.metaboanalyst.ca/)
- [XCMS documentation](https://www.bioconductor.org/packages/release/bioc/html/xcms.html)
- [COBRApy documentation](https://cobrapy.readthedocs.io/)
- [GNPS documentation](https://ccms-ucsd.github.io/GNPSDocumentation/)

## Related Skill
`metabolomics.md` *(planned)*

---

[← Previous: 3.35 Small RNA & ncRNA](../35_Small_RNA_and_ncRNA/) | [Course Overview](../../README.md) | [Next: 3.37 Virology Bioinformatics →](../37_Virology_Bioinformatics/)
