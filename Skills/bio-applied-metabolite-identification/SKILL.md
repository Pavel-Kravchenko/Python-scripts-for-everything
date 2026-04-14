---
name: bio-applied-metabolite-identification
description: "Metabolite identification from MS/MS spectra: spectral matching, molecular formula prediction, and database searching (HMDB, KEGG). Use when annotating unknown metabolites."
tool_type: python
primary_tool: RDKit
---

# Metabolite Identification and Annotation

- [HMDB metabolite database](https://hmdb.ca/)
- [MZmine 3 documentation](https://mzmine.github.io/mzmine_documentation/)
- [SIRIUS + CSI:FingerID](https://bio.informatik.uni-jena.de/software/sirius/)
- [MetaboAnalyst pathway analysis](https://www.metaboanalyst.ca/MetaboAnalyst/upload/PathUploadView.xhtml)

## MSI Identification Levels

| Level | Requirement | Typical method |
|-------|-------------|----------------|
| 1 | Exact mass + MS/MS + RT match to authentic standard | Reference standard |
| 2 | Spectral library match | GNPS, MassBank, NIST |
| 3 | Putative chemical class | Mass + isotope pattern |
| 4 | Uncharacterized | Unknown feature |

## Molecular Formula Assignment

- Exact mass matching: 5 ppm tolerance
- Adducts: [M+H]+, [M+Na]+, [M+NH4]+, [M-H]-
- Isotope pattern scoring, RDBE filter

## MS/MS Spectral Matching

- Cosine similarity to GNPS / MassBank / NIST libraries
- GNPS molecular networking for unknown metabolite families
- SIRIUS + CSI:FingerID for de novo structure prediction

## Differential Abundance

- Two groups: t-test / Mann-Whitney
- Multi-factor: limma linear models
- Volcano plot: fold change vs -log10(p)
- FDR correction (Benjamini-Hochberg)

## Pathway Enrichment (MSEA)

- ORA with KEGG metabolic pathways
- QEA using all metabolite ranks
- MetaboAnalyst pathway impact plots

## Pitfalls

- **Mass accuracy thresholds**: Use ppm, not absolute Da; 5 ppm at m/z 500 = 0.0025 Da
- **Schymanski confidence levels**: Level 1 requires reference standard match -- most identifications are Level 2-3
- **Ion mode bias**: ESI+ and ESI- detect different metabolite classes; always run both modes
