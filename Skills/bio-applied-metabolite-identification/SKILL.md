---
name: bio-applied-metabolite-identification
description: "Metabolite identification from MS/MS spectra: spectral matching, molecular formula prediction, and database searching (HMDB, KEGG). Use when annotating unknown metabolites."
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/36_Metabolomics/02_metabolite_identification.ipynb"
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: rdkit 2024.03+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Metabolite Identification and Annotation

*Source: Course notebook `Tier_3_Applied_Bioinformatics/36_Metabolomics/02_metabolite_identification.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 36 · Notebook 2**

*Prerequisites: Notebook 1 (LC-MS Preprocessing)*

---

**By the end of this notebook you will be able to:**
1. Apply exact mass and isotope pattern matching for molecular formula assignment
2. Use MS/MS fragmentation spectra for compound identification via spectral matching
3. Query HMDB, KEGG, and MetaboLights for metabolite annotation
4. Perform differential abundance testing with appropriate statistical tests
5. Conduct metabolite set enrichment analysis (MSEA) for pathway interpretation



**Key resources:**
- [HMDB metabolite database](https://hmdb.ca/)
- [MZmine 3 documentation](https://mzmine.github.io/mzmine_documentation/)
- [SIRIUS + CSI:FingerID](https://bio.informatik.uni-jena.de/software/sirius/)
- [MetaboAnalyst pathway analysis](https://www.metaboanalyst.ca/MetaboAnalyst/upload/PathUploadView.xhtml)

## 1. Levels of Metabolite Identification

> MSI reporting standards (Level 1–4). Level 1: exact mass + MS/MS + RT match to authentic standard. Level 2: spectral match to library. Level 3: putative based on chemical class. Level 4: uncharacterized.

## 2. Molecular Formula Assignment

> Exact mass matching with 5 ppm tolerance. Adduct consideration ([M+H]+, [M+Na]+, [M+NH4]+, [M-H]-). Isotope pattern scoring. RDBE (rings + double bond equivalents) filter.

```python
# Example: Exact mass matching
# import pyteomics.mass as mass
# from pyteomics import fasta
# # For metabolomics, use molmass or rdkit
# from rdkit.Chem import Descriptors
# from rdkit.Chem import MolFromSmiles
# mol = MolFromSmiles('C(C(=O)O)N')  # Alanine
# exact_mass = Descriptors.ExactMolWt(mol)
```python

## 3. MS/MS Spectral Matching

> Cosine similarity to GNPS / MassBank / NIST spectral libraries. GNPS molecular networking for unknown metabolite families. SIRIUS + CSI:FingerID for de novo structure prediction.

## 4. Differential Abundance Analysis

> t-test / Mann-Whitney for two groups. Linear models with limma for multi-factor designs. Volcano plot: fold change vs -log10 p-value. Multiple testing correction (Benjamini-Hochberg).

## 5. Pathway Enrichment (MSEA)

> Over-representation analysis (ORA) with KEGG metabolic pathways. Quantitative enrichment analysis (QEA) using all metabolite ranks. MetaboAnalyst pathway impact plots. iPath3 visualization.

## Common Pitfalls

- **Mass accuracy thresholds**: Use ppm, not absolute Da; 5 ppm at m/z 500 = 0.0025 Da
- **Schymanski confidence levels**: Level 1 requires reference standard match — most identifications are Level 2-3
- **Ion mode bias**: ESI+ and ESI- detect different metabolite classes; always run both modes
