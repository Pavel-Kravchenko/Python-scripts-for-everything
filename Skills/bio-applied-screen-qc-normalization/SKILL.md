---
name: bio-applied-screen-qc-normalization
description: "**Tier 3 — Applied Bioinformatics | Module 33 · Notebook 2**"
tool_type: reference
source_notebook: "Tier_3_Applied_Bioinformatics/33_CRISPR_Screen_Analysis/02_screen_qc_normalization.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# CRISPR Screen QC, Normalization, and Advanced Methods

*Source: Course notebook `Tier_3_Applied_Bioinformatics/33_CRISPR_Screen_Analysis/02_screen_qc_normalization.ipynb`*

# CRISPR Screen QC, Normalization, and Advanced Methods

**Tier 3 — Applied Bioinformatics | Module 33 · Notebook 2**

*Prerequisites: Notebook 1 (MAGeCK Gene Essentiality)*

---

**By the end of this notebook you will be able to:**
1. Assess screen quality with Gini index, read depth distribution, and sgRNA evenness
2. Apply CRISPRcleanR for copy-number-bias correction
3. Analyze CRISPRi / CRISPRa screens for regulatory element discovery
4. Use MAGeCK-VISPR for interactive result exploration
5. Compare screen hits to functional genomics databases (DepMap, BioGRID)



**Key resources:**
- [CRISPRcleanR paper (Iorio et al. 2018)](https://www.nature.com/articles/s41598-018-29pompier)
- [MAGeCK-VISPR](https://bitbucket.org/liulab/mageck-vispr)
- [DepMap 22Q4 data](https://depmap.org/portal/download/)

## 1. Screen Quality Metrics

> Read depth distribution. Gini index for library representation evenness. Non-targeting control guide behavior. Fraction of reads mapping to library. Replicate correlation.

## 2. Copy-Number Bias Correction with CRISPRcleanR

> How genomic amplifications inflate LFC scores. CRISPRcleanR genomic-segment normalization. Before/after scatter plots for amplified loci.

## 3. CRISPRi / CRISPRa Screens

> dCas9-KRAB (CRISPRi) for silencing vs dCas9-VPR (CRISPRa) for activation. TSS-targeting guide design. CASA algorithm for CRISPRi screens. Non-coding regulatory screen analysis.

## 4. Comparing Screen Hits to Databases

> Cross-reference hits with DepMap CERES scores. PPI network enrichment (STRING). Pathway analysis of hits (MSigDB). Drug target overlap (DGIdb).
