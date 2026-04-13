---
name: bio-applied-lc-ms-preprocessing
description: "**Tier 3 — Applied Bioinformatics | Module 36 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/36_Metabolomics/01_lc_ms_preprocessing.ipynb"
---

# LC-MS Metabolomics Data Preprocessing

*Source: Course notebook `Tier_3_Applied_Bioinformatics/36_Metabolomics/01_lc_ms_preprocessing.ipynb`*

# LC-MS Metabolomics Data Preprocessing

**Tier 3 — Applied Bioinformatics | Module 36 · Notebook 1**

*Prerequisites: Module 06 (Statistics for Bioinformatics), Module 18 (Proteomics & Structural Methods)*

---

**By the end of this notebook you will be able to:**
1. Describe LC-MS and GC-MS metabolomics experimental workflows
2. Process raw data with XCMS: peak picking, alignment, gap filling
3. Understand m/z, retention time, and adduct formation in mass spectra
4. Apply normalization strategies (PQN, IS-based, LOESS) for batch correction
5. Perform exploratory analysis: PCA, heatmaps, missing value imputation



**Key resources:**
- [XCMS documentation](https://www.bioconductor.org/packages/release/bioc/html/xcms.html)
- [MetaboAnalyst 6.0](https://www.metaboanalyst.ca/)
- [HMDB (Human Metabolome Database)](https://hmdb.ca/)
- [Galaxy Metabolomics training](https://training.galaxyproject.org/training-material/topics/metabolomics/)

## 1. Metabolomics Overview

> Metabolome definition: small molecules < 1500 Da. Untargeted vs targeted vs semi-targeted approaches. LC-MS (liquid chromatography) vs GC-MS (volatile compounds, requires derivatization) vs NMR. Sample types: plasma, urine, tissue, cells.

## 2. LC-MS Data Acquisition

> Chromatographic separation (RPLC, HILIC). ESI ionization: positive vs negative mode. Data-dependent acquisition (DDA) vs data-independent (SWATH/DIA). mzML open format.

## 3. XCMS Peak Picking and Alignment

> centWave algorithm for peak detection. Retention time alignment with obiwarp. Correspondence: match peaks across samples. Gap filling for missing values. Feature table: rows = features (m/z, RT), columns = samples.

```python
# Example: XCMS processing in R (via subprocess or rpy2)
# # Alternative: use ms-data-core-api or pyOpenMS in Python
# from pyopenms import *
# exp = MSExperiment()
# MzMLFile().load('sample.mzML', exp)
```

## 4. Normalization and Quality Control

> Probabilistic Quotient Normalization (PQN). Internal standard normalization. QC pooled sample injection for LOESS correction. Coefficient of variation (CV) threshold for feature filtering (< 30% in QC samples).

## 5. Exploratory Analysis

> PCA score plot (colored by condition, batch). OPLS-DA for supervised class separation. Coefficient of variation distribution. Missing value pattern (MCAR vs MNAR) and imputation strategies (kNN, min/2).
