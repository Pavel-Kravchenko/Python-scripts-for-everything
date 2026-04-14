---
name: bio-applied-lc-ms-preprocessing
description: "LC-MS metabolomics data preprocessing: peak picking, retention time alignment, gap filling, and adduct detection. Use when processing raw mass spectrometry data for metabolomics studies."
tool_type: python
primary_tool: Python
---

# LC-MS Metabolomics Data Preprocessing

- [XCMS documentation](https://www.bioconductor.org/packages/release/bioc/html/xcms.html)
- [MetaboAnalyst 6.0](https://www.metaboanalyst.ca/)
- [HMDB (Human Metabolome Database)](https://hmdb.ca/)
- [Galaxy Metabolomics training](https://training.galaxyproject.org/training-material/topics/metabolomics/)

## Metabolomics Overview

> Metabolome definition: small molecules < 1500 Da. Untargeted vs targeted vs semi-targeted approaches. LC-MS (liquid chromatography) vs GC-MS (volatile compounds, requires derivatization) vs NMR. Sample types: plasma, urine, tissue, cells.

## LC-MS Data Acquisition

> Chromatographic separation (RPLC, HILIC). ESI ionization: positive vs negative mode. Data-dependent acquisition (DDA) vs data-independent (SWATH/DIA). mzML open format.

## XCMS Peak Picking and Alignment

> centWave algorithm for peak detection. Retention time alignment with obiwarp. Correspondence: match peaks across samples. Gap filling for missing values. Feature table: rows = features (m/z, RT), columns = samples.


## Normalization and Quality Control

> Probabilistic Quotient Normalization (PQN). Internal standard normalization. QC pooled sample injection for LOESS correction. Coefficient of variation (CV) threshold for feature filtering (< 30% in QC samples).

## Exploratory Analysis

> PCA score plot (colored by condition, batch). OPLS-DA for supervised class separation. Coefficient of variation distribution. Missing value pattern (MCAR vs MNAR) and imputation strategies (kNN, min/2).

## Pitfalls

- **Mass accuracy drift**: Calibrate m/z values before peak picking; uncalibrated data produces false metabolite IDs
- **Adduct confusion**: The same metabolite produces [M+H]+, [M+Na]+, [M+K]+ ions — group adducts before statistical analysis
- **Missing value handling**: Zeros in LC-MS data may mean below detection limit, not absent — use min/2 or KNN imputation, not zero-fill
