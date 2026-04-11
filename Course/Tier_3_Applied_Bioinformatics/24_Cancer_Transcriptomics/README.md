# Module 24: Cancer Transcriptomics: Subtype Classification

**Tier 3 — Applied Bioinformatics | Module 24**

Classify melanoma tumors into molecular subtypes using variance filtering, hierarchical clustering, semi-supervised Random Forest, and Kaplan–Meier survival analysis on TCGA-SKCM-style expression data.

## What You'll Learn

- Loading and preprocessing cBioPortal RNA-seq expression matrices (log1p + z-score)
- Variance filtering to select high-information genes
- Unsupervised exploration: hierarchically clustered heatmaps (seaborn clustermap), PCA, t-SNE
- Semi-supervised classification: training a Random Forest on labeled samples to assign subtypes to unlabeled tumors (Tirosh framework)
- Hierarchical clustering with Ward linkage to define four Harbst subtypes from marker gene expression
- Kaplan–Meier survival curves using `lifelines` (with manual fallback)
- Quantifying agreement between classification schemes using normalized mutual information (NMI) and cross-tabulation

## Prerequisites

- Module 03 (RNA-seq Analysis) — expression quantification, count matrices
- Module 06 (Statistics for Bioinformatics) — normalization, variance, z-scores
- Module 07 (Machine Learning for Biology) — Random Forest, train/test split

## Data

Synthetic melanoma expression data (200 samples × 5,000 genes) generated to mirror the TCGA-SKCM cohort structure available on [cBioPortal](https://www.cbioportal.org/study/summary?id=skcm_tcga). Real data can be downloaded directly from cBioPortal in tab-delimited format.

## Related Skills

- `cancer-transcriptomics.md` — variance filtering, clustering, KM survival, subtype comparison patterns
- `ml-deep-learning-bio.md` — Random Forest and other classifiers for biological data
- `rnaseq-metagenomics.md` — upstream RNA-seq processing pipeline

---

[← Previous: 3.23 TF Footprinting](../23_TF_Footprinting/23_tf_footprinting.ipynb) | [Course Overview](../../README.md)
